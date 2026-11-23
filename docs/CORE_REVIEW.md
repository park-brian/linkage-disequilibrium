# Core review: running LDlink locally, WASM vs JS, and a smaller core

2026-11-22. Three things:
- LDlink is now runnable locally, so it can be compared with our code module by module.
- A fair measurement of what WASM buys over plain JS.
- A size and elegance review of `experiments/`, with the `src/` core it should become.

## 1. LDlink running locally (`D:\Projects\LDlink\ldlink-local`)

| piece | how |
|---|---|
| backend | LDlink's own `docker/backend.dockerfile`, built from the unmodified `nci-webtools-dceg-linkage` checkout (htslib 1.21, Python 3.13, mod_wsgi). Exposed on `127.0.0.1:8090` (`/LDlinkRest/...`). |
| code | A copy of `server/` mounted over `/opt/ldlink`, changed by `patch_local_vcf.py` in exactly **6 places / 11 lines**: `"s3://%s/%s" % (bucket, path)` → `"/data/%s" % path` (`LDcommon` ×1, `LDassoc`, `LDassoc_sub`, `LDexpress_ld_sub`, `LDproxy_sub`), and `checkS3File()` returns `True`. Nothing else differs (`diff -r`). |
| genotypes | `data/ldlink/1000G/<build>/`: NTFS hardlinks to the local 1000G VCFs and `.tbi`. They use no extra disk, and `db/` is untouched. |
| database | `mongo:7`. The 2023 mongodump is restored from S3 with `--maintainInsertionOrder`, so natural order = dump order, the order `find_one` sees. `gwas_catalog` comes from the 2025-10 backup (the same source as our release) instead of the 2023 dump. `api_log`, `*_tmp` and `*_bak` are skipped. |
| config | `server.env`: local-only placeholder credentials (a Mongo user that exists only in this container; dummy AWS keys, because `LDcommon.tabix()` always asks boto3 for credentials). |

**Why the patch is needed.** The dev bucket `s3://nci-cbiit-dceg-data-nonprod/ldlink/1000G/` has the VCFs but **no `.tbi` indexes** (0 for all three builds). LDlink's `tabix -fhD s3://…` cannot work against that bucket as it stands.

**Status:**
- The backend image is built.
- The restore is running. dbSNP's 22.9 GB at about 50 MB/min means about 8 h, plus index builds.

**Comparison plan (once the restore is done).** Drive LDlink's REST API with the same inputs as our code, and diff the results.
- **LDmatrix:** r² and D′ files for the test vectors, plus indels, multi-allelic sites, chrX males, monomorphic pairs, merged rsIDs and coordinate inputs.
- **LDpair, LDpop, LDhap:** the same inputs.
- **LDproxy:** proxies with r², D′, MAF, distance, FORGEdb, RegulomeDB and function.
- **LDtrait, LDexpress:** GWAS and GTEx hits.
- **SNPclip, SNPchip, LDassoc:** using the bundled example.
- **Also check:** warnings and errors verbatim, including the §15 cases where LDlink returns LD for a variant that isn't the rsID's.

This answers the open REQUIREMENTS §5 item ("not yet cross-checked against LDlink") without a token, and replaces our reimplemented ground truth with LDlink itself.

## 2. Is WASM an improvement over plain JS?

The earlier "5–35×" figure compared blocked WASM with the step-1 JS code: a 16-bit lookup-table popcount, unblocked and unpruned, with stats computed per pair. That baseline was never optimized.

The fair comparison is `wasm/ld_kernel_js.mjs`, a line-for-line plain-JS mirror of `band()` and `band_blocked()`. It uses the same packed `Uint32Array` from `prepare()`, SWAR popcount, the same r² bound pruning, the same 4-row register tiling and the same output order. Its output is **identical** to the WASM kernel's: same pairs, same order, same float32 values. See `wasm/bench_js_vs_wasm.mjs`.

**Throughput.** GRCh37 chr22 common, 2 Mb window, 500 kb band, 42.9 M pairs, single thread, best of 3 (the Phase B build was running on 3 of 32 threads):

| population | threshold | WASM band | **WASM blocked** | WASM SIMD blocked | JS band | **JS blocked** | blocked speedup |
|---|---|---:|---:|---:|---:|---:|---:|
| CEU (198 haps) | 0 | 98 | **140** | 120 | 36 | **40** | **3.5×** |
| EUR (1,006) | 0 | 65 | **82** | 71 | 15 | **18** | **4.6×** |
| ALL (5,008) | 0 | 22 | **24** | 25 | 3.9 | **4.5** | **5.3×** |
| CEU | r² ≥ 0.8 | 7.2 | 12.7 | 12.2 | 5.3 | 7.4 | 1.7× |
| EUR | r² ≥ 0.8 | 7.5 | 13.0 | 12.3 | 4.6 | 6.9 | 1.9× |
| ALL | r² ≥ 0.8 | 7.7 | 10.9 | 12.2 | 2.9 | 3.5 | 3.1× |

(Mpairs/s of evaluated pairs. In the pruned rows, both spend most of their time in the per-pair bound test.)

**Per request, at the sizes LDlink's modules use** (ALL population, all pairs, prepare + kernel):

| request | pairs | WASM | plain JS |
|---|---:|---:|---:|
| 30 SNPs | 435 | 0.07 ms | 0.37 ms |
| 300 SNPs (LDmatrix web limit) | 44,850 | 2 ms | 10 ms |
| 3,000 SNPs (API limit) | 4.5 M | 0.35 s | 1.0 s |

WASM costs a 6.4 KB binary and 0.6 ms to instantiate.

**Conclusions:**
- **WASM is a real but modest gain: 3.5–5.3× unthresholded, 1.7–3.1× pruned.** It is 5–35× only against unoptimized JS.
- **It matters only at pyramid scale** (multi-Mb windows, large populations). For every LDlink module's own limits, plain JS takes milliseconds.
- **SIMD never beats scalar** in V8 (except slightly for ALL with pruning). Drop the SIMD build.
- **C `band()` is `band_blocked(B = 1)`**, apart from output order, so one blocked entry point is enough.

**Recommendation:**
- Keep **one** kernel API with two interchangeable backends:
  - **WASM scalar blocked** (`band_blocked`, `band_masked`, `tile`, `pair_counts`);
  - the **JS mirror**, as the test oracle and the fallback where WASM is unavailable.
- Remove the SIMD build and the unblocked `band()`.
- If simplicity matters more than about 4× on large windows, plain JS alone is a defensible choice: no toolchain, about 90 lines, identical results.

## 3. Code review: elegance and size

An independent read-only review. I verified the two most consequential items: the browser blocker and the planner's region reads.

**Inventory:**

| category | lines | contents |
|---|---:|---|
| runtime core (browser) | ~1,590 | `lookup`, `resolve`, `planner`, `select_rows`, `pyfmt`, `wasm/ld_kernel.{mjs,c}`, `ldsc/{ldsc_regress,munge,ldscore_l2}` |
| data pipeline (offline) | ~1,600 | dbSNP, NCBI, genotypes (`build_genotypes_v2`, `rsid_label`, `build_all.sh`), annotations, LD scores, manifest |
| verification and benchmarks | ~2,300 | `verify_*`, `kernel_check`, `parity_inputs.py`, benches, `test/phaseA.test.mjs` |
| obsolete | ~550 + snapshots | `build_genotypes.mjs`, `build_dbsnp.mjs` and their benches, `prior-work/`, `out/phaseA_*` copies |

**Must fix for the browser.** `ldsc/ldscore_l2.mjs:117` reads `process.argv` at import time, so importing it in a browser throws. Move the CLI into its own file.

**Main simplifications:**
1. **`lookup.mjs`:** `sortedLookup` and `rangeScan` duplicate about 95 lines of page pruning, and return different shapes. Merge them into one `pageRanges(meta, column, predicate)` + `readRanges()`, with one return shape.
2. **`planner.query` re-implements `resolveInputs`** only to inject its hybrid rsID path. Instead, call `dbsnp.inputs(tokens, build, { rsids: hybridRsids })`. `planner.query` (80 lines, six jobs) then shrinks to cost planning plus execution.
3. **Use numbers, not BigInt.** Genome keys (< 2.6e10) and rsIDs (< 2^31) are exact as JS numbers. Convert INT64 once at read time and drop the BigInt arithmetic and the `String()` Map keys.
4. **One `genome.mjs` for chromosome and build tables.** There are 6 copies today; `planner.mjs:19` is unused.
5. **One rsID parser.** `resolve.mjs` and `planner.mjs` use different regexes.
6. **Region reads:** `planner.readRange` uses hyparquet's `filter`, which prunes row groups only. Use the page-pruned `rangeScan`, which `select_rows` already uses.
7. **Kernel:**
   - `prepare()` calls `reset()`, so only one prepared set can be live at a time, and the API doesn't say so.
   - It also decodes the Parquet sparse format (`rowBits`), which is a storage concern. Move that to `genotypes.decodeHaps`.
   - Measure whether the identity fast-copy path, where R1 lived, is worth keeping.
8. **ldsc:**
   - Move the numerics (`erfc`, `normIsf`, `median`, `npSum`, `quantile`) and the allele sets into `ldsc/stats.mjs`. Then `munge` stops importing the regression module and the l2 CLI.
   - Split `ldscoreL2` (decode, standardize, window, block sum), and decode the BED once instead of twice.
   - Delete the dead `Mtot`, `Math.max(x)` and an unreachable throw.
9. **Small duplicates:**
   - `minimalForm` (3 copies; `rsid_label.mjs` should be the only one)
   - `median` (2), `pop32`/`stats` (2)
   - `countingFile` (5 bench scripts), `readParquet` (2)
   - the Map-of-arrays idiom (8 sites; use a 3-line `groupBy`)

**Target `src/` (about 1,400 lines, down from about 1,710 including duplicates):**

| module | responsibility | ~lines |
|---|---|---:|
| `source.mjs` | file access and footer cache, manifest-aware `exists()` (the only test for "absent") | 30 |
| `parquet.mjs` | page-index reads: `lookupKeys`, `scanRange`, `pool` | 100 |
| `genome.mjs` | chromosome codes and names, keys, build paths, `firstToken`, `parseCoord`, `parseRsid`, `groupBy` | 50 |
| `dbsnp.mjs` | LDlink-parity resolution: `rsids`, `coords`, `inputs` (pluggable rsID strategy) | 85 |
| `genotypes.mjs` | tier files, row fetch, `decodeHaps`, `pickParity` / `pickAllele` (`select_rows` merged in) | 110 |
| `planner.mjs` | cost model plus hybrid rsID strategy only | 130 |
| `annotations.mjs` | track reads: `findOne`, `range`, `overlaps` (missing today; the logic is only in `verify_annotations`) | 55 |
| `ld/kernel.{c,mjs}` + `ld/kernel_js.mjs` | one API, WASM and JS backends | 165 + 130 + 90 |
| `pyfmt.mjs` | Python formatting | 29 |
| `ldsc/{stats,regress,munge,l2}.mjs` | LDscore family, no CLI code | 110 + 180 + 130 + 95 |

**Pipeline:**
- **Canonical scripts:** `build_genotypes_v2`, `build_dbsnp_v2`, `build_all.sh`, `build_annotations`, `convert_ldscore_ref`, `build_manifest`, `check_release`.
- **Retire:** v1 builders, the page-size sweep (`FINAL` becomes the only mode), and the `_v2` suffixes.
- **Entry point:** one `pipeline/build.sh <stage…>` with idempotent stages: `dbsnp`, `ncbi`, `genotypes`, `annotations`, `ldscore`, `release`.

## 4. Proposed order

1. **Phase B keeps running** (machine time). When it finishes: `check_release`, spot verification, and upload to `v3/`.
2. **The LDlink comparison**, when the restore is done:
   - Write `compare_ldlink.mjs`: the module × case matrix above, LDlink REST vs our code, results diffed.
   - Record every behaviour we don't reproduce, as a fix or as a documented deliberate difference.
3. **The `src/` refactor (§3)**, driven by the existing tests, verifiers and the LDlink comparison as the safety net:
   - `genome`, `source` and `parquet` first;
   - then `dbsnp`, `genotypes`, `planner`, the kernel API and `ldsc`.
   - Every verifier must pass unchanged after each step.
4. **Then Phase C** (the app) on `src/`, starting with the LDmatrix slice.

## 5. Decision: WASM + JS fallback (implemented 2026-11-22)

**One kernel API, `wasm/ld_kernel.mjs`, with two backends that give identical results:**
- `loadKernel(bytes)` returns the WASM backend. It falls back to JS when no bytes are given or WebAssembly is
  unavailable.
- `jsKernel()` gives the JS backend directly.

**What the backends share:** a front end with `prepare`/pack, output sizing, `strict`, `exact` and the pair bound.

**What differs:** only five engine primitives: `rowCounts`, `bandBlocked`, `bandMasked`, `tile`, `pairCounts`.

**Changes:**
- **C kernel:** scalar only, no SIMD. `band()` is gone; `band_blocked` covers it, and the default block is 128.
  `band_masked` shares `put()`. Rows are `ceil(nh/64)` words, no longer padded to 128-bit vectors. A single
  `ld_kernel.wasm` (5.7 KB) replaces the scalar and SIMD files.
- **`prepare()` no longer silently invalidates earlier sets.** WASM memory is a bump arena freed by an explicit
  `reset()`. Views are taken after any allocation, because `memory.grow` detaches them. The new tests caught exactly
  that mistake in `exact()`.
- **Retired:** `bench_kernel.mjs` and `bench_blocked.mjs` (SIMD/unblocked comparisons; the correctness part is
  superseded by `kernel_check` against LDlink).

**Verification:**
- `kernel_check` runs both backends against LDlink's math on 5 windows × 6 populations × every path: 0 failures, and
  WASM and JS identical in all 30 runs.
- `test/phaseA.test.mjs` runs the kernel tests on both backends: 23/23.
- `bench_js_vs_wasm` passes outputs identical in every case.

**Speed under load** (Phase B build + Mongo restore running): WASM is 5.1–6.6× faster unthresholded and 1.8–4.0×
pruned. An idle machine gave 3.5–5.3× and 1.7–3.1× earlier. Re-measure when idle.

## 6. First result against real LDlink: LDmatrix is identical (2026-11-22)

**Local LDlink, two fixes to the setup:**
- **Indexes:** dbSNP's index builds must finish before querying. The last one ran after the data load, and until it
  finished `find_one({"id"})` was a full collection scan.
- **Working directories:** LDlink runs `tabix` with its working directory set to `/data/1000G/<build>`, so those
  (empty) directories must exist.

After that, LDmatrix answers in under 1 s.

**Comparison:** `experiments/compare_ldmatrix.mjs` runs the same inputs through LDlink's REST API and through our
pipeline, and diffs the SNP list, the order and every cell.
- **Our pipeline:**
  1. `resolveInputs` (via the planner)
  2. parity row selection
  3. WASM `exact()`
  4. Python `round(·, 3)`
  5. LDmatrix's module rules: input gate, NA positions, not in 1000G, not biallelic, position order
- **Result: 30 / 30 matrices identical** (15 cases × r² and D′, about 13,700 cells each), with the same SNP lists and
  order in every case.
- **Cases** (`out/ldlink_compare/ldmatrix_cases.json`):
  - the test vectors on GRCh37, GRCh38 and High Coverage
  - a chr22 window in EUR, YRI, ALL and CEU+YRI
  - indel labels found at POS+1, plus the allele-mismatch rsIDs from §15, where LDlink computes LD on an
    overlapping variant (we reproduce that in parity mode)
  - coordinate inputs, including malformed and upper-case ones
  - chrX across PAR1 and outside it (haploid males)
  - GRCh38 chr7 and High Coverage chr2
  - chr1 around the 12-ALT site (dropped as not biallelic, as LDlink does)
  - chr11 with a merged rsID
- **What it settles:**
  - The rest of the chain matches the real service, not only our reading of the code: the input chain, `find_one`
    order, tabix-overlap row choice, the haploid/missing handling, and the 3-decimal output.
  - The REQUIREMENTS §5 cross-check, open until now, is closed for LDmatrix.
- **Facts learned on the way:**
  - LDmatrix recognises only the GT strings `0|0`, `0|1`, `1|0`, `1|1`, `0` and `1`; anything else makes both
    haplotypes missing. The only `.` calls in 1000G are haploid chrY, where both rules agree.
  - A `+` in a query string must be encoded (`CEU%2BYRI`), or LDlink sees a space and rejects the population.

**Next:** the same harness pattern for LDpair, LDpop, LDhap and LDproxy (the same engine, different outputs), then
LDtrait, LDexpress, SNPclip, SNPchip and LDassoc.

## 7. LDpair, LDpop and LDhap against real LDlink (2026-11-22)

**Ports.** `experiments/modules/` holds the module ports, on the core:
- `context.mjs`: planner, tier rows, samples, population haplotype indexes, dbSNP and annotation lookups.
- `common.mjs`: LDlink's shared pieces: `set_alleles`, the 2×2 statistics, `get_query_variant_c` and its
  warnings/errors.
- One file per module.

Each port returns LDlink's own output, byte for byte:
- **LDpair:** the API JSON.
- **LDpop:** the API's text table.
- **LDhap:** the full `calculate_hap` JSON (warnings included) and the API's text.

**Kernel addition:** `kernel.counts(p, I, J)` → `[n, alt_i, alt_j, both]` per pair, over the haplotypes called in both
rows. These are the 2×2 margins every module needs. `exact()` is now built on it.

**Harness:**
- `compare_ldpair.mjs <cases> <out> [ldpair|ldpop]`: REST.
- `compare_ldhap.mjs`: calls LDlink's Python functions directly in the backend container via
  `ldlink-local/drive.py`, so the comparison sees the warnings the API drops.
- Inputs:
  - every consecutive pair, or chunk of ≤ 30, from the 15 LDmatrix cases: indels, allele-mismatch rsIDs, coordinates,
    chrX PAR / non-PAR, a merged rsID, the 12-ALT site, three builds;
  - edge cases: invalid, missing and short inputs; not in dbSNP; bad and combined populations; different chromosomes;
    chrY on GRCh38 / GRCh37; over 30 inputs; blanks and duplicates; > 1 Mb range; monomorphic.

| module | cases | identical | other |
|---|---:|---:|---|
| LDpair | 429 pairs | **418** | 6 tie order only (LDlink is nondeterministic, below); 5 VCF ID (below) |
| LDpop | 537 pairs (134 with pop=ALL, 32 populations each) | **531** | 6 VCF ID |
| LDhap | 70 lists | **70** | — |

**LDlink behaviours the ports reproduce** (each is easy to get wrong):
- **LDpair:**
  - A..D follow the *sorted allele strings* (`sorted(hap)`), not REF/ALT.
  - χ² uses Python's exact integer numerator (BigInt in JS).
  - p is `2·(1 − ½(1 + erf(√χ²/√2)))`, formatted `str(round(·, 4))` or `"<0.0001"`.
  - Only the error is returned for a failed pair.
  - The "not a biallelic" error appends the *warning* text (sic).
- **LDpop:**
  - A haploid call (`"0"` → `"0|."`) counts as one REF and one ALT allele over two.
  - SNP2's frequencies are divided by SNP1's allele total.
  - Frequencies are keyed by the raw REF/ALT strings.
  - Rows are in `pop_order`; a super-population adds its sub-populations, and ALL adds all 31.
- **LDhap:**
  - Tokens are deduplicated before coordinates are replaced, so a coordinate and its rsID both survive and collapse
    later in `parse_vcf`.
  - The same-chromosome error is raised after `parse_vcf`.
  - A haplotype containing any haploid or uncalled allele is dropped whole (chrX males lose their second
    haplotype).
  - Output order is count descending, ties by haplotype string.
  - Frequencies are `round(·, 4)`, and the allele text is `round(round(f, 4), 3)`.
  - A multi-allelic SNP is a warning when other SNPs remain, else an error.

**LDlink is nondeterministic in LDpair's haplotype order.** `get_population` returns `list(set(ids))`. String
hashing is randomized per process, and mod_wsgi runs `(nproc+1)/2` processes, so the order of equal-count haplotypes
(`hap1..hap4`) follows a different random sample order per process. The same request run twice gave 10, then 6, pairs
with ties reordered. Our port breaks ties in key order; the harness reports these as `tie_order`.

**The one real gap: the raw VCF ID.** LDlink prints the VCF's own ID in several places:
- in "Genomic position does not match RS number at 1000G position (chr:pos = **rs…**)" (LDpair, LDpop, LDproxy,
  LDmatrix warnings);
- as LDproxy's RS_Number, and as the key of LDproxy's Function and FORGEdb lookups.

Our rows keep the ID only when it is the label (`rsid_source` 1). On GRCh37 chr22, 0.65% of rows (7,224 of 1.1 M) have
an rsID we cannot reproduce:
- a merged VCF rsID (4,863);
- an rsID absent at that position in the dump (2,361).

GRCh38 VCFs have no rsIDs: their IDs are `.` or `chr:pos:ref:alt`, so LDlink falls back to the first dbSNP doc at the
position, which we reproduce. **Fix:**
- `build_genotypes_v2.mjs` writes `vcf_id`, the raw ID only when it contains "rs" and is not `rs<rsid>`, otherwise
  null.
- `common.vcfId(row)` uses it.
- The GRCh37 rebuild is in `out/phaseB/geno_vcfid/`.

### LDproxy (`modules/ldproxy.mjs`, `compare_ldproxy.mjs`)

**Comparison:**
- **Ours:** the API text, line by line, plus the warning.
- **LDlink:** the API text, plus `warning` read from LDlink's own `proxy<ref>.json`. We pass `reference=<uuid>` and read
  the shared tmp directory.
- **Queries:** 57: the first three inputs of every LDmatrix case, plus `r2_d=d`, a 100 kb window, bad `r2_d`/window,
  not in dbSNP, a coordinate, a bad population, multi-allelic, chrY on GRCh37/GRCh38, ALL.

**Result on the old release: 38 / 57 identical**, including every GRCh38 and High Coverage query (1,200–2,200 lines
each). The other 19 are all GRCh37, and differ only on rows whose raw VCF rsID we did not store (§7 above): 1–41 lines
out of 90–1,700. With `vcf_id` (release_v4, below), all 57 are identical.

**LDproxy mechanics reproduced:**
- **The window is cut into `NUM_SUBPROCESSES` chunks** (`chunkWindow`, here 4), each read by tabix with overlap
  semantics:
  - a variant spanning a chunk boundary is listed twice;
  - the last chunk ends 3 bp past `pos + window`.
- **Rows:** every biallelic row with Ms ≠ 0.
- **Sort:** stable, |distance| ascending, then the statistic descending.
- **The first line is the top row:** RS, distance and statistics come from that row, but Coord and Alleles from the
  query.
- **Proxies:** R² > 0.01 and RS ≠ the query label, where the label becomes the VCF ID after an RS-mismatch warning.
  That warning is appended twice (sic).
- **RS_Number:** the VCF ID when it starts with "rs", else the first dbSNP doc at the VCF position.
- **Function:** `find_one` by that ID, or the position doc; the dump's `"NA"` is null in our index.
- **FORGEdb:** by RS, empty → `NA`.
- **RegulomeDB:** at the VCF position.
- **MAF:** the proxy's minor allele frequency over the haplotypes called in both.
- **Distance:** measured from the query's *VCF* position.
- **The REST layer validates `r2_d` ∈ {r2, d} and `window` ∈ [0, 1e6] first** ("Invalid r2_d parameter." /
  "Invalid window parameter."). The module's own `r2_d` fallback warning is unreachable through the API.

**Context additions:** `windowRows` (tabix-overlap rows across tiers in VCF order), `dbsnpAt` / `dbsnpById` (LDlink's
two `find_one` forms), `forge`, `regulome`.

**Local LDlink fix:** the image provides `python3` as a symlink in `/opt/ldlink/apache-bin`, which our mount of
`server/` over `/opt/ldlink` hides. LDproxy's subprocesses (`python3 LDproxy_sub.py`) then fail with "Internal server
error". `server/apache-bin/python3 → /usr/bin/python3.13` restores it.

### Result with `vcf_id`: every comparison identical (2026-11-22)

**The GRCh37 rebuild** (`OUT=out/phaseB/geno_vcfid ./build_all.sh 3 grch37`):
- 24 chromosomes; `check_release` reports 0 failures.
- 85.2 M rows, the same as before.
- `vcf_id` is set on 7,224 chr22 rows, exactly the gap measured above.
- Files grow by about 0.15%.

It took about 5 h because the machine was under memory and disk pressure (another disk scanner was running); Phase B
took 18 min for chr1.

**Release root:** `experiments/out/release_v4`:
- GRCh37 genotypes from the rebuild;
- everything else hardlinked to `out/release`;
- its own `manifest.json` (324 files).

| comparison (ROOT=out/release_v4) | identical |
|---|---|
| LDmatrix | **30 / 30** matrices |
| LDpair | **420 / 429** pairs; the other 9 differ only in tie order (LDlink randomizes it) |
| LDpop | **537 / 537** pairs |
| LDhap | **70 / 70** lists (JSON with warnings, and API text) |
| LDproxy | **57 / 57** queries, every line |

**Every number, label, warning and error the five modules print now matches LDlink.** The only difference left is
nondeterministic on LDlink's side.

## 8. SNPclip, SNPchip, LDtrait, LDexpress and LDassoc against real LDlink (2026-11-22)

**Same pattern as §7.** One port per module in `experiments/modules/` and one comparer per module
(`compare_<module>.mjs`), all on `out/release_v4`, each diffing LDlink's API text line by line, plus errors.
- **REST:** SNPclip, SNPchip, LDtrait, LDexpress.
- **LDassoc:** has no API (the route says programmatic access is not available), so its comparer calls
  `calculate_assoc` in the container through `ldlink-local/drive.py`, on the bundled prostate example.

**`context.mjs` gained:** `platforms`/`arrays` (SNPchip), `gwas` (LDtrait), `eqtl` (LDexpress), `gene` (LDassoc),
built on the same `sortedLookup`/`rangeScan` reads.

**`pyfmt.mjs` gained `cFormatG`** (C `%g`). Every one of the 143 M GTEx v8 `slope`/`pval_nominal` strings
round-trips through `%.6g`, so the doubles we store print exactly as LDlink prints Mongo's strings.

| module | cases | identical |
|---|---|---|
| SNPclip | 49: every LDmatrix list × 3 r²/MAF threshold settings, a 1,500-rsID dense list, edges | **49** |
| SNPchip | 66: every LDmatrix list × 4 platform sets, 2,000-rsID lists on two builds, edges | **66** |
| LDtrait | 55: chunks of every LDmatrix list, GWAS-hit lists by threshold / D′ / window / population / build, edges | **55** |
| LDexpress | 58: chunks of every LDmatrix list × 5 tissues, thresholds, all 54 tissues, build variants, edges | **58** (one case re-checked after a fix) |
| LDassoc | 21: variant / gene / region on three builds, origins and lowest-p origins, chrX, errors | **21** |

**SNPclip mechanics reproduced:**
- **Input filter:** the route keeps only lines starting `rs`/`chr`.
- **MAF:** counted from genotype strings (a haploid call is one allele); kept when MAF ≥ threshold, unrounded.
- **Greedy pruning in input order:** a kept variant removes every later variant with r² ≥ threshold.
- **Pairs:** r² pairs genotype strings individual by individual, so a chrX male who is diploid in PAR1 and haploid
  outside it is skipped for that pair. The port prepares a second set without the haploid males for exactly those
  pairs.
- **Multi-allelic records:** printed with their raw REF and ALT list (`A=NA, G,T=NA`).
- **`maf ≥ 0`:** a surviving monomorphic variant makes LDlink compare `None ≥ threshold` and fail with "Internal
  server error". The port does the same.

**SNPchip mechanics reproduced:**
- **Order:** sorted by (chromosome number, position); positions are written `22:123`.
- **Warnings** start with ". " because an empty warning is joined.
- **Platform matching:** by name or by `name.rstrip(' Array')`, which strips the *characters* " Aray"; the
  " Array" suffix is removed for printing.
- **An empty platform request crashes LDlink** (`append(...).removesuffix` on None) whenever any position has arrays.

**LDassoc mechanics reproduced:**
- **Query row:** the origin genotype is the *first* row overlapping its position (`snp_no_dups` line 0), while
  `get_query_variant_c` checks the *last*.
- **Tabix regions:** one per association position, in file order. A row is listed under every position it overlaps,
  and kept when its own POS is in the file (with that position's last p-value).
- **Sort:** |distance|, then p ascending.
- **Correlated alleles:** written proxy-first, with the allele indexes crossed.
- **Coordinates:** the region-mode check compares coordinate *strings*.
- **A coordinate origin raises KeyError.**
- **High Coverage lowest-p search:** it asks tabix for `chrchr8:…`, so it always reports the variants as not in
  1000G.

**LDtrait mechanics reproduced:**
- **GWAS records** come in index order: position, then insertion, which is our file's `(key, dump order)`.
- **Pair LD** comes from one tabix row per position (the last overlapping row). It is NA, and so never matched, when
  that row is multi-allelic or its POS differs from the queried position (an indel listed at POS+1).
- **R2 / D′** are printed unrounded (`str(float)`).
- **Columns are shifted against their headers:**
  - "Risk Allele" holds the RAF;
  - "Effect Size" holds `float(OR or BETA)`;
  - "Beta or OR" holds the CI with brackets, letters and spaces removed.

**LDexpress mechanics reproduced:**
- **Window:** LDproxy's chunked window, r² (or D′) ≥ threshold. The RS ID is the VCF rsID, or the first dbSNP doc at
  the VCF position.
- **eQTLs:** in index order, inner-joined to `gtex_tissues` and `gtex_genes`, kept when `pval_nominal` < threshold;
  values are printed `%.6g`.
- **A query that is not in 1000G deletes the request's population file** (`get_query_variant_c`'s `is_output`
  cleanup), so every later query's subprocesses fail and return nothing.
- **A multi-allelic query is computed on its GT digits** (`0|2`: only 0/1 haplotypes count). The port models it as
  ALT 1 with the haplotypes carrying ALT ≥ 2 uncalled.
- **Input handling:**
  - inputs are `"+".join(line.strip().lower())`;
  - blank entries count toward the limit of 10, then are dropped;
  - a 1–2 character entry fails the request ("Input list of RS numbers is empty").

**Every LDlink module is now reproduced on the core** (LDmatrix, LDpair, LDpop, LDhap, LDproxy, SNPclip, SNPchip,
LDtrait, LDexpress, LDassoc), on all three builds, against the real service. The only differences left are LDlink's
own nondeterminism (LDpair tie order).

**Not covered:**
- the LDscore family (§5 / Phase A ported it separately against ldsc);
- the plots (Bokeh);
- web-only paths behind the internal auth header, such as LDtrait's "too many pairs" prompt.

## 9. The `src/` refactor (2026-11-22)

**The runtime core now lives in `src/`** (layout and line counts in `src/README.md`, 2,728 lines including the ten
LDlink module ports and the C kernel). `experiments/` keeps the pipeline, verifiers, benchmarks and comparers, which
import from `../src`.

**The move:**
- `scripts/move_to_src.py` moved the files and rewrote every relative import across both trees, resolving each import
  to its pre-move target.
- `node_modules` and `package.json` moved to the repo root, so both trees resolve `hyparquet`.
- `ld/build.sh` rebuilds `kernel.wasm` byte-identical from its new place.
- Before the move, the code was snapshotted in `../backups/experiments_code_before_src_refactor.tgz`.

**§3 items done:**
- **One `genome.mjs`.** It holds chromosome codes and names, keys and the `BUILDS` table: genotype directory,
  position build, dbSNP key column, title, hg label, 1000G prefix. It replaces the copies in `dbsnp`, `planner`,
  `context`, `ldhap` and `ldassoc` (including the unused `CHROM_NAME` in the planner). `firstToken` and `parseCoord`
  live here; `isNotFound` moved to `source.mjs`.
- **`source.mjs`** is the one cached file opener; `dbsnp.createResolver` uses it.
- **The planner exposes its file cache** (`planner.file`, footer lengths from the manifest) and its resolver. The
  module context reuses both instead of keeping a second cache.
- **`genotypes.mjs`** gained `decodeHaps` (moved out of the kernel: storage, not LD) and the tier-row reads
  (`createGenotypes`).
- **`annotations.mjs`** holds the track lookups that were in the module context; the dbSNP `find_one` forms
  (`findAt`, `findById`) moved to `dbsnp.mjs`. `ldlink/context.mjs` is now 29 lines of composition.
- **ldsc:**
  - `ldsc/stats.mjs` holds `erfc`, `normIsf`, `normPdf`, `chi2sf1`, `median` and the allele sets, so `munge` no
    longer imports the regression module or `l2`.
  - `l2.mjs` has no CLI code, so it is browser-safe (the §3 "must fix"); the CLI is `experiments/ldsc/ldscore_l2.mjs`.
  - The duplicate `median` is gone.
- **The kernel module exports `KERNEL_WASM`**, the binary's URL next to it, for both Node and browser loading.
- **The LDmatrix logic** moved out of `compare_ldmatrix.mjs` into `src/ldlink/ldmatrix.mjs`, so all ten ports are in
  `src/`.

**Not done (deliberately, with reasons):**
- **Merging `sortedLookup` / `rangeScan` into one page-range primitive (§3.1).** It is the hottest path and the
  planner's byte predictions are calibrated against it. Worth doing with `verify_planner`'s measured bytes as the
  check, as its own change.
- **Numbers instead of BigInt for genome keys (§3.3).** It touches every key comparison, including `sortedLookup`'s
  key ordering, and is also a separate change.
- **Splitting `planner.query`'s six jobs (§3.2).**

**Safety net:** `experiments/run_regression.sh` runs the unit tests, the data verifiers and all ten LDlink
comparisons. Its result on the refactored tree is recorded below.

**Result on the refactored tree** (`./run_regression.sh out/release_v4`, about 1.5 h; logs in
`experiments/out/regression/`): every step passes.

| step | result |
|---|---|
| unit tests | 24 / 24 |
| verify_planner | hybrid == direct, select == select_rows, non-canonical inputs == resolver: 0 differences |
| verify_parity | 0 mismatches in all 12 checks (rsID, coordinate, test-vector and input-chain parity, both builds) |
| verify_annotations | 0 mismatches |
| kernel_check (GRCh37 chr1) | 0 failures, WASM and JS identical |
| verify_select (GRCh37 chr22, vcf_id rebuild) | parity 22,000 / 22,000, allele 22,000 / 22,000 |
| ldsc munge / regress | worst relative difference 3e-14 against ldsc |
| LDmatrix · LDpair · LDpop · LDhap · LDproxy | 30/30 · 419/429 + 10 tie order · 537/537 · 70/70 · 57/57 |
| SNPclip · SNPchip · LDtrait · LDexpress · LDassoc | 49/49 · 66/66 · 55/55 · 58/58 · 21/21 |

## 10. Kernel: row ranges on the masked path (2026-11-23)

**The change:** `band(p, { rows: [r0, r1) })` now works on masked sets (rows with haploid or uncalled haplotypes) as
well. `band_masked` (C) and `bandMasked` (JS) take `iStart`/`iEnd`, and the front end no longer rejects ranges there.
This lets the app's worker pool split chrX regions.

**Checks:**
- `kernel.wasm` rebuilt (5,684 bytes).
- New unit test K4: on both backends, any partition of the masked band into row ranges equals the full band, pair
  for pair. 25/25 tests pass.
- `kernel_check` on GRCh37 chr1 and chrX (42 masked rows): 0 failures.
