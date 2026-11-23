# Roadmap: from verified experiments to the hyparquet LDlink app

Status as of 2026-11-21 (Phase A done; see §3). Companion docs:
- [REQUIREMENTS.md](REQUIREMENTS.md): what each module must do.
- [PARQUET_LAYOUT.md](PARQUET_LAYOUT.md): data layout and measured results, §1–16.
- [PLAN_REVIEW.md](PLAN_REVIEW.md): design findings, steps 1–9.
- [DATA.md](DATA.md): source inventory.

This doc says what is left, in what order, and what "done" means for each piece. It starts with a code review of
everything under `experiments/`.

## 1. Where things stand

**Done and verified** (steps 1–9 of PLAN_REVIEW.md):

| area | state | evidence |
|---|---|---|
| Genotype layout | MAF-tiered bitset/sparse Parquet, ZSTD, page index. `vcf_order` + allele columns added (§16) | LD = reference math on 4,800 pairs; v3 columns identical to v2 |
| dbSNP | Full rsID index (33 partitions) + position indexes for both builds, from LDlink's 2023 dump | 0 parity mismatches / 13,400 cases |
| dbSNP vs NCBI | b155 reproduces the dump (chr22 exactly; genome-wide 99.995%); b157 drift measured | §15 |
| Allele-aware mode | Opt-in row selection (`select_rows.mjs`, planner `select`) | 4 × 50,000 rsIDs vs independent ground truth |
| LD kernel | WASM popcount, blocked, r²-bound pruning, workers | exact vs reference (but see finding R1) |
| Planner | Scan vs point-read cost model; hybrid rsID resolution | hybrid == direct; cheaper choice in every case |
| Network | Private S3 bucket, presigned Range requests | windows < 1 s cold (§11) |
| Annotations | 20 files, 3.2 GB, position-keyed | 0 mismatches / 6,950 cases |
| LDscore family | `--l2`, munge, h²/rg ported to JS | byte-identical on LDlink's examples; ~1e-14 on h²/rg |

**Not started:**
- The app itself. The root `index.html`/`app.js`/`utils.js`/`sources.js` are the old prototype: base64-deflate
  genotypes, and the known missing-haplotype zero-fill bug in `utils.js`.
- Every chromosome other than chr22 (3 builds) and GRCh37 chrX.

## 2. Code review findings (2026-11-21)

Three independent reviews covered:
- **reader/query path:** `lookup`, `resolve`, `planner`, `select_rows`, the WASM kernel, the verifiers
- **builders and pipeline:** genotypes, dbSNP, annotations, manifest, NCBI comparison
- **LDscore ports**

Findings marked **✔** were re-checked by reproducing them or by querying the built data. The rest are as reported:
read against the code, not re-run.

### 2.1 Must fix before building more data or shipping

| # | where | defect | impact | fix | |
|---|---|---|---|---|---|
| R1 | `wasm/ld_kernel.mjs:48-50,64-65` (`prepare`) | Identity-copy path copies whole 64-bit words, but masks only the one 32-bit word that contains bit `nh`. Later words keep other samples' bits; with `nh % 32 == 0` nothing is masked | Wrong allele counts, r² and D′ whenever the population's haplotypes are a leading prefix of the sample order and `nh < 2·samples` (e.g. 64 of 200 haps → pa = 128) | Zero every word from `nh >> 5` (partial) through `2W − 1`, or copy only `ceil(nh/8)` bytes | ✔ |
| B1 | `build_genotypes_v2.mjs:119-126,259` (`labelRsid`) | rsID labels are looked up at the raw VCF POS. The dump stores deletions and anchored delins/insertions at **POS+1** (§15 rule) | GRCh37 chr22: only 4 of 27,831 deletion rows carry their allele-matching rsID; ~⅔ of insertions are unlabelled or mislabelled. Parity resolution is unaffected (a label is accepted only if the dump has it at that POS), but labels, hybrid coverage and the §15 (a) table are wrong for indels | For non-SNV classes, also search POS+1 (the `compare_ncbi.sql` rule) and use allele matches as the primary label when available. Rebuild, then re-run §15 (a) and §16 | ✔ |
| B2 | `build_genotypes_v2.mjs:115,259` (`cls`) | Type class is computed on padded alleles. A split multi-allelic row REF=AT ALT=GT (really an SNV) is classed MNV | 423 GRCh37 chr22 rows get no label or a wrong one | Classify on `minimalForm(pos, ref, alt)` | ✔ |
| L1 | `ldsc/munge.mjs:114` (`printCols`) | Output columns come from the file's translated names; an N derived from N_CAS/N_CON is never printed | Every case/control sumstats file loses N, so h²/rg drop all rows | Print `SNP, N, Z, A1, A2` in `dat` column order, with derived N appended last (as pandas does) | ✔ |
| L2 | `ldsc/munge.mjs:90` | `!("N" in dat[0])` guard: an existing N column wins over N_CAS/N_CON. In the reference, the derived N overwrites it | Wrong N and `n_min` filter when a file has N, N_CAS and N_CON | Drop the guard | |
| L3 | `ldsc/ldsc_regress.mjs:74` (`getSeparators`) | Exact `floor(i·N/nb)` vs numpy's `floor(linspace)`, which computes `i·(N/nb)` in floating point | Block boundaries shift by 1 SNP for ~4% of N (e.g. N = 228, 452: 3 boundaries); SEs and p change up to ~1% | Compute `Math.floor(i * (N / nb))` exactly as numpy does (check `linspace`'s last-point rule) | ✔ |
| L4 | `ldsc/ldscore_l2.mjs:17-25` (`pyFixed3`) | The tie test `Number.isInteger(x*2000)` is true after float rounding even for non-ties | 0.0005 → "0.000" (Python 0.001), 31.6765 → "31.676"; hits munge N output for N with ≥ 4 decimals | Decide the tie on the exact binary value (e.g. compare `x` with the decimal midpoint via `toPrecision(20)`, or BigInt-scaled exact arithmetic) | ✔ |
| P1 | `planner.mjs:180,185` | Coordinates are parsed with `Number(t[1])`; NaN or a fraction reaches `BigInt()` in `coordBytes` | One malformed `chr22:abc` or `chr22:12_345` input throws for the whole query; LDlink returns "no match" | Reuse the resolver's parsed coordinate | ✔ (by reading) |

### 2.2 Should fix (wrong in edge cases, or a verifier blind spot)

**Query path**
- **P2** `planner.mjs:161`: label-resolved rsIDs carry no `ref_rsid`. Low impact today, because accepted labels are
  always current IDs, but output differs from the direct path.
- **P3** `resolve.mjs:78-80` + planner: coordinate inputs resolve to a doc that may be a merged rsID, with no
  `ref_rsid`. Allele mode then searches under the old ID and wrongly reports `alleles_not_in_1000g`. Carry `ref_rsid`
  in the position index (or look it up).
- **P4** `planner.mjs:192-199`: overlapping user regions are scanned separately, so rows come back twice. `g.chrom`
  isn't normalized (`"chr22"`, `"x"`, `23` match nothing). Merge regions per chromosome; normalize chromosome names.
- **P5** `resolve.mjs:39,69`: `.catch(() => ({}))` turns any I/O error into "not found". A network failure looks like
  "not in dbSNP". Only treat "file absent from the manifest" as not found; rethrow everything else.
- **P6** `lookup.mjs` `sortedLookup`: returns nothing if `columns` omits the key column. `rangeScan` adds it;
  `sortedLookup` should too.

**Kernel**
- **K1** `wasm/ld_kernel.c:84,88`: thresholds compare unrounded r²/D′. LDlink rounds to 3 decimals first
  (LDmatrix.py:313-318), so r² = 0.0996 passes LDlink's 0.1 threshold but not ours. Add a parity mode that applies
  `round(·, 3)` before comparing, and widen the pruning margin accordingly.
- **K2** `wasm/ld_kernel.mjs:71`: `band()` defaults to `maxOut = 1 << 24`, which allocates 4 × 64 MB per call (WASM
  memory never shrinks). Size it from the estimated pair count, or stream the output.
- **K3** `ld_kernel.c` `band_blocked`: `B <= 0` loops forever, and an empty worker split `[0, 0]` means "all rows",
  so that worker recomputes the whole band and duplicates output.
- **K4 (parity note):** LDlink builds its 2×2 table from `sorted()` concatenated allele strings. For some indel allele
  pairs (e.g. {A, AC} × {A, T}) that order differs from index order, and LDlink's table differs. Document it as a
  deliberate difference, or reproduce it in parity mode.

**Builders**
- **B3** `build_annotations.mjs:79` vs `annotations.sql`: genes are sorted by (key, txEnd, ord), but the check
  expects (key, ord). Both `genes.parquet` report `sorted: false`, so the check is ignored. Use one order and check
  that one. ✔
- **B4** `build_annotations.mjs:90`, `build_dbsnp_v2.mjs`: an unsorted input is only reported, never rejected.
  Parity depends on order; throw instead.
- **B5** `json_array_stream.mjs` / `dbsnp_bson2tsv.mjs`: numeric-vs-string BSON types collapse to strings, while
  LDlink's Mongo queries are type-sensitive (numeric `$gte` for GWAS, `str(pos)` for dbSNP). Assert the `nonString`
  counts are 0, and keep a type histogram per field.
- **B6** `compare_ldlink_parse.sql`: `JOIN USING` never matches NULL to NULL (`type`, `fset`), and
  `row_number() OVER ()` relies on scan order. This only inflates mismatch counts, so the reported 0 stands, but fix
  it before reusing the script.
- **B7** `build_genotypes_v2.mjs:238`: allele indexes ≥ 10 (`"10|0"`) are parsed as single characters. Not present in
  the current inputs; add a guard (`alts.length > 9` → throw) before building all chromosomes.
- **B8** `build_genotypes_v2.mjs:53`: an unknown dbSNP `type` becomes 255, which is "compatible" with every indel.
  Map it to an explicit "unknown" class.
- **B9** `allele_join.mjs:73`: still breaks ties by tier; switch to `vcf_order` so the §15 (b) numbers are exact.
- **B10 (design):** `minimalForm` trims but doesn't left-align, so the same indel written at different places in a
  repeat won't match. Both inputs are left-normalized today; left-align properly (needs reference sequence) if other
  VCF sources are added.

**LDscore**
- **L5** `ldscore_l2.mjs:70-71`: `--ld-wind-* 0` is accepted. ldsc rejects it ("Must specify exactly one --ld-wind
  option"), and LDlink allows `ldwindow=0` through. Reject 0.
- **L6** `ldsc_regress.mjs`: no "More blocks than data points" error when the step-1 count is below `n_blocks`, so
  JS returns numbers where ldsc errors. Add the check.
- **L7** `ldsc_regress.mjs:40`: `erfc(±Infinity)` is NaN, so se = 0 gives p = NaN (reference: 0). `chi2Isf1(5e-324)`
  is NaN (scipy 1480.9). Handle the limits.
- **L8** `munge.mjs`: smaller divergences:
  - INFO-NA bucketing for a lower-case `info` header (log only)
  - duplicate-column checks (Python checks only the first translated field, and errors on two columns mapping to
    the same name)
  - number parsing: `inf`, hex, BOM, extra fields

**Verifier gaps** (checks that could pass while code is wrong). Close these before trusting more results:
- **V1** `ld_check.mjs` compares a JS re-implementation with the reference, **not the WASM kernel**, and skips
  indels, duplicate positions and multi-allelic sites. R1 slipped through for this reason. Run the kernel (all paths:
  scalar, blocked, pruned, masked, threads) against `ldmatrix_reference.mjs`, with population subsets including
  prefix-shaped ones.
- **V2** `verify_planner.mjs` compares hybrid vs direct only on found/pos/chrom, and `same_rows` compares counts.
  Compare full variant objects and row contents; exercise `select` and non-canonical inputs.
- **V3** `verify_parity.mjs` tests only canonical inputs. Add LDlink's input normalization cases (case, whitespace,
  `rs`/`chr` stripping, `_` in integers) and warning text.
- **V4** `verify_select.mjs` counts the SV case as agreement, and barely tests `not_in_1000g`. Report SV positions
  as a separate class; add a sample of rsIDs that are not near 1000G rows.

## 3. Remaining work, in order

### Phase A: correctness fixes (≈ 2–3 days)
1. ✅ **Fix R1, B1, B2, L1–L4, P1** (done 2026-11-21), plus B7, B8 and P2 along the way.
   - **Tests:** `experiments/test/phaseA.test.mjs` (`node --test test/phaseA.test.mjs`), 13 tests.
     - Reference fixtures come from the real `munge_sumstats.py` and numpy (`test/fixtures/`, `test/numpy_separators.json`).
     - The same checks run against the pre-fix code (`out/phaseA_before/check_before.mjs`) reproduce every failure:
       - **R1:** 64/200 haplotypes counted as 128.
       - **L1:** no N column. **L2:** N = 3000, not 2700. **L4:** N 2000.000, not 2000.001.
       - **L3:** separators 57, 114, 171 instead of 56, 113, 170.
       - **P1:** RangeError.
   - **Code changes:**
     - **R1:** `wasm/ld_kernel.mjs` `prepare` clears all bits ≥ nh.
     - **B1/B2/B8:** new `rsid_label.mjs` (used by `build_genotypes_v2.mjs`).
       - Labels are looked up where LDlink stores the variant: POS+1 for deletions/delins, POS or POS+1 for insertions,
         minimal POS for SNV/MNV.
       - The class is computed on the minimal form. NCBI allele matches are a new label source (4).
       - Unknown dbSNP types are compatible with nothing.
       - New column `rsid_pos_delta` (dump position − POS), used by the planner's hybrid path.
     - **B7:** more than 9 ALTs → error.
     - **L1/L2:** `ldsc/munge.mjs` follows `process_n` and pandas column order, with numpy pairwise mean.
     - **L3:** `getSeparators` mirrors `np.linspace`; equal to numpy for every N in 200–5000.
     - **L4:** `pyFixed3` decides on the exact decimal expansion; equal to Python on 45,018 values.
     - **P1:** `resolve.mjs` exports `parseCoord`, and the planner uses it.
     - **P2:** label hits carry `ref_rsid`.
   - **Rebuild** (`out/v4_*`, installed in the release and uploaded):
     - All non-label columns are identical to v3 in all 12 files.
     - Indel labels: deletions labelled with their allele match went 4 → 27,638 of 27,831 (GRCh37 chr22);
       insertions 1,652 → 16,120 of 17,035; padded split SNVs 0 → 421 of 423.
     - Label ≠ allele match: GRCh37 10,733 → 731 (VCF IDs that disagree with NCBI), GRCh38 9,278 → 0,
       High Coverage 54,482 → 1, chrX 38,853 → 2,783.
   - **Re-verified:**
     - `verify_select`: 4 × 50,000, both modes.
     - `verify_planner`: hybrid == direct (0 differences), cheaper plan chosen.
     - `verify_parity`: 0 / 13,400.
     - `--l2` byte-identical; munge identical; h²/rg ≤ 3.3e-14.
   - **Superseded by the rebuild:** §15 (a)'s label table and §16's "allele ≠ rsID label" column describe the old
     labels.
2. ✅ **Close verifier gaps V1–V4** (done 2026-11-21). The new checks found two more bugs (**R2**, **P7**), both
   fixed.
   - **V1** `kernel_check.mjs` runs the WASM kernel itself against LDlink's math (`ldmatrix_reference.mjs` on the VCF).
     - **Coverage:**
       - every path: scalar + SIMD × band, unpruned, blocked, blocked split into worker ranges, tile, masked
       - six populations: CEU, EUR, ALL, and leading-prefix sets with nh = 40, 64, 200
       - 42 SNPs per window from all three tiers, 14 of them indels
       - five windows: GRCh37/GRCh38/High Coverage chr22, GRCh37 chrX across PAR1, chrX non-PAR
     - **Result:** 0 failures. One ALL pair is the known K4 allele-sort difference, reported separately.
   - **R2 (found by V1, fixed):** hyparquet decodes `BYTE_ARRAY` columns **without a UTF8 annotation as text by
     default**. That corrupted the sparse-tier `haps` (55–63% of rows) and every `missing` bitset on read.
     - **Effect:** wrong LD for any rare-tier variant (e.g. a monomorphic SNP looked variable).
     - **Why it was missed:** `ld_check.mjs` read only the common tier, whose fixed-length `haps` aren't affected.
       No file on disk was wrong.
     - **Fix:** `lookup.mjs` (`sortedLookup`, `rangeScan`), `planner.mjs` and the bench/check scripts read with
       `utf8: false`. Only genotype `haps`/`missing` lack the annotation (scan of all 126 release files); annotated
       strings stay strings.
     - **Test:** `test/phaseA.test.mjs`.
   - **V2** `verify_planner.mjs` now checks four things:
     - hybrid == direct on **every field** (rsid, ref_rsid, chrom, pos, warning)
     - scan vs points on row **contents** (pos, tier, sha1 of haps+missing)
     - `select` (parity and allele, auto and direct; 601 inputs incl. indels with `rsid_pos_delta` = 1) against
       `select_rows` on a whole-chromosome read
     - 45 non-canonical inputs against the resolver

     All show 0 differences. It caught a mistake in the first P2 fix: label hits now return `ref_rsid: null` like the
     direct path, because the index stores `ref_rsid` only for merged docs.
   - **V3** `parity_inputs.py` re-implements LDlink's input chain expression for expression and resolves against the
     dump with Mongo order:
     - `validsnp`: lower, strip, split, first token
     - `replace_coord_rsid`, including warning text
     - `get_coords`: `strip("rs")` + `int()`
     - **Cases:** 2,327 inputs per build: case, whitespace, `rs`/`chr` stripping, `+`, `_`, leading zeros, extra
       tokens, merged IDs, multi-doc coordinates, malformed coordinates.
     - **Result:** the resolver and the planner (auto and direct) match on every field, including 480 warning texts:
       0 mismatches.
   - **P7 (found by V3, fixed):** the resolver rejected inputs with a second token, where LDlink uses the first token.
     Coordinate inputs also kept the typed position, while LDlink takes the replaced rsID's `find_one` position.
     - **Fix:** `resolve.mjs` `firstToken` + `resolveInputs()` (the whole chain), and the planner uses the same chain.
     - **Module logic for Phase C:** LDmatrix's own gate (input must end in a digit) rejects 120 inputs that the
       resolver accepts (`rs123rs`).
   - **V4** `verify_select.mjs`:
     - SV picks are a separate class. 1 per 50,000 on GRCh37: chr22:24,198,505 and chrX:114,542,484, where LDlink's
       last overlapping record is a symbolic SV that the build skips.
     - The sample adds 5,000 rsIDs overlapping no 1000G row.
     - **Result:** parity agrees on all other samples, and allele mode on all 55,000 per set (about 4,990 correct
       `not_in_1000g` each).
3. ✅ **The §2.2 items** (done 2026-11-21). Each fix has a regression test in `test/phaseA.test.mjs` (23 tests); the
   reference-based ones use fixtures from ldsc, pandas, numpy and scipy.
   - **K1 (reframed by reading LDlink):** only LDmatrix rounds to 3 decimals. LDtrait/LDexpress (`LD_calcs`),
     LDproxy, LDassoc and SNPclip threshold the **unrounded double**.
     - **Bit-identical math:** the kernel's double r²/D′ equal Python's `delta**2/Ms` and `abs(delta/min(...))`,
       checked on 3,000 count tuples (Python fixtures, including near-ties). So `band()` thresholds fall exactly like
       LDlink's.
     - **New API:**
       - `strict: true` for modules that compare with `>`;
       - `exact(p, I, J)` returns double values, because band/tile outputs are float32;
       - `pyRound(x, 3)` gives LDmatrix's display rounding. It lives in the new `pyfmt.mjs`, which `pyFixed3` also
         uses now.
   - **P3:** can't occur with this data. Not one dump position (984 M GRCh37, 1.0 B GRCh38) has only merged docs,
     so a coordinate always resolves to a current rsID. P7's chain would carry `ref_rsid` anyway.
   - **P5:** `resolve.mjs` treats only a missing file (manifest `exists()`, ENOENT, 404) as "not found". Network/403
     errors now throw.
   - **P4:** planner regions get chromosome names normalized (`chr22`, `x`, `23`), are merged when they overlap or
     touch, and invalid ones are reported in `plan.invalidRegions`.
   - **P6:** `sortedLookup` always reads its key column.
   - **K2:** `band()` sizes its output from the exact in-window pair bound (cap 4 M) and re-runs once when more pairs
     pass. A new C export `written()` makes truncation detectable; an explicit `maxOut` sets `truncated`.
   - **K3 (C fix, WASM rebuilt):**
     - `band_blocked` clamps `B` to at least 1, so it no longer loops forever.
     - An empty row range is now empty instead of "all rows".
     - Row ranges also work without `block`.
     - `kernel_check` re-run: 0 failures.
   - **B3/B4:** the genes order check now matches the SQL order (key, txEnd, dump order; LDlink's refGene query has no
     sort), and the gene files were rebuilt with identical rows. Annotation and dbSNP builders now **throw** on
     out-of-order input. A partial `build_annotations.mjs` run no longer overwrites the full build report.
   - **B5:** the BSON and JSON converters count numeric-looking strings as `numeric_string`. For position fields that
     LDlink queries numerically, they write them as `s:<value>`, so they stay un-keyed just as LDlink's numeric
     query never matches them. `dbsnp_bson2tsv.mjs` fails on non-string query fields (the 2023 dump has none).
     - **Found:** in all 7 sources, every string-typed position is `"NA"`.
     - **Correction:** the "one GWAS numeric string" reported earlier was an artefact of the old converter, which
       labelled extended-JSON numbers as strings. Re-fetched with the fixed converter (credentials refreshed), GWAS
       GRCh37 has 835,536 numbers, 96 `"NA"` strings and 0 numeric strings.
     - **The record in question:** it stores `{"$numberLong": "3411805015"}`, a real number but an impossible position
       (3.4 Gb). LDlink's window queries can't reach it, and `mkkey` doesn't key it (positions must be < 1e9).
     - **Result:** the regenerated TSV is byte-identical, so no release file changes.
       `fetch_annotations.sh [collection ...]` can now fetch single collections.
   - **B6:** `compare_ldlink_parse.sql` joins are NULL-safe and its order is deterministic. Results are unchanged
     (14,499,527 identical; 139,322 merged copies).
   - **B9:** `allele_join.mjs` orders by `vcf_order`. Exact §15 (b) numbers:
     - **GRCh37 chr22:** 1,087,038 same row, 3,307 different row, 151,831 alleles not on LDlink's row.
     - **GRCh38:** 1,018,866 same row, 6,372 different row, 177,129 alleles not on LDlink's row.
     - **Labels after the Phase A rebuild:** GRCh37 1,107,484 agree, 731 disagree, 66 unlabelled with an allele match.
   - **L5:** `--ld-wind` uses Python truthiness, so 0 or None means unset and exactly one window must be set.
   - **L6:** "More blocks than data points." (jackknife.py:104).
   - **L7:** `erfc(±∞)`. `chi2Isf1` now does Newton in log space with an asymptotic tail, so it matches scipy's `%.3f`
     Z for all 4,870 normal-range P values down to 2.2e-308 (it was NaN at 5e-324). Deep-subnormal P (≤ 4e-323)
     differs from scipy, whose cephes code loses precision there.
   - **L8:** `munge.mjs` now follows pandas and numpy:
     - Python header splitting (a BOM is not whitespace)
     - only the first mapped field is checked for duplicates; the last mapped value for "different ... columns"
     - repeated names: the first copy is read
     - INFO exemption by raw name
     - pandas float parsing (`inf`, no hex; signed-column read error)
     - numpy `_lerp` quantile (NaN for `[x, inf]`)
     - Python `repr` in the median warning
     - extra fields ignored

     10 edge-case files give the same output or error as the reference (`test/fixtures/l8`).
   - **Re-verified after all of the above:**
     - `--l2` byte-identical; BBJ munge identical; h²/rg ≤ 3.3e-14
     - kernel check: 0 failures
     - `verify_parity` and `verify_planner`: see §5
   - **Still open:**
     - **K4:** LDlink's sorted-string allele table, documented as a deliberate difference; 1 pair in the ALL check.
     - **B10:** left-alignment, a latent design limit.
     - **SV gap:** 1 per 50,000 picks.
4. ✅ **Rebuild the four genotype sets** (done with item 1). Still to do: re-run `allele_join.mjs` for the §15 (a)
   label table, after B9.

**Done when** every fix has a test that failed before it, and all verifiers pass on the rebuilt sets.

### Phase B: a complete data release (≈ 1 week, mostly machine time)

**Progress (2026-11-22):**
- **Inputs:** per-chromosome dbSNP and NCBI slices (`phaseB_slices.sql`), `build_all.sh` (resumable, 3 parallel) and
  `check_release.mjs` are written. The first full run then exposed two problems; both are fixed and the run was
  restarted from scratch.
- **12-ALT records:** the B7 guard stopped GRCh37 chr1 at a real **12-ALT** record (1:207,237,233). Genotypes are now
  parsed once per record by `vcf_gt.mjs` `gtCodes()`: a character fast path for indexes ≤ 9, an exact split
  otherwise.
  - GRCh37 chr22 rebuilt byte-identical to v4 in every column, and faster (230 → 152 s).
  - chrX genotypes were identical too.
  - Tested in `test/phaseA.test.mjs`.
- **Slice determinism:** the dbSNP slices were sorted by `(key, rsid)` only. DuckDB's sort is unstable, so the dump's
  duplicates of one rsID at one position (e.g. rs59096789 both current and merged into rs1556846434) came out in
  arbitrary order, and the labeller takes the first.
  - It showed up as 3 label differences in a chrX rebuild.
  - All 47 slices changed order in 34–11,898 rows when re-sorted.
  - Slices now sort by `(key, rsid, dump ord)`, LDlink's `find_one` order, which is a total and deterministic order.
  - Every chromosome was rebuilt, including the four chr22/chrX sets in the current release, which came from the
    same kind of unordered slices.
- ✅ **Done 2026-11-22: all 70 chromosomes** (GRCh37 1–22, X, Y; GRCh38 and High Coverage 1–22, X), 0 build
  failures.
  - **`check_release`:** 0 problems (counts, columns, metadata, order, every source row exactly once, tier bounds,
    labels).

    | build | chromosomes | rows | size |
    |---|---:|---:|---:|
    | GRCh37 | 24 | 85,210,348 | 3.26 GB |
    | GRCh38 | 23 | 78,229,218 | 3.05 GB |
    | High Coverage | 23 | 73,626,409 | 4.34 GB |

  - **Spot checks on new chromosomes:**
    - `kernel_check` vs LDlink math (GRCh37 chr1, GRCh38 chr7, HC chr2): 0 failures, both backends identical.
    - `verify_select` (GRCh37 chr1, HC chr2; 55,000 each): parity agrees apart from 1 known SV case; allele mode
      55,000/55,000.
    - The chr1 12-ALT site: 12 consecutive rows whose `ac` equal the VCF's `AC` exactly.
    - `verify_planner`: 15/15. Its whole-chromosome reference now reads each variant's own chromosome; it had
      assumed chr22, which surfaced once chr11 existed.
  - **chrY** (GRCh37): 62,042 records, 110 SVs skipped, every row haploid.
  - **Published:** `s3://eaap-ac-test-bucket/ldlink-parquet/v3/`, 333 objects, 34.8 GB. Genotypes were uploaded;
    dbSNP, annotations and LDscore were copied server-side from `v2/`, which is kept. The manifest covers 324 files
    (7 MB: per-build manifests would keep browser start-up small).
1. **Genotypes for every chromosome.**
   - Coverage: GRCh37 chr1–22, X, Y; GRCh38 chr1–22, X; High Coverage chr1–22, X.
   - Each build needs its dbSNP position slice and NCBI allele slice (one DuckDB export per chromosome, from
     `out/ncbi_dbsnp/b155/vcf3{7,8}.parquet`).
   - Write one driver script (`build_all.sh`), with the B7 guard.
   - Chromosome size: chr1 is ~4× chr22 (watch memory: the chrX build needed ~30 GB heap with both slices).
2. **Decide the open data questions** (§4 below) before the full build.
3. **Release integrity check:** `build_manifest.mjs`, plus a script that checks every file's row count, sort order,
   `vcf_order` monotonicity and metadata keys.
4. **Re-measure sizes** (PLAN_REVIEW finding 6 asks for this): the expected total is ~3–4 GB per build.
5. **Upload a versioned prefix** (e.g. `ldlink-parquet/v3/`) and keep `v2/` until the app moves over.

**Done when** the manifest covers every chromosome and build, and integrity plus parity checks pass on a sample
from each.

### Phase C: the app (the main remaining work)

**Progress (2026-11-22): the first end-to-end slice runs in a browser and matches LDlink.**
- **`app/`** (no build step: ES modules plus an import map for `hyparquet`) runs every module port from `src/` in the
  page. The release is read over HTTP Range requests and the LD comes from the WASM kernel. It has:
  - one form per module;
  - outputs rendered as LDlink prints them (tables, downloads);
  - an LDmatrix canvas plot with D′ below the diagonal and r² above, the connector, a Mb rug, the refGene track and
    hover values.
- **`app/check.html`** is the browser parity check. It runs 18 cases in the page and fetches the same requests from
  the local LDlink (which allows cross-origin requests), then diffs them: **18 / 18 identical**. That covers LDmatrix on
  three builds, LDpair, LDpop, LDproxy, SNPclip, SNPchip, LDtrait and LDexpress. LDassoc was checked in the app on the
  bundled example.
- **`scripts/serve.mjs`** is a local dev server with Range support. It mounts the repo, the release (`/release/`),
  the population panels (`/panels/`) and LDlink's LDassoc examples (`/examples/`). The launch config is
  `D:\Projects\LDlink\.claude\launch.json`, name `app`.
- **Speed in the browser (cold):**
  - LDmatrix 0.25–0.4 s, LDpair 0.1 s, SNPclip / SNPchip 0.25–0.3 s, LDexpress 0.5 s.
  - LDproxy: 8.7 s → 1.3 s. Function and FORGEdb lookups now go by position where that is exact; 11,912 reads → 874.
  - LDtrait: 33 s → 0.74 s. One window read and one kernel call per query, instead of a read and a call per GWAS
    record.
- **Core additions:**
  - `annotations.genes`, `annotations.recomb`, `annotations.forgeAt`;
  - `dbsnp.findAllAt`, `dbsnp.duplicated`;
  - `ldlink/common.rsAnnotations`, the positional fast path, with the `find_one` fallback for duplicated rsIDs.
- **Region view ("Pyramid" tab).** Any region up to 5 Mb, or a variant ± width:
  - **LD:** the band kernel runs in a module Web Worker (`app/pyramid_worker.mjs`; the kernel's import chain is now
    dependency-free via `src/haps.mjs`).
  - **Pixels:** each pair is placed at (midpoint, half distance) and each pixel shows the strongest pair there,
    splatted when variants are sparse.
  - **Tracks on the same axis:** GWAS Catalog (−log10 p, genome-wide line), GTEx eQTL for chosen tissues (colour =
    effect sign), recombination rate, RefSeq genes (transcripts collapsed per gene, packed with label room). Hover
    shows the pair / record.
  - **Inputs:** MAF tiers, max pair distance, threshold.
  - **Speed:**
    - 1 Mb EUR common: 4,675 variants, 2.8 M pairs, 0.44 s total (data 0.30 s, LD 0.06 s).
    - 3 Mb ALL common + low: 34,816 variants, 228 M pairs, 7.5 s (LD 7.0 s on one WASM thread; the page stays
      responsive).
  - **Check:** a hovered pair (rs9625886 × rs79955051, EUR) reads r² 0.350; LDlink's LDpair gives 0.3497.
  - **Core additions:** `windowRows(…, tiers)`, `annotations.eqtlRange`.
- **Worker pool** (`app/pyramid_pool.mjs`):
  - **Split:** the band's first-variant rows are cut into shares of equal pair counts (two pointers), up to
    min(cores − 1, 8) workers. A share needs at least 2 M pairs, so small regions stay on one worker.
  - **Transfer:** each worker receives only its share plus the rows within `maxDist` after it.
  - **Merge:** the page merges the pixel grids by max. Ties go to the smallest (i, j) in the worker and in the merge,
    so the result does not depend on the worker count.
  - **Check:** chr22 3 Mb ALL common + low (34,816 variants, 228 M pairs): 1 worker vs 8 workers, 0 value
    differences and 0 pair differences.
  - **Speed:** 7.0 s → 1.04 s (6.8×) on the 32-core dev machine.
  - **Masked path:** rows with haploid or uncalled haplotypes (chrX males) use the kernel's masked path, which now
    takes row ranges too (`band_masked` in C, `bandMasked` in JS). Masked regions split the same way:
    - chrX non-PAR 2 Mb EUR (6,938 variants): 1.18 s → 0.14 s in the app;
    - chrX non-PAR 3 Mb ALL common + low (26,581 variants, 192 M pairs): 23.9 s → 3.9 s.
    - 1 worker vs 8: 0 differences in both.
    - **Unit test K4:** on both backends, any split of the masked band equals the full band, pair for pair.
    - **kernel_check (GRCh37 chr1 and chrX):** 0 failures. The masked path runs about 4× slower per pair than the
      unmasked one (four popcounts per pair, no pruning).
- **Still to do:**
  - a faster masked path: prune with a bound, or group samples by mask;
  - FORGEdb and RegulomeDB as point tracks;
  - a click-to-LDproxy link from the pyramid;
  - per-module plots (LDproxy / LDassoc locus plots);
  - the LDscore views;
  - population panels in the release (served from `db/` today);
  - a worker for long windows;
  - replacing the root prototype.

Target: one LD module (REQUIREMENTS.md §1) with annotation tracks on a matrix/pyramid view; then each LDlink module
as a thin view over it.
1. **Core library** (`src/core/`, promoted from `experiments/`):
   - **Pipeline:** resolve (`resolve.mjs`), plan (`planner.mjs`), select rows (`select_rows.mjs`; parity default,
     allele opt-in), decode haplotypes, then the LD kernel (WASM + worker pool).
   - **Output contract** (finding 8): per-pair stats, tile schema and zoom levels for the pyramid, thresholded-pair
     export.
   - **Parity details:** LDlink's input normalization, warnings and errors (REQUIREMENTS §1.2), 3-decimal rounding
     (K1), and the deliberate differences list (PLAN_REVIEW finding 9).
2. **Annotation track API** (REQUIREMENTS §3): `(build, chrom, start, end)` or rsID list → rows, using `rangeScan`
   (the planner should use it for scans as well, per §11's reader lesson).
3. **Views:**
   - LDmatrix: matrix plus pyramid, genes/recombination/FORGEdb/RegulomeDB/GWAS/eQTL tracks.
   - LDpair, LDpop, LDproxy, LDhap, LDtrait, LDexpress.
   - LDassoc (user file), SNPclip, SNPchip.
   - LDscore/LDherit/LDcorrelation (JS ports, after Phase A's L fixes).
4. **Exports:** TSV/text outputs with LDlink's columns (REQUIREMENTS §4: same values to 3 decimals, not
   byte-identical), plus PNG/SVG from the browser.
5. **Replace the old prototype** at the repo root. Its `sources.js` example SNPs are the test vectors.

**Done when** each module reproduces LDlink's numbers on its reference cases (Phase E), in a browser, against the
hosted release.

### Phase D: hosting
1. **CloudFront (HTTP/2) in front of the bucket, with CORS.** S3's REST endpoint is HTTP/1.1, which allows about 6
   browser connections (§11). This is outward-facing, so it needs your go-ahead and a decision on the domain/account.
2. **Browser network benchmark:** 1/10/100 Mb windows and 300/3000-SNP lists, cold and warm. Recalibrate the
   planner's cost constants on real HTTP (§10 open item).
3. **Cache headers and versioned prefixes** so file updates are atomic.

### Phase E: parity against the live service

**Progress (2026-11-22):** LDlink now runs locally (`D:\Projects\LDlink\ldlink-local`, CORE_REVIEW.md §1 and §6). **LDmatrix: 30 / 30 matrices identical** to LDlink across all three builds, six populations, indels, multi-allelic sites, coordinates, chrX and merged rsIDs (`experiments/compare_ldmatrix.mjs`). No token is needed.

**Module ports compared** (`experiments/modules/`, CORE_REVIEW.md §7), on `experiments/out/release_v4`: GRCh37
rebuilt with the new `vcf_id` column (the raw VCF ID LDlink prints when it is a merged or unlisted rsID).

| module | identical |
|---|---|
| LDmatrix | 30 / 30 matrices |
| LDpair | 420 / 429 pairs; 9 differ only in tie order, which LDlink randomizes per process |
| LDpop | 537 / 537 pairs |
| LDhap | 70 / 70 lists |
| LDproxy | 57 / 57 queries |

**Uploaded (2026-11-22):** `s3://eaap-ac-test-bucket/ldlink-parquet/v4/`, 333 objects, 34.8 GB.
- **How:**
  1. In-bucket copy of `v3/`.
  2. `aws s3 sync` of `out/release_v4`, which uploaded the 72 rebuilt GRCh37 genotype files, `manifest.json` and three
     annotation files with newer local copies.
- **Verified:** every object's size equals the local file's; sampled objects are byte-identical.
- **Access:** the bucket is still private (all four public-access blocks on). `v2/` and `v3/` are kept.
- **Pending:** make `release_v4` the local `out/release`.

**The other five modules** (CORE_REVIEW.md §8):

| module | identical |
|---|---|
| SNPclip | 49 / 49 |
| SNPchip | 66 / 66 |
| LDtrait | 55 / 55 |
| LDexpress | 58 / 58 |
| LDassoc | 21 / 21 |

Each was compared against LDlink's REST text, or against the Python function for LDassoc, which has no API.

**Phase E is done for every module.**

**The `src/` refactor is done** (CORE_REVIEW §9, `src/README.md`): the runtime core, the LD kernel, ldsc and all ten
module ports are in `src/`, and `experiments/run_regression.sh` passes in full on the refactored tree.

**Next:** Phase C (the app) on `src/`. Optional core follow-ups from §9:
- merge `sortedLookup` / `rangeScan`;
- use numbers for genome keys;
- split `planner.query`.
1. **Cross-check the test vectors on ldlink.nih.gov** (REQUIREMENTS §5 says this is not done yet). The web endpoint
   needs a browser session; the API needs a token, which you would supply.
2. **Extend the vectors:** indels, a multi-allelic site, chrX males, an NA (monomorphic) pair, a High Coverage
   duplicate position, a merged rsID, a coordinate input, and one case per §15 failure class (alleles not in 1000G,
   different row).
3. **One reference case per module,** with expected outputs saved under `scripts/profile/results/`.

### Phase F: performance (after the app works)
- Balance worker ranges by pair count.
- chrX mask-grouped fast path.
- WebGPU for ALL-population tiles over 100 Mb windows (§9 open items).
- A slimmer rsID index (drop `chrom`; move `ref_rsid` to its own file): about 45 KB per scattered rsID cold (§10).
- Optionally, a separate rare-variant tier file (REQUIREMENTS §6 Q4).

### Phase G: housekeeping
- **Move code:** reusable code from `experiments/` goes to `src/` with tests; experiments stay as measurement scripts.
- **Disk** (all regenerable; your call to delete):

| path | size | regenerate with |
|---|---:|---|
| `experiments/out/ncbi_dbsnp/b15{5,7}/vcf3{7,8}.tsv` | 232 GB | `ncbi_vcf2tsv.mjs` (~40 min each); the Parquet copies hold the same data |
| `experiments/out/ncbi_dbsnp/b15{5,7}/compare.duckdb` | 242 GB | `compare_ncbi.sql` (~1 h each); the summary CSVs are kept |
| `experiments/out/dbsnp_dump/` intermediates (`work.duckdb`, sorted dirs) | ~170 GB | `dbsnp_sort.sql`. Keep `dbsnp_dump.tsv` + `work.duckdb` while parity checks still use them |
| `experiments/out/v2_*`, `grch37_chr22/` | ~1 GB | superseded by `v3_*` |

## 4. Decisions needed from you

1. **GRCh38 chrX** (2019 release: ~107k records, no haploid calls): keep, replace or hide? *(REQUIREMENTS §6 Q1)*
2. **GWAS Catalog:** the 2025-10 backup, or a fresh EBI pull at build time? GTEx v8 or newer? *(PLAN_REVIEW
   finding 10)*
3. **dbSNP for allele mode:** stay on b155 (matches the dump, so parity and allele mode agree on rsIDs), or offer
   b157 (+120.6 M rsIDs)? Recommendation: b155 now; b157 as a later "current dbSNP" mode.
4. **Hosting:** CloudFront distribution + CORS on `eaap-ac-test-bucket`, or another account/domain? (Phase D)
5. **LDlink API token** for live cross-checks. (Phase E)
6. **Disk cleanup:** which of the Phase G intermediates to delete.

## 5. Suggested next step

Phase A is complete.

**Next is Phase B:** build every chromosome once with the fixed labeller. Before starting it, decide §4 Q1 (GRCh38
chrX) and Q2 (GWAS source / GTEx version).
