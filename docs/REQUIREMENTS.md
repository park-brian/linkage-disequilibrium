# Requirements: a browser-only LDlink

Goal: replace LDlink's modules with **one LD module**. You give it a variant set and it computes pairwise
statistics over it. Every annotation can be computed on demand and drawn on any result, so the separate LDlink tools
become presets and views of the same engine instead of separate code paths:

- The **variant set** can be any of these, or a mix:
  - rsIDs or `chr:pos` values, **with no count limit**
  - one or more regions, typically **1–10 Mb**, sometimes up to 100 Mb
- Populations, genome build, and filters (MAF tier, r²/D′ thresholds, a distance band for the pyramid) apply to that
  set.
- **Outputs:**
  - streamed tiles of pairwise stats (r², D′, haplotype counts, correlated alleles)
  - matrix/pyramid views
  - "proxies of X" (one row of the matrix)
  - pruning (SNPclip)
  - haplotype tables (LDhap)
  - per-population comparison (LDpop)
- The **annotation tracks** (§3) attach by position or rsID to any of these: genes, recombination, RegulomeDB,
  FORGEdb, GWAS, GTEx eQTL, arrays.

The LDlink limits in §2 are reference behaviour to compare against, not limits we keep. The data layout and fetch
strategy for large windows are in [DATA.md](DATA.md) (the release contract) and [history/PARQUET_LAYOUT.md](history/PARQUET_LAYOUT.md) (the measurements).

**Where it stands (2026-11-26).** Built toward the one module:
- the engine: one LD kernel (WASM and JS, identical results) with banded and tiled computation, on Parquet reads that
  return every matching row (`src/core`, docs/ARCHITECTURE.md);
- LDlink itself: all thirteen modules reproduced byte for byte, uploads included (`src/parity`, [PARITY.md](PARITY.md));
- `src/core/ld.js`: pair, matrix and proxy over any variant list with no count limit, numbers and typed errors, the
  allele-aware row by default;
- the region view: an LD pyramid up to 5 Mb under the annotation tracks, on a pool of Web Workers;
- the set page (`app/set.html`, the visual work in progress): one variant set, the view chosen by what is held, on
  the parity layer's context and ports.

Open: the views on `core/ld.js` rather than the ports (so LDlink's limits stop applying); annotations as core reads
(they are LDlink's queries in `src/parity/annotations.js`); regions and streamed tiles in the core API; pruning,
haplotype tables and per-population comparison outside the ports.

- Runtime: browser only. Static Parquet (DATA.md) read with **hyparquet** over HTTP Range; the region view's LD in
  Web Workers. No query server.
- Build: Parquet produced from the sources in DATA.md with **hyparquet-writer**.
- Reference implementation: `../nci-webtools-dceg-linkage` (paths below are relative to its `server/`).

## 1. Core LD engine (shared by every module)

### 1.1 Inputs common to all modules
| input | values | reference |
|---|---|---|
| genome build | `grch37`, `grch38`, `grch38_high_coverage` (HC shares GRCh38 coordinates and refGene) | `LDcommon.py:24-65` |
| populations | `ALL`, `AFR AMR EAS EUR SAS`, 26 populations; union of selected panels, de-duplicated | `LDcommon.py:460-481` |
| variants | `rs123` or `chr7:12345` (case-insensitive, de-duplicated, max per module below) | `LDcommon.py:341-376` |
| statistic | `r2` or `d` (D′) for thresholds/sorting/colors | `LDlink.py` param validation ~599-736 |

### 1.2 Variant resolution
1. `rsN` → the first dbSNP record with that id, in dump order (`find_one`; `LDcommon.py:378-442`).
2. `chr:pos` → the dbSNP records at that position for the build → the first with `id == ref_id` (warn "Multiple
   rsIDs" when several), else the first record (`LDcommon.py:220`).
3. Errors/warnings to reproduce: position `NA` for the build; chrY with a GRCh38 build; not all on one chromosome
   (LDmatrix, LDhap, SNPclip); span > 1 Mb → "switch rate errors" warning (LDmatrix).
4. Genotype row lookup is **by position** (VCF IDs are absent in GRCh38). Reference behaviour with several rows:
   the last row overlapping the position in VCF order (`LDcommon.py:590-605`; its POS check is dead code). Warnings:
   "missing from 1000G", "Genomic position does not match RS number at 1000G position", "variant is not a
   biallelic", "monoallelic in the chosen population(s)" (`LDcommon.py:504-571`).
   *New behaviour (built: `pickAllele`, the default of `core/ld.js`)*: the row whose alleles match the rsID's NCBI
   alleles; several matches or none are reported instead of silently picking.

### 1.3 Alleles and haplotypes
- Allele normalization `set_alleles` (`LDcommon.py:486-499`): SNV unchanged; insertion → `-`/`ALT[1:]`; deletion →
  `REF[1:]`/`-`; both longer → strip first base of each.
- GT mapping (`LDmatrix.py:221-244`): `0|0 0|1 1|0 1|1` → two haplotypes; haploid `0`/`1` → one allele + missing;
  anything else → both missing. Missing haplotypes are excluded per pair.
- Only samples in the selected panels count; haplotypes = 2 × samples (minus missing).

### 1.4 Pairwise statistics (`LDmatrix.py:268-372`, `LDpair.py:190-225`, `LDcommon.py:655-678`)
For variants i, j, count the four 2-locus haplotypes; zero-fill absent classes. Reference orders the classes by
**lexicographic sort of the concatenated allele strings** → A, B, C, D. Then:

```
Δ  = A·D − B·C
Ms = (A+C)(B+D)(A+B)(C+D)
D′ = |Δ| / min((A+C)(A+B), (B+D)(C+D))   if Δ < 0
     |Δ| / min((A+C)(C+D), (A+B)(B+D))   otherwise
r² = Δ² / Ms
Ms = 0 → D′ = r² = "NA"
χ² = N·Δ²/Ms,  p = 2·(1 − Φ(√χ²))          (LDpair, LDpop)
correlated alleles (only if r² > 0.1): if A·D / max(B·C, 0.01) > 1 → pair class[0], class[3] else class[1], class[2]
MAF (LDproxy) = min(allele freq) per locus (LDproxy_sub.py:55-56)
rounding: 3 decimals (Python round = banker's rounding on exact ties)
```

Implementation requirement: count by **allele index (0/1)**, not by string concatenation. The reference's string
keys can collide for normalized indels (e.g. `"T"+"TA"` vs `"TT"+"A"`) or MNPs that normalize to identical alleles;
document any such divergence instead of copying it. D′ and r² are invariant to allele relabeling, so index order is
safe for every case where the reference is well-defined.

### 1.5 Performance (browser)
The targets were set from history/PARQUET_LAYOUT.md; what was built and measured is in [STATUS.md](STATUS.md) ("Speed",
"Requests") and [decisions.md](decisions.md):
- **Any SNP list:** resolved to regions first when that is cheaper (the planner's hybrid rsID strategy), so cost grows
  with the regions, not the list length.
- **Pairwise math:** a WASM kernel with scalar popcount (SIMD was never faster in V8), banded and tiled, with pyramid
  pixels aggregated in Web Workers (a 3 Mb pyramid over ALL, 228 M pairs: 1.0 s on 8 workers); never a full N×N
  matrix.

## 2. Modules

These are LDlink's own limits and outputs. Use them for parity tests; the new engine has no SNP-count limit.

| module | inputs & limits | output | reference |
|---|---|---|---|
| **LDmatrix** | 2–300 variants (API: 3000), one chromosome; pops; `r2`/`d` | N×N D′ and r² (tab text: header `RS_number` + rsIDs); plot: D′ lower triangle (blue, alpha=D′), r² upper (red, alpha=r²), diagonal purple, NA gray; connector from evenly spaced columns to real positions; position rug (Mb); refGene track (all or collapsed transcripts, 80 kb packing gap, ±2.5% buffer); optional FORGEdb score on diagonal | `LDmatrix.py:40-44, 139-158, 268-389, 394-1229` |
| **LDpair** | up to 10 pairs | per pair: allele counts/freqs, 2×2 table, 4 haplotypes (count, freq), D′, r², χ², p, correlated alleles | `LDpair.py:33, 190-225` |
| **LDpop** | 1 pair; pops (super-pops expand to sub-pops) | per population: N, allele freqs of both, r², D′, χ², p; map (lat/long table) | `LDpop.py:120-158, 343-518` |
| **LDproxy** | 1 variant; window 0–1,000,000 (default 500,000); `r2`/`d`; collapse transcripts; annotate FORGE/Regulome | proxies with r² > 0.01: RS, coord, alleles, MAF, distance, D′, r², correlated alleles, FORGEdb, RegulomeDB, dbSNP function; sorted by \|distance\| then stat; UCSC bed/bedGraph track (r² bins 0.8–1 … 0–0.2); plot with recombination rate + genes | `LDproxy.py:80-84`, `LDproxy_sub.py` |
| **LDhap** | ≤ 30 variants, one chromosome | haplotype strings with count & frequency (desc), per-SNP alleles/frequencies | `LDhap.py:26` |
| **LDtrait** | ≤ 50 variants; window ≤ 1 Mb (default 500 kb); `r2`/`d` threshold (default 0.1); web warns > 5,500 pairs | GWAS Catalog hits in LD: query, trait, PMID, RS, position, alleles, r², D′, risk allele, effect size (95% CI), beta/OR, p | `LDtrait.py:192, 329-341` |
| **LDexpress** | ≤ 10 variants; window ≤ 1 Mb; r²/D′ threshold 0.1; p threshold 0.1; tissue list | GTEx v8 eQTLs in LD: query, RS, position, r², D′, gene symbol, Gencode ID, tissue, non-effect/effect allele freq, effect size, p | `LDexpress.py:59`, `LDexpress_tissues_sub.py:29-73` |
| **LDassoc** | user GWAS file (chr, pos, p columns); mode variant / gene / region; windows 500 kb / 100 kb / 0; region ≤ 3,000,000 bp | LD vs index variant (given or lowest p) for every file variant: same columns as LDproxy + p; Manhattan-style plot with LD colouring, recombination, genes; example dataset available | `LDassoc.py:342-344` |
| **SNPclip** | ≤ 5000 variants, one chromosome; r² threshold 0.1; MAF threshold 0.01 | greedy LD pruning in input order: kept list + per-variant reason (MAF, monoallelic, in LD with X, not in 1000G) | `SNPclip.py:20, 143-171` |
| **SNPchip** | ≤ 5000 variants; platform list | per variant: position + arrays containing it | `SNPchip.py:66` |
| **LDscore** | PLINK bed/bim/fam (the example or an upload); cM or kb window | ldsc `--l2` log from "Summary of LD Scores" | `LDlink.py` ldscore route, `ldsc_utils.run_ldsc_command` |
| **LDherit** | sumstats (example or upload, normalized as LDlink does); LD score population; observed or liability scale | ldsc `--h2` log | `LDlink.py` ldherit route, `sumstats_normalizer.py` |
| **LDcorrelation** | two sumstats; population; scale | ldsc `--rg` log | `LDlink.py` ldcorrelation route |

## 3. Annotation track contract

Every annotation is a source that, given `(build, chrom, start, end)` or a list of rsIDs, returns rows the plots can
draw. Tracks are independent of the LD computation and can be attached to any view.

| track | key | fields a view needs | source (DATA.md, "annotations/") |
|---|---|---|---|
| genes | region | symbol, transcript, strand, tx/cds bounds, exons; collapse-by-symbol option | `refGene_grch37/38` |
| recombination | region | pos, rate (cM/Mb) | `recomb` |
| RegulomeDB | position | score category | `regulome` |
| FORGEdb | rsID | score 0–10 | `forge_score` |
| dbSNP function | rsID/position | function list | `dbsnp.function` |
| GWAS Catalog | region | trait, PMID, risk allele, effect, CI, p | `gwas_catalog` (2025 backup newer) |
| GTEx eQTL | region (+ tissue filter) | gene, tissue, effect size, p, maf | `gtex_tissue_eqtl` + `gtex_genes` + `gtex_tissues` |
| genotyping arrays | position | platform codes | `snp_col` + `platforms` |

## 4. Out of scope / deferred
- API tokens, user registration, rate limiting, email (`ApiAccess.py`), server-side SVG/PDF export (`bokehExport.py`,
  phantomjs) — browser renders and exports directly.
- ldsc's `validBfile` schema messages for malformed uploads (they carry LDlink's server paths; [PARITY.md](PARITY.md)).

The ports go further than this document first asked: their outputs are byte-identical to LDlink's, checked against a
local LDlink (`npm run verify:ldlink`) and on every push by the parity snapshot.

## 5. Reference test vectors
Generated by `verify/reference/ldmatrix_reference.js` — a line-for-line port of the reference LDmatrix math run on
the local 1000G VCFs (same files as S3). Results: `verify/reference/results/tv_*.json` (include raw haplotype counts).

Example SNPs (`sources.js` `exampleSnps`), all chr2, population **CEU (99 samples)**. Ordered by position
(GRCh37 / GRCh38): rs12465506 217894234/217029511 G/A · rs13412666 217905779/217041056 G/A ·
rs13387042 217905832/217041109 A/G · rs13426489 217906246/217041523 T/G · rs12620095 217907514/217042791 C/G ·
rs6721996 217909463/217044740 G/A · rs12613955 217916052/217051329 T/G · rs4442975 217920769/217056046 G/T ·
rs13409592 217926328/217061605 T/C · rs12622536 217936015/217071292 T/C

Full symmetric matrices; values are identical for GRCh37, GRCh38 and GRCh38 High Coverage:

```
r2          rs12465506 rs13412666 rs13387042 rs13426489 rs12620095 rs6721996 rs12613955 rs4442975 rs13409592 rs12622536
rs12465506  1          0.92       0.92       0.92       0.792      0.881     0.775      0.843     0.153      0.458
rs13412666  0.92       1          1          1          0.864      0.96      0.846      0.919     0.205      0.511
rs13387042  0.92       1          1          1          0.864      0.96      0.846      0.919     0.205      0.511
rs13426489  0.92       1          1          1          0.864      0.96      0.846      0.919     0.205      0.511
rs12620095  0.792      0.864      0.864      0.864      1          0.829     0.979      0.826     0.176      0.604
rs6721996   0.881      0.96       0.96       0.96       0.829      1         0.812      0.96      0.214      0.487
rs12613955  0.775      0.846      0.846      0.846      0.979      0.812     1          0.846     0.193      0.618
rs4442975   0.843      0.919      0.919      0.919      0.826      0.96      0.846      1         0.228      0.511
rs13409592  0.153      0.205      0.205      0.205      0.176      0.214     0.193      0.228     1          0.162
rs12622536  0.458      0.511      0.511      0.511      0.604      0.487     0.618      0.511     0.162      1

D'          rs12465506 rs13412666 rs13387042 rs13426489 rs12620095 rs6721996 rs12613955 rs4442975 rs13409592 rs12622536
rs12465506  1          0.979      0.979      0.979      0.977      0.938     0.977      0.937     0.802      0.819
rs13412666  0.979      1          1          1          1          1         1          0.959     0.949      0.848
rs13387042  0.979      1          1          1          1          1         1          0.959     0.949      0.848
rs13426489  0.979      1          1          1          1          1         1          0.959     0.949      0.848
rs12620095  0.977      1          1          1          1          1         1          0.978     0.945      0.857
rs6721996   0.938      1          1          1          1          1         1          1         0.951      0.845
rs12613955  0.977      1          1          1          1          1         1          1         1          0.858
rs4442975   0.937      0.959      0.959      0.959      0.978      1         1          1         1          0.848
rs13409592  0.802      0.949      0.949      0.949      0.945      0.951     1          1         1          1
rs12622536  0.819      0.848      0.848      0.848      0.857      0.845     0.858      0.848     1          1
```

Status: these vectors only exercise common biallelic SNVs. The LDmatrix port is compared against a local LDlink (the
same code and data as ldlink.nih.gov) on 15 variant lists (30 cases: both statistics), which include indel labels,
the 12-ALT multi-allelic site on chr1, chrX inside and outside PAR1, coordinates, a merged rsID, GRCh38 and High
Coverage (`verify/ldlink/cases/ldmatrix_cases.json`); `test/ld.test.js` checks `core/ld.js` against the port on the
example SNPs.

## 6. Open questions
1. GRCh38 (2019) chrX has only ~107k records and no haploid calls — keep, replace with another release, or hide?
2. ~~GWAS Catalog source~~ **Decided:** the 2025-10 backup in LDlink's S3 (`pipeline/annotations/fetch.sh`), as LDlink
   serves it.
3. ~~dbSNP source~~ **Decided:** keep the 2023 b155 mongodump. It is LDlink's own data and is needed to reproduce LDlink's exact rsID resolution and LD values. See PLAN_REVIEW.md finding 4 for the join rule and parity mode.
4. ~~Rare variants in a separate file~~ **Decided:** the genotypes are split into MAF tiers (`common`, `low`, `rare`
   files per chromosome; DATA.md), so common queries read only the common tier.
