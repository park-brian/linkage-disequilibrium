# Parquet layout: small and fast for 1–100 Mb LD queries

Goal: keep LDlink's data as small as possible and fast to query from the browser (hyparquet over HTTP Range),
for **one LD engine** that takes any number of variants. Typical windows are **1–10 Mb**, sometimes up to 100 Mb.

All numbers are measured on GRCh37 chr22 (35 Mb, 1.1M variants, 2504 samples) with hyparquet 1.33.1 /
hyparquet-writer 0.18.0. Scripts are in [`experiments/`](../experiments/), raw results are in
[`experiments/results/`](../experiments/results/). "Requests" and "bytes" come from an instrumented AsyncBuffer
over a local file. Timings are local decode only, so add network time on top.

## 1. Scale of the problem

| window | all 1000G variants | super-pop MAF ≥ 5% | dbSNP rows |
|---|---|---|---|
| 1 Mb | ~33k | ~5.9k | ~380k |
| 10 Mb | ~300k | ~45k | ~3.7M |
| 100 Mb | ~3M | ~450k | ~37M |

- Most variants are rare: 73% have max super-population MAF < 1%.
- All-pairs LD over 300k variants is 4.5e10 pairs. Even a 500 kb band over 10 Mb of common variants is 103M pairs.
  So the engine must work band-limited and tiled, and it must not build an N×N matrix.

## 2. Genotypes: what we measured

Each row is one biallelic allele (multi-allelic sites split), sorted by `pos`. The candidate layouts differ only in
the `haps` column.

Full chr22 file sizes (`results/genotypes_build_grch37_chr22.json`):

| layout | size | B/variant |
|---|---|---|
| VCF.gz (source) | 208 MB | 188 |
| bitset, SNAPPY, 64 KB pages | 99.9 MB | 90.1 |
| bitset, GZIP, 64 KB | 62.7 MB | 56.5 |
| bitset, ZSTD-3, 16 KB / 64 KB / 256 KB / 1 MB | 74.0 / 58.3 / 49.0 / 44.0 MB | 66.7 / 52.5 / 44.1 / 39.6 |
| sparse minor-allele varint deltas, ZSTD-3, 64 KB | 48.7 MB | 43.9 |
| **tiered**, ZSTD-3, 256 KB: common (≥5%) bitset + low (1–5%) bitset + rare (<1%) sparse | 26.3 + 5.6 + 7.5 = **39.4 MB** | 155 / 43 / 9 |

Region scans (`results/genotypes_bench_grch37_chr22.json`):

| query | layout | rows | requests | bytes | decode |
|---|---|---|---|---|---|
| 1 Mb | bitset ZSTD 256 KB, all variants | 33,451 | 3 | 2.6 MB | 77 ms |
| 1 Mb | common tier only | 5,890 | 3 | 1.7 MB | 20 ms |
| 10 Mb | bitset ZSTD 1 MB, all variants | 302,633 | 15 | 11.7 MB | 0.48 s |
| 10 Mb | common tier only | 45,369 | 15 | **7.6 MB** | 0.15 s |
| 10 Mb | common + low + rare(sparse) | 302,633 | 45 | 11.5 MB | 0.44 s |

What the scans show:
- **Requests stay low.** hyparquet merges adjacent pages into single reads of up to 2 MB, so a 10 Mb window takes
  about 15 requests per file. They arrive in about 3 sequential phases: footer, page index, data.
- **Tiering is the biggest saving.** It works because common variants barely compress (155 B each) while rare ones
  are almost free (9 B sparse). Most LD use only needs common variants.
- **Bigger pages compress better.** The trade-off is that a point lookup reads more data than it needs.
- **A SNP list inside a window costs the same as scanning the window.** For example, 300 or 3000 SNPs spread over
  10 Mb touch essentially every page.
- **Don't use hyparquet's `$in` filter for lists.** It is slow (1.6–8 s) because it builds row objects. Use the
  page-index lookup in `experiments/lookup.mjs`, or scan the window and select by `pos`.

LD math, banded at 500 kb over 10 Mb of common variants (103M pairs), plain JS 32-bit popcount:

| population | haplotypes | pairs/s | time |
|---|---|---|---|
| CEU | 198 | 26M | 3.9 s |
| EUR | 1006 | 9.3M | 11 s |
| ALL | 5008 | 2.2M | 46 s |

Before packing, the haplotypes are gathered down to the selected population (`packPopulation`). This takes 0.05–1 s.
**For 10 Mb windows the math, not I/O, is the bottleneck.** We need WASM SIMD popcount (expected 10–20×) or
WebGPU, plus computing only what is drawn: aggregated tiles for the pyramid, and full detail on zoom.

## 3. dbSNP: two indexes plus a hybrid lookup

The composite key is `chrom_code × 1e9 + pos`, with chrom codes 1–22, X=23, Y=24, MT=25. It is a single INT64 that
gives one total order across the genome and reads well (e.g. `22025000000`). There is one key per build:
`key37` and `key38`.

**Measured fact:** rsIDs from one region do not cluster in rsID space. One 10 Mb chr22 window contains 3.77M dbSNP
rows that spread across about 405k distinct 1k-row pages of the rsID-sorted file, covering nearly the whole ID
range. So an rsID-sorted lookup costs about one page per rsID, whatever the list's locality.

| index | sort | columns | measured size | full-genome estimate |
|---|---|---|---|---|
| rsID index | `rsid` (INT64, delta) | `rsid`, `key37`, `key38`; merged IDs (`ref_rsid`) in a separate small file | 10.0 B/row | ~10.7 GB (all 1.07B) or ~1–1.5 GB if limited to 1000G rsIDs + merges |
| position index | `key` (INT64, delta) | `key`, `rsid`, `type`, `function` | 4.5 B/row | ~4.8 GB per build |

Lookup costs (`results/dbsnp_bench*.json*`), with page size 8 KB, 1M-row row groups, and a cached AsyncBuffer:

| query | requests | bytes |
|---|---|---|
| 1 rsID (cold: footer + page index + page) | 7 | 0.3–0.8 MB |
| 16 rsIDs | ~74 | 1.9 MB |
| 300 rsIDs | ~920 | 7.7 MB (26 KB/rsID) |
| 3000 rsIDs | ~7,200 | 52 MB (17 KB/rsID) |
| position window 1 Mb, all dbSNP | 3 | 4.8 MB |
| position window 10 Mb, all dbSNP | 11 | 22 MB |

### Recommended hybrid
1. **Direct path, ≤ ~32 rsIDs.** Do batched page-index lookups in the rsID index, in parallel. This takes about
   2 requests and 17–25 KB per rsID.
2. **Region path, larger lists:**
   1. Resolve a sample of about 32 input rsIDs through the direct path.
   2. Cluster their positions into regions, merging gaps under about 2 Mb.
   3. Scan the **genotype tiers** for those regions. They carry an `rsid` column joined from dbSNP at build time,
      and the LD computation needs these rows anyway.
   4. Match every input rsID in memory, so the cost depends on window size, not list length.
   5. Send the few unmatched IDs (merged rsIDs, variants not in 1000G, outliers) back through the direct path.
3. **Window queries (1–100 Mb).** No dbSNP lookup is needed. The genotype tiers already provide pos → rsid for every
   variant that LD can use.
4. **The full position-sorted dbSNP index** is only for annotation (dbSNP `function`) and for `chr:pos` inputs that
   are not in 1000G. It is too heavy to scan routinely (22 MB per 10 Mb), so read it page-wise like the rsID index,
   or limit it to annotated rows.
5. **Partition the rsID index by rsID range** (e.g. 64M IDs per file, about 33 files). The client computes the file
   name from the rsID, so each footer stays a few KB and nothing has to be fetched just to locate a partition.

## 4. Recommended layout (v1)

> Revised by [PLAN_REVIEW.md](PLAN_REVIEW.md): tier on max MAF over the **26 sub-populations** (super-population
> tiers miss 20% of sub-population-common variants), join rsIDs by position + type from the 2023 dump (kept for LDlink parity), keep the full rsID index, and use a
> ploidy mask for chrX.

```
genotypes/<build>/chr<N>.common.parquet   max sub-pop MAF ≥ 5%,    bitset haps, ZSTD-3, 64 KB pages
genotypes/<build>/chr<N>.low.parquet      1% ≤ maf_max < 5%,       bitset haps, ZSTD-3, 64 KB pages
genotypes/<build>/chr<N>.rare.parquet     maf_max < 1%,            sparse haps, ZSTD-3, 16 KB pages
   columns: pos INT32 (delta, column index) · rsid INT64 · rsid_candidates INT8 · ref · alt · multiallelic_source
   BOOLEAN · ac INT32 · haps · missing (null unless missing calls) · maf_max FLOAT; rows in VCF order within a
   position (parity); 50k-row row groups; sample order + ploidy mask in key-value metadata
dbsnp/rsid/part-<rsid>>26>.parquet        full 2023 dump (1.07B, 32 files): rsid INT64 (delta, column index) · ref_rsid
                                          (null unless merged) · chrom · key37 · key38; rows in (rsid, dump order); 8 KB pages
dbsnp/<grch37|grch38>/chr<N>.parquet      key INT64 (delta, column index) · rsid · is_ref · type · function;
                                          rows in (key, dump order); 16 KB pages
annotations/<build>/...                   per DATA.md §5, sorted by key/pos
manifest.json                             footer ranges, row-group key ranges, page-index locations (§7)
populations.json                          per build: sample order, population → sample indexes
```

**Estimated size per build:** about 3.6 GB of genotypes. That is chr22's v2 total of 46.6 MB (§6) scaled by
variant count, about 77× for the genome. The full rsID index (needed for LDlink parity) adds ~10.7 GB, shared by all builds, against 14–37 GB of VCF.gz per build.

**Codec:** ZSTD needs `hyparquet-compressors` in the browser (fzstd). SNAPPY would avoid that
dependency but doubles the size, so ZSTD is worth it.

## 5. Open items / next measurements
- **WASM SIMD popcount kernel.** Target ≥ 200M pairs/s for ALL. Also tile aggregation for the pyramid view.
- **Re-check the GRCh38 and High Coverage builds.** High Coverage has 3202 samples and 6% duplicate positions on
  chr22, so tier sizes will shift.
- **Real network test** against S3/CloudFront: HTTP/2 concurrency, measuring wall time for 10 Mb and 100 Mb windows
  with progressive tier loading.
- **Build the full rsID index** from the 2023 dump (all 1.07B rows, plus merges) and confirm the partition footer sizes.

## 6. Step 1 results: v2 build (corrected tiers, rsID labels)

Built with `experiments/build_genotypes_v2.mjs`, benchmarked with `experiments/bench_v2.mjs`. Results are in
`experiments/results/v2_*.json`. Data: GRCh37 chr22, cached AsyncBuffer, local files.

**Build:** 1,109,368 rows in 2 min 49 s.
- **Tiers on sub-population `maf_max`:** common 218,038 rows, low 274,858, rare 616,472. With super-population tiers
  the common tier had 170,030.
- **rsID labels:**

| how the row was labelled | rows |
|---|---|
| VCF rsID | 500,435 |
| VCF rsID, merged → current | 2,190 |
| position + type | 573,155 |
| no candidate (null) | 33,588 (3.0%) |
| ambiguous (`rsid_candidates` > 1) | 3,703 |

- **rsID check:** where the VCF rsID is present at the position, position + type agrees on 502,299 of 502,625
  rows (99.94%).

**Size per tier and page size:**

| tier | 16 KB pages | 64 KB pages | 256 KB pages |
|---|---|---|---|
| common (bitset) | 41.8 MB | **32.1 MB** | 29.0 MB |
| low (bitset) | — | **8.0 MB** | 7.5 MB |
| rare (sparse) | **6.5 MB** | 6.4 MB | 6.4 MB |

The chosen pages (bold) total **46.6 MB for chr22**, including the new `rsid`, `function` and flag columns. v1 at
256 KB was 39.4 MB.

**Window scans against scattered point reads:**

| file | 10 Mb window | 30 scattered points | 300 scattered points |
|---|---|---|---|
| common 16 KB | 12.3 MB | 48 KB/pt | 13 KB/pt |
| **common 64 KB** | 9.1 MB | **36 KB/pt** | **17 KB/pt** |
| common 256 KB | **8.0 MB** | 79 KB/pt | 42 KB/pt |
| **low 64 KB** | 2.2 MB | 32 KB/pt | 7 KB/pt |
| low 256 KB | 2.1 MB | 55 KB/pt | 12 KB/pt |
| **rare 16 KB** | 1.7 MB | 46–58 KB/pt | 18 KB/pt |
| rare 256 KB | 1.7 MB | 151 KB/pt | 20 KB/pt |

**Choice:** 64 KB pages for common and low, 16 KB for rare.
- 64 KB costs +14% on common-tier windows and makes scattered points 2–2.5× cheaper.
- 16 KB rare pages cost +1% on windows and make sparse points 3× cheaper.
- A 10 Mb window is then about **9.1 MB for common only, or 13 MB for all tiers**.

**Planner rule of thumb:** scanning a window costs about 0.9 MB per Mb of common tier. A point costs 17–36 KB. So
scan a region once it holds more than about 25–50 requested variants per Mb; below that, read pages.

**Lessons for the reader:**
- **Pages don't line up across columns.** Parquet sizes pages per column in bytes, so a `pos` page (4 B/row) spans
  about 4,000 rows while a common `haps` page spans 104. Pruning on the `pos` column index alone reads about 550 KB
  per point. Point lookups must get exact row numbers from `pos`, then map them to the **`haps` column's own offset
  index**. `experiments/lookup.mjs` does this through its `payloadColumn` option.
- **No offset index on single-page chunks.** hyparquet-writer only writes an offset index when a column chunk has
  more than one page. A small chunk (such as rare-tier `haps` at 256 KB pages) then has no page index and is read
  whole.
- **Row groups are cut once per flush.** The builder starts a row group per 50k-row source chunk in every tier, so
  common row groups hold about 10k rows. That is fine for pruning; revisit it if footer size matters genome-wide.

## 7. Step 2 results: dbSNP indexes from LDlink's dump (parity)

Pipeline:
1. `dbsnp_bson2tsv.mjs` streams the 22.9 GB `dbsnp.bson.gz` from S3 in 45 minutes and keeps each document's dump
   order.
2. `dbsnp_sort.sql` (DuckDB, intermediates only) sorts it in 12.5 minutes.
3. `build_dbsnp_v2.mjs` writes the final files with hyparquet-writer in about 19 minutes, all three index sets in
   parallel.

Results are in `experiments/results/v2_dbsnp_*`.

**What the dump contains:**
- **Docs:** 1,073,486,098, all fields stored as strings, chromosomes 1–22, X and Y only.
- **rsIDs:** 1,073,106,380 distinct. **296,308 rsIDs have more than one doc** (379,718 extra docs), so "first doc in
  dump order" matters.
- **Merged docs** (`id ≠ ref_id`): 20,546,366.
- **Missing positions:** 21,961,126 docs have no GRCh37 position and 834,161 have no GRCh38 position.
- **9,208 anomalies:** impossible positions ≥ 1e9 (9,202 GRCh37, 6 GRCh38, e.g. chr11:69,716,399,231). They are
  kept in `anomalies.parquet`, and their keys are null in the indexes.
  - **Known parity difference:** for these rsIDs LDlink returns the bogus position and then fails to find the variant
    in 1000G, while we report "no position for this build". Either way no LD is computed; only the message differs.

**Index files:**

| index | files | rows | size | bytes/row |
|---|---|---|---|---|
| rsID (`part-<rsid>>26>`, rsid · ref_rsid · chrom · key37 · key38) | 32 | 1,073,486,098 (= all docs) | 9.98 GB | 9.3 |
| GRCh37 position (key · rsid · is_ref · type · function) | 24 | 1,051,524,972 (= docs − no GRCh37 position) | 4.73 GB | 4.5 |
| GRCh38 position | 24 | 1,072,651,937 (= docs − no GRCh38 position) | 4.82 GB | 4.5 |

Every file is checked for (key, dump order) ordering and for the exact row count of its intermediate.

**Manifest:** 83 files (dbSNP plus chr22 genotypes) in 414 KB of JSON, 80 KB gzip or 46 KB brotli. Per file it holds
the exact footer range, each row group's key range and its page-index locations.

**Parity check** (`verify_parity.mjs`): `parity_expected.sql` computes the expected answers straight from the dump
with Mongo semantics. The cases are sampled to stress the hard paths.

| check | cases | mismatches |
|---|---|---|
| rsID → position, GRCh37 / GRCh38 (`find_one`: first doc in dump order). Includes 1,123 merged, 523 duplicated and ~300 with no position per build | 3,700 each | **0 / 0** |
| `chr:pos` → rsID, GRCh37 / GRCh38 (`replace_coord_rsid`). About 1,370 have several docs and 1,235 must warn | 3,000 each | **0 / 0** |
| reference test-vector SNPs → same positions as the reference LD run | 10 each | **0 / 0** |

- The resolver (`resolve.mjs`) also reproduces LDlink's input handling: `validsnp` lower-cases and strips, rsIDs go
  through `strip("rs")`, coordinates through `strip("chr")`, and positions parse like Python's `int()`. Only
  chromosomes 1–22, X and Y are accepted.
- "Several docs at a position but none current" never occurs in this dump, so that branch of `replace_coord_rsid`
  can't be exercised.

**Reader lessons** (the fixes are in `experiments/lookup.mjs`):
- **`onChunk` can return rows outside `[rowStart, rowEnd)`.** Whole pages come back, and chunks can overlap. Place
  values by absolute row number. The first index build wrote duplicate rows at every 1M-row chunk boundary (chr1 had
  91.8M rows instead of 81.9M); the order and row-count checks caught it.
- **Fetch the column index and the offset index separately.** Writers put all column indexes, then all offset
  indexes, at the end of the file. One span covering both read about 1 MB per key.
- **Bound concurrency** to 16 reads in flight. Each concurrent `parquetRead` allocates row-group-sized decode buffers,
  and 3,700 unbounded lookups peaked around 5 GB of heap.
- **Cold rsID lookups cost about 60 KB per key** when the keys are spread across row groups. Most of that is the
  5 columns' offset indexes, fetched once per row group touched. To cut it, drop `chrom` (it can be derived from the
  keys) or keep it in a separate file.

## 8. Step 3 results: GRCh38, High Coverage and chrX builds

These use the same v2 builder with the chosen layouts (`FINAL=1`: common 64 KB, low 64 KB, rare 16 KB). Inputs are
the local copies of the S3 VCFs plus dbSNP slices exported from the dump. The builder now loads dbSNP into typed
arrays, placing values by row number, which is needed for chrX's 42M dbSNP rows. Results are in
`experiments/results/v2_{grch38_chr22,hc_chr22,grch37_chrX}_*` and `results/ld_check/`.

| build | rows | size (common / low / rare) | no rsID | ambiguous | build time |
|---|---|---|---|---|---|
| GRCh37 chr22 (§6) | 1,109,368 | 46.6 MB (32.1 / 8.0 / 6.5) | 3.0% | 0.33% | 2.8 min |
| GRCh38 chr22 | 1,059,079 | **44.7 MB** (30.7 / 7.8 / 6.2) | 6.0% | 0.25% | 3.1 min |
| High Coverage chr22 | 1,070,401 | **66.4 MB** (47.8 / 13.3 / 5.3) | 8.6% | 2.4% | 3.5 min |
| GRCh37 chrX | 3,501,213 | **112.0 MB** (69.5 / 24.6 / 17.9) | 8.6% | 0.34% | 7.5 min |

- **rsID sources:**
  - GRCh38 and High Coverage have no VCF rsIDs, so every label comes from position + type.
  - On chrX, VCF rsIDs (1,393,810 direct + 11,501 merged → current) agree with position + type 99.88% of the time
    (1,403,564 of 1,405,311).
- **High Coverage is about 50% larger** than GRCh38 for the same chromosome. It has 3202 samples (801-byte bitsets)
  and more low-frequency variants, and its common tier has 265k rows against 213k.
- **GRCh38 VCF samples:** the VCF has 2548 samples, while the panels list 2504 (NA18498 is missing from the VCF).
  The 45 extra samples are stored but never selected. Dropping them would save about 2%.
- **chrX ploidy is exact:**
  - `male_haploid_samples` lists 1,233 samples, the same number of males as the chrY VCF.
  - 3,389,180 rows are flagged `male_haploid`, all outside PAR1 (60,001–2,699,520) and PAR2
    (154,931,044–155,260,560), and every PAR row is diploid. That's 0 violations either way.
- **Window and point costs:**

| file set | 1 Mb window, common | 10 Mb window, all tiers | 300 points |
|---|---|---|---|
| GRCh38 chr22 | 2.9 MB | 9.1 + 2.3 + 1.7 = 13.1 MB | 7–18 KB/pt |
| High Coverage chr22 | 4.8 MB | 15.4 + 3.7 + 1.4 = 20.6 MB | 12–17 KB/pt |
| GRCh37 chrX | 0.9 MB | 5.3 + 1.6 + 1.1 = 8.0 MB | 21–38 KB/pt |

**End-to-end LD check:** `ld_check.mjs` computes r² and D′ for 25 common SNVs (25×24 off-diagonal pairs) two ways.
One uses the Parquet bitsets, including the haploid-male mask. The other uses the reference LDlink math
(`scripts/profile/ldmatrix_reference.mjs`) on the VCF.

| build | population | region | mismatches |
|---|---|---|---|
| GRCh37 chr22 | CEU | 30–32 Mb | 0 / 600 |
| GRCh38 chr22 | CEU, EUR | 30–32, 40–42 Mb | 0 / 600 each |
| High Coverage chr22 | CEU, AFR | 30–32, 40–42 Mb | 0 / 600 each |
| GRCh37 chrX | CEU, EUR | 30–32 Mb (haploid males), 1.0–2.9 Mb (spans the PAR1 boundary), 100–102 Mb | 0 / 600 each |

That's 4,800 pairs in total, all identical to 3 decimals. This confirms the bitset encoding, the per-build sample
order, the population masks and the chrX ploidy handling.

**Genome-wide size estimates:** chr22's Parquet size scaled by each build's VCF.gz size ratio.
- GRCh37: about 3.8 GB
- GRCh38: about 3.4 GB
- High Coverage: about 4.8 GB

Add the shared dbSNP indexes (19.5 GB) to that. chrX compresses better than the autosomes because males are
haploid: 112 MB from a 1.92 GB VCF.

## 9. Step 4 results: the LD kernel (WebAssembly)

Source: `experiments/wasm/ld_kernel.c`, built by `wasm/build.sh` with emsdk clang into standalone modules (no libc,
no JS glue):
- **`ld_kernel_scalar.wasm`** (6.3 KB): uses `i64.popcnt`.
- **`ld_kernel_simd.wasm`** (7.2 KB): uses `v128.and` + `i8x16.popcnt`.

The JS wrapper is `wasm/ld_kernel.mjs`. Results are in `experiments/results/wasm/`.

**What the kernel does:** the wrapper packs the selected population's haplotypes from the v2 rows (bitset or sparse
tier, missing calls, haploid-male mask) into W 64-bit words per variant. Per-variant alt counts are computed once.
For each pair only `x = popcount(alt_i ∧ alt_j)` is needed, then r² and D′ come from LDlink's formulas in f64,
counted by allele index. When valid masks differ (chrX across a PAR boundary, `.` calls), `band_masked` uses four
popcounts per pair.

**Entry points:**
- `band`: distance-banded pairs above an r² or D′ threshold, as sparse output (i, j, r², D′).
- `band_blocked`: the same result, cache-blocked with 4-row register tiling and a row range so workers can split
  the work.
- `tile`: dense r²/D′ blocks quantized to bytes (0–250, 255 = NA) for the matrix and pyramid views.
- `pair_counts`: exact counts for a list of pairs.

**r² bound pruning:** with minor allele frequencies m₁ ≤ m₂, r² can't exceed m₁(1−m₂)/(m₂(1−m₁)), so a pair is
skipped when that bound is below the threshold. The bound equals r² exactly for perfectly coupled pairs. Without a
1e-9 relative margin, rounding dropped 23 pairs sitting exactly at r² = 0.2; with the margin, pruning is lossless.

**Correctness** (`bench_kernel.mjs`, `bench_blocked.mjs`, `bench_threads.mjs`):
- **Counts and statistics:** 0 count mismatches and 0 r²/D′ mismatches against an independent JS reference that
  counts straight from the Parquet bitsets. That's 3,000 random in-band pairs × 3 datasets × 2 kernels: GRCh37
  chr22 common with CEU; all three tiers merged, including sparse rows, with EUR; and chrX across PAR1 with EUR.
- **Pruning:** results with and without pruning are identical in every configuration.
- **Blocking:** blocked and unblocked bands give identical (i, j, r²) sets across 18 configurations (3 populations ×
  2 kernels × 3 block sizes).
- **Threads:** the union of worker results equals the single-thread result (same pair count, kept count and r² sum).

**Speed:** single thread, 10 Mb of GRCh37 chr22 common variants (59,372), 500 kb band, 175M pairs. Timings on this
machine vary by up to about 2× between runs.

| kernel | CEU (198 haps) | EUR (1,006) | ALL (5,008) |
|---|---|---|---|
| JS, 32-bit table popcount (step 1 code) | 35 Mpairs/s | 11 | 1–2 |
| WASM scalar, unblocked | 103–143 | 49–55 | 16–18 |
| **WASM scalar, blocked (B = 128)** | **178–224** | **120** | **34–37** |
| WASM SIMD, blocked (B = 128) | 105–179 | 96 | 30 |
| scalar blocked + r² ≥ 0.8 pruning (whole 10 Mb band) | 0.5 s | 0.6 s | **1.1 s** |

- **Scalar beats SIMD.** V8 lowers `i8x16.popcnt` to a multi-instruction byte-shuffle sequence, while x86 has a native
  64-bit `popcnt`. Ship the scalar kernel and keep SIMD as an option to re-measure on other CPUs.
- **Pruning skips 91–96% of pairs** at r² ≥ 0.8, and 62–81% at r² ≥ 0.2.
- **Tiles:** a 1000 × 1000 tile takes 10–70 ms, so matrix and pyramid views can render interactively.
- **Workers** (`worker_threads`, the same pattern as Web Workers), scalar blocked kernel, 10 Mb, 500 kb band:

| threads | EUR kernel time | ALL kernel time |
|---|---|---|
| 1 | 3.2 s | 12.4 s |
| 4 | 1.2 s | 4.3 s |
| 8 | **0.61 s** | 2.4 s |
| 16 | 0.77 s | **1.65 s** |

These single-thread times came from a slower run of this machine (about 2× the blocked benchmark). The scaling
ratios are the reliable part: about 5× at 8 threads and 7.5× at 16. Rows are split evenly by count, not by pair
count, which leaves some load imbalance.

**Takeaways for the engine:**
- Use the scalar blocked kernel in Web Workers, with `navigator.hardwareConcurrency` workers capped at about 8.
- Prune whenever there is a threshold.
- Use `tile` for views, and `band` (sparse pairs above a threshold) for tables and export.
- A 10 Mb window of common variants with a 500 kb band then takes about 1 s for EUR and 2–3 s for ALL unthresholded
  across 8 workers, or about 1 s for ALL at r² ≥ 0.8 on one thread.
- **Open items:**
  - balance worker ranges by pair count
  - a mask-grouped fast path for chrX, so that only pairs straddling the PAR boundary take the masked path
  - WebGPU, if ALL-population tiles over 100 Mb windows are needed

## 10. Step 5 results: query planner

`experiments/planner.mjs` takes a variant set (rsIDs, `chr:pos`, regions) and does three things:
1. Resolves it with LDlink parity.
2. Groups the positions per chromosome into clusters (gaps under 2 Mb).
3. For each cluster and tier file, chooses a **region scan** or **page-level point reads**, whichever has fewer
   estimated bytes.

It costs everything from `manifest.json` alone, without opening a file. The manifest now carries compressed bytes
per row group (and per column for genotype files): 627 KB, or 129 KB gzipped. Validation is `verify_planner.mjs`,
with results in `experiments/results/v2_planner_report.json`. Every scenario runs cold and measures every byte read,
broken down by file group.

**rsID resolution, exact and hybrid:**
- **Direct:** used for 32 rsIDs or fewer, for the 296,308 rsIDs with several dump documents (listed in
  `dbsnp/duplicated_rsids.parquet`, 374 KB), and whenever the hybrid is estimated to cost more.
- **Hybrid:** resolve a sample of 32 rsIDs directly to find their regions. Then scan the genotype files' (`pos`,
  `rsid`, `rsid_source`) columns over those regions and match the rest in memory. Only labels that guarantee a dbSNP
  document for that rsID at that position count: `rsid_source` 1 (VCF) or 3 (position + type). Leftovers go direct.
- **Results:**

| inputs | strategy | differences vs direct | bytes (vs direct) |
|---|---|---|---|
| 325 rsIDs clustered in 2 Mb, with duplicated, merged and absent rsIDs | hybrid | **0** | **6.6 MB** (19.1 MB) |
| 3,025 rsIDs clustered in 2 Mb, with the same extras | hybrid | **0** | **7.4 MB** (81.3 MB) |
| 300 rsIDs scattered over chr22 + chrX | direct (chosen) | 0 | 26.8 MB (26.6 MB) |

**Genotype fetch:**

| scenario | chosen | predicted | measured | alternative (measured) |
|---|---|---|---|---|
| window 1 Mb, all tiers | scan | 1.48 MB | 2.19 MB | — |
| window 10 Mb, common | scan | 8.12 MB | 9.28 MB | — |
| window 10 Mb, all tiers | scan | 11.6 MB | 13.3 MB | — |
| 30 coordinates scattered on chr22 | points | 2.93 MB | 2.68 MB | scan: 33.7 MB |
| 300 coordinates in 2 Mb | points | 3.01 MB | 3.17 MB | scan: 5.65 MB |
| 3,000 coordinates in 2 Mb | points | 4.06 MB | 5.43 MB | scan: 7.0 MB |

The coordinate rows include LDlink-parity resolution through the dbSNP position index (0.7–3.3 MB). In every case
checked, the planner's choice was the cheaper one when measured, and both strategies returned the same rows.

**How the cost model is built:**
- **Scans:** column bytes prorated by position span, plus a fixed per-file cost (exact footer length plus offset
  indexes of the touched row groups).
- **Points:** expected distinct pages hit, P·(1 − (1 − 1/P)ᵏ), times the average `haps` page, plus other columns and
  the fixed cost.
- **Direct rsIDs:** footer per partition touched, plus offset indexes per row group touched, plus about 8 KB per key.
- **Coordinates:** the same structure over the dbSNP position index.
- **Accuracy:** genotype predictions are 0.6–1.5× of measured, rsID predictions 0.7–0.85×, and coordinate predictions
  0.56–0.69× (consistently under). That's good enough to choose correctly, which is all the planner needs.

**Open items:**
- **Duplicate positions across tiers.** When one position has several rows in different tiers (split multi-allelic
  sites in GRCh38 and High Coverage), parity's "last row in VCF order" needs a `vcf_order` column. Today order is
  kept only within a tier. ✅ Resolved: `vcf_order` is in the genotype files (§16).
- **Large scattered rsID lists cost about 45 KB per rsID** cold, mostly per-row-group offset indexes. A slimmer rsID
  index (drop `chrom`; move `ref_rsid` to its own file) would cut that.
- **Calibrate on real HTTP:** these are local-file byte counts. Request counts and latency will shift the
  scan-vs-points break-even.

## 11. Step 6 results: over the network (S3, HTTPS Range requests)

**Setup:**
- **Bucket:** `s3://eaap-ac-test-bucket/ldlink-parquet/v2/`, AWS account 285397596138 (profile `eaap`), us-east-1.
  - It is **private**: all four block-public-access settings are on, with SSE-S3 (AES-256) default encryption.
  - The upload was the whole release, 99 objects and 19.8 GB, in 2 min 52 s.
- **Client:** Node `fetch` from this machine through presigned GET URLs (`presign.sh`, 12 h expiry). Each file is an
  `AsyncBuffer` making Range requests. Byte lengths come from the manifest, so no HEAD requests are needed.
- **Method:** every scenario runs cold. Baseline round trip is about 51 ms median (33–220 ms). `MAX_INFLIGHT` caps
  concurrent requests to model browser connection limits.
- **Results:** `experiments/results/v2_net_bench*.json`.

**Query time** (planner, cold, including the manifest fetch of about 50 ms):

| query (GRCh37) | requests | MB | unlimited | ≤ 100 in flight | ≤ 6 in flight |
|---|---|---|---|---|---|
| window 1 Mb, all tiers | 12 | 2.8 | 0.58 s | 0.29 s | 0.30 s |
| window 10 Mb, common | 24 | 10.4 | 0.41 s | 0.42 s | 0.52 s |
| window 10 Mb, all tiers | 61 | 13.9 | 0.86 s | 0.97 s | 0.90 s |
| 300 rsIDs clustered in 2 Mb (hybrid) | 441 | 6.4 | 1.5 s | 1.5 s | 4.3 s |
| 3,000 rsIDs clustered in 2 Mb (hybrid) | 353 | 7.3 | 1.2 s | 1.2 s | 3.7 s |
| 300 rsIDs scattered over chr22 + chrX (direct) | 2,843 | 27.1 | 4.4 s | 4.3 s | **19.5 s** |
| 30 coordinates scattered on chr22 | 305 | 3.1 | 1.0 s | 0.8 s | 2.2 s |
| 300 coordinates in 2 Mb | 292 | 3.8 | 1.2 s | 1.1 s | 2.1 s |

- **Window queries** read 16–26 MB/s with few requests, and are interactive even with only 6 connections.
- **End to end:** fetching a 10 Mb window of common variants then computing LD (EUR, 500 kb band, r² ≥ 0.8, pruned,
  one thread) is 0.4 s + 1.5 s. That's 8.0M pairs evaluated after pruning 167M, keeping 713,124 pairs.
- **Byte predictions** hold over HTTP: measured is 0.8–1.0× of predicted for every query except the 1 Mb window
  (1.9×).

**What this means for the deployment:**
1. **Put CloudFront (HTTP/2 or HTTP/3) in front of the bucket.**
   - S3's REST endpoint only speaks HTTP/1.1, so a browser gets about 6 connections per host. At 6 in flight, list
     queries slow by 3–4.5×, and 300 scattered rsIDs take 19.5 s.
   - HTTP/2 multiplexing (~100 streams) brings them back to 0.8–4.3 s.
   - CloudFront also gives the browser a CORS-enabled, cacheable origin. A test bucket used directly from a browser
     would also need a CORS rule allowing `GET` with the `Range` header.
2. **Cut requests per lookup.** Direct rsID lookups cost about 9.5 requests each: one per column page plus offset
   indexes. Options:
   - a slimmer rsID index (drop `chrom`, move `ref_rsid` to its own file)
   - coalescing adjacent column ranges within a row group in the reader
   - shipping per-partition offset indexes in the manifest, or a sidecar, so cold lookups skip them
3. **Keep window queries as the main path.** They're already within budget: under 1 s for 10 Mb, all tiers, cold.

## 12. Step 7 results: annotation tracks (position-keyed, LDlink parity)

Sources: LDlink's 2023 mongodump (`refGene_grch37/38`, `recomb`, `regulome`, `forge_score`, `genes_name_coords`,
`gtex_tissue_eqtl`, `gtex_tissues`, `gtex_genes`, `snp_col`, `platforms`) and the **2025-10-24 backup of the live
`gwas_catalog`** (835,632 docs; newer than the dump's 462,053).

Pipeline:
1. `fetch_annotations.sh` streams each collection from S3 in 11 minutes. It uses `bson_stream.mjs`, plus
   `json_array_stream.mjs` for the GWAS backup, keeps dump order, and keeps only the fields LDlink reads.
2. `annotations.sql` (DuckDB) builds the key and sorts by (key, dump order).
3. `build_annotations.mjs` writes the files with hyparquet-writer. Every file is checked for key order and exact row
   count.

Results are in `experiments/results/v2_annotations_*`.

**Files.** Rows not keyable (position `NA`, chromosome 0/XY/M, impossible positions) are counted in
`v2_annotation_counts.csv`.

| file (per build) | rows (GRCh37 / GRCh38) | size | key / LDlink query |
|---|---|---|---|
| `annotations/<b>/recomb.parquet` | 3,283,362 / 3,282,596 | 28 MB | key → rate; `getRecomb` range |
| `annotations/<b>/regulome.parquet` | 16,357,032 / 16,334,982 | 24 MB | key → score; `get_regDB` find_one (first in dump order) |
| `annotations/<b>/forge.parquet` | 37,953,042 / 37,949,750 | 211 MB | key → rsid, score, placed at each rsID's LDlink position |
| `annotations/<b>/genes.parquet` | 78,105 / 78,362 | 4.3 MB | key = chrom·1e9 + txStart; `getRefGene` overlap |
| `annotations/<b>/gwas.parquet` | 835,534 / 835,562 | 28 MB | key → LDtrait fields + study metadata; window query |
| `annotations/<b>/gtex_eqtl.parquet` | 70,955,246 / 71,478,479 | 1.28 GB | key → gene_id, tissue, maf, pval, slope, slope_se, tss_distance |
| `annotations/<b>/arrays.parquet` | 10,393,379 / 10,383,296 | 79 MB | key → platform list (pipe-joined, doc order); SNPchip |
| `annotations/forge_rsid.parquet` | 37,962,826 | 25 MB | rsid → score; `get_forgeDB` find_one (exact parity) |
| `annotations/gene_names.parquet`, `gtex_tissues`, `gtex_genes`, `platforms` | 26,195 / 54 / 56,200 / 83 | < 1 MB each | lookups |

In total that's 3.2 GB for both builds, and GTEx is 80% of it. The manifest now covers 112 files: 677 KB, or 140 KB
gzipped.

**Parity** (`parity_annotations.sql` + `verify_annotations.mjs`): the expected answers are computed from the **raw
stored strings** with Mongo/LDlink query semantics, independently of the key-building code.

| check (both builds) | cases | mismatches |
|---|---|---|
| RegulomeDB `find_one({chr, pos})`, including multi-doc positions and absent positions | 1,900 × 2 | **0** |
| recombination range: count, positions in order, rate sum | 100 × 2 | **0** |
| refGene overlap: transcript set (LDlink's query has no sort, so it's a set comparison) | 100 × 2 | **0** |
| GWAS window ±500 kb: records in index order | 100 × 2 | **0** |
| GTEx: tissue set + position, inner-joined to tissues and genes (`$unwind`) | 150 × 2 | **0** |
| SNPchip: every doc's platform list at a position, in dump order | 1,100 × 2 | **0** |
| FORGEdb `find_one({snp_id})`, including duplicated and absent rsIDs | 2,050 | **0** |

The first run reported one GWAS "mismatch". It was an artifact of the check: SQL's `'a' || NULL` is NULL, and one
record has an empty trait. The fixed check uses `coalesce`.

**Window cost** (GRCh37 chr22, cold, local byte counts, using `lookup.mjs` `rangeScan`):

| track | 1 Mb window | 10 Mb window | rows in 10 Mb |
|---|---|---|---|
| recombination | 107 KB | 209 KB | 13,737 |
| RegulomeDB | 106 KB | 183 KB | 59,987 |
| FORGEdb | 187 KB | 826 KB | 134,169 |
| genes (overlap; scan starts 2.5 Mb left) | 557 KB | 694 KB | 432 |
| GWAS | 1.6 MB | 1.6 MB | 3,055 |
| GTEx eQTL | 523 KB | 2.7 MB | 256,563 |
| arrays | 5.0 MB | 5.2 MB | 46,920 |

Each costs 6–20 requests and 5–200 ms of decoding.

**Reader lesson: hyparquet's `filter` option prunes row groups, not pages.** With 1M-row row groups it read whole
column chunks: the same 1.5–13 MB for 1 Mb and 10 Mb windows. `rangeScan` picks pages from the key column's page index
and reads them through the offset index, cutting bytes 10–25×. The planner should use it for scans as well.

The remaining fixed costs are GWAS dictionary pages (one 835k-row group) and arrays' long platform-list strings. A
per-row platform bitset would shrink arrays, but SNPchip needs each doc's platform order, so the strings stay for now.

## 13. LDscore (inventory; follow-on work)

LDlink's LDscore family is three modules:
- **LDscore** runs `ldsc --l2` on a user-uploaded PLINK `bed/bim/fam`, with a cM or kb window.
- **LDherit** runs ldsc's h² regression on summary statistics.
- **LDcorrelation** runs ldsc's rg regression on summary statistics.

All three call the external `ldsc` package (`CBIIT/ldsc@ldsc313`, `ldscore.ldsc_utils`). Their reference data lives
only on LDlink's EFS volume, `/data/ldscore/`. It was inventoried and copied with a temporary ECS task, tagged
`ldlink-ldscore-probe`, which was stopped afterwards:
- **Copied as one tarball:** `s3://eaap-ac-test-bucket/ldlink-source/ldscore_efs.tar`, 1,818,624,000 bytes, sha256
  `708f46a4…`, verified locally after download.
- **Contents (1.7 GB):**
  - **Reference LD scores** per ancestry group: `afr`, `amr`, `asj`, `eas`, `est`, `fin`, `nfe`, `nwe`, `seu`,
    plus `eur`, `eas_fak`, `eas_copy_eur`, `seu_all` and an older nested `ldscore/` copy.
    - Each group has `<chr>.l2.ldscore.gz` (columns `CHR SNP BP L2`, e.g. 493,515 SNPs for amr chr12) and
      `<chr>.l2.M_5_50`.
    - `eur` is the standard Broad `eur_w_ld_chr` set: 1000G EUR, r² summed over HapMap3 SNPs, HLA handled. Most of the
      others are gnomAD-style ancestry LD scores, **not computable from the 1000G panel**, so they must be converted
      as data.
  - `w_hm3.snplist` (1,217,312 HapMap3 SNPs with A1/A2).
  - LDscore example inputs: `22.bed/bim/fam` (379 individuals, 19,156 SNPs) with its ldsc 3.0.1 output
    (`--l2 --ld-wind-cm 1.0`: `22.l2.ldscore.gz`, `.M`, `.M_5_50`, `22.log`). These are ready-made test vectors.
  - LDherit/LDcorrelation examples: `BBJ_HDLC22.txt`, `BBJ_LDLC22.txt`, and the full `BBJ_HDLC.txt`.

**What browser LDscore needs** (done, see §14):
1. **Convert the reference LD scores to Parquet:** per group, a key plus SNP, L2 and M_5_50 (in metadata).
2. **A genotype LD score kernel for user `bed` files.** ldsc uses *unphased genotype* correlation with its N-adjusted
   r², not our phased haplotype r², so it needs a second kernel (2-bit genotype planes plus popcount) and a JS `.bed`
   reader. Verify it against the `22.*` example output.
3. **Port ldsc's h² and rg regressions** (weighted least squares with block jackknife) to JS. Verify them against
   LDlink's results for the BBJ examples.


## 14. Step 9 results: LDscore / LDherit / LDcorrelation in JS (ldsc parity)

All three modules now run in plain JS, with no Python or ldsc. The code is in `experiments/ldsc/` and is a
line-for-line port of `CBIIT/ldsc@ldsc313`. It is checked against the reference ldsc, run locally in a Python 3.13
venv at `D:\Projects\LDlink\ldsc-ref` with `ref_dump.py` dumping full-precision JSON, and against LDlink's stored
example outputs.

**Reference panels → Parquet** (`convert_ldscore_ref.mjs`; `release/ldscore/`, 1.30 GB, 14 files):

| file | rows | size | M (Σ M_5_50) |
|---|---:|---:|---:|
| afr | 17,447,760 | 199.5 MB | 2,060,476 |
| amr / asj / est / nfe / nwe | 9.9–11.2 M | 112–131 MB | 0.87–1.55 M |
| eas | 9,222,652 | 106.1 MB | 1,475,758 |
| fin / seu / seu_all | 8.9–9.6 M | 103–107 MB | 0.38–0.76 M |
| eur, eas_fak, eas_copy_eur | 1,290,028 | 27.4 MB each | 1,173,569 |
| w_hm3.snplist | 1,217,311 | 5.7 MB | — |

- **Columns:** `CHR` INT32, `SNP`, `BP` INT32, `L2` DOUBLE, plus `CM`/`MAF` where the source has them.
- **Row order:** rows are sorted stably by (CHR, BP), which keeps the source order. Duplicate SNPs are kept, and
  readers apply ldsc's `drop_duplicates(SNP)` (keep first).
- **Key-value metadata:**
  - `M_5_50` and `M` per chromosome, as JSON
  - `chromosomes`
  - `source`
  - `ldsc_order`
- **Manifest:** keyed by `CHR`, so a client picks chromosomes by row group.
- **Note:** `eur`, `eas_fak` and `eas_copy_eur` are byte-identical copies in LDlink. The release keeps them as
  aliases, and a client could share one file.

**`--l2` for user PLINK files** (`ldscore_l2.mjs`):
- The port covers:
  - the `.bed` reader
  - the MAF and missingness filter
  - standardization (missing → mean, ddof = 0)
  - `getBlockLefts`
  - the chunked `__corSumVarBlocks__` with c = 50
  - `l2_unbiased`
  - Python `%.3f` output: ties go to even, and negative zero prints as `-0.000`
- **Results:**
  - LDlink's example (`22.bed`, 379 × 19,156, `--ld-wind-cm 1`): `22.l2.ldscore.gz`, `.M` and `.M_5_50` are
    **byte-identical** to LDlink's stored output. Runtime is 2.9 s in Node.
  - Also identical to the reference for `--ld-wind-kb 500`, `--ld-wind-cm 2`, and an input with missing genotypes
    plus monomorphic SNPs.
- **This is not the popcount kernel proposed in §13.** It is the float64 dot-product algorithm, which is what exact
  parity needs. A 2-bit popcount kernel is an optional speedup, to be checked against this port.

**`munge_sumstats`** (`munge.mjs`), with LDlink's flags `--merge-alleles w_hm3.snplist`:
- The port covers:
  - header cleaning and column-name mapping
  - exactly one signed statistic
  - the pandas NA strings
  - all row filters, in order: NA, merge-alleles, INFO, FRQ, P, alleles
  - `drop_duplicates`
  - `N ≥ quantile(0.9)/1.5`
  - the median sign check
  - Z = √chi2.isf(P, 1), computed with a log-space Newton step on Cody `erfc`
  - the allele merge in w_hm3 order
  - `%.3f` output
- **Results:** for `BBJ_HDLC22.txt` and `BBJ_LDLC22.txt`, the decompressed `.sumstats` text is **identical** to the
  reference (ignoring Python's CRLF line endings). That is 1,217,311 rows each with 14,270 non-missing, and both
  runs take about 1 s.

**h² / rg regressions** (`ldsc_regress.mjs`):
- The port covers:
  - IRWLS with 2 reweighting rounds
  - the two-step estimator (h² only)
  - the fast block jackknife with 200 blocks
  - the ratio jackknife for rg
  - the allele match and flip
  - `h2_obs_to_liab`
- **Results** (the full pipeline, raw `BBJ_*22.txt` → JS munge → JS regression, eas panel, against ldsc):

| quantity | JS | ldsc | rel. diff |
|---|---:|---:|---:|
| h² (HDL) | 0.0816302013237094 | 0.0816302013237097 | 4e-15 |
| h² se | 0.0348314748009528 | 0.0348314748009528 | 4e-16 |
| rg (HDL, LDL) | 0.396291483820782 | 0.396291483820785 | 6e-15 |
| rg se | 0.255960815765345 | 0.255960815765347 | 7e-15 |
| rg p | 0.121561975672461 | 0.121561975672461 | 4e-15 |

- Also matching to ≤ 3.3e-14:
  - h² on the eur panel
  - the liability-scale factor (P = 0.3, K = 0.1 → 1.252332750492982)
  - every `hsq1`, `hsq2` and `gencov` field
- **Reproduced ldsc quirk:** `estimate_rg` never applies the default `--two-step 30`. When it checks
  `args.intercept_h2 is None`, that value has already been turned into `[None, None]`. The JS rg therefore
  defaults to `twostep: null`. Fixing this would change LDlink's numbers, so it is a deliberate parity choice.
- **Timing:** rg takes about 6.5 s in Node (three regressions over about 14k SNPs plus the jackknife).

**Published:** `s3://eaap-ac-test-bucket/ldlink-parquet/v2/ldscore/` holds 15 objects, 1.30 GB. The manifest has 126
files, and the release totals 24.4 GB.

## 15. Step 8 results: NCBI dbSNP rebuild vs. LDlink's 2023 dump (b155 exact; b157 drift)

**Sources** (in `experiments/out/ncbi_dbsnp/`; `fetch_ncbi_dbsnp.sh`):
- b155 VCFs: `GCF_000001405.25.gz` (GRCh37, 23.2 GB) and `.39.gz` (GRCh38, 24.0 GB), plus `refsnp-merged.json.bz2`
  (784 MB) and `refsnp-chr22.json.bz2` (5.5 GB). NCBI publishes no md5 files for b155, so each file was checked
  against the server's Content-Length; the gzip/bzip2 streams also decoded cleanly.
- b157 (2025-01, the latest release): `.25.gz` (28.2 GB) and `.40.gz` (29.6 GB) with `.tbi`, plus the merge history.
  All passed md5.

**How LDlink built the collection** (`scripts/dbSNP_mongo_update`, b155, July 2021). It reads the per-chromosome
refSNP **JSON**, not the VCFs. For each RefSNP it writes one doc per rsID in `[refsnp_id] + dbsnp1_merges`, with
`ref_id = refsnp_id`:
- **chromosome:** taken from the file name.
- **position:** the first number in the reference allele's HGVS on the GRCh37.p13 / GRCh38.p13 chromosome
  placement. If there are several placements, the last one wins.
- **type:** `variant_type`.
- **function:** the protein sequence-ontology term of `genes[0].rnas[0]` for each allele, with `_variant` removed,
  joined in Python **set order**, which is hash-randomized. Compare it as a set.
- **Merged rsIDs:** `rsjson_process_merges.py` then adds, for each merged rsID not already present, a copy of
  `find_one(target)` with `ref_id = target`.

**Exact reproduction on chr22.** `compare_ldlink_parse.sql` re-runs LDlink's unmodified parser on the b155 chr22 JSON
(112 min) and compares the result with the dump:

| | records |
|---|---:|
| dump chr22 docs | 14,638,849 |
| parser output, field-identical to the dump (function as a set) | **14,499,527** (all of them) |
| merged-rsID copies (`id ≠ ref_id`) that are valid `refsnp-merged` pairs and equal the target's primary doc | 139,134 |
| …where the target is itself a merged copy (so `find_one` returned that copy; positions equal) | 131 |
| unexplained | **0** |

The one thing NCBI data cannot reproduce is **dump order.** The dump holds docs in Mongo storage order: chr22 docs
are spread over 662 M dump positions, and only 4% of consecutive parser records stay adjacent. Parity mode's
"first doc" rule therefore still needs the dump, which confirms the decision to keep it.

**Genome-wide, b155 VCF vs. dump** (`ncbi_vcf2tsv.mjs` then `compare_ncbi.sql`; primary docs, one per rsID):
- **rsIDs:**
  - All **1,052,939,732** primary dump rsIDs are in the b155 VCFs; 0 are dump-only.
  - 87,820 rsIDs are only in NCBI (mostly SNVs, clustered, e.g. at chr1:142–143 Mb). They are probably in
    `refsnp-other`/`unsupported` JSON, which LDlink never parsed.
- **Type:** identical for all of them.
- **`NA` positions:** an `NA` position in the dump means exactly "absent from that assembly's VCF"; there are no
  one-sided cases.
- **Positions.** The rule is **POS for SNV/MNV and VC=INS; POS+1 for every other anchored record** (deletions and
  delins, including insertions inside repeats, because the HGVS of the REF allele covers the repeat):

| | equal | differ | of which multi-placement rsIDs |
|---|---:|---:|---:|
| GRCh37 | 1,030,991,578 | 55,571 (0.005%) | 231 |
| GRCh38 | 1,052,083,827 | 18,126 (0.002%) | 17,819 |

- **Remaining GRCh37 differences** are multi-allelic repeat indels: NCBI's REF-allele HGVS starts somewhere other
  than the VCF anchor. They also include a **LDlink parser bug.** `re.sub(r"\D", "", …)` turns HGVS repeat notation
  such as `g.16907955AT[11]` into `1690795511`. At least **4,830 rsIDs have positions beyond the end of any
  chromosome** (4,822 GRCh37, 12 GRCh38; e.g. rs1569175796 at chr22:1,690,795,511). Corrupt values that happen to
  land inside a chromosome are not counted. Parity mode reproduces these as stored.
- **Chromosome:** 1,320,234 rsIDs differ, and all have several placements:
  - 1.29 M are X/Y PAR: LDlink says Y because it reads the chromosome from the chrY file name, while the VCF lists
    X first.
  - The rest are GRCh38 multi-chromosome placements.
- **function vs. VCF flags,** restricted to the 4 terms both can express (missense, synonymous, stop_gained,
  frameshift): 1,050,782,553 equal. In 2,059,379 + 97,199 cases the VCF flags more, because they cover any
  transcript while LDlink keeps only the first gene/RNA. Only 601 are dump-only.
- **Merges:**
  - 18,914,049 of the dump's merged pairs are `refsnp-merged` pairs.
  - In 1,629,716 the dump's target is another rsID in the merge chain (it came from `dbsnp1_merges` of the
    current record).
  - 0 merged pairs are unexplained.

**Allele-aware join on chr22** (`allele_join.mjs`). It compares minimal-form REF/ALT from NCBI b155 with the 1000G
release rows. There are two checks.

*(a) The release's rsID labels* (finding 4's rule: VCF rsID, else position + type):

| GRCh37 rows (1,109,368) | agree | disagree | no allele match | unlabelled, allele match exists |
|---|---:|---:|---:|---:|
| label from VCF rsID (source 1/2) | 502,201 | 400 | 24 | — |
| label from position + type (source 3) | 562,718 | 10,333 (1.8%) | 104 | — |
| no label (source 0) | — | — | 959 | 32,629 |

GRCh38 (1,059,079 rows, all source 3 or 0): 982,896 agree, **9,278 disagree (0.9%)**, 3,414 have no allele match,
and 48,451 unlabelled rows have an allele-matching rsID.

*(b) LD impact.* For each b155 rsID on chr22, LDlink resolves it to the dump position, then takes the **last tabix
row overlapping that position** (`LDcommon.py:590-603`). The table compares that row with the row whose alleles
match the rsID:

| rsIDs that LDlink maps to a 1000G row | GRCh37: 1,242,176 | GRCh38: 1,202,367 |
|---|---:|---:|
| same row as allele-aware | 1,087,531 (87.5%) | 1,020,075 (84.8%) |
| a different 1000G row carries the rsID's alleles (LDlink picks the wrong one) | 2,814 | 5,163 |
| **no 1000G row has the rsID's alleles,** yet LDlink returns LD for the overlapping variant | **151,831 (12.2%)** | **177,129 (14.7%)** |
| …same position, different alleles | 93,440 | 103,317 |
| …inside another variant's REF span (overlap only) | 58,391 | 73,812 |

- **Warnings:** LDlink warns only when the 1000G row's ID starts with `rs` and differs from the query
  (`LDcommon.py:546`). GRCh38 1000G rows have no IDs, so **all GRCh38 cases are silent**. In GRCh37, rows without an
  rs ID are also silent.
- **Order ambiguity:** 8,174 GRCh37 and 2,833 GRCh38 picks fall on positions where rows from different MAF tiers tie.
  Exact parity there needs the `vcf_order` column (§10; now added, §16).

**b157 drift vs. the dump** (same comparison, b157 VCFs):
- **New and retired rsIDs:**
  - **+120,579,738 new rsIDs** (+11%).
  - **1,424,372 dump rsIDs are gone** from b157 (56,273 merged, 1,368,099 removed with no merge record).
  - +624,266 new merges.
- **Shared rsIDs** (1,051,515,360) are very stable:
  - GRCh37 positions: 69,058 differ, 21,051 lost their placement, 4,272 gained one.
  - GRCh38: 17,913 differ, 3,522 lost, 5 gained.
  - type: 138 changed.
- **Implication:** moving LDlink to b157 would mostly *add* resolvable rsIDs. Existing rsID → position results would
  barely change.

**Decision input** (PLAN_REVIEW finding 4, "allele-aware as an opt-in mode"):
- The positional rule is right for about 87% of the rsIDs LDlink can resolve to a 1000G row.
- For 12–15% it silently, or with a GRCh37-only warning, reports LD for a *different* variant that occupies or
  spans the position.
- An allele-aware mode needs the NCBI alleles only at build time. Each genotype row can carry the list of
  allele-matching rsIDs, and a query then picks the row whose list contains the rsID (or reports "alleles not in
  1000G"). The dbSNP index files stay unchanged.

## 16. Allele-aware mode (opt-in) and exact VCF order in the genotype files

Building on §15: the genotype files now carry what both row-selection modes need. Parity stays the default.

**New columns** (`build_genotypes_v2.mjs`; `ALLELES=<NCBI slice>` enables the allele columns):

| column | type | meaning |
|---|---|---|
| `vcf_order` | INT32, delta-encoded | row index in the source VCF (multi-allelic records split in ALT order). Exact VCF order across MAF tiers; closes the §10 `vcf_order` item |
| `allele_is_label` | BOOLEAN | the row's allele matches are exactly `[rsid]`, the common case (about 96% of rows) |
| `allele_rsids` | UTF8, optional | otherwise: the matching current rsIDs, comma list, ascending; null means no match |

- **New key-value metadata:**
  - `max_ref_len`: rows overlapping position p start in `[p − max_ref_len + 1, p]`.
  - `vcf_order_rule` and `allele_rule`.
- **How a match is defined:** an NCBI b155 ALT whose minimal form (trim the common suffix, then the common prefix)
  equals the row's.
- **Inputs:** one chromosome of the b155 VCF Parquet from §15 (`out/ncbi_dbsnp/b155/chr22_37.parquet` etc.). b155 is
  used so that rsIDs agree with the dump that resolves them.
- **Build cost:** about 3 min per chr22 set, 7 min for chrX.
- **Every pre-existing column is unchanged.** All 12 files are identical row for row to v2 in `pos`, `rsid`, `haps`,
  `missing` and the rest (DuckDB `EXCEPT ALL`, 0 rows), so the earlier LD checks still hold.

| set | rows | no allele match | 1 match | 2+ | allele ≠ rsID label | size (v2 → v3) |
|---|---:|---:|---:|---:|---:|---|
| GRCh37 chr22 | 1,109,368 | 1,087 | 1,107,259 | 1,022 | 10,733 | 46.6 → 47.4 MB |
| GRCh38 chr22 | 1,059,079 | 18,454 | 1,040,506 | 119 | 9,278 | 44.7 → 45.7 MB |
| High Coverage chr22 | 1,070,401 | 26,063 | 1,044,176 | 162 | 54,482 | 66.4 → 67.7 MB |
| GRCh37 chrX | 3,501,213 | 1,069 | 3,498,549 | 1,595 | 38,853 | 112.0 → 115.7 MB |

**Row selection** (`select_rows.mjs`; planner: `query({ …, select: "parity" | "allele" })` sets
`variants[i].selected`). Both modes read the overlap window `[pos − max_ref_len + 1, pos]` across all tiers with
page-level `rangeScan`, because LDlink's tabix query returns overlapping records, not only rows at exactly POS.
- **`parity`:** the overlapping row with the highest `vcf_order` (LDlink keeps the last tabix row). A row split from
  a multi-allelic record gives `not_biallelic`.
- **`allele`:** the row whose matches contain the query's **current** rsID (`ref_id` for merged inputs). Statuses:
  - `ok`
  - `ok_split_multiallelic`
  - `multiallelic_in_1000g`: the matches are the split ALTs of one record (LDlink: "not biallelic").
  - `ambiguous`: several distinct rows; all are returned.
  - `alleles_not_in_1000g`: rows overlap the position, but none carry the rsID's alleles. The parity row is returned
    alongside for comparison.
  - `not_in_1000g`.

**Verification** (`verify_select.mjs`). 50,000 rsIDs per set, all of them rsIDs whose dump position overlaps a 1000G
row. Ground truth is independent of the new columns:
- **Parity:** one streaming pass over the source VCF, taking the last record that overlaps each position.
- **Allele:** allele matches recomputed from the NCBI slice.

| set | parity agree | allele agree | ok | alleles not in 1000G | multi-allelic | ambiguous |
|---|---:|---:|---:|---:|---:|---:|
| GRCh37 chr22 | 50,000 / 50,000 | 50,000 / 50,000 | 43,742 | **6,065 (12.1%)** | 190 | 3 |
| GRCh38 chr22 | 50,000 / 50,000 | 50,000 / 50,000 | 42,704 | **7,296 (14.6%)** | 0 | 0 |
| High Coverage chr22 | 50,000 / 50,000 | 50,000 / 50,000 | 34,218 | **14,811 (29.6%)** | 0 | 971 |
| GRCh37 chrX | 50,000 / 50,000 | 50,000 / 50,000 | 43,739 | **6,000 (12.0%)** | 261 (incl. ambiguous) | — |

(The chrX multi-allelic and ambiguous counts come from the run before the `multiallelic_in_1000g` split, so they are
combined there.)

- **Known parity gap:** the build skips symbolic SV records (`<…>` ALTs). In 1 of 50,000 GRCh37 positions (chr22 and
  chrX each), LDlink's last overlapping record is such an SV.
- **Smoke test (planner, GRCh37):** for rs370150588 (`16062095:GT>G`), parity uses the T>G SNV rs587600953 at
  16,062,096, while allele mode returns the GT>G deletion. For rs1430893789, allele mode reports
  `alleles_not_in_1000g` where LDlink silently uses rs62224611's row.
- **Planner check** (`verify_planner.mjs`) on the new files: hybrid rsID resolution is identical to direct (0
  differences), and the cheaper strategy is still chosen in every measured case.
- **Published:** the 12 genotype files and the manifest are in `s3://eaap-ac-test-bucket/ldlink-parquet/v2/` (release
  now 24.45 GB).

**Update 2026-11-21 (ROADMAP Phase A):** the genotype files were rebuilt with the fixed labeller (`rsid_label.mjs`).
- **How labels are found now:** the labeller searches the dump where LDlink stores each variant:
  - POS+1 for deletions/delins
  - POS or POS+1 for insertions
  - the minimal POS for SNVs/MNVs, including padded split rows
- **Source order:** NCBI allele match (source 4) ranks below the VCF rsID and above position + type.
- **New column `rsid_pos_delta`:** dump position − POS. Usually 0 or 1; up to 40 for padded split rows whose minimal
  variant sits further right.
- **Results (GRCh37 chr22):**
  - Deletions labelled with their own allele match: 4 → 27,638 of 27,831.
  - Insertions: 1,652 → 16,120 of 17,035.
  - Unlabelled rows: 33,588 → 941 (GRCh38: 63,491 → 14,801).
- **Supersedes:** the label-quality numbers in §15 (a) and the "allele ≠ rsID label" column above, which describe
  the old labels.
- **Unchanged:** all other columns (12/12 files identical), parity and allele selection (4 × 50,000), and hybrid ==
  direct.
