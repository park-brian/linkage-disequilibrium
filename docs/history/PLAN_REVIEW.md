# Plan review (2026-11-19)

This review covers [REQUIREMENTS.md](REQUIREMENTS.md), [PARQUET_LAYOUT.md](PARQUET_LAYOUT.md) and
[DATA.md](DATA.md), checked against the goals:
- one LD module
- any number of SNPs
- windows of 1–10 Mb (up to 100 Mb)
- data as small and as fast as possible
- two dbSNP indexes with a hybrid lookup

Findings are ordered by impact. Items marked **measured** come from checks run during this review.

## Findings

### 1. The MAF tiers must use sub-population MAF, not super-population MAF (**measured**)
- **Check:** `scripts/profile/tier_check.py` over the first 150k GRCh37 chr22 records.
- **Result at 5%:** of the 29,764 variants with MAF ≥ 5% in at least one of the 26 populations, 5,977 (**20%**) are
  below 5% in every super-population.
- **Result at 1%:** 26,022 of 67,302 (**39%**) would be missed.
- **Consequence:** a "common only" query for CEU, FIN, PEL and so on would silently drop population-specific
  common variants.
- **Fix:** define `maf_max` over the 26 sub-populations. The common tier grows by about 25% in rows (23.8k → 29.8k
  of 150k), so expect about +20% on common-tier bytes.
- **How it's used:** the tier is only a storage partition. The engine still filters by the selected population's
  actual MAF after loading.

### 2. Scattered, genome-wide SNP lists are not covered yet (✅ addressed by the step 5 planner)
- **The gap:** the hybrid assumes the inputs cluster into a few regions. A list like 3000 GWAS hits across the
  genome does not.
- **Cost of each isolated SNP:**
  - 1 rsID-index page: about 17 KB
  - 1 genotype page: a 256 KB uncompressed common-tier page holds about 418 rows, roughly 65 KB compressed
- **At 3000 SNPs** that is about 240 MB, which is not acceptable.
- **Fix:** a cost-based planner for each chromosome cluster. It picks the cheaper of:
  - scanning the region: bytes known from the manifest
  - reading the distinct pages: page sizes known from the offset index
- **Page size:** the common tier's page size is a real trade-off. **Measured in step 1:** 64 KB pages are +14% on
  10 Mb scans and make scattered point reads 2–2.5× cheaper (36 vs 79 KB/point) than 256 KB pages. See
  PARQUET_LAYOUT.md §6.
- **Pairs:** LD is computed only within a chromosome, and optionally within a distance band. So for scattered lists
  the pair count is small, even when fetching is not.

### 3. rsID-keyed annotations must be re-keyed by position at build time
- **The problem:** FORGEdb (38M rows) and dbSNP `function` are keyed by rsID. Drawing them over a 10 Mb window
  (300k variants) would mean 300k rsID lookups.
- **Fix:** join them into the genotype rows, or into position-sorted annotation files, when building.
- **General rule:** every annotation must be reachable by `(build, chrom, pos)` range scans.

### 4. The rsID ↔ genotype join: keep the dump and use position + type (**measured**; decision: keep the dump)
- **Ambiguity:** 5.3% of chr22 dbSNP positions have more than one rsID; the maximum is 186 at one position. The dump
  has no alleles.
- **Why we keep it anyway:** it is LDlink's own data. It is needed to **reproduce LDlink's exact rsID → position
  resolution, and therefore its exact LD values,** for parity testing.
- **Measured on GRCh37 chr22** (1,109,368 rows; 527,823 have a single VCF rsID to use as ground truth):
  - **Position alone:** 5.6% of rows are ambiguous among current rsIDs. 1.7% have no dbSNP record at that position.
  - **Position + compatible `type`** (snv ↔ snv; indel classes match each other): only 0.33% of rows stay ambiguous.
  - **Accuracy:** with ties broken by lowest rsID, the rule picks the VCF's rsID for **500,113 of 500,435 (99.94%)**
    rows whose rsID is present at that position.
- **Rule for building rsID labels:**
  1. Use the VCF's rsID if it exists in dbSNP, mapping merged IDs to `ref_rsid`.
  2. Otherwise use position + type among current rsIDs, with ties broken by lowest rsID. Store `rsid_candidates`
     (the count) so ambiguity is visible.
  3. Otherwise leave it null.
  - GRCh38 and High Coverage have no VCF rsIDs, so they use step 2 with `key38`.
- **Parity mode (default when comparing with LDlink):**
  - Resolve input rsIDs → position from the dump, exactly as LDlink does (`id`, preferring `id == ref_id`).
  - Pick the genotype row by position: the row whose POS matches, otherwise the last row in VCF order.
  - Report sites that were multi-allelic in the source as "not biallelic".
  - **Build requirements:**
    - Preserve VCF order within a position.
    - Keep a `multiallelic_source` flag (or `n_alt`) on split rows.
    - Keep the **full** rsID index (all 1.07B rows, ~10.7 GB). This is what lets parity mode tell "not in dbSNP"
      apart from "not in 1000G", as LDlink does.
- **Later:** rebuild dbSNP from a current NCBI release (with alleles) and **compare it with the dump**. Measure:
  - rsID → position changes
  - how often the position + type rule disagrees with an allele-aware join
  - how many LD results change

  Then decide whether allele-aware matching becomes an opt-in mode.

### 5. Haploid chrX should not use a per-row missing bitset (✅ fixed in v2; verified in step 3)
- **The problem:** males are haploid on non-PAR chrX, so the current design writes a 626-byte `missing` bitset on
  nearly every chrX row.
- **Fix:** store a per-chromosome ploidy mask in the file metadata (male sample indexes, PAR boundaries). Keep
  per-row `missing` only for real `.` calls.
- **Status:** not measured yet. chrX is not in the benchmark.

### 6. The size estimates are extrapolations from one chromosome of one build
- **Genotypes:** "about 3 GB per build" is chr22 × 77, for GRCh37 only, and before fix 1 (+~20% common tier) and the
  INT64 `rsid` column.
  - The `rsid` column could add about 4 B/row. On the 9 B/row rare tier, that is close to +50%.
  - High Coverage has 801-byte bitsets, more indels and 6% duplicate positions, so expect it to be larger.
- **rsID index:** parity needs the full index, about 10.7 GB from 10 B/row × 1.07B. A 1000G-limited copy would be
  about 1–1.5 GB, not "under 1 GB", but it can't tell "not in dbSNP" apart from "not in 1000G".
- **To do:** re-measure after the fixes.

### 7. Round trips: add a per-build manifest (✅ done in steps 2 and 5)
- **Today's cost:** every file costs a footer read, then a page-index read, then data. hyparquet's default first
  footer read is 512 KB.
- **With 3 tiers × several chromosomes,** those metadata phases dominate small queries.
- **Fix:** a manifest file (one per build) holding, for each file:
  - footer offset and length
  - each row group's position range
  - page-index offsets
- **Effect:** the client can skip the footer phase or fetch exact byte ranges, and the planner (finding 2) gets its
  cost inputs for free.

### 8. The LD engine needs an output contract and algorithmic pruning, not only faster popcount
- **Speed estimate (measured in step 4):** blocked scalar WASM is 5–35× faster than the JS code, and SIMD is not
  faster than scalar `i64.popcnt` in V8. See PARQUET_LAYOUT.md §9.
- **r² bound pruning:** the r² bound from allele frequencies, `r² ≤ (p_min·(1−p_max)) / (p_max·(1−p_min))`, lets
  thresholded queries (proxies, pruning, LDtrait/LDexpress-style joins) skip most pairs without counting them.
- **Tile computation:** a tile is a bit-matrix product (Gram matrix of haplotype vectors), which suits WebGPU and
  blocked WASM.
- **"Any number of SNPs" needs a defined output:**
  - which statistics are kept per pair
  - the tile schema and zoom levels for the pyramid (aggregate, e.g. max r² per pixel bin)
  - the export format for thresholded pairs, e.g. a Parquet/TSV stream instead of an N×N text matrix

### 9. List the deliberate differences from LDlink
These change results compared with LDlink, and parity tests need to know about them:
- multi-allelic sites are split (LDlink rejects them)
- haplotypes are counted by allele index (LDlink uses string keys)
- duplicate positions are resolved by alleles (LDlink keeps the last row)
- High Coverage uses only the 2504 panel samples (same as LDlink, but worth saying)
- MAF tiers can hide rare variants unless the user requests them

### 10. Data versions are decisions, not defaults
- **dbSNP b155** (2023 dump): **kept** as the parity source; see finding 4.
- **GWAS Catalog:** 2025 backup or fresh pull?
- **GTEx v8:** newer releases exist.
- **GRCh38 chrX (2019 release):** only 107k records and no haploid calls, so it is suspect.

## What is solid
- **Bitset genotypes + ZSTD + page index + coalesced range reads:** 10 Mb of common variants is 7.6 MB in 15
  requests; all of chr22 is 39 MB against a 208 MB VCF.
- **Two dbSNP indexes with the composite key `chrom × 1e9 + pos`.** The measurement that rsIDs don't cluster
  justifies the hybrid. Small lists use direct rsID-page lookups; everything else resolves through region scans of
  the genotype files.
- **Cost driven by window size, not list size, for clustered inputs.** That matches "any number of SNPs".
- **The reference LD math is pinned down** with test vectors in REQUIREMENTS.md §5.

## Revised next steps
1. ✅ **Done 2026-11-19** (PARQUET_LAYOUT.md §6): chr22 v2 is 46.6 MB, using 64/64/16 KB pages. The rsID labels
   are 99.94% consistent with the VCF, and point lookups go through the payload column's offset index.
   **Rebuild chr22 with the corrected tiers.** Use sub-population `maf_max`, the INT64 `rsid` column, and a ploidy
   mask. Measure page sizes of 16/64/256 KB per tier, including compressed bytes per page for point reads.
2. ✅ **Done 2026-11-19** (PARQUET_LAYOUT.md §7). The indexes were built from the dump and the manifest is written.
   The parity resolver has 0 mismatches over 13,400 sampled cases and the test vectors.
   **Build the dbSNP indexes from the 2023 dump.** That means the full rsID index (partitioned), the merged-ID map,
   the position + type rsID labels in the genotype rows, and the per-build manifest. Then add a parity-mode
   resolver and verify it against the reference test vectors.
3. ✅ **Done 2026-11-19** (PARQUET_LAYOUT.md §8): GRCh38 chr22 is 44.7 MB, High Coverage chr22 66.4 MB and GRCh37
   chrX 112 MB. chrX ploidy is exact, and Parquet LD matches the reference on 4,800 pairs.
   **Extend the builds and benchmarks** to GRCh38 and High Coverage chr22, plus GRCh37 chrX.
4. ✅ **Done 2026-11-19** (PARQUET_LAYOUT.md §9): the WASM kernel is exact against the reference and pruning is
   lossless. The scalar blocked kernel runs at 34–224 Mpairs/s on one thread, 5–7.5× more with workers.
   **LD kernel:** WASM SIMD, blocked tiles, and r² bound pruning. Benchmark CEU/EUR/ALL on 1 and 10 Mb.
5. ✅ **Done 2026-11-19** (PARQUET_LAYOUT.md §10): hybrid rsID resolution equals direct (0 differences, 3–11× fewer
   bytes for clustered lists). Scan-vs-points choices were cheaper in every measured case.
   **Query planner + manifest.** Cover clustered lists, scattered genome-wide lists and windows, and check the cost
   model against measured bytes.
6. ✅ **Done 2026-11-20** (PARQUET_LAYOUT.md §11): a private test bucket with presigned Range requests. Windows take
   under 1 s cold even at 6 connections. List queries need HTTP/2 (CloudFront) and fewer requests per lookup.
   **Network test from S3/CloudFront.** This needs a place to publish test files, which is your call.
7. ✅ **Done 2026-11-20** (PARQUET_LAYOUT.md §12): 20 annotation files (3.2 GB), 0 parity mismatches over 6,950
   cases, and 1–10 Mb windows at 0.1–5 MB per track.
   7. **Annotation files keyed by position,** joined at build where rsID-keyed: FORGEdb, function, GWAS, GTEx,
   RegulomeDB, recombination, genes, arrays.
8. ✅ **Done 2026-11-20** (PARQUET_LAYOUT.md §15): NCBI dbSNP rebuild compared with the 2023 dump.
   - **chr22, exact:** LDlink's own parser on b155 JSON reproduces all 14,638,849 docs. Only the dump's Mongo order
     can't be rebuilt, so the dump is still needed.
   - **Genome-wide (b155 VCF):** all 1.05 B rsIDs present; type identical; positions match LDlink's HGVS rule for
     99.995% (GRCh37) and 99.998% (GRCh38).
   - **LDlink parser bug:** at least 4,830 rsIDs have impossible positions (HGVS repeat notation was concatenated).
   - **chr22 allele-aware join:** LDlink returns LD for a 1000G variant whose alleles don't match the rsID in 12% of
     rsID → 1000G resolutions (GRCh37) and 15% (GRCh38, never warned).
   - **b157 drift:** +120.6 M new rsIDs and 1.42 M removed; shared rsIDs barely move.
   - ✅ **Allele-aware opt-in mode built** (PARQUET_LAYOUT.md §16):
     - The genotype files gain `vcf_order`, `allele_is_label` and `allele_rsids` (+1–3% size).
     - `select_rows.mjs` and the planner's `select` option choose rows in `parity` or `allele` mode.
     - Both modes agree with independent ground truth on 4 × 50,000 rsIDs.
9. ✅ **Done 2026-11-20** (PARQUET_LAYOUT.md §14): LDscore family in JS.
   - 14 reference-panel Parquet files (1.30 GB).
   - `--l2` byte-identical to LDlink's `22.*` example.
   - `munge_sumstats` output identical to the reference.
   - h² and rg match ldsc to about 1e-14, end to end from the raw BBJ files.
   - The ldsc rg quirk (no two-step) is reproduced on purpose.
