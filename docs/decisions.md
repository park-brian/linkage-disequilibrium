# Decisions

Short records of choices that shape the code. The newest are at the bottom. Each says what was decided, why, and
what would change it.

## Browser-only, static Parquet, hyparquet (2026-11)

**What:** the app reads Parquet files with Range requests. There is no API server.

**Why:** LDlink's cost is a fleet of servers plus a Mongo cluster; static files on a CDN cost almost nothing to
host. hyparquet reads page indexes, so a query touches megabytes, not gigabytes.

**Revisit if:** a module needs data too large to stream per query.

## Parity first, opt-in improvements (2026-11)

**What:** every port reproduces LDlink's output exactly, including LDlink's bugs ([PARITY.md](PARITY.md)). The
allele-aware row choice is a separate mode (`select: "allele"`).

**Why:** a byte-for-byte oracle makes every change checkable. A new implementation of an existing tool is only
trusted when it can show it agrees.

**Revisit if:** we stop needing LDlink as the reference. Then the bug-compatibility rules become fixes.

## One LD kernel API, WASM + plain-JS backends (2026-11-22)

**What:** `src/core/ld/kernel.js` exposes one API; `kernel.wasm` (C, scalar popcount, blocked) and `kernel_js.js` give
identical results.

**Why:** WASM is 3.5–5.3× faster on large unthresholded bands and 1.7–3.1× with pruning; at module sizes both take
milliseconds. JS is the test oracle and the fallback. SIMD was dropped: it was never faster in V8.

**Revisit if:** WebGPU becomes worth it for genome-wide views.

## The rsID label rule (2026-11)

**What:** each 1000G row carries the rsID LDlink would show:
1. the VCF rsID (merged → current);
2. else an NCBI allele match;
3. else position + variant type.

The label is searched where LDlink's dump stores the variant (POS + 1 for deletions).

**Why:** this reproduces LDlink's labels without its Mongo lookups at query time.

## `vcf_id` column (2026-11-22)

**What:** GRCh37 rows store the raw VCF ID when it contains "rs" and is not the label.

**Why:** LDlink prints it in warnings and in LDproxy's RS_Number; 0.65% of rows differed without it. The cost is
about 0.15% in size.

## Positional annotation fast path (2026-11-23)

**What:** dbSNP function and FORGEdb score are looked up at the rsID's known position (POS + `rsid_pos_delta`). An
rsID with several dump docs (`duplicated_rsids`) falls back to the exact `find_one` path.

**Why:** LDproxy went from 7.4 s / 11,912 reads to 1.2 s / 874 reads, with identical output.

## Data outside the repo (2026-11-23)

**What:** generated data lives under `$LD_DATA` (default `../data`). Every machine path comes from
`tools/config.js` and `config.sh`, and scripts run with the working directory at `$LD_DATA`.

**Why:** 993 GB of outputs inside the tree made search and backup dangerous; there were about 25 hardcoded paths.

## Our page-index lookups, not hyparquet's (2026-11-23)

**What:** `src/core/parquet.js` `sortedLookup` stays. hyparquet 1.33's own pruning (`parquetScan` with a `$in`
`pruningFilter` and `usePageIndex: true`) returned identical rows in every case of bench/pruning.js (since deleted; in
git history), but it did not
win on every measure (cold reads, local release; `$LD_DATA/bench/results/pruning.json`):

| case | ours: requests / bytes / ms | hyparquet: requests / bytes / ms |
|---|---|---|
| rsID index, 32 keys | 264 / 3.0 MB / 218 | 193 / 1.8 MB / 357 |
| rsID index, 512 keys | 3,620 / 38.1 MB / 2,009 | 3,109 / 8.6 MB / 4,809 |
| dbSNP by position, 200 keys | 1,154 / 4.8 MB / 2,203 | 762 / 6.8 MB / 7,331 |
| genotypes, 50 points (`haps` payload) | 243 / 2.9 MB / 216 | 233 / 2.4 MB / 264 |

**Why:** the plan's rule was to adopt only if neither requests nor bytes regress; bytes regress on the position file
and CPU is 1.2–3.3× higher.

**What it showed about ours:** without a payload column, `sortedLookup` reads every requested column over the key
pages' whole row spans; hyparquet reads each column only on its own pages, which is why 512 rsIDs cost 8.6 MB there
and 38 MB here. Mapping every column to its own pages, as the `payloadColumn` path already does for one column, is
the follow-up.

**Revisit if:** hyparquet batches row reads per page (its per-row `readColumn` decodes a page once per row), or the
follow-up above lands and a re-run favours either side.

**Follow-up (2026-11-23, cleanup slice 6):** done. `sortedLookup` and `rangeScan` share one read path
(`src/core/parquet.js` `matchRows`): the key column alone over its candidate pages gives the exact matched rows, and every
other column is read on its own pages (its offset index) holding those rows. A range keeps one read of all columns
over the key pages, because a range's rows are contiguous and every column's pages there are needed anyway
(per-column reads measured 3 → 11 reads for a 1 Mb genotype window). Measured as the app reads, one cached
AsyncBuffer per file, before → after (reads / MB):

| case | before | after |
|---|---|---|
| rsID index, 32 keys | 240 / 2.40 | 300 / 1.82 |
| rsID index, 512 keys | 2,634 / 13.8 | 2,644 / 8.38 (hyparquet, uncached: 3,109 / 8.6) |
| dbSNP by position, 200 keys | 542 / 2.27 | 473 / 2.17 |
| genotypes, 50 points | 149 / 1.77 | 137 / 1.74 |
| genotypes 1 Mb window (rangeScan) | 3 / 1.16 | 3 / 1.16 |
| recombination 3 Mb (rangeScan) | 5 / 0.04 | 4 / 0.08 (the column and offset indexes in one span) |

Single-column reads take hyparquet's decoded arrays (`onChunk`), not row objects: building objects for whole page
spans was 3x slower than before and ran the parity check out of memory. With arrays, CPU falls too (local, 200
position keys: 1.18 → 0.62 s; 512 rsIDs: 1.49 → 1.17 s). The `payloadColumn` option is gone (every column now gets
that treatment).
