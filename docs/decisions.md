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

**What:** `src/ld/kernel.mjs` exposes one API; `kernel.wasm` (C, scalar popcount, blocked) and `kernel_js.mjs` give
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
`tools/config.{mjs,sh}`, and scripts run with the working directory at `$LD_DATA`.

**Why:** 993 GB of outputs inside the tree made search and backup dangerous; there were about 25 hardcoded paths.
