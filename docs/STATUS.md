# Status

**Updated:** 2026-11-23 · **branch:** `restructure` (the live site's `origin/main` plus the cleanup's commits, not pushed) · **release:** `release_v4` · **live:** https://park-brian.github.io/linkage-disequilibrium/ (data: https://d2x2yo644ey89o.cloudfront.net)

## What works

- **Every LDlink module runs on the core and reproduces LDlink.** LDmatrix, LDpair, LDpop, LDhap, LDproxy, SNPclip,
  SNPchip, LDtrait, LDexpress and LDassoc on GRCh37, GRCh38 and GRCh38 High Coverage, and LDherit, LDcorrelation and
  LDscore (ldsc 3.0.1), byte for byte against a real local LDlink. Bugs are reproduced too ([PARITY.md](PARITY.md)).
- **The browser app** (`npm run serve` → `/app/`; live on GitHub Pages, deployed by `.github/workflows/pages.yml` on
  every push to `main`, with the release URL from the repository variable `LD_DATA_URL`):
  - all thirteen modules, plus the region view;
  - LDlink's LDproxy and LDassoc figures; an LDmatrix plot with a gene track;
  - a region view: an LD pyramid up to 5 Mb under GWAS, eQTL, FORGEdb / RegulomeDB, recombination and gene tracks,
    with LD on a pool of up to 8 Web Workers;
  - links: every form reads and writes its inputs in the URL; a pyramid pair opens LDpair / LDproxy, a locus point
    LDpair.
- **The release:** 326 Parquet files + manifest + `examples/`, about 34.8 GB, in `$LD_RELEASE` (self-contained: panels
  and example inputs inside), and **public at `https://d2x2yo644ey89o.cloudfront.net`** (CloudFront over the private
  bucket's `ldlink-parquet/v5/`; DATA.md). The built site runs every view against it.

## Verified

| gate | command | result | last run |
|---|---|---|---|
| G1 | `npm test` | 37 / 37 (docs and usage lines, line length, layers, K5, D1 panels, F1 numpy formatting, the core LD API against the ports) | 2026-11-23 |
| G2 | `npm run verify` | planner, parity (12 checks), annotations (incl. 34,805 range-track rows): 0 differences · kernel vs LDlink math, chr1 and chrX (masked): 0 failures · row selection 22,000 / 22,000 (parity and allele) · ldsc ≤ 3e-14 relative | 2026-11-23 (final code) |
| G3 | `npm run verify:ldlink` | LDmatrix 30/30 · LDpair 419/429 + 10 tie order (LDlink randomizes ties) · LDpop 537/537 · LDhap 70/70 · LDproxy 57/57 · SNPclip 49/49 · SNPchip 66/66 · LDtrait 55/55 · LDexpress 58/58 · LDassoc 21/21 · LDherit 29/29 · LDcorrelation 19/20 + 1 degenerate_z (rg = 1, se ≈ 1e-7) · LDscore 19/19 | 2026-11-23 (cleanup slice 10, `25c749b`; slice 9 after it only inserts line breaks) |
| G3s | `npm run verify:snapshot` | 1,440 / 1,440 hashes after every slice of the cleanup; over HTTP (`--http`, dev server): LDmatrix, LDassoc, LDscore, LDherit 99 / 99 | 2026-11-23 (final code) |
| G4 | `npm run verify:app` | 14 / 14 views in headless Chromium, no page or console errors, every figure's pixels unchanged through the cleanup, a run abandoned by a module switch shows nothing, back returns to the inputs; coverage `app/` 94.9%, `src/` 87% (snapshot: 92%); the built site ran LDproxy on the local release | 2026-11-23 (final code) |

## Speed (cold, browser, local dev server)

| query | time |
|---|---|
| LDmatrix (8 SNPs) | 0.25–0.4 s |
| LDpair | 0.1 s |
| LDproxy (±500 kb) | 1.3 s |
| LDtrait | 0.74 s |
| LDexpress | 0.5 s |
| LDherit / LDcorrelation (chr22 example, EUR) | 3–4 s |
| LDscore (chr22 example, 1 cM) | 3.5 s |

Over CloudFront (cold edge, from a home connection): LDmatrix 2.7 s, LDpair 1.1 s, LDproxy 14 s (505 requests: the
`sortedLookup` follow-up below), LDtrait 1.9 s, LDassoc 3.0 s, LDherit 4.1 s, Pyramid 0.9 s.
| Pyramid 1 Mb EUR, common | 0.22 s on 4 workers |
| Pyramid 3 Mb ALL, common + low (228 M pairs) | 1.0 s on 8 workers (7.0 s on 1) |
| Pyramid chrX 2 Mb EUR, common | 0.1 s on 8 workers (0.17–0.24 s before the same-mask path) |
| Pyramid chrX 3 Mb EUR, common + low (226 M pairs) | 0.64 s on 8 workers (1.45 s before) |

**Requests (cold, browser, local server; `bench/app_requests.js`):**

| | before the restructure | now |
|---|---|---|
| startup (manifest) | 2 requests, 7.15 MB | 1 request, 1.92 MB (0.61 MB gzipped: hosting should compress it) |
| LDmatrix (8 SNPs) | 97 requests, 2.70 MB | 62 requests, 2.75 MB |
| LDproxy rs3 ±500 kb | 886 requests, 12.6 MB | 505 requests, 13.4 MB |

The cuts: no HEAD per file (sizes from the manifest), one exact footer read per file (the resolver had its own 64 KB
read), reads issued together share Range requests (`coalesce`, 16 KB gap: −24% requests for +9% bytes on LDproxy),
and rsIDs with an unknown doc position try the VCF POS before the rsID index (LDproxy: 88 → 74 index lookups).

## Done: the cleanup ([history/CLEANUP_PLAN.md](history/CLEANUP_PLAN.md), outcome at its end; the review it came from: [history/CLEANUP_REVIEW.md](history/CLEANUP_REVIEW.md))

Eleven slices, committed on `restructure`, not pushed. What changed for a reader of the code:
- **Two layers:** `src/core` (reads return every match in file order; the LD kernel; `core/ld.js`, an LD API with
  numbers and typed errors) and `src/parity` (LDlink's choices, messages and formatting, the 13 ports). A test keeps
  core from importing parity.
- **Gates that fail**, a 13-minute parity snapshot (`verify:snapshot`) for most changes, the app driven in Chromium
  (`verify:app`), coverage (`npm run coverage`, `coverage:snapshot`).
- **One read path** for lookups and ranges (each column on its own pages): 512 rsIDs 13.8 → 8.4 MB.
- **The bugs of the review fixed:** kernel arena reuse, poisoned caches, missing files, stale results, pointer
  offsets, quadratic overlap scans, the server and the Pages workflow.
- **A pipeline runbook** (`pipeline/build.sh`) with every producer scripted; chr22 rebuilt byte-identical (GRCh37).
- **About 6,000 lines deleted**; code lines within 140 characters except 63 verbatim LDlink strings.

## Earlier: the restructure and Phase C

[history/RESTRUCTURE.md](history/RESTRUCTURE.md) (outcome in §12) and [history/PLAN.md](history/PLAN.md) (outcome in §9).

## Follow-ups

- **LDproxy over CloudFront** took 14 s (505 requests) before the cleanup's read path; not re-measured over the CDN.
- **Annotations** moved to the parity layer whole: the core has no annotation reads of its own yet (`core/ld.js`
  has no `annotate`).
- **LDscore family uploads** beyond the examples are untested against LDlink (its upload route); `sumstats_normalizer`
  and the upload compatibility checks are not ported.

## Open decisions (yours)

- **Push `restructure` to `main`** (it redeploys the site): the cleanup is local.
- **The published verifier reports** under the bucket's `ldlink-parquet/v5/` (slice 1 moved them out of the
  release locally): delete them from the bucket, or keep.
- **Scratch output of the slice 8 gate** in `$LD_DATA`: `scratch_gate_release/`, `phaseB/geno_gate/` (about 200 MB).

- **The test bucket's older prefixes** (`ldlink-parquet/v2/`–`v4/`, 94 GB: 24.5 + 34.8 + 34.8): keep or delete. `v4/` still has
  manifest version 1.
- **About 870 GB of regenerable intermediates** in `$LD_DATA` (`ncbi_dbsnp`, `dbsnp_dump`, `phaseB` build copies):
  keep or delete.
