# Status

**Updated:** 2026-11-23 · **branch:** `restructure` (from tag `parity-complete`) · **release:** `release_v4`

## What works

- **Every LDlink module runs on the core and reproduces LDlink.** LDmatrix, LDpair, LDpop, LDhap, LDproxy, SNPclip,
  SNPchip, LDtrait, LDexpress and LDassoc on GRCh37, GRCh38 and GRCh38 High Coverage, and LDherit, LDcorrelation and
  LDscore (ldsc 3.0.1), byte for byte against a real local LDlink. Bugs are reproduced too ([PARITY.md](PARITY.md)).
- **The browser app** (`npm run serve` → `/app/`; `npm run build:site` for GitHub Pages):
  - all thirteen modules, plus the region view;
  - LDlink's LDproxy and LDassoc figures; an LDmatrix plot with a gene track;
  - a region view: an LD pyramid up to 5 Mb under GWAS, eQTL, FORGEdb / RegulomeDB, recombination and gene tracks,
    with LD on a pool of up to 8 Web Workers;
  - links: every form reads and writes its inputs in the URL; a pyramid pair opens LDpair / LDproxy, a locus point
    LDpair.
- **The release:** 326 Parquet files + manifest + `examples/`, about 34.8 GB, in `$LD_RELEASE` (self-contained: panels
  and example inputs inside). The private test bucket's `v4/` predates the panels, examples, counts and manifest v2.

## Verified

| gate | command | result | last run |
|---|---|---|---|
| G1 | `npm test` | 29 / 29 (docs check, K5, D1 panels, F1 numpy formatting) | 2026-11-23 |
| G2 | `npm run verify` | planner, parity (12 checks), annotations (incl. 34,805 range-track rows): 0 differences · kernel vs LDlink math, chr1 and chrX (masked): 0 failures · row selection 22,000 / 22,000 (parity and allele) · ldsc ≤ 3e-14 relative | 2026-11-23 (final code) |
| G3 | `npm run verify:ldlink` | LDmatrix 30/30 · LDpair 416/429 + 13 tie order (LDlink randomizes ties) · LDpop 537/537 · LDhap 70/70 · LDproxy 57/57 · SNPclip 49/49 · SNPchip 66/66 · LDtrait 55/55 · LDexpress 58/58 · LDassoc 21/21 · LDherit 29/29 · LDcorrelation 19/20 + 1 degenerate_z (rg = 1, se ≈ 1e-7) · LDscore 19/19 | 2026-11-23 (final code) |
| G4 | `/app/check.html` + the app | 18/18 identical on the dev server and on the built site with the release on another origin; all 14 views run, locus plots on GRCh37 / GRCh38 / High Coverage, links open and run, region track hashes unchanged | 2026-11-23 (final code) |

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

## Done: the restructure ([history/RESTRUCTURE.md](history/RESTRUCTURE.md); outcome in §12)

Slices 0–6 are done, one commit each, on branch `restructure`; not pushed. Where the work differed from the plan,
§12 says so: `genotypes.calls` and `pyFloat` were not extracted, hyparquet's native page pruning was not evaluated,
and genome keys stay BigInt where they meet INT64 columns.

## Done: Phase C and hosting prep ([history/PLAN.md](history/PLAN.md); outcome in §9)

All slices are committed on `restructure`, not pushed. The static site builds (`dist/`, 0.7 MB) and runs against a
release on another origin; `.github/workflows/pages.yml` deploys it once Pages uses Actions.

## Follow-ups

- `sortedLookup` reads every requested column over whole key-page spans; mapping each column to its own pages (as
  the payload path does) cut 512 rsID lookups from 38 MB to 8.6 MB in hyparquet's version (`decisions.md`).
- LDscore family uploads beyond the examples are untested against LDlink (its upload route); `sumstats_normalizer` and
  the upload compatibility checks are not ported.

## Open decisions (yours)

- **Pushing** `restructure` to `main` (a fast-forward). Switch Pages to the Actions workflow first, or in the same
  step: the legacy build would serve `app/` with an import map into the untracked `node_modules/`.
- **A public data host** for the release (Range + CORS): a public S3 bucket behind CloudFront, a Cloudflare R2
  bucket, or a Hugging Face dataset; then set the repository variable `LD_DATA_URL` (history/PLAN.md §6).
- **Uploads:** the release to that host (copy `examples/` too: the manifest lists Parquet only), and the version 2
  manifest to the test bucket's `v4/`.
- **About 870 GB of regenerable intermediates** in `$LD_DATA` (`ncbi_dbsnp`, `dbsnp_dump`, `phaseB` build copies):
  keep or delete.
