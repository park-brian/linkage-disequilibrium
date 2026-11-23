# Status

**Updated:** 2026-11-23 · **branch:** `restructure` (from tag `parity-complete`) · **release:** `release_v4`

## What works

- **Every LDlink module runs on the core and reproduces LDlink.** LDmatrix, LDpair, LDpop, LDhap, LDproxy, SNPclip,
  SNPchip, LDtrait, LDexpress and LDassoc, on GRCh37, GRCh38 and GRCh38 High Coverage, byte for byte against a real
  local LDlink. Bugs are reproduced too ([PARITY.md](PARITY.md)).
- **The browser app** (`npm run serve` → `/app/`):
  - all ten modules;
  - an LDmatrix plot with a gene track;
  - a region view: an LD pyramid up to 5 Mb under GWAS, eQTL, recombination and gene tracks, with LD on a pool of
    up to 8 Web Workers.
- **The release:** 324 Parquet files + manifest, 34.8 GB, in `$LD_RELEASE` and the private test bucket (`v4/`).

## Verified

| gate | command | result | last run |
|---|---|---|---|
| G1 | `npm test` | 27 / 27 (with the docs check and K5) | 2026-11-23 |
| G2 | `npm run verify` | planner, parity (12 checks), annotations: 0 differences · kernel vs LDlink math, chr1 and chrX (masked): 0 failures · row selection 22,000 / 22,000 (parity and allele) · ldsc ≤ 3e-14 relative | 2026-11-23 (final code) |
| G3 | `npm run verify:ldlink` | LDmatrix 30/30 · LDpair 420/429 + 9 tie order (LDlink randomizes ties) · LDpop 537/537 · LDhap 70/70 · LDproxy 57/57 · SNPclip 49/49 · SNPchip 66/66 · LDtrait 55/55 · LDexpress 58/58 · LDassoc 21/21 | 2026-11-23 (final code) |
| G4 | `/app/check.html` + the app | 18/18 identical in the browser; all 11 views run and render (SNPchip with no array chosen shows LDlink's own error) | 2026-11-23 (final code) |

## Speed (cold, browser, local dev server)

| query | time |
|---|---|
| LDmatrix (8 SNPs) | 0.25–0.4 s |
| LDpair | 0.1 s |
| LDproxy (±500 kb) | 1.3 s |
| LDtrait | 0.74 s |
| LDexpress | 0.5 s |
| Pyramid 1 Mb EUR, common | 0.22 s on 4 workers |
| Pyramid 3 Mb ALL, common + low (228 M pairs) | 1.0 s on 8 workers (7.0 s on 1) |
| Pyramid chrX 2 Mb EUR, common | 0.1 s on 8 workers (0.17–0.24 s before the same-mask path) |
| Pyramid chrX 3 Mb EUR, common + low (226 M pairs) | 0.64 s on 8 workers (1.45 s before) |

**Requests (cold, browser, local server; `bench/app_requests.js`):**

| | before the restructure | now |
|---|---|---|
| startup (manifest) | 2 requests, 7.15 MB | 1 request, 7.15 MB (1.49 MB gzipped: hosting should compress it) |
| LDmatrix (8 SNPs) | 97 requests, 2.70 MB | 62 requests, 2.75 MB |
| LDproxy rs3 ±500 kb | 886 requests, 12.6 MB | 505 requests, 13.4 MB |

The cuts: no HEAD per file (sizes from the manifest), one exact footer read per file (the resolver had its own 64 KB
read), reads issued together share Range requests (`coalesce`, 16 KB gap: −24% requests for +9% bytes on LDproxy),
and rsIDs with an unknown doc position try the VCF POS before the rsID index (LDproxy: 88 → 74 index lookups).

## Done: the restructure ([history/RESTRUCTURE.md](history/RESTRUCTURE.md); outcome in §12)

Slices 0–6 are done, one commit each, on branch `restructure`; not pushed. Where the work differed from the plan,
§12 says so: `genotypes.calls` and `pyFloat` were not extracted, hyparquet's native page pruning was not evaluated,
and genome keys stay BigInt where they meet INT64 columns.

## After that (Phase C, the app)

- FORGEdb / RegulomeDB tracks;
- LDproxy and LDassoc locus plots;
- a link from a pyramid pair to LDproxy;
- the LDscore views;
- population panels inside the release;
- hosting (Phase D: CloudFront + HTTP/2, which needs your decision).

## Open decisions (yours)

- **Pushing** the `restructure` branch to `github.com/park-brian/linkage-disequilibrium`.
- **About 870 GB of regenerable intermediates** in `$LD_DATA` (`ncbi_dbsnp`, `dbsnp_dump`, `phaseB` build copies):
  keep or delete.
- **Hosting** (Phase D).
- **A slimmer manifest** (a release format change): 93% of its 7.15 MB is the genotype files' per-row-group
  `colBytes` for 18 columns, while the planner prices 5. Rewriting only `manifest.json` would cut the app's startup
  download several-fold, but the bucket's `v4/` would then differ from the local release.
