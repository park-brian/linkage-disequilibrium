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
| G1 | `npm test` | 25 / 25 | 2026-11-23 |
| G2 | `npm run verify` | planner, parity (12 checks), annotations: 0 differences · kernel vs LDlink math, chr1 and chrX (masked): 0 failures · row selection 22,000 / 22,000 (parity and allele) · ldsc ≤ 3e-14 relative | 2026-11-23 |
| G3 | `npm run verify:ldlink` | LDmatrix 30/30 · LDpair 419/429 + 10 tie order (LDlink randomizes ties) · LDpop 537/537 · LDhap 70/70 · LDproxy 57/57 · SNPclip 49/49 · SNPchip 66/66 · LDtrait 55/55 · LDexpress 58/58 · LDassoc 21/21 | 2026-11-23 (before slice 1; the moves don't touch module code) |
| G4 | `/app/check.html` + the app | 18/18 identical in the browser; Pyramid chr22 / chrX and LDassoc render | 2026-11-23 |

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
| Pyramid chrX 2 Mb EUR | 0.2 s on 6 workers |

## In progress: the restructure ([RESTRUCTURE.md](RESTRUCTURE.md) §11)

| slice | state |
|---|---|
| 0: baseline commit + fixes | ✅ |
| 1: layout, data root, config | ✅ |
| 2: docs | ✅ this file |
| 3: one LDlink comparer runner | next |
| 4: narrow dedup in `src/ldlink`; the app module table into its own forms module | |
| 5: fewer reads, numeric genome keys | |
| 6: annotation fast path, masked-kernel same-mask path | |

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
