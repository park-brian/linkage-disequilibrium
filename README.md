# linkage-disequilibrium

NCI's [LDlink](https://ldlink.nih.gov), rebuilt as a **browser-only** app.
- **Data:** 1000 Genomes haplotypes, LDlink's dbSNP dump and its annotation collections, as static Parquet files.
- **Reading:** the files are read with HTTP Range requests ([hyparquet](https://github.com/hyparam/hyparquet)), with
  LD computed in the page by a WASM kernel (plain-JS fallback).
- **Parity:** all thirteen LDlink modules (the ten LD tools plus LDherit, LDcorrelation and LDscore) reproduce LDlink's
  output byte for byte, checked against a real local LDlink.

**Where things stand:** [docs/STATUS.md](docs/STATUS.md).

## Layout

```
src/        the runtime: browser + Node, no build step   → docs/ARCHITECTURE.md
  core/       reads that return every match in file order (source, parquet, dbsnp, genotypes, genome, collect), haps,
              the LD kernel (ld/: kernel.c → kernel.wasm, kernel_js.js, one API in kernel.js), numerics (stats)
  parity/     everything LDlink-shaped: resolution and row choice (resolve, select, input), the planner, annotation
              queries, ldsc (ldsc/), Python number formatting (pyfmt), the thirteen ports + context (ldlink/)
              → docs/PARITY.md
app/        the browser UI (forms, render, region view, LDmatrix plot, pyramid + tracks, worker pool); app/set.html +
            app/set/ is a design prototype: one variant set, one context bar, nine lenses on the same ports
pipeline/   offline release build: dbsnp, ncbi, genotypes, annotations, ldscore, release
verify/     data verifiers (data/), ldsc checks (ldsc/), LDlink comparisons (ldlink/), run.sh
test/       unit tests and fixtures
bench/      benchmarks and recorded results
tools/      config.js / config.sh (every machine path), serve.js, build_site.js, run_sql.py, check_docs.js, coverage.js, ldlink-local/ (oracle kit)
docs/       STATUS · ARCHITECTURE · DATA (release contract) · PARITY · decisions · REQUIREMENTS · history/
```

Generated data (releases, build intermediates, reports) lives **outside** the repo in `$LD_DATA` (default `../data`).
See `tools/config.js` for every path and its environment variable.

## Commands

| | |
|---|---|
| `npm run serve` | dev server with Range support → `http://localhost:8095/app/` (the repo and the release, from config) |
| `npm run build:site -- --data <url>` | the static site for GitHub Pages in `dist/` (the release is read from `<url>`; `.github/workflows/pages.yml` runs it) |
| `npm test` | unit tests (seconds; needs the release unless `LD_SKIP_DATA=1`) |
| `npm run verify` | data verifiers: planner, parity, annotations, kernel (chr1, chrX), selection, ldsc (about 20 min) |
| `npm run verify:snapshot` | the parity snapshot: every LDlink comparison case through our ports, hashed against `verify/ldlink/snapshot.json` (about 13 min, no Docker; `-- --http <url>` reads the release over HTTP; CI runs it against CloudFront) |
| `npm run verify:app` | every view in headless Chromium against the dev server (about 1 min): fails on page or console errors, error results, missing figures, a stale run shown after a switch; prints `app/` and `src/` coverage (Playwright; once: `npx playwright install chromium`) |
| `npm run verify:set` | the set page prototype (`app/set.html`) in headless Chromium: every lens, a region, one variant, mixed inputs (about 30 s) |
| `npm run coverage` | line coverage of `src/` from the unit tests (`coverage:snapshot`: from the parity snapshot, about 92%; `tools/coverage.js`) |
| `npm run verify:ldlink` | the thirteen module comparisons against a local LDlink (about 1.5 h; `tools/ldlink-local/README.md`) |
| `npm run bench` | WASM vs JS kernel benchmark (`bench/` has the others) |
| `npm run build:kernel` | rebuild `src/core/ld/kernel.wasm` (`CLANG=<emsdk>/upstream/bin/clang`) |

The pipeline stages are run per directory (`pipeline/<stage>/`). Each script's header has its usage, and every
script runs with its working directory at `$LD_DATA`.

## Requirements

- Node 22+;
- for the pipeline: DuckDB and Python 3;
- for the oracle: Docker.

The browser app needs only a static file host with Range support.
