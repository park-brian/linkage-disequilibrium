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
src/        the runtime core: browser + Node, no build step   → docs/ARCHITECTURE.md
  ld/         LD kernel: kernel.c → kernel.wasm, kernel_js.js, one API in kernel.js
  ldlink/     the thirteen LDlink module ports + shared context → docs/PARITY.md
  ldsc/       LD score regression, munging, LD scores
app/        the browser UI (forms, render, region view, LDmatrix plot, pyramid + tracks, worker pool, parity check page)
pipeline/   offline release build: dbsnp, ncbi, genotypes, annotations, ldscore, release
verify/     data verifiers (data/), ldsc checks (ldsc/), LDlink comparisons (ldlink/), run.sh
test/       unit tests and fixtures
bench/      benchmarks and recorded results
tools/      config.{mjs,sh} (every machine path), serve.js, run_sql.py, move.py, ldlink-local/ (oracle kit)
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
| `npm run verify:ldlink` | the thirteen module comparisons against a local LDlink (about 1.5 h; `tools/ldlink-local/README.md`) |
| `npm run bench` | WASM vs JS kernel benchmark (`bench/` has the others) |
| `npm run build:kernel` | rebuild `src/ld/kernel.wasm` (`CLANG=<emsdk>/upstream/bin/clang`) |

The pipeline stages are run per directory (`pipeline/<stage>/`). Each script's header has its usage, and every
script runs with its working directory at `$LD_DATA`.

## Requirements

- Node 22+;
- for the pipeline: DuckDB and Python 3;
- for the oracle: Docker.

The browser app needs only a static file host with Range support.
