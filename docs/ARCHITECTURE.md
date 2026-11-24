# Architecture

**The whole system:** static Parquet files, read in the browser with HTTP Range requests, and LD computed on the
page.
- There is no server-side computation.
- The only server is a static file host.
- The offline pipeline produces the files.
- A local copy of the real LDlink is the oracle for parity.

```
            offline (pipeline/)                                    runtime (src/ + app/, browser or Node)
 1000G VCFs ─┐                                              ┌──────────────────────────────────────────────┐
 LDlink mongodump ─┼─► pipeline/<stage> ─► release/ ─HTTP──►│ source → parquet → dbsnp · genotypes ·        │
 NCBI dbSNP ─┤      (dbsnp, ncbi, genotypes,   (Parquet +    │           annotations → planner                │
 ldsc refs ──┘       annotations, ldscore,      manifest)    │ ld/kernel (WASM | JS, worker-safe)             │
                     release)                                │ ldlink/<module> ports · ldsc/ ports            │
                                                             │ app/ UI · app/workers (LD pyramid pool)        │
                                                             └──────────────────────────────────────────────┘
 local LDlink (tools/ldlink-local) ◄── verify/ldlink compares every module, byte for byte
```

## src/: the runtime (browser-safe, no build step)

Two layers (REVIEW.md §4): `core/` reads and computes without LDlink's choices; `parity/` holds everything
LDlink-shaped (the ports, their messages and Python's number formatting). The files still at the top of `src/` are
core reads with LDlink's choices inside them (first or last match, coordinate parsing); they split in the next step.
Dependencies flow one way, core never importing parity:

```
core/source ─► core/parquet ─► dbsnp, genotypes, annotations ─► planner ─► parity/ldlink/* (context, common, modules)
genome ─────────────────────────┘           core/haps ─► core/ld/kernel (no other imports: workers load it by path)
core/stats ─► ldsc/regress, ldsc/munge; ldsc/l2 ◄─ parity/pyfmt, parity/ldsc/{log,pandas}
```

| module | one line |
|---|---|
| `genome.js` | chromosome codes and names, genome keys, the `BUILDS` table, LDlink's `firstToken` / `parseCoord` |
| `core/source.js` | the one cached `{ file, metadata }` per path (sizes and footer lengths from the manifest: no size probe, one exact footer read); `isNotFound`, `kv` |
| `core/parquet.js` | page-index reads: `sortedLookup` (point keys), `rangeScan` (key range), `pool`; `coalesce` (reads issued together share Range requests) |
| `core/collect.js` | `groupBy`, `firstRows` (the first row per key: the row a Mongo `find_one` returns) |
| `dbsnp.js` | LDlink-parity resolution: `resolveRsids`, `resolveCoords`, `resolveInputs`; the `find_one` forms `findAt`, `findAllAt`, `findById`; `duplicated` |
| `core/haps.js` | `decodeHaps`: bitset or sparse → bitset (no dependencies, so the kernel stays worker-safe) |
| `genotypes.js` | parity / allele-aware row selection (`pickParity`, `pickAllele`, `candidateRows`); tier-row reads (`row`, `windowRows`, `samples`) |
| `annotations.js` | genes, recombination, RegulomeDB, FORGEdb (by rsID and by position), arrays, GWAS, GTEx eQTL, gene coordinates |
| `planner.js` | the query planner: resolve inputs (hybrid rsID strategy), price scan vs point reads from the manifest, fetch rows; owns the shared file cache (`planner.file`) |
| `parity/pyfmt.js` | Python-exact formatting: `round`, `'%.nf'`, `str(float)`, `nextafter`; C `%g` |
| `core/ld/kernel.js` | one LD API: `prepare`, `band` (row ranges, masked or unmasked), `tile`, `counts`, `exact`. Backends: `kernel.wasm` (C, scalar popcount) and `kernel_js.js`, with identical results. Masked rows carry a mask class: pairs with identical masks take one popcount |
| `ldsc/` | LD score regression (`regress`), sumstats munging (`munge`), LD score estimation (`l2`), the release's panels (`panels`); the shared numerics are `core/stats.js`, ldsc's log text and pandas' layouts `parity/ldsc/log.js` and `pandas.js` |
| `parity/ldlink/context.js` | composes planner, genotypes, annotations, dbSNP, panels (from the release unless the caller passes `readPanel`) and kernel for the module ports |
| `parity/ldlink/common.js` | LDlink's shared pieces: populations, the input `gate`, `validsnp`'s `dedupeLists`, the chrY message, `set_alleles`, kernel counts → 2×2 (`cells`, `table2x2`, `ldStats`), the chunked window (`chunkWindow`, `chunkedWindow`, `countsAgainst`), `get_query_variant_c`, positional annotation lookups |
| `parity/ldlink/<module>.js` | LDmatrix, LDpair, LDpop, LDhap, LDproxy, SNPclip, SNPchip, LDtrait, LDexpress, LDassoc, LDherit, LDcorrelation, LDscore: each returns LDlink's own output (text / JSON / error), byte for byte, including LDlink's bugs ([PARITY.md](PARITY.md)) |

**Rules:**
- **No Node imports in `src/`.** The only bare imports are `hyparquet` and `hyparquet-compressors`.
- **The kernel chain** (`core/ld/kernel.js`, `kernel_js.js`, `core/haps.js`) has no package imports, so
  `app/workers/` can load it without an import map.
- **Data access goes through a caller-supplied `open(path, byteLength?)`** (local files in Node, Range requests in the
  browser; `byteLength` comes from the manifest, so the browser skips its HEAD request). The browser also passes the
  fetched `manifest` to `createContext`.

## app/: the UI

- **`index.html`:** an import map for `hyparquet*` (no bundler).
- **`site.js`:** where the release is: `?data=<url>`, else the URL `tools/build_site.js` wrote, else `../release`.
- **`main.js`:** the core in the page (kernel and release) and routing: `#<Module>` or `#<Module>?<inputs>[&run=1]`.
- **`forms.js`:** the module table: each module's fields, the port it calls, how its result is shown. A form fills
  from the hash query and writes its inputs back when it runs, so results can be bookmarked and linked.
- **`render.js`:** outputs as LDlink prints them (text blocks → tables, LDpair, the matrix, pyramid and locus
  cards), and the links between views (a pyramid pair → LDpair / LDproxy; a locus point → LDpair).
- **`region.js`:** the region view's data (window rows and tracks), LD + pixels on the worker pool.
- **`views/`:**
  - `tracks.js`: annotation tracks on one genomic axis (GWAS, eQTL, FORGEdb / RegulomeDB, recombination, genes
    collapsed or per transcript, the axis), shared by the region view and the locus plots;
  - `matrix.js`: the LDmatrix plot;
  - `pyramid.js`: the region view (tracks + pyramid);
  - `locus.js`: LDlink's LDproxy and LDassoc figures (points, recombination, labels, rug, genes).
- **`workers/`:**
  - `pyramid.js`: band + pixel aggregation in a worker;
  - `pool.js`: shares of equal pair counts, merged by pixel max, deterministic ties.

## pipeline/: building a release (offline)

`pipeline/build.sh <release> [stage ...]` runs every stage below in order and skips a stage whose output exists
(its header lists each stage's output). SQL files run through `tools/run_sql.py`, whose `@for` blocks expand the
per-chromosome and per-partition statements.

| stage | scripts | produces |
|---|---|---|
| `dbsnp/` | `bson_stream`, `bson2tsv`, `sort.sql`, `dup_ids.sql`, `build`, `dup_rsids` | `dbsnp/` |
| `ncbi/` | `fetch.sh`, `vcf2tsv`, `vcf2parquet.sql`, `merged2tsv.py` | NCBI b155 allele slices (the `allele_rsids` input); the one-off comparisons are in `docs/history/ncbi_analysis/` |
| `genotypes/` | `slices.sql`, `build_all.sh` → `build.js` (uses `rsid_label`, `vcf_gt`), `to_release.sh` | `genotypes/` |
| `annotations/` | `fetch.sh`, `annotations.sql`, `build.js` | `annotations/` |
| `ldscore/` | `convert_ref`, `l2_cli` | `ldscore/` |
| `release/` | `manifest`, `extras`, `check` | `manifest.json`, panels, examples, ldsc counts; the integrity report |

**Working directory:** every script runs with it set to `$LD_DATA` (`tools/config.js` `useDataCwd()`, or
`tools/config.sh` + `cd`). Data paths in JS, SQL and Python are therefore data-relative, and repo files resolve from
the script's own location.

## verify/, test/, bench/, tools/

- **`test/`:** unit tests (`npm test`).
- **`verify/data/`:** data verifiers against the dump and the VCFs (planner, parity, annotations, the kernel against
  LDlink's math, row selection).
- **`verify/ldsc/`:** the ldsc ports against ldsc.
- **`verify/ldlink/`:** the module ports against a local LDlink.
- **`verify/run.sh`:** runs them (`npm run verify`, `npm run verify:ldlink`).
- **`bench/`:** benchmarks; they write their results to `$LD_DATA/bench/results/`.
- **`tools/`:**
  - `config.js`, `config.sh`: every machine path;
  - `serve.js`: the dev server; `build_site.js`: the static site;
  - `check_docs.js`: every path the docs name exists (in `npm test`);
  - `run_sql.py`: the DuckDB statement runner;
  - `ldlink-local/`: the oracle kit.
