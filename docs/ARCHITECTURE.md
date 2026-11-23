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

## src/: the core (browser-safe, no build step)

Dependencies flow one way:

```
genome ─► source ─► parquet ─► dbsnp, genotypes, annotations ─► planner ─► ldlink/* (context, common, modules)
haps (no deps) ─► ld/kernel  ◄─ pyfmt          ldsc/* (stats ─► regress, munge; l2)
```

| module | one line |
|---|---|
| `genome.mjs` | chromosome codes and names, genome keys, the `BUILDS` table, LDlink's `firstToken` / `parseCoord` |
| `source.mjs` | the one cached `{ file, metadata }` per path (sizes and footer lengths from the manifest: no size probe, one exact footer read); `isNotFound`, `kv` |
| `parquet.mjs` | page-index reads: `sortedLookup` (point keys), `rangeScan` (key range), `pool`; `coalesce` (reads issued together share Range requests) |
| `dbsnp.mjs` | LDlink-parity resolution: `resolveRsids`, `resolveCoords`, `resolveInputs`; the `find_one` forms `findAt`, `findAllAt`, `findById`; `duplicated` |
| `haps.mjs` | `decodeHaps`: bitset or sparse → bitset (no dependencies, so the kernel stays worker-safe) |
| `genotypes.mjs` | parity / allele-aware row selection (`pickParity`, `pickAllele`, `candidateRows`); tier-row reads (`row`, `windowRows`, `samples`) |
| `annotations.mjs` | genes, recombination, RegulomeDB, FORGEdb (by rsID and by position), arrays, GWAS, GTEx eQTL, gene coordinates |
| `planner.mjs` | the query planner: resolve inputs (hybrid rsID strategy), price scan vs point reads from the manifest, fetch rows; owns the shared file cache (`planner.file`) |
| `pyfmt.mjs` | Python-exact formatting: `round`, `'%.nf'`, `str(float)`, `nextafter`; C `%g` |
| `ld/kernel.mjs` | one LD API: `prepare`, `band` (row ranges, masked or unmasked), `tile`, `counts`, `exact`. Backends: `kernel.wasm` (C, scalar popcount) and `kernel_js.mjs`, with identical results. Masked rows carry a mask class: pairs with identical masks take one popcount |
| `ldsc/` | LD score regression (`regress`), sumstats munging (`munge`), LD score estimation (`l2`), shared numerics (`stats`) |
| `ldlink/context.mjs` | composes planner, genotypes, annotations, dbSNP, panels (from the release unless the caller passes `readPanel`) and kernel for the module ports |
| `ldlink/common.mjs` | LDlink's shared pieces: populations, the input `gate`, `validsnp`'s `dedupeLists`, the chrY message, `set_alleles`, kernel counts → 2×2 (`cells`, `table2x2`, `ldStats`), the chunked window (`chunkWindow`, `chunkedWindow`, `countsAgainst`), `get_query_variant_c`, positional annotation lookups |
| `ldlink/<module>.mjs` | LDmatrix, LDpair, LDpop, LDhap, LDproxy, SNPclip, SNPchip, LDtrait, LDexpress, LDassoc: each returns LDlink's own output (text / JSON / error), byte for byte, including LDlink's bugs ([PARITY.md](PARITY.md)) |

**Rules:**
- **No Node imports in `src/`.** The only bare imports are `hyparquet` and `hyparquet-compressors`.
- **The kernel chain** (`ld/kernel.mjs`, `kernel_js.mjs`, `haps.mjs`, `pyfmt.mjs`) has no package imports, so
  `app/workers/` can load it without an import map.
- **Data access goes through a caller-supplied `open(path, byteLength?)`** (local files in Node, Range requests in the
  browser; `byteLength` comes from the manifest, so the browser skips its HEAD request). The browser also passes the
  fetched `manifest` to `createContext`.

## app/: the UI

- **`index.html`:** an import map for `hyparquet*` (no bundler).
- **`site.mjs`:** where the release is: `?data=<url>`, else the URL `tools/build_site.mjs` wrote, else `../release`.
- **`main.mjs`:** the core in the page (kernel and release) and navigation.
- **`forms.mjs`:** the module table: each module's fields, the port it calls, how its result is shown.
- **`render.mjs`:** outputs as LDlink prints them (text blocks → tables, LDpair, the matrix and pyramid cards).
- **`region.mjs`:** the region view's data (window rows and tracks), LD + pixels on the worker pool.
- **`views/`:**
  - `matrix.mjs`: the LDmatrix plot;
  - `pyramid.mjs`: the region tracks + pyramid.
- **`workers/`:**
  - `pyramid.mjs`: band + pixel aggregation in a worker;
  - `pool.mjs`: shares of equal pair counts, merged by pixel max, deterministic ties.
- **`check.html`:** runs 18 cases in the page against the local LDlink.

## pipeline/: building a release (offline)

| stage | scripts | produces |
|---|---|---|
| `dbsnp/` | `bson_stream`, `bson2tsv`, `sort.sql`, `build`, `dup_rsids` | `dbsnp/` |
| `ncbi/` | `fetch.sh`, `vcf2tsv`, `merged2tsv.py`, `compare.sql`, `allele_join` | NCBI b155 allele slices (the `allele_rsids` input) |
| `genotypes/` | `slices.sql`, `build_all.sh` → `build.mjs` (uses `rsid_label`, `vcf_gt`) | `genotypes/` |
| `annotations/` | `fetch.sh`, `annotations.sql`, `build.mjs` | `annotations/` |
| `ldscore/` | `convert_ref`, `l2_cli` | `ldscore/` |
| `release/` | `manifest`, `check`, `presign.sh` | `manifest.json`; the integrity report |
| `profile/` | one-off profiling scripts and their results | — |

**Working directory:** every script runs with it set to `$LD_DATA` (`tools/config.mjs` `useDataCwd()`, or
`tools/config.sh` + `cd`). Data paths in JS, SQL and Python are therefore data-relative, and repo files resolve from
the script's own location.

## verify/, test/, bench/, tools/

- **`test/`:** unit tests (`npm test`).
- **`verify/data/`:** data verifiers against the dump and the VCFs (planner, parity, annotations, the kernel against
  LDlink's math, row selection).
- **`verify/ldsc/`:** the ldsc ports against ldsc.
- **`verify/ldlink/`:** the module ports against a local LDlink.
- **`verify/run.sh`:** runs them (`npm run verify`, `npm run verify:ldlink`).
- **`bench/`:** benchmarks and their recorded results.
- **`tools/`:**
  - `config.{mjs,sh}`: every machine path;
  - `serve.mjs`: the dev server;
  - `run_sql.py`: the DuckDB statement runner;
  - `move.py`: scripted moves with import rewriting;
  - `ldlink-local/`: the oracle kit.
