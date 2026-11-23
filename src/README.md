# src/: the browser-ready core

Everything here runs in a browser or in Node, with no build step. It reads a release (Parquet files + `manifest.json`)
through a caller-supplied `open(path) → AsyncBuffer` (local files, or HTTP Range requests). The only dependencies are
`hyparquet` and `hyparquet-compressors`. The offline pipeline, verifiers and benchmarks live in `experiments/`.

| module | responsibility | lines |
|---|---|---:|
| `genome.mjs` | chromosome codes and names, genome keys, the build table (`BUILDS`, `TITLE`, `tierPath`); LDlink's `firstToken` / `parseCoord` | 39 |
| `source.mjs` | cached `{ file, metadata }` per path, `isNotFound` (the only test for "absent") | 26 |
| `parquet.mjs` | page-index reads: `sortedLookup` (point keys), `rangeScan` (key range), `pool` | 161 |
| `dbsnp.mjs` | LDlink-parity resolution (`resolveRsids`, `resolveCoords`, `resolveInputs`) and the two `find_one` forms (`findAt`, `findById`) | 113 |
| `genotypes.mjs` | parity / allele-aware row selection, `decodeHaps`, tier-row reads (`row`, `windowRows`, `samples`) | 104 |
| `planner.mjs` | cost model (scan vs point reads, hybrid rsID strategy) and query execution; owns the shared file cache | 256 |
| `annotations.mjs` | FORGEdb, RegulomeDB, arrays/platforms, GWAS Catalog, GTEx eQTL (+ tissues, genes), gene coordinates | 75 |
| `pyfmt.mjs` | Python-exact formatting (`round`, `'%.nf'`, `str(float)`, `nextafter`) and C `%g` | 48 |
| `ld/kernel.{mjs,c}`, `ld/kernel_js.mjs` | one LD kernel API, WASM (`kernel.wasm`, built by `ld/build.sh`) and plain-JS backends with identical results | 174 + 145 + 116 |
| `ldsc/{stats,regress,munge,l2}.mjs` | the LDscore family (ldsc ports); no CLI code (the CLI is `experiments/ldsc/ldscore_l2.mjs`) | 59 + 211 + 181 + 114 |
| `ldlink/context.mjs` | composes planner + genotypes + annotations + dbSNP + panels + kernel for the module ports | 29 |
| `ldlink/common.mjs` | LDlink's shared pieces: populations, `set_alleles`, 2×2 statistics, `get_query_variant_c` | 45 |
| `ldlink/<module>.mjs` | LDmatrix, LDpair, LDpop, LDhap, LDproxy, SNPclip, SNPchip, LDtrait, LDexpress, LDassoc: each returns LDlink's own output, byte for byte | 32–164 |

**Dependencies flow one way:**

genome → source → parquet → dbsnp / genotypes / annotations → planner → ldlink/*

`haps.mjs` (the haplotype decoder) has no dependencies, so `ld/kernel` and its imports run in a Web Worker without an import map.

**Using a module port:**

```js
import { loadKernel, KERNEL_WASM } from "./ld/kernel.mjs";
import { createContext } from "./ldlink/context.mjs";
import { ldpair } from "./ldlink/ldpair.mjs";
const kernel = await loadKernel(new Uint8Array(await (await fetch(KERNEL_WASM)).arrayBuffer()));
const ctx = await createContext({ root, open, readPanel, kernel }); // readPanel(pop) → sample IDs
const out = await ldpair(ctx, { var1: "rs3", var2: "rs4", pop: "CEU", build: "grch37" });
```

**Regression suite:** `experiments/run_regression.sh` runs:
- the unit tests;
- the data verifiers (planner, parity, annotations, kernel against LDlink's math, row selection, ldsc);
- all ten module comparisons against a local LDlink (docs/CORE_REVIEW.md §6–8).
