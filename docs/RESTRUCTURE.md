# Restructure plan: from working prototype to a codebase you can hold in your head

2026-11-23. The parity work is done: every LDlink module reproduces LDlink byte for byte, in Node and in a browser.
The code that got us there is spread over `experiments/` (60 entries), with six overlapping docs. This plan moves it
to a layout where:
- the directory tree says what each part is;
- one page says where we are;
- an agent or a new reader can get the right context in two or three reads.

## 1. Where we are (measured)

| area | size | state |
|---|---|---|
| `src/` runtime core | 2,900 lines, 30 files | clean layering. Duplication remains in `ldlink/` (below) |
| `app/` browser UI | 585 lines | one 200-line `main.mjs` holds forms, rendering and region logic |
| `experiments/` | 4,400 lines in about 60 files | pipeline, verifiers, comparers, benchmarks, v1 dead code and one-offs, all flat, plus 993 GB of outputs in `out/` |
| `docs/` | 2,470 lines in 6 files | ROADMAP (529) and CORE_REVIEW (514) are dated narratives; facts about the current state are scattered across them |
| repo root | old prototype (`app.js`, `index.html`, `sources.js` 114 KB, `utils.js`) | the only tracked files. **None of the current work is committed.** `.gitignore` names `experiments/node_modules`, but the packages now sit at the root |
| machine paths | about 25 hard-coded `D:/Projects/LDlink/...` and `"out/release_v4"` strings | in comparers, the server, build scripts |
| outside the repo | `../ldlink-local` (docker-compose, `patch_local_vcf.py`, `drive.py`, `server.env`) | the parity oracle; not versioned |

**Duplication in `src/ldlink/`** (each item is one concept written several times):

| concept | copies | target |
|---|---:|---|
| LDlink's list input chain: lines → lower → token lists → dedupe → replace coordinates → gate → not-in-dbSNP / chrY-on-GRCh38 warnings | 6 (ldmatrix, ldhap, snpclip, snpchip, ldtrait, ldexpress) | `common.parseList(ctx, text, build, rules)` |
| 2×2 statistics (A..D → D′, r², MAF, correlated alleles) inline | 4 (ldproxy, ldexpress, ldassoc, ldtrait), besides `common.ldStats` | `ldStats` everywhere |
| chunked LDproxy-style window (`chunkWindow` → rows per chunk → counts) | 2 (ldproxy, ldexpress), plus a variant in ldassoc | `common.windowLD(ctx, build, q, window, chunks)` |
| per-sample genotype classes (diploid / haploid / uncalled) from `decodeHaps` | 6 | `genotypes.calls(row, sampleIdx, haploidMales)` |
| Python `float()` acceptance regex | 2 | `pyfmt.pyFloat` |
| chrY-on-GRCh38 message | 9 | `common.MSG.chrY38(v)` (and the other recurring messages) |
| comparer setup (`createContext` with node file access + panels) | 9 | `verify/lib/node_context.mjs` |

## 2. Principles

1. **One concept, one place.** If a rule from LDlink's code appears twice, it is a function in `ldlink/common.mjs`.
2. **The tree is the table of contents.** Top-level directories are roles (runtime, UI, pipeline, verification,
   benchmarks, tools); file names are nouns from the domain. No `_v2` suffixes, no `experiments/`.
3. **Code and data never mix.** Generated files live under one data root outside the repo, found through one
   config module. No machine paths in code.
4. **Every claim has a command.** Status numbers come from `npm run verify` output; a doc quotes the command that
   produced each number.
5. **Docs split into "now" and "history".** A small set of living docs describe the current system. Dated narratives
   are frozen in `docs/history/`, and decisions are short records.
6. **The browser constraint shapes `src/`:**
   - no Node imports;
   - no bare imports outside `hyparquet*`;
   - the kernel chain stays dependency-free, so it is worker-safe.

## 3. Target layout

```
README.md                 what this is · 60-second map · commands · status line (links docs/STATUS.md)
AGENTS.md                 context for agents/new contributors: invariants, commands, where things live, don'ts
package.json              scripts: serve · test · verify · verify:ldlink · bench · build:<stage>
src/                      runtime core, browser + node, no build step (unchanged layering)
  genome  source  parquet  haps  genotypes  dbsnp  annotations  planner  pyfmt
  ld/       kernel.{mjs,c,wasm} kernel_js.mjs build.sh
  ldsc/     stats regress munge l2
  ldlink/   context common  ldmatrix ldpair ldpop ldhap ldproxy snpclip snpchip ldtrait ldexpress ldassoc
            registry.mjs   ← one record per module: name, params (types, defaults, limits), run, format
app/                      UI only
  index.html  main.mjs (shell, ≈60 lines)  forms.mjs (from the registry)  render.mjs
  views/      matrix.mjs  pyramid.mjs  tables.mjs
  workers/    pyramid.mjs  pool.mjs
  check.html  check.mjs   (browser parity page)
pipeline/                 offline data build (was experiments/build_*, *.sql, fetch_*)
  build.sh <stage…>       stages: fetch · dbsnp · ncbi · genotypes · annotations · ldscore · release
  dbsnp/      bson_stream  json_array_stream  bson2tsv  sort.sql  build.mjs  dup_rsids.mjs
  ncbi/       fetch.sh  vcf2tsv  merged2tsv.py  compare.sql
  genotypes/  build.mjs (was build_genotypes_v2)  rsid_label  vcf_gt  slices.sql  build_all.sh
  annotations/ fetch.sh  annotations.sql  build.mjs
  ldscore/    convert_ref.mjs  l2_cli.mjs
  release/    manifest.mjs  check.mjs  upload.sh (presign/sync)
test/                     unit tests + fixtures (node --test test/)
verify/                   data verifiers and the LDlink oracle
  run.sh                  the regression suite (was run_regression.sh) → $LD_DATA/verify/summary.txt
  data/       planner  parity (+ parity_inputs.py, parity_*.sql)  annotations  select  kernel
  ldlink/     compare.mjs <module>  (one runner) + cases/*.json + adapters/<module>.mjs (fetch LDlink, normalize)
  ldsc/       munge  regress  ref_dump.py
  lib/        node_context.mjs  diff.mjs
tools/
  config.mjs              paths: LD_DATA (data root), LD_RELEASE, LD_PANELS, LD_VCF, LDLINK_URL — env with defaults
  serve.mjs               dev server (Range; mounts from config)
  ldlink-local/           docker-compose.yml  patch_local_vcf.py  drive.py  server.env.example  README (how to stand it up)
bench/                    kept benchmarks (js_vs_wasm, ldmatrix_size, threads, net, planner cost)
docs/
  STATUS.md               one page: what works (numbers + the command), what's next, open decisions
  ARCHITECTURE.md         data flow, module map, dependency rule, release layout, where parity rules live
  DATA.md                 the release contract: every file, schema, key, sort order, metadata (from PARQUET_LAYOUT)
  PARITY.md               LDlink behaviours reproduced, per module (from CORE_REVIEW §6–8; each with LDlink file:line)
  REQUIREMENTS.md         product scope (kept)
  decisions/              0001-wasm-plus-js.md · 0002-parity-vs-allele.md · 0003-vcf-id.md · … (≤ 30 lines each)
  history/                ROADMAP.md · CORE_REVIEW.md · PLAN_REVIEW.md · PARQUET_LAYOUT.md (frozen, dated)
archive/ (or deleted, with git history)  v1 builders, prior-work, superseded benches, old root prototype
```

**Data root (outside the repo; `LD_DATA`, default `../data`):**

```
release/<version>/        what the app reads (now experiments/out/release_v4)
build/                    pipeline intermediates (dbsnp_dump, ncbi, phaseB, annot_src, …)
verify/                   reports (ldlink_compare/*, regression logs, summary)
```

**What moves where** (about 60 `experiments/` entries):
- **pipeline (24):** `build_genotypes_v2`, `rsid_label`, `vcf_gt`, `build_all.sh`, `phaseB_slices.sql`,
  `build_dbsnp_v2`, `bson_stream`, `json_array_stream`, `dbsnp_bson2tsv`, `dbsnp_sort.sql`, `build_dup_rsids`,
  `ncbi_vcf2tsv`, `ncbi_merged2tsv.py`, `fetch_ncbi_dbsnp.sh`, `compare_ncbi.sql`, `allele_join`,
  `build_annotations`, `annotations.sql`, `fetch_annotations.sh`, `ldsc/convert_ldscore_ref`, `ldsc/ldscore_l2`,
  `build_manifest`, `check_release`, `presign.sh`, `run_sql_statements.py`.
- **verify (19):** `verify_*` ×4, `kernel_check`, `parity_inputs.py`, `parity_*.sql`, `compare_ldlink_parse.sql`,
  `compare_*` ×9 (→ one runner + adapters), `ldsc/verify_*`, `ldsc/ref_dump.py`, `run_regression.sh`.
- **bench (6):** `wasm/bench_*` ×3, `net_bench`, `bench_v2`, `bench_annotations`.
- **archive / delete (≈10):**
  - v1 builders: `build_genotypes`, `build_dbsnp`;
  - v1 benches: `bench_dbsnp*`, `bench_genotypes`;
  - `ld_check` (superseded by `kernel_check`);
  - `prior-work/`, `results/` (old step-1/2 JSON, which the docs quote);
  - `out_grch37_chr22.log`, the empty `modules/`, and the root prototype.

## 4. Simpler code (target: about −400 lines, same behaviour)

| change | where | ≈ lines | safety net |
|---|---|---:|---|
| `common.parseList`, `common.MSG`, `ldStats` everywhere, `windowLD`, `genotypes.calls`, `pyfmt.pyFloat` | `src/ldlink`, `src/genotypes` | −200 | the 10 comparisons (all identical today) |
| **module registry:** each port exports `{ name, params, run, format }`. The app's forms, the comparers' argument handling and the docs' parameter tables derive from it | `src/ldlink/registry.mjs` | −120 in `app` + `verify` | browser check page, comparers |
| **one comparer runner** plus 10 adapters (LDlink request + normalization) and case files | `verify/ldlink` | −250 | same reports, same totals |
| `app/main.mjs` split into shell / forms / render / views; drop the `window.__lastPyramid` hook (or gate it on `?debug`) | `app` | ±0, clearer | manual + check page |
| genome keys as numbers, not BigInt (< 2.6e10 is exact; one conversion at read time) | `genome`, `dbsnp`, `annotations`, `parquet` key compare | −30 (46 BigInt sites) | verify_planner, verify_parity, verify_annotations |
| merge `sortedLookup` / `rangeScan` into `pageRanges` + `readRanges` | `parquet.mjs` | −60 | verify_planner (bytes predicted vs measured) |

## 5. Faster (ordered by gain per effort)

| change | effect | measure with |
|---|---|---|
| pass `byteLength` from the manifest to `asyncBufferFromUrl` | removes one HEAD request per file opened (dozens per cold query) | browser check page timings |
| coalesce adjacent page ranges into one Range request (gap ≤ 64 KB) in `parquet.mjs` | LDproxy: 874 reads → about 100; matters more over real HTTP/S3 | `bench/net` |
| `rsid_dup` boolean column on genotype rows (at build) | drops the whole-file `duplicated_rsids` read from the annotation fast path | LDproxy bytes (19.5 MB today) |
| masked kernel: a pruning bound per pair (valid counts) or grouping samples by mask | chrX band about 4× slower per pair than autosomes today | `bench/threads`, K4 |
| numbers for keys (above) | fewer BigInt ops in every lookup | verify_planner timings |
| HTTP/2 hosting (Phase D, your decision) | S3 HTTP/1.1 caps browsers at about 6 connections | `bench/net` |

## 6. Easier to maintain

- **Commit the current state first**, with a correct `.gitignore` (`node_modules/`, `.claude/`, `*.log`, `data/`,
  `server.env`), then tag it `parity-complete`. Today a disk failure would lose everything since the prototype. After
  that, one commit per step below.
- **`package.json` scripts:**
  - `npm run serve`;
  - `npm test` (unit, < 2 s);
  - `npm run verify` (data, about 15 min);
  - `npm run verify:ldlink` (needs local LDlink, about 1.5 h);
  - `npm run bench`.
- **`tools/config.mjs`** is the only place that knows paths. Every script takes its paths from it (env overrides).
  `tools/ldlink-local/README.md` documents standing up the oracle: the 6-line patch, restore, indexes, the
  `apache-bin/python3` link, `PATH` for `docker exec`.
- **Optional checking:** `// @ts-check` + JSDoc types on `src/` public functions, checked by `tsc --noEmit
  --allowJs`, with no build step. This catches shape errors in the module ports, where most of the subtle code is.
- **Style:** keep dense, but cap lines at about 140 characters, and keep a one-paragraph header per file saying what
  LDlink behaviour it reproduces (most already have one).

## 7. Getting context fast (for people and agents)

1. **`README.md`** (≤ 80 lines): what it is, the tree with one line per directory, the five commands, and the latest
   status line.
2. **`AGENTS.md`** (≤ 100 lines):
   - **invariants:** parity first; `src/` stays browser-safe; the kernel chain stays dependency-free;
   - **how to check a change:** which verify step covers which directory;
   - **where the LDlink source is:** `../nci-webtools-dceg-linkage/server`, with the `file:line` habit;
   - **the don'ts:** no machine paths, no outputs in the repo, never push or upload without asking.
3. **`docs/STATUS.md`:** one page, regenerated in part by `verify/run.sh` (it writes the totals table). The
   "what next" list is maintained by hand.
4. **`docs/ARCHITECTURE.md`:** one diagram's worth of text: release → source → parquet → dbsnp / genotypes /
   annotations → planner → ldlink ports → app; plus the worker and kernel path.
5. Everything else is reached from those four.

## 8. Order of work (each step ends green on `npm test` + `npm run verify`; LDlink comparisons at steps 4 and 6)

| step | what | risk | size |
|---|---|---|---|
| 0 | `.gitignore`, commit everything as-is, tag `parity-complete` | none | minutes |
| 1 | data root out of the tree (`LD_DATA`); `tools/config.mjs`; replace hard-coded paths | low (paths only) | ½ day |
| 2 | directory moves (scripted import rewrite, like `scripts/move_to_src.py`); archive dead code; `package.json` scripts | low (mechanical) | ½ day |
| 3 | docs: README, AGENTS, STATUS, ARCHITECTURE, DATA, PARITY, decisions; move the narratives to `history/` | none | ½ day |
| 4 | `src/ldlink` dedup (`parseList`, `ldStats`, `windowLD`, `calls`, `MSG`) + registry; one comparer runner | medium (behaviour); fully covered by the comparisons | 1 day |
| 5 | speed quick wins: manifest `byteLength`, range coalescing | low | ½ day |
| 6 | numbers for keys; `pageRanges` / `readRanges` | medium; covered by verify_planner / parity | 1 day |
| 7 | `rsid_dup` column (GRCh37 + GRCh38 rebuild, about 1–5 h machine time) and the masked kernel bound | low / medium | 1 day + machine time |
| 8 | back to Phase C features (tracks, locus plots, LDscore views) on the cleaner base | — | — |

## 9. Decisions for you

1. **Commit and push (step 0).** Commit locally now; pushing to `github.com/park-brian/linkage-disequilibrium` is
   outward-facing, so it needs your go-ahead. The root prototype would be replaced in the same push.
2. **Data root location**, e.g. `D:\Projects\LDlink\data` (contains `release/`, `build/`, `verify/`), and whether to
   move `ldlink-local` into the repo as `tools/ldlink-local` (code only; its 200 GB of Mongo and dump data stays out).
3. **Archive or delete the dead code.** Deleting is cleaner once committed, since git keeps it.
4. **Disk:** about 870 GB of regenerable intermediates (`ncbi_dbsnp` 632 GB, `dbsnp_dump` 192 GB, the `phaseB` build
   copies). Keep or delete after step 1.
5. **Type-checking** (`@ts-check` + `tsc`): adopt, or skip.

---

## 10. Implementation plan (approved 2026-11-23: "plan carefully, then fully implement")

**Decisions taken by default.** These are reversible; each slice states the default it uses.

| § 9 question | default | why |
|---|---|---|
| commit / push | commit on branch `restructure`; tag the baseline `parity-complete`; **no push** | pushing is outward-facing; you push when ready |
| data root | `D:/Projects/LDlink/data` (`LD_DATA`), moved from `experiments/out` by rename | same volume: instant, and the hardlinks are kept |
| `ldlink-local` | its code is versioned in `tools/ldlink-local/` (source of truth); the running stack in `../ldlink-local` is left in place and keeps working | moving the live directory would mean restarting Docker volumes for nothing |
| dead code | deleted after the baseline commit (git keeps it) | an archive directory is clutter |
| regenerable intermediates (about 870 GB) | **kept**; nothing is deleted | deleting is irreversible and was not asked for |
| `@ts-check` + `tsc` | **not adopted**: it needs a TypeScript download | it can be added later without layout changes |
| rebuilds and uploads | no bucket upload. A release rebuild only if a slice needs it (slice 7 avoids that) | the upload needs your OK |

**Acceptance gates.**

| gate | command | runtime | covers |
|---|---|---|---|
| G1 | `npm test` | < 5 s | unit tests (25), docs path check |
| G2 | `npm run verify` | about 20 min | planner, parity, annotations, kernel_check, select, ldsc |
| G3 | `npm run verify:ldlink` | about 1.5 h; needs `../ldlink-local` running | the 10 module comparisons |
| G4 | browser | — | `app/check.html` 18/18; the app's Pyramid and LDmatrix views render |

Every slice must pass G1 and G2 before its commit. Slices that change module behaviour (4, 6, 7) also pass G3. Slices
that touch the app (2, 4, 5) also pass G4.

### Slice 1: data root and one config module

- **Move:**
  - `experiments/out` → `../data` (rename);
  - `experiments/results` → `../data/results` (historical measurement JSON).
- **New `tools/config.mjs`.** Exports:
  - `DATA` (`LD_DATA`, default `<repo>/../data`);
  - `RELEASE` (`LD_RELEASE`, default `DATA/release_v4`);
  - `PANELS` (`LD_PANELS`, default `<repo>/../db/data/1000G_population_samples`);
  - `VCF` (`LD_VCF`, default `<repo>/../db/1000G`);
  - `LDLINK_URL` (default `http://127.0.0.1:8090`);
  - `LDLINK_LOCAL` (default `<repo>/../ldlink-local`).

  Shell scripts source `tools/config.sh`, which exports the same names.
- **Every script:** `out/…` and `D:/Projects/…` literals become config paths. Scripts that took a root argument keep
  it, defaulting to `RELEASE`.
- **Accept:** a grep for machine paths and `"out/` in code finds nothing outside docs and comments; G1, G2.

### Slice 2: directory moves, dead code out, npm scripts

- **Scripted move:** `tools/move.py`, a mapping table that rewrites relative imports *and*
  `new URL("…", import.meta.url)` strings. Layout per § 3:
  - `experiments/*` → `pipeline/`, `verify/`, `bench/`;
  - `experiments/test` → `test/`;
  - `scripts/serve.mjs` → `tools/serve.mjs`;
  - `scripts/profile` → `pipeline/profile`;
  - `app/pyramid_*` → `app/workers/`, `app/*_view` → `app/views/`.
- **Delete** (after a grep confirms no references):
  - v1 `build_genotypes` and `build_dbsnp`;
  - `bench_dbsnp*`, `bench_genotypes`, `ld_check`;
  - `prior-work/`, empty `modules/`, the stray log;
  - the root prototype (`app.js`, `index.html`, `sources.js`, `utils.js`); its test vectors already live in the app;
  - `scripts/move_to_src.py`.
- **`tools/ldlink-local/`:** `docker-compose.yml`, `patch_local_vcf.py`, `drive.py`, `server.env.example`
  (placeholders only), and a README with the stand-up steps.
- **`package.json` scripts:** `serve`, `test`, `verify`, `verify:ldlink`, `bench`, `build` (→ `pipeline/build.sh
  <stage…>`, dispatching to the existing stage scripts without rewriting them).
- **Accept:** G1, G2, G4; `node --check` on every moved `.mjs`.

### Slice 3: docs

- **New, all short:**
  - `README.md` (≤ 80 lines);
  - `AGENTS.md` (≤ 100);
  - `docs/STATUS.md` (one page);
  - `docs/ARCHITECTURE.md`;
  - `docs/DATA.md` (the release contract: every file, schema, key, sort, metadata);
  - `docs/PARITY.md` (per module, the reproduced LDlink behaviours with LDlink file:line);
  - `docs/decisions.md` (≤ 15 lines per decision).
- **Moved unchanged to `docs/history/`:** ROADMAP, CORE_REVIEW, PLAN_REVIEW, PARQUET_LAYOUT, the old DATA,
  RESTRUCTURE (at the end). REQUIREMENTS stays.
- **Rule:** new docs summarize and link; nothing is deleted.
- **Accept:** every number in STATUS cites its command; `tools/check_docs.mjs` (part of `npm test`) confirms that
  every repo path named in README / AGENTS / ARCHITECTURE / STATUS exists.

### Slice 4: `src/ldlink` dedup, registry, one comparer runner

**`common.mjs` gains:**
- **`parseList`:** LDlink's list chain (route pre-filter, `validsnp`, `replace_coords_rsid_list`, gate, not in
  dbSNP, chrY on GRCh38). The per-module message variants are passed in, so the exact texts stay as they are.
- **`windowLD`:** LDproxy's chunked window, shared by ldproxy and ldexpress.
- **`ldStats`** used everywhere instead of the inline 2×2 code, returning `maf_q` / `maf_p` too. Correlated-allele
  orientation is passed in.

**Elsewhere:**
- **`genotypes.calls`:** per-sample called / haploid / alleles.
- **`pyfmt.pyFloat`:** Python `float()` acceptance and value.
- **`src/ldlink/registry.mjs`:** `{ module: { params, run } }`. Params are plain `{ name, type, default, label }`.
  `app/forms.mjs` builds its forms from it; escape hatches stay in the modules.

**Comparers:** `verify/ldlink/compare.mjs <module|all>`, one runner over per-module adapters (LDlink request +
normalization) and case generators. Reports go to `DATA/verify/ldlink/`. The 9 `compare_*.mjs` files go.

**Accept:** G1, G2, G4, and G3 at today's totals: LDmatrix 30/30, LDpair 429 (identical + tie order), LDpop 537,
LDhap 70, LDproxy 57, SNPclip 49, SNPchip 66, LDtrait 55, LDexpress 58, LDassoc 21.

### Slice 5: fewer requests

- **Manifest sizes:** pass `byteLength` from the manifest when opening files (no HEAD per file in the browser).
- **Range coalescing** in `parquet.mjs`: page ranges within 64 KB of each other become one read.
- **Accept:** G1, G2, G4; LDproxy read counts and browser cold timings recorded before and after in STATUS.

### Slice 6: numeric keys, unified page reads

- **Numeric genome keys:** `genome.key` returns a Number; INT64 is converted once at read time.
- **Unified page reads:** `sortedLookup` and `rangeScan` share one page-range primitive.
- **Accept:** G1, G2 (verify_planner byte ratios unchanged), G3.

### Slice 7: annotation fast path and masked kernel

- **Duplicate-rsID test without a rebuild.** Measure LDproxy's bytes first. If the whole-file
  `duplicated_rsids` read dominates, replace it with a page-pruned lookup of the window's rsIDs. The `rsid_dup`
  column waits for the next rebuild (noted in STATUS).
- **Masked kernel fast path.** Pairs whose two rows share the same valid mask use the unmasked formula with the shared
  valid count. Masks are grouped by hash at prepare time. In C and JS, checked by K4, a new K5 (fast path == masked
  path) and `kernel_check` on chrX.
- **Accept:** G1, G2, G3 (chrX cases), G4 (chrX pyramid timing before and after).

### Final

- A full G1–G4 run;
- STATUS updated with the final numbers;
- this file moved to `docs/history/`;
- a final commit.

No push.

---

## 11. Review outcome and revised slices (2026-11-23)

An independent read-only review of §10 against the code found 18 concrete problems. The top three were verified
directly:
- the `parity-complete` commit lacks the `ref_*.log` test fixtures, so it fails its own tests from a fresh clone;
- `.gitignore`'s `data/` also hides the planned `verify/data/`;
- G1 *skips* the release-dependent tests, which read `../out/release`, when the data is absent.

**The revisions:**

| # | finding | change |
|---|---|---|
| 1 | `.gitignore` drops the `ref_*.log` fixtures; `data/` matches `verify/data/` | `/data/`, `!**/fixtures/**/*.log`; commit the fixtures now |
| 2 | G1 skips silently; the tests read `out/release` (an older release than `release_v4`) | tests take the release from `tools/config.mjs`; they **fail** when it is configured but missing |
| 3 | path literals are mostly in SQL (single-quoted), shell (unquoted), Python and template strings, which the slice 1 grep missed | **SQL / `run_sql_statements.py` / duckdb run with `cwd = LD_DATA`** and the `out/` prefix stripped, not templated; shell sources `tools/config.sh`; JS and Python read the config |
| 4 | uncovered paths: the `ldlink-local` tmp exchange directory, the LDassoc examples, the 1000G VCFs, the emsdk clang, `kernel_check`'s reference script with a sibling repo's `node_modules` | config gains `LDLINK_TMP`, `LDASSOC_EXAMPLES`, `VCF`, `CLANG`, `REF_NODE_MODULES` |
| 5 | test inputs under `out/` (`ldlink_compare/ldmatrix_cases.json`); `experiments/results` is tracked and written by benchmarks | the cases move to `verify/ldlink/cases/`; results stay in git as `bench/results/` |
| 6 | `scripts/profile` holds the kernel reference and test-vector expectations used by verifiers | → `verify/reference/` (ldmatrix_reference, results/); the pure profiling scripts → `pipeline/profile/` |
| 7 | every script assumes `cwd = experiments/` | resolve sibling files from `import.meta.url` / `$(dirname "$0")`. **Slices 1 + 2 become one scripted pass**, so about 60 files are edited once |
| 8 | `tools/ldlink-local` compose mounts relative dirs, builds from a sibling repo, and has a literal (local throwaway) password | the kit is versioned as a **template**: compose takes `LDLINK_HOME` / `LDLINK_SRC` / `MONGO_ROOT_PASSWORD` from `.env` (gitignored; `.env.example` committed). The comparers pipe the repo's `drive.py` (one source of truth). The running stack is untouched |
| 9 | unplaced files: `experiments/README.md`, the `wasm/*.json` reports, `run_sql_statements.py`, `allele_join` | README → folded into the docs; reports → `bench/results/`; `run_sql_statements.py` → `tools/`; `allele_join` → `pipeline/ncbi/` |
| 10 | deletes are safe (only docs mention them) | unchanged |
| 11 | G4 misses LDhap / LDassoc and a rendered pyramid from the moved worker | G4 = check page 18/18 **+** the app runs Pyramid (chr22 and chrX) and LDassoc (example) and renders them |
| 12 | `parseList` would be one function with 6 behaviours (split, limit order, labels, chrY messages, warning joins, select) | **dropped.** Extract only `gate`, `dedupeLists` and recurring message builders |
| 13 | `ldStats` everywhere would drag BigInt χ² and erfc into window loops and change per-module filter order | **narrower:** `table2x2(n, a1, a2, both)` → A–D, Ms, delta, r², D′; modules keep their own order and allele strings |
| 14 | `windowLD` is shared by ldproxy and ldexpress only (ldexpress inserts `multiQuery`; ldassoc differs) | **narrower:** `windowCounts` used by those two |
| 15 | a cross-cutting registry won't pay (app converts types, comparers send strings, ports validate in their own order) | **dropped.** The app's `MODULES` table moves to `app/forms.mjs` as is |
| 16 | the comparer runner must land separately from the dedup, or a G3 difference can't be attributed | **new slice 3:** runner first, at today's totals; dedup after |
| 17 | `asyncBufferFromUrl` takes `byteLength`, but `openFile(path)` has one argument everywhere; the manifest fetch still costs a HEAD; `source.mjs` re-reads dbSNP footers (its 64 KB key ≠ the planner's exact one); "range coalescing" can't work in `parquet.mjs` (row ranges, pool(16); hyparquet already coalesces inside a read); hyparquet 1.33 does page-index pruning itself | slice 5: `open(path, size)` signature everywhere, the manifest via `fetch`, the footer double-read fixed, coalescing via planned reads (`maxGapRows`) or none. Evaluate hyparquet's native pruning before writing any primitive. **Slices 5 + 6 merge** so `parquet.mjs` changes once |
| 18 | only 5 BigInt sites are genome keys; the rest reproduce Python's unbounded `int()`. The masked same-mask fast path is sound only per pair inside `band_masked` (grouping changes output order); G2/G3 don't exercise `band()` on masked rows; module chrX speed is `js.maskedCounts` | numbers for **keys only**; the fast path is checked per pair (hash, then full-mask compare, per-row valid counts) in C, JS and `maskedCounts`. G2 gains `kernel_check` on chrX; K5 tests fast == masked |

**Revised slices** (each slice ends green on G1 + G2, and G3 / G4 where marked; one commit per slice):

0. **Baseline fix:** `.gitignore`; commit the fixtures and this section.
1. **Layout + paths (one scripted pass):**
   - `tools/config.{mjs,sh}`;
   - move the data root;
   - the directory moves (§3 as amended above), with `import.meta.url` resolution;
   - SQL with `cwd = LD_DATA`;
   - dead code deleted;
   - npm scripts and the `tools/ldlink-local` template;
   - tests read the config and fail when the release is missing.

   Gates: G1, G2, G4.
2. **Docs:** README, AGENTS, STATUS, ARCHITECTURE, DATA, PARITY, decisions; history; `tools/check_docs.mjs` in G1.
3. **Comparer runner:** `verify/ldlink/compare.mjs` + adapters + cases, at today's totals. Gate: G3.
4. **Narrow dedup:** `gate`, `dedupeLists`, messages, `table2x2`, `windowCounts`, `genotypes.calls`,
   `pyfmt.pyFloat`; `app/forms.mjs`. Gates: G3, G4.
5. **Reads:** `open(path, size)`, the manifest via `fetch`, the footer double-read, numeric genome keys, an
   evaluation of hyparquet's native page pruning (adopted only if `verify_planner`'s bytes do not regress), and
   coalescing via planned reads. Gates: G2, G3, G4; before / after numbers in STATUS.
6. **Annotation fast path** (page-pruned duplicate lookup) **and the masked kernel same-mask fast path**:
   C + JS + `maskedCounts`, K5, `kernel_check` on chrX. Gates: G3, G4 (chrX pyramid timing).
7. **Final:** full G1–G4, STATUS numbers, this file → `docs/history/`.
