# Working in this repo (for agents and new contributors)

**Read first:**
1. `README.md` (the layout);
2. `docs/STATUS.md` (where we are);
3. `docs/ARCHITECTURE.md` (how it fits).

Then the doc for the area you touch:
- `docs/DATA.md` for release files and schemas;
- `docs/PARITY.md` for module behaviour;
- `docs/decisions.md` for why things are as they are.

Dated narratives in `docs/history/` are background only. Don't treat them as the current state.

## Invariants

1. **Parity first.** `src/parity/ldlink/*` must keep producing LDlink's exact output, bugs included. Change module behaviour
   only with the LDlink source open (`../nci-webtools-dceg-linkage/server`; cite `file:line` in comments), then run
   that module's comparison.
2. **`src/` is browser code.**
   - No Node imports; the only bare imports are `hyparquet` and `hyparquet-compressors`.
   - Data access goes through a caller-supplied `open(path, byteLength?)`.
3. **The kernel chain stays dependency-free.** `src/core/ld/kernel.js`, `kernel_js.js` and `src/core/haps.js`
   import nothing else, because `app/workers/` loads them without an import map.
4. **WASM and JS backends give identical results.**
   - Change `kernel.c` and `kernel_js.js` together.
   - Rebuild with `npm run build:kernel`; commit `kernel.wasm`.
5. **No machine paths, no outputs in the repo.**
   - Every path comes from `tools/config.js` and `config.sh`.
   - Pipeline and verify scripts call `useDataCwd()` (or `cd "$LD_DATA"`) and use data-relative paths.
   - Repo files resolve from `import.meta.url`.
6. **Dense, commented code.** Each file opens with a paragraph saying what it does and which LDlink behaviour it
   reproduces. Keep lines ≤ 140 characters. The lint (`eslint.config.js`, in `npm test`) checks correctness only
   (undefined names, unused values, unreachable code); add no style rules.

## Which check covers what

| you changed | run |
|---|---|
| anything | `npm test` |
| `src/core/*` (reads), `src/parity/{planner,resolve,select,annotations}.js` | `npm run verify` (planner, parity, annotations, select) |
| `src/core/ld/*` | `npm test` (K tests, both backends) + `npm run verify` (kernel chr1 + chrX) |
| `src/parity/ldsc/*` | `npm run verify` (ldsc steps) |
| `src/parity/ldlink/*` (and anything above) | `npm run verify:snapshot` first; then `npm run verify:ldlink` (or one module: `node verify/ldlink/compare.js <module>`) |
| `app/*` | `npm run verify:app` (every view in Chromium), and a look at the views you touched |
| `pipeline/*` | the stage on one chromosome, then `node pipeline/release/check.js <root>` |

Reports land in `$LD_DATA/verify/`; `summary.txt` has one line per step.

CI (`.github/workflows/checks.yml`) runs `npm test` without data, the snapshot over the hosted release and
`verify:app` against it on every pull request and push, and Pages deploys only after them. It cannot run
`verify` or `verify:ldlink` (the local release and oracle): run those yourself where the table says so.

## Don'ts

- **Don't** push, upload to buckets, or call external services without the user's OK.
- **Don't** delete generated data. It is expensive to rebuild (hours to days of machine time).
- **Don't** scan `$LD_DATA` or `node_modules` with broad searches; `$LD_DATA` is about 1 TB.
- **Don't** print credentials (`../.env`, the oracle kit's `server.env`), presigned URLs or keys.
