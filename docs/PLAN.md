# Plan: smaller, checked on every push, easier to keep correct

**Written:** 2026-11-23, after the cleanup (`docs/history/CLEANUP_PLAN.md`), on `main` = `c83ad9a`. Each slice ends
green on its gates and is one commit; this file moves to `docs/history/` with an outcome when done.

**Where the size is** (lines, without data files): `src/core` 1,137 · `src/parity` 2,689 · `app` 1,032 ·
`pipeline` 1,729 · `verify` 2,033 · `test` 1,751 (with fixtures) · `tools` 492 · living docs 699.

**Principles:** measure before deleting; parity first (LDlink's output may not change by a byte); every slice
leaves fewer concepts, not just fewer lines; a check that can't run on every push is a check that rots.

## Slices

1. **Measure the read path as users get it** (no code changes).
   - The question: the planner's byte-cost model (`scanBytes`, `pointBytes`, `coordBytes`, `directRsidBytes`, the
     hybrid rsID strategy, `COST`, `predict`: about 150 of `planner.js`'s 320 lines) was tuned before slice 6 changed
     how reads work. Does it still beat simple rules: point reads for positions, range scans for regions, every rsID
     resolved directly?
   - How: `bench/app_requests.js` per module (LDmatrix 8 and 300 SNPs, LDproxy, SNPclip 1,000 and 5,000, LDtrait,
     LDassoc, the region view), over the dev server and over CloudFront: requests, bytes and wall time, with the
     planner as is and forced to the simple rules (`force: "points"`, `forceRsid: "direct"` exist already).
   - Also re-measures what STATUS lists as unknown (LDproxy over CloudFront was 505 requests / 14 s before slice 6).
   - Output: a table in `decisions.md` and a yes / no for slice 2.
   - **Done** (`decisions.md`, "The planner fetches no genotype rows for the ports"): the simple rules lose (resolving
     big lists directly costs 3–7× the requests; always scanning reads up to 1.8× the bytes on sparse lists), but
     `query` fetches genotype rows that no port uses. Without that fetch every list module needs 5–37% fewer
     requests and up to 29% fewer bytes, with identical output.
2. **The planner without its row fetch.**
   - `query` resolves and selects; it no longer fetches genotype rows. Gone with the fetch: `fetchRows`, `readRange`
     / `readPoints`, the scan-or-points pricing (`pointBytes`, the genotype side of `scanBytes`), `plan.fetch`,
     `predicted.genotypes`, the `force` and `regions` options (the region view reads with `windowRows`).
   - Stays: the hybrid rsID strategy and the label-scan pricing it decides with (`scanBytes` over the label
     columns), `select`.
   - The manifest keeps per-column sizes only for the label columns the hybrid strategy prices, so the app's first
     download (1.9 MB, 0.6 MB gzipped) shrinks; `PRICED_COLUMNS` moves next to it.
   - `verify/data/planner.js` keeps its resolution checks; its fetch-plan checks go. `bench/planner_modes.js` keeps
     the rsID modes.
   - Gates: G1, G3s, then full G2 + G3; `verify:app`; the slice 1 bench (no regression). Publishing the smaller
     manifest is an upload to the bucket (your OK); until then the current one still works.
   - **Done** (`c62a8f3`): `planner.js` 320 → 214 lines; manifest v3 0.91 MB (0.36 MB gzipped), down from 1.92 MB
     (0.61 MB), carrying only each file's size and footer and the row groups the hybrid strategy prices (the rest had
     no reader). Gates: G1 36/36, G3s 1,440/1,440 (and 459/459 over CloudFront on its v2 manifest), G2 all pass, G3
     13 modules 0 differ, `verify:app` 14/14 with figures unchanged, the bench reproduces slice 1's numbers.
     Published 2026-11-24, after the code went live: `v5/manifest.json` is v3 (CloudFront invalidated; the v2 file
     backed up in `$LD_DATA/manifest_backup/`); the live app and 149 snapshot cases over CloudFront run on it.
     Differences from the plan: `PRICED_COLUMNS` stays in the planner, its one reader (the manifest builder imports
     it); the manifest also lost what nothing read, beyond the priced columns.
3. **The parity checks in CI.** A workflow on every push and pull request:
   - `npm test` with the release (not `LD_SKIP_DATA`) is impossible in CI (34.8 GB), but the snapshot runs over HTTP:
     `npm run verify:snapshot -- --http <CloudFront URL>`, and `npm run verify:app` with the app pointed at
     CloudFront (Playwright's Chromium, cached).
   - First measure what one snapshot run reads over HTTP (bytes × CloudFront's price per GB). If a full run is
     too costly, run a fixed subset on pushes and the full set weekly.
   - The Pages deploy then waits for it. Gates: a green run, and a deliberately broken port turning it red.
   - **Done** (`5536f66`): `checks.yml` (unit, parity, app jobs), called by `pages.yml` before every deploy. A full
     snapshot over CloudFront reads 1.1 GB in 7,400 requests (about $0.09 at list price; CloudFront's free tier is
     1 TB a month), so the full set runs on every push. On GitHub: unit under a minute, parity 20 minutes
     1,440 / 1,440, app 14 / 14; a one-word change to an LDpair message failed parity on its one case. Also:
     Node reads over HTTP cap at 64 in flight (`verify/ldlink/lib.js` `httpFetch`); SNPclip's and SNPchip's dense
     cases read through the context, so the snapshot needs no local data; `verify:app --data <url>`. Figure hashes
     stay a local check (Linux Chromium draws text differently).
4. **A linter for what tests miss** (adds `eslint` as a dev dependency: your OK). A flat config with only
   correctness rules (`no-undef`, `no-unused-vars`, `no-unreachable`, `no-dupe-keys`, `no-self-assign`,
   `no-constant-condition`), no style rules (the code's density is deliberate), in `npm test`. The cleanup found
   two bugs of exactly this kind by hand (an undefined `labels`, an unused leftover `join`). Gate: G1 clean, every
   finding fixed or justified in the config.
   - **Done:** `eslint` 10 with 29 correctness rules (`eslint.config.js`), run by `test/lint.test.js`. First run: 58
     findings, none a live bug. 11 unused imports and values were removed (`TITLE` in LDexpress, `MASTHEAD` in LDherit,
     `normPdf` / `erfc` / `VALID` in the ldsc ports, bench, tools and verifier leftovers). 12 were `verify/app/run.js`'s
     `page.evaluate` callbacks (browser globals now). 35 were `src/core/stats.js`'s erfc coefficients, kept as
     published (rewriting them to the shortest form gave the same 118 doubles; the rule is off for that file, with
     the reason). Gates: G1 37 / 37 (a probe with an undefined name fails it), G3s, `verify:app` 14 / 14.
5. **The one parity gap left:** LDherit, LDcorrelation and LDscore on uploaded files (only LDlink's examples are
   compared). Add G3 cases that upload through LDlink's upload route (a modified copy of an example, plus files
   that trip `sumstats_normalizer` and the compatibility checks), port what differs, refresh the snapshot after a
   full G3. Gate: G3 with the new cases, 0 differ.
6. **Docs that match the code.** `docs/REQUIREMENTS.md` still describes the target of the first phase ("one LD module");
   mark what `core/ld.js` and the app now cover, and keep the rest as the short list of open feature work, or move
   it to history and keep that list in STATUS. Close the open decisions in STATUS (below). Gate: G1 (docs check).

## Later, only if wanted (new features, not tidying)

- `annotate` and regions in `core/ld.js` (core annotation reads: the parity layer's queries return every match).
- The app's own views (region view, LDmatrix plot) on `core/ld.js` instead of the ports' outputs.
- `// @ts-check` on `src/core` (JSDoc types; no build step).

## Waits for your OK

- Slice 2: uploading a new manifest to the bucket. Slice 4: adding `eslint`. Slice 3: CloudFront egress for CI runs.
- Still open from the cleanup: the published verifier reports under `ldlink-parquet/v5/` (delete or keep); the
  bucket's old prefixes `v2/`–`v4/` (94 GB); the gate's scratch output in `$LD_DATA` (about 200 MB); about 870 GB of
  regenerable intermediates in `$LD_DATA`.
