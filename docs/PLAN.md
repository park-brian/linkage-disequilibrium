# Plan: Phase C (the app) and hosting

**Written:** 2026-11-23, on branch `restructure`. Implemented in slices; each slice ends green on its gates and is one
commit. When done, this file moves to `docs/history/` with an outcome section, as RESTRUCTURE.md did.

## 1. Where we start

- **Parity:** all ten modules reproduce a local LDlink byte for byte (STATUS.md).
- **The app** runs every module, an LDmatrix plot, and a region pyramid under GWAS, eQTL, recombination and gene
  tracks.
- **Outside the release, the app still reads two things:** the population panels (`/panels/`, from `$LD_PANELS`)
  and the LDassoc examples (`/examples/`, from the oracle kit). Both are dev-server mounts.
- **Hosting:** the repo is public and GitHub Pages is on: a legacy build serving `main`'s root at
  `park-brian.github.io/linkage-disequilibrium/`. Remote `main` (`64bdf76`) is an ancestor of `restructure`.

## 2. Decisions this plan makes

1. **The release becomes self-contained.** Panels go in as one Parquet file; the example inputs go in as plain files
   under `examples/`. After that the app needs one data URL and nothing else.
2. **App and data are hosted apart.** GitHub Pages limits a site to 1 GB and each file to 100 MB; the release is
   34.8 GB. So Pages serves the app (a static build), and the release lives on a host that serves byte ranges with
   CORS. The app takes the data URL from `?data=` or a build-time default.
3. **LDlink's plots are the spec.** For the LDproxy and LDassoc locus plots, colours, sizes, axes, labels and
   hovers follow `server/LDcommon.py` `ldproxy_figure` and `server/LDassoc.py`. The region view's FORGEdb and
   RegulomeDB track is our addition (LDlink shows those scores only as point labels and hovers), and is labelled as
   such.
4. **The LDscore family joins the parity suite.** The oracle container has `ldsc.py` and LDlink's
   `ldscore.ldsc_utils` but no `/data/ldscore`. The data is local (`$LD_DATA/ldscore_src/ldscore`, 3.4 GB), so
   mounting it there makes LDherit, LDcorrelation and LDscore comparable byte for byte, like the other modules.
5. **Measure before replacing our page-index code.** hyparquet 1.33 prunes pages itself (`filter` over column +
   offset indexes). It replaces `sortedLookup` only if requests and bytes do not regress.

## 3. Gates

| gate | what | added in this plan |
|---|---|---|
| G1 | `npm test` | tests for panels, the data URL, range tracks, LDscore text formatting |
| G2 | `npm run verify` | range FORGEdb / RegulomeDB lookups equal point lookups over sample windows |
| G3 | `npm run verify:ldlink` | adapters for `ldherit`, `ldcorrelation`, `ldscore` |
| G4 | browser: `/app/check.html` + the app | the built site, cross-origin; locus plots; cross-links; LDscore views |

A slice that changes a module's code reruns that module's G3 comparison; G3 runs in full at the end.

## 4. Slices

### 1. A self-contained release

- `pipeline/release/extras.mjs <release>`:
  - writes `panels.parquet` (`pop`, `sample`, file order; 32 panels, about 140 KB of text);
  - copies the example inputs into `examples/`:
    - LDassoc's `prostate_example_grch37.txt` and `prostate_example_grch38.txt`;
    - LDscore's `22.bed`, `22.bim`, `22.fam`;
    - LDherit / LDcorrelation's `BBJ_HDLC22.txt` and `BBJ_LDLC22.txt`;
  - then rebuilds the manifest. `panels.parquet` is listed like any Parquet file. `examples/` is listed with
    `byteLength` only, under a new top-level `extras` key.
- `createContext`: `readPanel` becomes optional. The default reads `panels.parquet` through the planner's file
  cache, so the app, check page and comparers stop passing one.
- Callers:
  - the app and check page fetch examples from `<data>/examples/`;
  - the dev server drops `/panels/` and `/examples/`;
  - `$LD_PANELS` stays for the pipeline and the kernel verifier, which read the original lists.
- Gates: G1; G3 for `ldmatrix`, `ldpair` and `ldassoc` (panels and the example path); G4.

### 2. The data URL and a static site

- `app/main.mjs` and `app/check.mjs` take the release URL from `?data=<url>`, else from `app/site.json`
  `{ "data": … }`, else `../release`. Examples resolve against that URL. A wrong URL fails with a message naming the
  URL.
- `npm run build:site` (`tools/build_site.mjs`) writes `dist/`:
  - `app/`, the `src/` files the app imports, and the four browser dependencies copied to `dist/vendor/`, with the
    import maps rewritten;
  - `app/site.json` from `--data`;
  - `.nojekyll`, and `index.html` redirecting to `app/`.

  Nothing else from the repo is copied.
- `tools/serve.mjs --cors` sends `Access-Control-Allow-Origin: *` and exposes the range headers. That lets a local
  run serve `dist/` and the release from different ports, which reproduces the hosted setup.
- `.github/workflows/pages.yml` builds `dist/` and deploys it with `actions/deploy-pages`. It is written and run
  locally through the build script only. Switching Pages to it is a step that needs your OK (§6).
- Gate: `dist/` served on `:8096` with the release on `:8095` (cross-origin); `check.html` 18/18, and every view
  runs.

### 3. One tracks module

- Move the region view's tracks out of `app/views/pyramid.mjs` into `app/views/tracks.mjs`: the axis, genes
  (collapsed), recombination, GWAS and eQTL. Each track has the same shape: `{ height, draw(g, x, top), hover(pos, y) }`.
  `renderPyramid` composes them. Behaviour must not change.
- Gate: a canvas hash of a fixed chr22 region, taken before and after, is identical; G4.

### 4. Locus plots for LDproxy and LDassoc

- **Ports:**
  - LDproxy and LDassoc already build their rows before formatting text. Each now also returns
    `{ plot: { query, rows, window } }` with those same rows: pos, rsid, R², D′, MAF, function, FORGEdb, RegulomeDB.
  - LDassoc's rows add p.
  - The text is unchanged: G3 for both modules proves it.
- **Plots** (`app/views/locus.mjs`), following LDlink:
  - **Points:**
    - x: position (Mb, ±10 kb padding).
    - LDproxy y: R² or D′ on [−0.03, 1.03], drawn only when R² > 0.01. Query point blue; a variant with a function
      class red; the rest orange. Size `9 + 14·MAF`.
    - LDassoc y: −log₁₀ p. Colour on a white-to-red ramp by R² (D′ when chosen); index point blue; variants not in
      1000G hollow gray. Size `9 + 14·MAF`; a 5e−8 line.
  - **Recombination:** a line on a right-hand 0–100 cM/Mb axis.
  - **Below:** a rug strip, then the gene track from slice 3.
  - **Labels:** FORGEdb (default), RegulomeDB, or none (LDlink's `annotate`).
  - **Hover:** LDlink's fields.
- **Form:** an "Annotate" select is added to both forms.
- Gates: G3 `ldproxy`, `ldassoc`; G4 (both plots render for the default inputs and on GRCh38).

### 5. A FORGEdb / RegulomeDB track on the region view

- `src/annotations.mjs` gains `forgeRange` and `regulomeRange`, which are `rangeScan`s on `key`.
- The region view adds a scores track:
  - FORGEdb score (0–10) as bars;
  - RegulomeDB rank as a categorical row;
  - titled "FORGEdb / RegulomeDB (not in LDlink's region views)".
- `verify/data/annotations.mjs` gains a check that the range reads equal the existing point reads on sample windows.
- Gates: G1, G2 (annotations), G4.

### 6. Links between views

- Each form fills its fields from the hash: `#LDpair?var1=rs3&var2=rs4&pop=CEU&build=grch37`. Running a form writes
  its inputs back to the hash, so a result can be bookmarked.
- Clicking a pyramid pixel opens a small menu: LDpair for the pair, and LDproxy for either variant, with the
  pyramid's population and build.
- Clicking a locus-plot point opens LDpair (query × point).
- Gate: G4 (each link opens the module with the right inputs and runs).

### 7. The LDscore family

- **7a. Oracle and LDherit.**
  - Mount `$LD_DATA/ldscore_src/ldscore` read-only at `/data/ldscore` in `tools/ldlink-local/docker-compose.yml`.
    Recreate only the backend container; Mongo is untouched.
  - Port `server/sumstats_normalizer.py` (LDSC-ready, PLINK, REGENIE and SAIGE inputs) to
    `src/ldlink/sumstats.mjs`.
  - Port the `ldherit` route + `ldsc_utils.run_herit_command` text (ldsc's log from "Total Observed scale", with
    LDlink's provenance line) to `src/ldlink/ldherit.mjs`, over `src/ldsc` and the release's LD score panels.
  - Add an `ldherit` adapter (the example, each population, liability scale, malformed inputs).
- **7b. LDcorrelation:** the same, for `--rg` (`ldcorrelation.mjs`, adapter).
- **7c. LDscore:** `ldscore.mjs` over `src/ldsc/l2.mjs` (PLINK bed/bim/fam uploads, cM or kb window). Output is
  ldsc's log text from "Summary of LD Scores"; adapter on the `22.*` example and window variants.
- **App:** three forms with file inputs and "use the example" buttons. Results are shown as text, as LDlink does.
  The LD score panels are read by column (`SNP`, `L2`) and by the chromosomes the sumstats cover.
- Gates: G1 (formatting units), G3 for each new module, G4.

### 8. SNPchip selects every array by default

LDlink's form starts with all platforms selected (`client/src/app/snpchip/form.tsx`). The app's SNPchip form will
too. Gate: G4 (the default run shows a table).

### 9. hyparquet's page pruning

- Prototype `sortedLookup` and `rangeScan` over `parquetRead({ filter })`.
- Measure requests and bytes on `verify/data/planner.mjs`'s cases and on `bench/app_requests.js`.
- Adopt it, and delete our page-index code, only if neither number regresses and G2 / G3 stay green. Otherwise record
  the numbers in `decisions.md` and keep ours.

### 10. Final

Full G1–G4, STATUS (gates, speed, what works), this file → `docs/history/` with its outcome; one commit.

## 5. Order and dependencies

- 1 → 2: the site needs one data URL.
- 3 → 4 and 3 → 5: shared tracks.
- 4 → 6: links into the locus plots.
- 7 depends only on 1 (the examples).
- 8 and 9 are independent.

The order above is also the commit order.

## 6. Steps that wait for your OK (outward-facing)

1. **Push** `restructure` to `main`. It is a fast-forward. Note that today's legacy Pages build serves `main`'s root,
   so pushing changes the live site even before step 2.
2. **Switch Pages** to the Actions workflow.
3. **Choose a public data host.** It must serve Range and CORS, and the private test bucket stays private:
   - **a public S3 bucket behind CloudFront in your AWS account** (HTTP/2 and compression, which shrinks the manifest
     to about 0.6 MB on the wire);
   - **a Cloudflare R2 public bucket** (no egress fees);
   - **a Hugging Face dataset repository** (free for public data; Range and CORS to be confirmed before choosing).

   GitHub Release assets are probably not an option: as far as I know their downloads send no CORS headers, so a
   browser can't read them cross-origin (to confirm).
4. **Upload** the release (34.8 GB plus extras) to the chosen host, and the version 2 manifest to the test bucket's
   `v4/`.

## 7. Not in this plan

- Deleting the ~870 GB of regenerable intermediates (your call; STATUS).
- The parquet-writer swap (planned in that repository; its slice 11 lands in `pipeline/`).
- LDlink's API-token system and server-side exports.

## 8. Review outcome and revisions (2026-11-23)

An independent read of this plan against the code found 11 problems. The revisions below replace the matching parts of
§4 and §6.

**Slices 1–2** (done, or folded in as they landed):
- The manifest `extras` key is dropped: examples are fetched whole with `fetch()`, which needs no size. Anything that
  copies the release (the §6 upload) must copy `examples/` explicitly; the manifest lists only Parquet files.
- `panels.parquet` joins the manifest builder's keyless small tables.
- `tools/serve.mjs` gains `--root <dir>`, `--no-release` and `--cors`. With `--cors` it answers `OPTIONS` (`Range`
  allowed) and sends `Timing-Allow-Origin`, so `bench/app_requests.js` (now matching the data URL, not `/release/`)
  measures cross-origin bytes. The site server mounts no `/release/`, so a broken data URL fails instead of quietly
  falling back to same-origin data.
- The build copies the four packages' whole directories (dynamic `import("hyparquet")` in `src/dbsnp.mjs`;
  `kernel.wasm` loads through `new URL`), keeps `app/` beside `src/`, leaves `check.html` out unless `--dev`, and
  CI installs with `npm ci`.

**Slice 2b (new, at the user's request): plain `.js`, not `.mjs`.** `package.json` already has `"type": "module"`.
Rename every `.mjs` in the repo to `.js` with `git mv` (`tools/move.py`, which rewrites relative import, `import()`
and `new URL()` specifiers), then fix the string-built paths that the move tool cannot see (`./adapters/${m}.mjs`,
`node --test` globs, npm scripts, shell scripts, the launch configs, doc paths). Gates: G1, G2 (planner, parity),
G3 (all ten), G4 (dev server and the built site).

**Slice 4 (locus plots), corrected against `LDcommon.py` / `LDassoc.py`:**
- LDproxy plots the full sorted list filtered only by R² > 0.01, including rows the text drops (RS equal to the
  query). The port returns the rows *before* printing. Each row carries the query's MAF too (LDproxy's hover shows
  "MAF (Query, Proxy)"). Points have alpha 0.5. Recombination is rate/100 on the main axis, with a right axis from
  −3 to 103.
- LDassoc's ramp runs `#FFCCCC` → `#FF0000` (101 steps), and every row at the index coordinate is blue, alpha 0.7.
  Association positions not in 1000G are drawn white with a gray outline, size 12.5, with their own hover, so the
  port returns them with their p-values. The y axis runs −0.03 → 1.03·max(y), and recombination is rate/100·max(y).
- The gene track gains LDlink's expanded mode (all transcripts) next to collapsed (`collapseTranscript`).

**Slice 6 (links):** routing changes first. The hash is `#<Module>` or `#<Module>?<query>`; `main.mjs` reads the
module before `?`, and a module button keeps the query only when it is that module's own.

**Slice 7 (the LDscore family), corrected:**
- The comparer calls `/LDlinkRest/`, where LDherit returns the text from "Total Observed scale" onward. With
  liability scale that marker is absent, so LDherit returns everything: the munge log, `\n--------\n` and the full ldsc
  log. LDcorrelation falls back to "Total Liability scale". There is no provenance line on this route. The exact
  text per route and scale is fixed from oracle fixtures before porting.
- Some lines carry the time ("Analysis finished at", "Total time elapsed", "Conversion finished at"). The
  adapters mask exactly those lines; everything else must match.
- None of the log text exists yet in `src/ldsc` (it returns numbers). The port includes ldsc's log lines, the
  pandas `describe()` / `corr()` layout of the LDscore summary, and the `to_string` rg table.
- LD score panels are read whole, not by chromosome: ldsc sorts by (CHR, BP) and then drops duplicate SNP IDs
  across the whole panel (afr has 613,138), and logs the count after that. A later release can store de-duplicated
  panels with their counts to make this cheaper.
- The REST routes validate inputs through `ldscore_compatibility.py` (`validate_sumstats_preanalysis`,
  `validate_bfile_compatibility`, `validate_ldscore_output`); these are ported, and they produce the bad-input
  errors. `sumstats_normalizer.py` is reached only through `/LDlinkRestWeb/validate_sumstats`; it is ported with an
  adapter on that route, or left out.
- The oracle mount uses `${LD_DATA:?}`, not a machine path. Read-only is fine: ldsc writes under `/data/tmp`.

**Slice 9 (hyparquet's pruning):** build the prototype on `parquetScan({ pruningFilter })`, which returns candidate
ranges and `readColumn(column, rowStart, rowEnd)`, so the payload column is still read only on its own pages.
`parquetRead({ filter })` is the wrong tool: it needs `usePageIndex: true` and reads every selected column over the
key's page spans. `$in` is O(rows × keys), and INT64 keys must be BigInt. Measure CPU as well as requests and bytes.

**§6, corrected order:** switch Pages to the workflow before pushing, or in the same step. Today's legacy build would
serve `app/index.html`, whose import map points at the untracked `node_modules/`. The data host's CORS rules must
allow the `Range` header.
