# Plan: the cleanup (from CLEANUP_REVIEW.md)

**Written:** 2026-11-23, on `restructure` (= `main` plus local commits). Input: [CLEANUP_REVIEW.md](CLEANUP_REVIEW.md). Each slice ends
green on its gates and is one commit. When done, this file and REVIEW.md move to `docs/history/` with an outcome.

## Gates

| gate | what | cost |
|---|---|---|
| G1 | `npm test` | seconds |
| G2 | `npm run verify` (data) | ~20 min |
| G3s | **new** `npm run verify:snapshot`: every G3 case through our ports, hashed and compared with the committed hashes | ~5 min, no Docker |
| G3 | `npm run verify:ldlink` against the local LDlink | ~1.5 h |
| G4 | `npm run verify:app`: every view in headless Chromium (and a look at the built site) | ~1 min |

G3s is valid because its hashes are taken from code that G3 has just passed. A slice that changes `src/` runs G1 and
G3s; G2 and G3 run in full at slices 2, 6, 10 and the end. A G3s difference is investigated with G3 on that module.

## Slices

1. **Gates that can fail.**
   - Every verifier sets `process.exitCode = 1` on a difference; `verify/run.sh` exits non-zero if any step did.
   - Reports go to `$LD_DATA/verify/`, not into the release.
   - `verify_select` reads the release, not a `phaseB` build copy.
   - The generated reports already in `release_v4` (and published under `v5/`) are removed: locally now; in the bucket
     only with your OK, since that deletes public objects.
2. **The parity snapshot (G3s).**
   - `verify/ldlink/snapshot.js`: runs each adapter's `cases()` and `ours()`, hashes each output (the LDscore family
     with its time lines masked), and compares with `verify/ldlink/snapshot.json` (`--update` rewrites it).
   - `nodeContext({ http })` reads the release over HTTP with `coalesce`, so G3s can also run against the dev server
     or CloudFront: that covers Range, CORS and coalescing from Node.
   - `app/check.html` and `check.js` are deleted (G3s over HTTP supersedes them).
   - Taken right after a full G3, which must pass first.
3. **Deletions** (about −6,000 lines, ~900 of them code):
   - `bench/results/` → `$LD_DATA/bench/` (benches write there from now on);
   - `tools/move.py`, `tools/moves_slice1.tsv`;
   - `pipeline/profile/`;
   - `bench/net.js`, `pipeline/release/presign.sh`, `bench/pruning.js` (its numbers stay in decisions.md);
   - `pipeline/dbsnp/bson2tsv.js` (use `bson_stream.js`), the b157 fetch in `ncbi/fetch.sh`, the non-final page-size
     sweep in `genotypes/build.js`;
   - `ncbi/compare.sql` and `ncbi/allele_join.js` → `docs/history/` (they are analyses, not build steps);
   - dead exports (REVIEW §4 of the `src/` read) and stale comments;
   - `hyparquet-writer` becomes a devDependency.
4. **Bug fixes** (REVIEW §3):
   - kernel use serialized per context (no shared-arena reuse);
   - failed reads evicted from every cache;
   - absent files reported, not destructured;
   - a run token in the app (stale results ignored), `res.ok` on example fetches, an option-list failure does not
     block a linked run, no stack traces shown;
   - pointer coordinates scaled for CSS-scaled canvases;
   - LDassoc's per-position overlap scans replaced by two pointers (likewise LDtrait, the planner);
   - `build_site.js` filters `check.js`; the Pages workflow runs `LD_SKIP_DATA=1 npm test` and requires `LD_DATA_URL`;
   - `serve.js`: bad `%` escapes and stream errors answer 400 / close instead of crashing;
   - `aria-live` on `#status` only; nav buttons become links (consistent history).
5. **Port helpers** (about −180 lines): resolve-and-classify a list; genotype calls and allele counts; the annotation
   pass; validation (build, population, `r2_d`, window, list limits); warning accumulation; Python number parsers;
   the sorted-allele 2×2 and correlated-allele test. Gate: G1, G3s.
6. **Core helpers** (about −120 lines): one page-range planner under `sortedLookup` and `rangeScan`; `sortedLookup`
   maps every column to its own pages (the over-read REVIEW §3.8 found); `planner.query` split into resolve /
   normalize / plan / select; one `firstBy`, `once`, `groupBy`, overlap test. Gates: G1, G3s, then full G2 + G3
   (the read path changed); before/after requests and bytes in STATUS.
7. **App de-duplication** (about −80 lines): one `niceStep`; LDmatrix on the shared tracks; one canvas + tooltip
   helper; one locus show function; number fields coerced in `value()`; one field per line in the module table.
   Gate: G4 (canvas hashes of the pyramid unchanged).
8. **Pipeline runbook.** `pipeline/build.sh` runs every stage in order with resume checks; the missing producers
   (`dbsnp_src/dup_ids`, the NCBI Parquet step, the copy into `release/genotypes`) become scripts; every `usage:`
   line names its real path; `tools/check_docs.js` also checks `usage:` lines and brace globs; the 178 unrolled
   `COPY` statements are generated. Gate: one chromosome through `build.sh` into a scratch release, then
   `pipeline/release/check.js` on it.
9. **Line length.** A check in `npm test` (140 characters, code files), then every over-long line rewrapped.
   Gate: G1, G3s.
10. **The parity layer** (REVIEW §4), in three commits:
    1. moves: a core directory and a parity directory (input parsing, resolution choices, selection, messages,
       formatting emulation, the ports); `nextUp` into the kernel first;
    2. the core returns every match in file order and the parity layer chooses (first / last / warnings);
    3. a clean API on the core (`ld.pair`, `ld.proxy`, `ld.matrix`, `annotate`) with tests.
    Gates: G1 and G3s after each; full G2 + G3 after the third.
11. **Final:** G1–G4 in full, STATUS, the plan and review to history, one commit.

## Waits for your OK

- Deleting the published reports under `ldlink-parquet/v5/` (slice 1).
- Pushing to `main` (it redeploys the site), at the end.

## Review outcome (2026-11-23)

An independent read of this plan against the code. Changes:
- **Slice 1:** only `differ` fails G3; `tie_order`, `ldlink_exception`, `ldlink_error` and `degenerate_z` are
  documented outcomes.
- **Slice 2:** the snapshot skips the LDlink side (`ldlink()`). The example inputs go through `lib.example()`, which
  reads from disk or over HTTP, so an HTTP run doesn't turn every LDscore-family case into "missing example". The
  time mask covers LDherit, LDcorrelation and LDscore. The snapshot runs sequentially until slice 4 fixes the arena.
  Case inputs (e.g. SNPclip's dense IDs) still come from the local release.
- **Slice 3:** `bench/pruning.js` stays until slice 6 has used it. With `check.html` gone, `build_site.js`'s
  check filter is removed rather than fixed.
- **Slice 4:** the arena fix is `ctx.withPrepared(build, rows, pop, use)`. It awaits the samples and panel, then
  prepares, runs the synchronous `use`, and resets in `finally`. The ports' own `reset()` calls go. No lock is needed,
  and SNPclip's early return is covered.
- **Slice 6:** per-column pages must keep file order (rows sorted by row number), every duplicate row, `String`
  key matching, the requested columns and BigInt keys. Fetch all offset indexes together, and measure wall time and
  requests over HTTP as well as bytes. `planner.query` is cut along the parity boundary (resolve / select apart).
- **Slice 9** (line length) runs after slice 10, which rewrites many of the same lines. Its check lands with a
  baseline now.
- **Slice 10:** step 1 moves only files that are cleanly core or cleanly parity. The mixed ones (`dbsnp.js`,
  `genotypes.js`, `annotations.js`, `ldsc/`) split in step 2.
- **Added:** document "the release's file order = LDlink's dump order" in DATA.md, and make LDassoc's 17 MB default
  example an explicit choice in the form rather than the default fetch.

## Outcome (2026-11-23)

All eleven slices are committed on `restructure` (local, not pushed), one commit per slice except where noted. Final gates:
G1 37/37 tests; G3s 1,440/1,440; G4 `verify:app` 14/14 with every figure hash unchanged, and the built site ran
LDproxy against the local release; G2 and G3 in full (STATUS.md has the counts).

| slice | commit | what landed | differences from the plan |
|---|---|---|---|
| 1 | `1d1e2b2` | every verifier fails on a difference; reports in `$LD_DATA/verify/` | `run.sh` crashed instead of reporting a failed step (`set -u` on an unset `FAILED`); found and fixed in slice 6 |
| 2–3 | `9577131`, `162024f` | G3s (`verify:snapshot`, also over HTTP); `check.html` gone; ~6,000 lines deleted; dead exports | slices 2 and 3 share one commit (they edited the same docs); `pipeline/dbsnp/bson2tsv.js` stays (it is the live dump step) |
| 4 | `b399bee` | the bug fixes of REVIEW §3 | — |
| 5 | `648ce4d` | `check`, `screen`, `resolveOne` … in `parity/ldlink/common.js` | net +28 lines, not −180: the helpers' doc comments cost what the copies saved |
| 6 | `526cd36` | one read path (`matchRows`), per-column pages, planner split, `collect.js` | ranges keep one all-column read (per-column reads cost 3 → 11 reads); single-column reads use decoded arrays (the first version ran the parity check out of memory) |
| 7 | `35edef0` | shared canvas / tooltip / locus helpers; figure hashes in `verify:app` | net −1 line, not −80; LDmatrix keeps its own gene track (the shared one would change the figure) |
| 8 | `7681ebc` | `pipeline/build.sh`, the three missing producers, `@for` blocks, usage-line check | the gate rebuilt chr22: GRCh37 byte-identical; GRCh38 identical rows plus an all-null `vcf_id` column (release_v4's GRCh38 predates it) |
| 10 | `9cac7ef`, `25c749b`, `ef66a83` | `src/core` and `src/parity`; core returns every match, parity chooses; `core/ld.js` | steps 2 and 3 share a commit; `ef66a83` adds pipeline import paths missed from `25c749b`; annotations moved to parity whole (their reads are LDlink's queries) rather than split |
| 9 | `73d8708` | 700 → 63 over-long lines | the 63 left are LDlink's message strings and data lists, kept verbatim |
| also | `3301a7c`, `07fe70b` | `tools/coverage.js`; `verify/app/run.js` (Playwright) | added on request during the cleanup |

Not done: the published verifier reports under `ldlink-parquet/v5/` are still in the bucket (deleting public objects
waits for an OK).
