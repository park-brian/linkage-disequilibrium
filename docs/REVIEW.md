# Code review (2026-11-23)

A full read of the repository after Phase C, by two independent reviewers (`src/`; everything else), with the
claims below checked against the code before they were written here. Nothing in this file has been fixed yet; it is
the input to the next plan.

## 1. Cost

**Now:** about **$3 a month**, almost all of it S3 storage (Standard, $0.023/GB-month):
- the live release, `ldlink-parquet/v5/`: 34.9 GB, $0.80;
- the older prefixes `v2/`–`v4/`: 94 GB, $2.16 (deleting them leaves about $0.85);
- `ldlink-source/ldscore_efs.tar`: 1.8 GB, $0.04.

GitHub Pages is free. CloudFront's always-free tier covers 1 TB of transfer and 10 million requests a month. Beyond
that it is $0.085/GB in price class 100 (first 10 TB) and $0.01 per 10,000 HTTPS requests. Cache misses also cost
S3 GETs ($0.0004 per 1,000); transfer from S3 to CloudFront is free.

**Per visit** (first visit, measured; the browser then caches for a day): the manifest is 0.6 MB compressed.

| view | transfer |
|---|---|
| LDpair, LDexpress, LDscore | 1–2.5 MB |
| LDmatrix, Pyramid, LDtrait, SNPchip | 2.7–4.7 MB |
| LDproxy | 12.4 MB in 484 requests |
| LDassoc (default example) | 20.6 MB, 17 MB of it the example file |
| LDherit / LDcorrelation | 25 MB with the EUR panel; up to 200 MB with AFR |

A typical session of 10–20 MB means the free tier carries **roughly 50,000–100,000 sessions a month**. Past that:

| sessions / month (15 MB each) | transfer | CloudFront, approx. |
|---|---|---|
| 100,000 | 1.5 TB | $40 |
| 1,000,000 | 15 TB | $1,200 |

Requests are cheap by comparison: an LDproxy run's 484 requests cost $0.0005 beyond the free tier.

The levers, in order:
1. Cut LDproxy's over-reads (§3, `sortedLookup`).
2. Drop LDassoc's 17 MB default fetch.
3. Read LD score panels by column chunks for chromosome-limited uploads.

## 2. How the code is

**Sound.** Comments cite their LDlink sources, the WASM and JS kernels are cleanly paired, the release layout and
planner are well reasoned, the module ports are isolated, and G1 is fast and meaningful.

**The weaknesses are structural, not stylistic:**
- **No gate can fail.** `verify/ldlink/compare.js` exits non-zero only on a usage error, never on a difference.
  The data verifiers never set an exit code, and `verify/run.sh` always exits 0. "Green" has meant "read the
  summaries".
- **Reports are written into the release.** `verify/data/planner.js`, `parity.js` and `annotations.js` write to
  `${root}/…`, so `parity_report.json` and siblings went public with the v5 upload. They hold test results, not
  secrets.
- **Dense lines.** About 700 lines across the repo exceed AGENTS.md's 140-character limit; about 260 of them are in
  `src/` and 143 in `app/`.
- **Duplication** across the 13 ports, about 180 lines, and in the core, about 120.
- **The pipeline isn't reproducible from documented commands.** Some producers are missing (`dbsnp_src/dup_ids`,
  the NCBI Parquet step, the copy into `release/genotypes`). Many `usage:` lines still name pre-restructure scripts.
- **About 6,000 lines are leftovers:** bench results (including byte-identical copies), one-off restructure tools,
  live-Mongo profiling, private-bucket benches.

## 3. Bugs and risks, ranked

1. **Shared WASM arena.**
   - `ctx.prepare` allocates, and every port calls `kernel.reset()` when it finishes.
   - A second run started while the first awaits can reset and reuse the first run's buffers in the one-microtask
     gap. The result would be wrong, silently.
   - Reachable from the UI: Run disables only its own form.
   - Fix: serialize kernel use per context.
2. **Rejected reads are cached forever.** In `source.js` and every memo in `annotations.js`, `dbsnp.js` and
   `context.js`, one transient network error poisons a file for the session.
3. **Stale results.** Switching modules mid-run renders the old result into the new module (`app/forms.js`). Fix
   with a run token.
4. **Missing files crash instead of reporting.**
   - `source.file` returns `null` for an absent file, and these callers destructure it anyway:
     `genotypes.js`, `annotations.js`, `ldsc/panels.js`, `planner.js`.
   - The app's example fetches never check `res.ok`.
5. **Quadratic loops on uploads.** In `ldassoc.js`, `at(p)` filters every window row per association position, so
   large uploads grow quadratically. Smaller cases are in `ldtrait.js` and `planner.js`.
6. **Hover and click are offset on narrow screens.** The canvases are CSS-scaled (`max-width: 100%`), but pointer
   math uses unscaled pixels.
7. **The site build.**
   - `tools/build_site.js` still filters `check\.(html|mjs)`, so `app/check.js` ships to the live site.
   - The Pages workflow deploys without running tests.
   - An unset `LD_DATA_URL` silently deploys a site pointed at `../release`.
8. **`sortedLookup` over-reads.** Without a payload column it reads every requested column over the key pages' whole
   row spans. In the comparison with hyparquet's version, 512 rsID lookups read 38 MB where 8.6 MB would do.

## 4. A parity layer

Yes, and it is the right long-term shape. Today LDlink's rules sit inside the core:
- **Resolution:** `dbsnp.js` is `find_one` = first document in dump order, with LDlink's warning text, including a
  bug kept verbatim.
- **Parsing:** `genome.js` parses coordinates as Python's `strip("chr")` and `int()` do.
- **Row selection:** `genotypes.pickParity` picks the last overlapping row, and the planner owns the
  `select: "parity"` mode.
- **Annotations:** `annotations.js` mixes `find_one` and `$lookup` semantics.
- **Release columns:** some (`vcf_id`, `rsid_source`, `rsid_pos_delta`) exist only for LDlink's labels.

**Target:**
- **A core directory (proposed `src/core`)** returns *every* match in file order: variants by ID or position, genotype rows overlapping a
  position, annotations, the LD kernel, the generic numerics.
- **A parity directory (proposed `src/parity`)** chooses from those lists (first, last, "Multiple rsIDs") and holds everything LDlink-shaped:
  - input parsing and Python `int()` / `float()` acceptors;
  - messages;
  - `pyfmt`, and the ldsc `log` and `pandas` emulation;
  - the 13 ports.
- **The contract:** the release's file order equals LDlink's dump order. That is a property of the data,
  documented in DATA.md, not code.
- **A clean API on the core:** `ld.pair`, `ld.proxy`, `ld.matrix` and `annotate`. These return numbers and typed
  errors, default to allele-aware selection, and leave rendering to the caller.

**Cost:** about a week, in three steps, with `npm run verify:ldlink` after each:
1. File moves (a day).
2. The core returns lists and parity wraps them (2–3 days; about 200 lines change).
3. The clean API (new code, about 300 lines).

It breaks imports in `verify/`, `app/`, `bench/` and the manifest builder. `nextUp` must move into the kernel first,
because workers load the kernel chain by path.

## 5. Making it as simple as possible

Ranked by value for effort:
1. **Gates that can fail.**
   - Exit non-zero on any difference.
   - Write reports to `$LD_DATA/verify/`, not the release.
   - Point `verify_select` at the release, not a build copy.
2. **A cheap parity snapshot.** Hash each G3 case's output and check the hashes in `npm test`, over HTTP with
   `coalesce`. Most changes then skip the 1.5-hour oracle run. The same mode replaces `check.html`.
3. **Delete about 6,000 lines.**
   - Committed bench results, which should live in `$LD_DATA/bench/`.
   - `tools/move.py` and its table.
   - `pipeline/profile/`.
   - `bench/net.js`, `presign.sh` and `bench/pruning.js` (its numbers are in `decisions.md`).
   - `dbsnp/bson2tsv.js`, the b157 fetch, and the page-size sweep in `genotypes/build.js`.
4. **Fix the bugs in §3 (1–7).**
5. **Shared port helpers** (about −180 lines):
   - resolve-and-classify a variant list;
   - genotype calls and allele counts;
   - the annotation pass;
   - validation (build, population, `r2_d`, window, list limits);
   - warning accumulation;
   - Python number parsers.
6. **Core helpers** (about −120 lines):
   - one page-range planner shared by `sortedLookup` and `rangeScan`;
   - split `planner.query`;
   - one `firstBy`, `once`, `groupBy` and overlap test.
7. **App de-duplication** (about −80 lines):
   - one `niceStep`;
   - LDmatrix on the shared tracks;
   - one canvas and tooltip helper;
   - merge `showProxy` / `showAssoc`;
   - number fields coerced in one place.
8. **Pipeline runbook.** One ordered pipeline build script (proposed: pipeline/build.sh) with resume checks, the missing producers, and corrected
   `usage:` lines. Generate the 178 hand-unrolled `COPY` statements.
9. **A line-length check in `npm test`**, then reflow the over-long lines.
10. **The parity layer (§4)**, on top of the helpers from 5 and 6.
