// Line coverage of repo files from raw V8 coverage (NODE_V8_COVERAGE=<dir> node <script>, or Chromium's
// page.coverage / CDP Profiler.takePreciseCoverage, which use the same format). Per file, the block ranges are
// painted outer to inner (an inner range's count overrides its enclosing one), and a line counts as missed when
// none of its non-blank characters ran. Prints one line per file, the total, and with --missed the missed lines.
// usage: node tools/coverage.js [--include src/] [--missed] (<coverage dir> | -- <node script> [args...])
//   -- runs the script with NODE_V8_COVERAGE set to a fresh temporary directory (child processes included), then reports
//   npm run coverage            the unit tests (node --test)      npm run coverage:snapshot   the parity snapshot (G3s)
import { readFileSync, readdirSync, mkdtempSync, rmSync } from "node:fs";
import { spawnSync } from "node:child_process";
import { tmpdir } from "node:os";
import { join, relative } from "node:path";
import { fileURLToPath } from "node:url";
import { REPO } from "./config.js";

const argv = process.argv.slice(2), sep = argv.indexOf("--");
const args = sep < 0 ? argv : argv.slice(0, sep), run = sep < 0 ? null : argv.slice(sep + 1);
const missedToo = args.includes("--missed"), include = args.includes("--include") ? args[args.indexOf("--include") + 1] : "src/";
const dir = run ? mkdtempSync(join(tmpdir(), "ldcov-")) : args.find((a, i) => !a.startsWith("--") && args[i - 1] !== "--include");
const USAGE = "usage: node tools/coverage.js [--include src/] [--missed] (<coverage dir> | -- <node script> [args...])";
if (!dir || (run && !run.length)) throw new Error(USAGE);
if (run) {
  const r = spawnSync(process.execPath, run, { stdio: "inherit", env: { ...process.env, NODE_V8_COVERAGE: dir } });
  if (r.status !== 0) process.exitCode = r.status ?? 1; // the run's failure is reported; coverage is still printed
}

// every process's report, per file: one list of ranges per script instance (painted apart, then OR-ed)
const scripts = new Map();
for (const f of readdirSync(dir).filter(f => f.endsWith(".json"))) {
  for (const s of JSON.parse(readFileSync(join(dir, f), "utf8")).result ?? []) {
    if (!s.url.startsWith("file:")) continue;
    const rel = relative(REPO, fileURLToPath(s.url)).replaceAll("\\", "/");
    if (!rel.startsWith(include) || rel.includes("node_modules/")) continue;
    if (!scripts.has(rel)) scripts.set(rel, []);
    scripts.get(rel).push(s.functions.flatMap(fn => fn.ranges));
  }
}

let hit = 0, all = 0;
for (const [rel, runsOf] of [...scripts].sort()) {
  const src = readFileSync(join(REPO, rel), "utf8"), ran = new Uint8Array(src.length), seen = new Uint8Array(src.length);
  for (const rs of runsOf) { // widest first, so a nested block's count overrides its parent's
    const here = new Int8Array(src.length).fill(-1);
    for (const r of [...rs].sort((a, b) => a.startOffset - b.startOffset || b.endOffset - a.endOffset))
      here.fill(r.count > 0 ? 1 : 0, r.startOffset, Math.min(r.endOffset, src.length));
    for (let k = 0; k < src.length; k++) if (here[k] >= 0) { seen[k] = 1; if (here[k]) ran[k] = 1; }
  }
  let off = 0, covered = 0, counted = 0;
  const missed = [];
  src.split("\n").forEach((line, i) => {
    let code = false, hit1 = false;
    for (let k = 0; k < line.length; k++) if (!/\s/.test(line[k])) { code = true; if (!seen[off + k] || ran[off + k]) hit1 = true; }
    if (code) { counted++; if (hit1) covered++; else missed.push(i + 1); }
    off += line.length + 1;
  });
  hit += covered; all += counted;
  const spans = missed.reduce((a, n) => { const l = a.at(-1); if (l && n === l[1] + 1) l[1] = n; else a.push([n, n]); return a; }, [])
    .map(([a, b]) => (a === b ? `${a}` : `${a}-${b}`));
  const tail = missedToo && spans.length ? `  missed: ${spans.join(" ")}` : "";
  console.log(`${(100 * covered / counted).toFixed(1).padStart(5)}%  ${rel} (${covered}/${counted})${tail}`);
}
if (all) console.log(`${(100 * hit / all).toFixed(1).padStart(5)}%  total (${hit}/${all} lines with code)`);
else { console.log(`no coverage for files under ${include}`); process.exitCode ||= 1; }
if (run) rmSync(dir, { recursive: true, force: true });
