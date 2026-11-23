// Line length: code lines stay within 140 characters. test/long_lines.json is the baseline (lines over the limit per
// file, from before the rule); a file may not exceed its count, and a file not in the baseline has none.
// Refresh the baseline after rewrapping: node test/lines.test.js --update
import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync, readdirSync, writeFileSync } from "node:fs";
import { join, relative, extname } from "node:path";
import { REPO } from "../tools/config.js";

const LIMIT = 140, DIRS = ["app", "bench", "pipeline", "src", "test", "tools", "verify"];
const CODE = new Set([".js", ".sh", ".py", ".sql", ".html", ".css"]), BASE = join(REPO, "test/long_lines.json");
const SKIP = /(^|\/)(node_modules|dist|__pycache__)\//;
function counts() {
  const out = {};
  for (const d of DIRS) for (const f of readdirSync(join(REPO, d), { recursive: true })) {
    const rel = relative(REPO, join(REPO, d, f)).replaceAll("\\", "/");
    if (!CODE.has(extname(rel)) || SKIP.test(rel)) continue;
    const n = readFileSync(join(REPO, rel), "utf8").split(/\r?\n/).filter(l => l.length > LIMIT).length;
    if (n) out[rel] = n;
  }
  return out;
}

if (process.argv.includes("--update")) {
  const c = counts();
  writeFileSync(BASE, `${JSON.stringify(Object.fromEntries(Object.entries(c).sort()), null, 1)}\n`);
  console.log(`baseline: ${Object.values(c).reduce((a, b) => a + b, 0)} lines in ${Object.keys(c).length} files`);
} else test(`code lines within ${LIMIT} characters (baseline test/long_lines.json)`, () => {
  const base = JSON.parse(readFileSync(BASE, "utf8")), over = Object.entries(counts()).filter(([f, n]) => n > (base[f] ?? 0));
  assert.deepEqual(over.map(([f, n]) => `${f}: ${n} > ${base[f] ?? 0}`), []);
});
