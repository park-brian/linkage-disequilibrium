// LDassoc: LDlink's calculate_assoc called in the local backend container (ldlink-local/drive.py; the web route is
// LDassoc's only entry point) on the bundled prostate example, vs modules/ldassoc.mjs: the table text, error, warning.
// usage (from experiments/): ROOT=out/release_v4 node compare_ldassoc.mjs <out.json>
import { readFileSync, writeFileSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { asyncBufferFromFile, cachedAsyncBuffer } from "hyparquet";
import { loadKernel, KERNEL_WASM } from "../src/ld/kernel.mjs";
import { createContext } from "../src/ldlink/context.mjs";
import { ldassoc } from "../src/ldlink/ldassoc.mjs";

const [outPath] = process.argv.slice(2);
const LOCAL = "D:/Projects/LDlink/ldlink-local", ROOT = process.env.ROOT ?? "out/release_v4", PANELS = "D:/Projects/LDlink/db/data/1000G_population_samples";
const ctx = await createContext({ root: ROOT, open: async p => cachedAsyncBuffer(await asyncBufferFromFile(p)),
  readPanel: async p => readFileSync(`${PANELS}/${p}.txt`, "utf8").split(/\s+/).filter(Boolean), kernel: await loadKernel(readFileSync(KERNEL_WASM)) });
const FILE = { grch37: "prostate_example_grch37.txt", grch38: "prostate_example_grch38.txt", grch38_high_coverage: "prostate_example_grch38.txt" };
const text = Object.fromEntries(Object.entries(FILE).map(([b, f]) => [b, readFileSync(`${LOCAL}/data/ldassoc_example/${f}`, "utf8")]));
const cols = ["chr", "pos", "p"], cases = [];
const add = (label, c) => cases.push({ label, columns: cols, build: "grch37", pop: "EUR", ...c });
for (const build of ["grch37", "grch38", "grch38_high_coverage"]) {
  add(`variant rs7837688 ±100 kb EUR ${build}`, { build, region: "variant", origin: "rs7837688", window: 100000 });
  add(`gene MYC ±100 kb (lowest p origin) ${build}`, { build, region: "gene", name: "MYC", window: 100000 });
  add(`region chr8:128.0–128.7 Mb YRI ${build}`, { build, region: "region", start: "chr8:128000000", end: "chr8:128700000", pop: "YRI" });
}
add("variant rs6983267 ±500 kb ALL", { region: "variant", origin: "rs6983267", window: 500000, pop: "ALL" });
add("gene MYC with origin rs6983267", { region: "gene", name: "MYC", window: 100000, origin: "rs6983267" });
add("region with origin, CEU+YRI", { region: "region", start: "chr8:128400000", end: "chr8:128600000", origin: "rs7837688", pop: "CEU+YRI" });
add("region chrX", { region: "region", start: "chrX:2700000", end: "chrX:3500000" });
add("region chr22 (no association rows)", { region: "region", start: "chr22:30000000", end: "chr22:30100000" });
add("region > 3 Mb", { region: "region", start: "chr8:120000000", end: "chr8:124000000" });
add("region start >= end (string compare)", { region: "region", start: "chr8:9000000", end: "chr8:10000000" });
add("coordinate origin", { region: "variant", origin: "chr8:128539360", window: 100000 });
add("origin not in dbSNP", { region: "variant", origin: "rs9999999999", window: 100000 });
add("gene unknown", { region: "gene", name: "NOTAGENE", window: 100000 });
add("origin missing p-value", { region: "variant", origin: "rs3", window: 100000 });
add("bad population", { region: "gene", name: "MYC", window: 100000, pop: "FOO" });

writeFileSync(`${LOCAL}/data/tmp/ldassoc_cases.json`, JSON.stringify(cases.map(c => ({ module: "ldassoc", file: FILE[c.build], ...c }))));
execFileSync("docker", ["compose", "exec", "-T", "backend", "python3.13", "-", "/data/tmp/ldassoc_cases.json", "/data/tmp/ldassoc_cases.out.jsonl"],
  { cwd: LOCAL, input: readFileSync(`${LOCAL}/drive.py`), env: { ...process.env, MSYS_NO_PATHCONV: "1" }, maxBuffer: 1 << 28, stdio: ["pipe", "ignore", "ignore"] });
const got = readFileSync(`${LOCAL}/data/tmp/ldassoc_cases.out.jsonl`, "utf8").trim().split("\n").map(l => JSON.parse(l));
const report = { cases: [], totals: { cases: 0, identical: 0, differ: 0 } };
for (const [i, c] of cases.entries()) {
  const r = got[i], l = r.exception ? { error: "Internal server error. Please contact LDlink admin." }
    : { ...(r.text != null ? { text: r.text } : {}), ...(r.json.error !== undefined ? { error: r.json.error } : {}), ...(r.json.warning !== undefined ? { warning: r.json.warning } : {}) };
  const o = await ldassoc(ctx, { ...c, text: text[c.build] });
  const norm = x => JSON.stringify({ text: x.text ?? null, error: x.error ?? null, warning: x.error ? null : x.warning ?? null });
  const same = norm(l) === norm(o);
  let detail = {};
  if (!same && l.text && o.text) {
    const a = l.text.split("\n"), b = o.text.split("\n"), d = [];
    for (let k = 0; k < Math.max(a.length, b.length); k++) if (a[k] !== b[k]) d.push({ line: k, ldlink: a[k], ours: b[k] });
    detail = { lines: [a.length, b.length], differing_lines: d.length, examples: d.slice(0, 5), warning: [l.warning, o.warning] };
  } else if (!same) detail = { ldlink: { ...l, text: l.text?.slice(0, 200) }, ours: { ...o, text: o.text?.slice(0, 200) }, exception: r.exception };
  report.totals.cases++; report.totals[same ? "identical" : "differ"]++;
  report.cases.push({ name: c.label, verdict: same ? "identical" : "differ", lines: l.text ? l.text.split("\n").length : 0, ...detail });
  console.log(JSON.stringify(report.cases.at(-1)).slice(0, 1500));
}
writeFileSync(outPath, JSON.stringify(report, null, 2));
console.log(JSON.stringify(report.totals));
