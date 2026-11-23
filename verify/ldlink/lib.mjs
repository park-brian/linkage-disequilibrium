// Shared plumbing for the LDlink comparisons: a Node context over the release, LDlink's REST API and its in-container
// Python functions (tools/ldlink-local/drive.py), and output normalization. Paths come from tools/config.mjs.
import { readFileSync, writeFileSync, existsSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { asyncBufferFromFile, cachedAsyncBuffer } from "hyparquet";
import { loadKernel, KERNEL_WASM } from "../../src/ld/kernel.mjs";
import { createContext } from "../../src/ldlink/context.mjs";
import * as CONFIG from "../../tools/config.mjs";

export { CONFIG };
export const BASE = CONFIG.LDLINK_URL;

/** the module context the ports run on: local release files, panels from config, the WASM kernel */
export async function nodeContext(root = process.env.ROOT || CONFIG.RELEASE) {
  return createContext({ root, open: async p => cachedAsyncBuffer(await asyncBufferFromFile(p)),
    readPanel: async p => readFileSync(`${CONFIG.PANELS}/${p}.txt`, "utf8").split(/\s+/).filter(Boolean), kernel: await loadKernel(readFileSync(KERNEL_WASM)) });
}

/** the LDmatrix cases every comparison draws inputs from (verify/ldlink/cases/) */
export const matrixCases = () => JSON.parse(readFileSync(CONFIG.repo("verify/ldlink/cases/ldmatrix_cases.json"), "utf8"));
export const clean = s => s.replace(/\r$/, "");

// at most LDLINK_CONCURRENCY requests in flight (default 4): adapters fan out with Promise.all, and a long LDlink queue
// makes the later requests outlive fetch's 300 s header timeout (LDexpress over all 54 tissues)
const LIMIT = +(process.env.LDLINK_CONCURRENCY || 4), waiting = [];
let active = 0;
async function limited(f) {
  if (active >= LIMIT) await new Promise(r => waiting.push(r)); else active++;
  try { return await f(); } finally { const next = waiting.shift(); if (next) next(); else active--; }
}
export const get = (path, query) => limited(async () => (await fetch(`${BASE}/LDlinkRest/${path}?${query}`)).text());
export const post = (path, body) => limited(async () => (await fetch(`${BASE}/LDlinkRest/${path}`,
  { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) })).text());
/** LDlink's API answer: a JSON error, or the text body */
export const asText = body => (body.trimStart().startsWith("{") ? { error: JSON.parse(body).error } : { text: body });

/** run LDlink's Python functions in the backend container: drive.py over a case file → one result per case */
export function drive(name, cases) {
  const tmp = CONFIG.LDLINK_TMP, inFile = `${name}_cases.json`, outFile = `${name}_cases.out.jsonl`;
  writeFileSync(`${tmp}/${inFile}`, JSON.stringify(cases));
  execFileSync("docker", ["compose", "exec", "-T", "backend", "python3.13", "-", `/data/tmp/${inFile}`, `/data/tmp/${outFile}`],
    { cwd: CONFIG.LDLINK_HOME, input: readFileSync(CONFIG.repo("tools/ldlink-local/drive.py")), env: { ...process.env, MSYS_NO_PATHCONV: "1" }, maxBuffer: 1 << 28, stdio: ["pipe", "ignore", "ignore"] });
  return readFileSync(`${tmp}/${outFile}`, "utf8").trim().split("\n").map(l => JSON.parse(l));
}
/** LDproxy's JSON side file for a request reference (its warning is not in the API text) */
export const proxyJson = ref => { const p = `${CONFIG.LDLINK_TMP}/proxy${ref}.json`; try { return existsSync(p) ? JSON.parse(readFileSync(p, "utf8")) : {}; } catch { return {}; } };

/** json.dumps(sort_keys=True): key order for comparing JSON */
export const canon = o => (Array.isArray(o) ? o.map(canon) : o && typeof o === "object" ? Object.fromEntries(Object.keys(o).sort().map(k => [k, canon(o[k])])) : o);
/** where two texts differ, line by line (report detail) */
export function textDiff(a, b) {
  const x = a.split("\n"), y = b.split("\n"), d = [];
  for (let i = 0; i < Math.max(x.length, y.length); i++) if (x[i] !== y[i]) d.push({ line: i, ldlink: x[i], ours: y[i] });
  return { lines: [x.length, y.length], differing_lines: d.length, examples: d.slice(0, 5) };
}
