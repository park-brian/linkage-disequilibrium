// Shared plumbing for the LDlink comparisons: a Node context over the release, LDlink's REST API and its in-container
// Python functions (tools/ldlink-local/drive.py), and output normalization. Paths come from tools/config.js.
import { readFileSync, writeFileSync, existsSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { asyncBufferFromFile, asyncBufferFromUrl, cachedAsyncBuffer } from "hyparquet";
import { cached, coalesce } from "../../src/core/parquet.js";
import { loadKernel, KERNEL_WASM } from "../../src/core/ld/kernel.js";
import { createContext } from "../../src/parity/ldlink/context.js";
import * as CONFIG from "../../tools/config.js";

export { CONFIG };
/** the comparisons, in run order (compare.js, snapshot.js) */
export const MODULES = ["ldmatrix", "ldpair", "ldpop", "ldhap", "ldproxy", "snpclip", "snpchip", "ldtrait", "ldexpress", "ldassoc", "ldherit",
  "ldcorrelation", "ldscore"];
export const BASE = CONFIG.LDLINK_URL;

/**
 * the module context the ports run on, with the WASM kernel: local release files (panels included), or with
 * http = a release URL, the app's read path (Range requests with sizes from the manifest, coalesced, cached)
 */
export async function nodeContext(root = process.env.ROOT || CONFIG.RELEASE, { http = process.env.LD_HTTP } = {}) {
  const kernel = await loadKernel(readFileSync(KERNEL_WASM));
  if (!http) return createContext({ root, open: async p => cachedAsyncBuffer(await asyncBufferFromFile(p)), kernel });
  const url = http.replace(/\/$/, ""), res = await fetch(`${url}/manifest.json`);
  if (!res.ok) throw new Error(`no release at ${url}: manifest.json returned HTTP ${res.status}`);
  return createContext({ root: url, manifest: await res.json(), kernel,
    open: async (u, byteLength) => cached(coalesce(await asyncBufferFromUrl({ url: u, byteLength }))) });
}

/** an example input from the release's examples/ (local file or over HTTP): text, or bytes when binary; null if absent */
const examples = new Map();
export function example(ctx, name, { binary = false } = {}) {
  const key = `${ctx.root}|${name}|${binary}`;
  if (!examples.has(key)) examples.set(key, (async () => {
    const at = `${ctx.root}/examples/${name}`;
    if (/^https?:/.test(ctx.root)) { const r = await fetch(at); return r.ok ? (binary ? new Uint8Array(await r.arrayBuffer()) : r.text()) : null; }
    return existsSync(at) ? (binary ? new Uint8Array(readFileSync(at)) : readFileSync(at, "utf8")) : null;
  })());
  return examples.get(key);
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
    { cwd: CONFIG.LDLINK_HOME, input: readFileSync(CONFIG.repo("tools/ldlink-local/drive.py")), env: { ...process.env, MSYS_NO_PATHCONV: "1", LD_DATA: CONFIG.DATA } /* the kit's compose file mounts ${LD_DATA}/ldscore_src */, maxBuffer: 1 << 28, stdio: ["pipe", "ignore", "ignore"] });
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
