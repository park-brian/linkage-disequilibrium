// Does the planner's rsID strategy still pay? Each case runs a module port cold (a fresh context: empty caches) over
// HTTP with the app's read path (verify/ldlink/lib.js nodeContext: manifest sizes, coalesced Range requests, cached),
// once per planner mode, counting every request and byte and the wall time:
//   auto    the planner as is (rsID lists > 32 resolve "hybrid" when the manifest prices the label scan lower)
//   direct  forceRsid: "direct" (every rsID looked up in the rsID index)
// Slice 1 of docs/history/PLAN_2026-11-24.md also measured the genotype-row fetch the planner no longer does
// (decisions.md: "The planner fetches no genotype rows for the ports").
// The outputs of the modes must be identical (checked by hash). The manifest fetch is counted apart.
// usage: node bench/planner_modes.js <release URL> [repeat=1] [case...]   → $LD_DATA/bench/results/planner_modes_<host>.json
import { createHash } from "node:crypto";
import { mkdirSync, writeFileSync } from "node:fs";
import { nodeContext, example, canon } from "../verify/ldlink/lib.js";
import { ldmatrix } from "../src/parity/ldlink/ldmatrix.js";
import { ldproxy } from "../src/parity/ldlink/ldproxy.js";
import { snpclip } from "../src/parity/ldlink/snpclip.js";
import { ldtrait } from "../src/parity/ldlink/ldtrait.js";
import { ldhap } from "../src/parity/ldlink/ldhap.js";
import { ldassoc } from "../src/parity/ldlink/ldassoc.js";
import { snpchip } from "../src/parity/ldlink/snpchip.js";
import { useDataCwd } from "../tools/config.js";

const [url, rep = "1", ...pick] = process.argv.slice(2);
if (!url) throw new Error("usage: node bench/planner_modes.js <release URL> [repeat] [case...]");
const MODES = { auto: {}, direct: { forceRsid: "direct" } };

// every fetch counted: requests and bytes (Range responses carry Content-Length)
const stats = { requests: 0, bytes: 0 }, realFetch = globalThis.fetch;
// at most 64 in flight, as a browser queues them (CloudFront refuses HTTP/2 streams beyond its limit)
let inFlight = 0; const queue = [];
globalThis.fetch = async (...a) => {
  if (inFlight >= 64) await new Promise(ok => queue.push(ok));
  inFlight++;
  try {
    const r = await realFetch(...a), body = await r.arrayBuffer(); // read here, so the slot frees when the bytes are in
    stats.requests++; stats.bytes += body.byteLength;
    return new Response(body, { status: r.status, statusText: r.statusText, headers: r.headers });
  } finally { inFlight--; queue.shift()?.(); }
};

/** a context whose planner queries run in the given mode */
async function contextIn(mode) {
  const ctx = await nodeContext(undefined, { http: url }), q = ctx.planner.query;
  ctx.planner.query = o => q({ ...o, ...MODES[mode] });
  return ctx;
}

// inputs: LDlink's example SNPs, and real rsIDs from chr22 (labels of 1000G rows), drawn evenly
const EX = ["rs13426489", "rs4442975", "rs13387042", "rs13412666", "rs6721996", "rs12465506", "rs12613955", "rs12620095"];
const probe = await nodeContext(undefined, { http: url });
const chr22 = (await probe.windowRows("grch37", "22", 20_000_000, 40_000_000, ["common"])).filter(r => r.rsid_source === 1);
const draw = (n, lo, hi) => { const w = chr22.filter(r => r.pos >= lo && r.pos <= hi), step = w.length / n;
  return Array.from({ length: n }, (_, k) => `rs${w[Math.floor(k * step)].rsid}`); };
const L300 = draw(300, 30_000_000, 31_000_000), L1000 = draw(1000, 25_000_000, 35_000_000), L5000 = draw(5000, 20_000_000, 40_000_000);
const assoc = await example(probe, "prostate_example_grch37.txt");
// sparse: 200 rsIDs spread over all of chr22; genome-wide: 250 rsIDs from each of four chromosomes
const S200 = draw(200, 16_000_000, 51_000_000);
const wide = [];
for (const [c, lo] of [["1", 100_000_000], ["6", 30_000_000], ["12", 50_000_000], ["22", 30_000_000]]) {
  const rs = (await probe.windowRows("grch37", c, lo, lo + 20_000_000, ["common"])).filter(r => r.rsid_source === 1);
  for (let k = 0; k < 250; k++) wide.push(`rs${rs[Math.floor((k * rs.length) / 250)].rsid}`);
}

const CASES = {
  "LDmatrix 8": c => ldmatrix(c, { inputs: EX, pop: "CEU", build: "grch37" }),
  "LDmatrix 300": c => ldmatrix(c, { inputs: L300, pop: "CEU", build: "grch37" }),
  "LDhap 8": c => ldhap(c, { snps: EX, pop: "CEU", build: "grch37" }),
  "LDproxy rs3": c => ldproxy(c, { var: "rs3", pop: "CEU", build: "grch37" }),
  "LDtrait 8": c => ldtrait(c, { snps: EX.join("\n"), pop: "CEU", build: "grch37" }),
  "SNPclip 1000": c => snpclip(c, { snps: L1000.join("\n"), pop: "CEU", build: "grch37" }),
  "SNPclip 5000": c => snpclip(c, { snps: L5000.join("\n"), pop: "CEU", build: "grch37" }),
  "SNPclip 200 sparse": c => snpclip(c, { snps: S200.join("\n"), pop: "CEU", build: "grch37" }),
  "SNPchip 1000 genome-wide": c => snpchip(c, { snps: wide.join("\n"), platforms: "all", build: "grch37" }),
  "LDassoc example": c => ldassoc(c, { text: assoc, columns: ["chr", "pos", "p"], region: "variant", origin: "rs7837688", window: 100000,
    pop: "EUR", build: "grch37" }),
};
const hash = o => createHash("sha256").update(JSON.stringify(canon(o), (k,
  v) => (typeof v === "bigint" ? `${v}n` : v instanceof Uint8Array ? `bytes:${v.length}` : v))).digest("hex").slice(0, 12);

const report = { url, repeat: +rep, cases: {} };
for (const [name, run] of Object.entries(CASES)) {
  if (pick.length && !pick.includes(name)) continue;
  const row = {};
  for (const mode of Object.keys(MODES)) {
    const runs = [];
    for (let k = 0; k < +rep; k++) {
      const ctx = await contextIn(mode); // the manifest fetch happens here, before the counters start
      const s0 = { ...stats }, t0 = performance.now(), out = await run(ctx);
      runs.push({ ms: Math.round(performance.now() - t0), requests: stats.requests - s0.requests, bytes: stats.bytes - s0.bytes,
        hash: hash(out) });
    }
    const med = key => runs.map(r => r[key]).sort((a, b) => a - b)[Math.floor(runs.length / 2)];
    row[mode] = { ms: med("ms"), requests: med("requests"), MB: +(med("bytes") / 1e6).toFixed(2), hash: runs[0].hash };
  }
  row.same = new Set(Object.values(row).map(m => m.hash)).size === 1;
  report.cases[name] = row;
  console.log(JSON.stringify({ case: name, ...Object.fromEntries(Object.entries(row).map(([m, v]) => [m,
    typeof v === "object" ? `${v.requests} req ${v.MB} MB ${v.ms} ms` : v])) }));
}
useDataCwd(); mkdirSync("bench/results", { recursive: true });
const file = `bench/results/planner_modes_${new URL(url).hostname}.json`;
writeFileSync(file, `${JSON.stringify(report, null, 1)}\n`); console.log(`wrote ${file}`);
