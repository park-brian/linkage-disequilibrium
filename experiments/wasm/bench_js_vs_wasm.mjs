// WASM vs JS backend of the same kernel API (kernel.mjs): identical inputs and calls; outputs must be identical
// (same pairs, order and float32 values). Reports Mpairs/s (evaluated pairs / kernel time), best of R, single thread.
// usage (from experiments/): node wasm/bench_js_vs_wasm.mjs <chrN.common.parquet> <panelDir> <lo> <hi> [maxDist] [R] [out.json]
import { readFileSync, writeFileSync } from "node:fs";
import { asyncBufferFromFile, parquetMetadataAsync, parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { loadKernel, jsKernel, KERNEL_WASM } from "../../src/ld/kernel.mjs";

const [path, panelDir, loArg, hiArg, distArg, rArg, outPath] = process.argv.slice(2);
const lo = +loArg, hi = +hiArg, maxDist = +(distArg ?? 500000), R = +(rArg ?? 3);
const file = await asyncBufferFromFile(path), metadata = await parquetMetadataAsync(file);
const samples = Object.fromEntries(metadata.key_value_metadata.map(x => [x.key, x.value])).samples.split(",");
let rows = [];
await parquetRead({ file, metadata, compressors, utf8: false, rowFormat: "object", columns: ["pos", "haps"], filter: { pos: { $gte: lo, $lte: hi } }, onComplete: r => { rows = r; } });
rows = rows.filter(r => r.pos >= lo && r.pos <= hi);
const backends = { wasm: await loadKernel(readFileSync(KERNEL_WASM)), js: jsKernel() };
const hapIdx = pop => pop === "ALL" ? samples.flatMap((_, i) => [2 * i, 2 * i + 1])
  : (s => samples.flatMap((x, i) => (s.has(x) ? [2 * i, 2 * i + 1] : [])))(new Set(readFileSync(`${panelDir}/${pop}.txt`, "utf8").split(/\s+/).filter(Boolean)));
const best = f => { let ms = Infinity, r; for (let k = 0; k < R; k++) { const t0 = performance.now(); r = f(); ms = Math.min(ms, performance.now() - t0); } return { r, ms }; };
const same = (a, b) => { if (a.count !== b.count) return false; for (let k = 0; k < a.count; k++) if (a.i[k] !== b.i[k] || a.j[k] !== b.j[k] || a.r2[k] !== b.r2[k] || a.dprime[k] !== b.dprime[k]) return false; return true; };
const results = { file: path.split(/[\\/]/).slice(-3).join("/"), window: [lo, hi], rows: rows.length, maxDist, runs: R, node: process.version, cases: [] };
for (const pop of ["CEU", "EUR", "ALL"]) for (const threshold of [0, 0.8]) {
  const opts = { maxDist, threshold }, t = {};
  for (const [name, k] of Object.entries(backends)) { const p = k.prepare(rows, { samples: samples.length, hapIndexes: hapIdx(pop) }); t[name] = best(() => k.band(p, opts)); k.reset(); }
  const ev = t.wasm.r.evaluated, c = { pop, threshold, evaluated: ev, kept: t.wasm.r.count, pruned: t.wasm.r.pruned, identical: same(t.wasm.r, t.js.r),
    wasm_mpairs_s: +(ev / t.wasm.ms / 1e3).toFixed(1), js_mpairs_s: +(ev / t.js.ms / 1e3).toFixed(1), wasm_ms: +t.wasm.ms.toFixed(1), js_ms: +t.js.ms.toFixed(1), speedup: +(t.js.ms / t.wasm.ms).toFixed(2) };
  results.cases.push(c); console.log(JSON.stringify(c));
}
if (outPath) writeFileSync(outPath, JSON.stringify(results, null, 2));
