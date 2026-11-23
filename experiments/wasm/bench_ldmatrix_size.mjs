// LDmatrix-sized requests (all pairs among N SNPs, ALL population): WASM vs JS backend, end to end per request
// (prepare + band), plus the one-off WASM compile/instantiate cost. usage (from experiments/): node wasm/bench_ldmatrix_size.mjs
import { readFileSync } from "node:fs";
import { asyncBufferFromFile, parquetMetadataAsync, parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { loadKernel, jsKernel, KERNEL_WASM } from "../../src/ld/kernel.mjs";

const file = await asyncBufferFromFile("out/release/genotypes/grch37/chr22.common.parquet"), metadata = await parquetMetadataAsync(file);
const samples = Object.fromEntries(metadata.key_value_metadata.map(x => [x.key, x.value])).samples.split(",");
let rows = []; await parquetRead({ file, metadata, compressors, utf8: false, rowFormat: "object", columns: ["pos", "haps"], filter: { pos: { $gte: 30e6, $lte: 31e6 } }, onComplete: r => { rows = r; } });
const bytes = readFileSync(KERNEL_WASM);
const t0 = performance.now(), wasm = await loadKernel(bytes), instMs = performance.now() - t0, js = jsKernel();
const all = samples.flatMap((_, i) => [2 * i, 2 * i + 1]);
const time = f => { f(); let best = Infinity; for (let r = 0; r < 5; r++) { const t = performance.now(); f(); best = Math.min(best, performance.now() - t); } return +best.toFixed(2); };
console.log(JSON.stringify({ wasm_instantiate_ms: +instMs.toFixed(2), wasm_bytes: bytes.length }));
for (const N of [30, 300, 3000]) {
  const rs = rows.slice(0, N), run = k => () => { const p = k.prepare(rs, { samples: samples.length, hapIndexes: all }); k.band(p, { maxDist: 1e9, prune: false }); k.reset(); };
  console.log(JSON.stringify({ snps: N, pairs: N * (N - 1) / 2, population: "ALL", wasm_ms: time(run(wasm)), js_ms: time(run(js)) }));
}
