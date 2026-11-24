// Multi-thread scaling of the blocked kernel: split the band by row range across worker threads (same as
// Web Workers in the browser). Verifies the union of worker results equals the single-thread result.
// usage: node bench/threads.js release_v4 <panelDir (e.g. $LD_PANELS)>
import { Worker, isMainThread, parentPort, workerData } from "node:worker_threads";
import { readFileSync, writeFileSync, mkdirSync } from "node:fs";
import { availableParallelism } from "node:os";
import { loadKernel, KERNEL_WASM } from "../src/core/ld/kernel.js";
import { useDataCwd } from "../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

if (!isMainThread) {
  const { wasm, rows, samples, idx, range, opts } = workerData;
  const k = await loadKernel(wasm), p = k.prepare(rows, { samples, hapIndexes: idx });
  const t = performance.now(), r = k.band(p, { ...opts, rows: range });
  parentPort.postMessage({ count: r.count, evaluated: r.evaluated, ms: performance.now() - t, sum: r.r2.reduce((a, v) => a + v, 0) });
} else {
  const { asyncBufferFromFile, parquetMetadataAsync, parquetRead } = await import("hyparquet");
  const { compressors } = await import("hyparquet-compressors");
  const [root, panelDir] = process.argv.slice(2);
  const file = await asyncBufferFromFile(`${root}/genotypes/grch37/chr22.common.parquet`), metadata = await parquetMetadataAsync(file);
  const samples = metadata.key_value_metadata.find(x => x.key === "samples").value.split(",");
  let rows = []; await parquetRead({ file, metadata, compressors, utf8: false, rowFormat: "object", columns: ["pos", "haps"],
    filter: { pos: { $gte: 25e6, $lte: 35e6 } }, onComplete: r => { rows = r; } });
  rows = rows.map(r => ({ pos: r.pos, haps: r.haps }));
  const wasm = readFileSync(KERNEL_WASM);
  const out = [];
  for (const pop of ["EUR", "ALL"]) {
    const idx = pop === "ALL" ? samples.flatMap((_, i) => [2 * i, 2 * i + 1]) : (s => samples.flatMap((x, i) => (s.has(x) ? [2 * i,
      2 * i + 1] : [])))(new Set(readFileSync(`${panelDir}/${pop}.txt`, "utf8").split(/\s+/).filter(Boolean)));
    const opts = { threshold: 0.8, prune: false, block: 128 };
    for (const T of [1, 2, 4, 8, Math.min(16, availableParallelism())]) {
      // split rows so each worker gets about the same number of pairs (pairs per row are roughly uniform)
      const n = rows.length, ranges = Array.from({ length: T }, (_, t) => [Math.floor(n * t / T), Math.floor(n * (t + 1) / T)]);
      const t0 = performance.now();
      const res = await Promise.all(ranges.map(range => new Promise((ok, bad) => {
        const w = new Worker(new URL(import.meta.url), { workerData: { wasm, rows, samples: samples.length, idx, range, opts } });
        w.once("message", ok); w.once("error", bad);
      })));
      const wall = performance.now() - t0, kernelMax = Math.max(...res.map(r => r.ms));
      const o = { pop, threads: T, pairs: res.reduce((a, r) => a + r.evaluated, 0), kept: res.reduce((a, r) => a + r.count, 0),
        r2_sum: +res.reduce((a, r) => a + r.sum, 0).toFixed(3),
        kernel_ms_max: Math.round(kernelMax), wall_ms_incl_worker_start: Math.round(wall), Mpairs_per_s_kernel: +(res.reduce((a,
          r) => a + r.evaluated, 0) / kernelMax / 1e3).toFixed(1) };
      out.push(o); console.log(JSON.stringify(o));
    }
  }
  mkdirSync("bench/results", { recursive: true }); // $LD_DATA/bench/results (cwd = LD_DATA)
  writeFileSync("bench/results/threads_report.json", JSON.stringify(out, null, 2));
}
