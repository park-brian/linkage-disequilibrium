// Pyramid worker: LD for every pair within maxDist (the kernel's band path), aggregated to screen pixels; with `range`,
// only pairs whose first variant is in [r0, r1) (one share of a pool, merged by pixel max on the page).
// Pair (i, j) sits at x = midpoint, y = half the distance (a 45° pyramid on the genomic axis); a pixel shows the
// strongest pair that lands on it, splatted to `radius` pixels when variants are sparse so cells touch.
// Imports only the kernel (dependency-free chain), so it needs no import map.
import { loadKernel, KERNEL_WASM } from "../src/ld/kernel.mjs";

let kernel;
self.onmessage = async ({ data: q }) => {
  try {
    kernel ??= await loadKernel(new Uint8Array(await (await fetch(KERNEL_WASM)).arrayBuffer()));
    const t0 = performance.now();
    const p = kernel.prepare(q.rows, { samples: q.samples, haploidMales: new Set(q.haploidMales), hapIndexes: q.hapIndexes });
    // rows [r0, r1): this worker's share of the band (pairs i in the range, j > i anywhere), on either kernel path
    const b = kernel.band(p, { maxDist: q.maxDist, threshold: q.threshold, stat: q.stat, prune: true, ...(q.range ? { rows: q.range } : {}) });
    const t1 = performance.now();
    const { W, lo, hi } = q, span = Math.max(1, hi - lo), H = Math.ceil(((q.maxDist / 2) / span) * W) + 1;
    const grid = new Float32Array(W * H).fill(-1), bi = new Int32Array(W * H), bj = new Int32Array(W * H);
    const pos = q.rows.map(r => r.pos), val = q.stat === "r2" ? b.r2 : b.dprime, R = q.radius, off = q.offset ?? 0;
    for (let k = 0; k < b.count; k++) {
      const li = b.i[k], lj = b.j[k], i = li + off, j = lj + off, v = val[k], a = pos[li], c = pos[lj]; // i, j: window indexes
      const x = Math.round((((a + c) / 2 - lo) / span) * W), y = Math.round((((c - a) / 2) / span) * W);
      for (let dy = -R; dy <= R; dy++) for (let dx = -R; dx <= R; dx++) {
        const xx = x + dx, yy = y + dy; if (xx < 0 || xx >= W || yy < 0 || yy >= H) continue;
        const o = yy * W + xx; // ties: the smallest (i, j), so the result does not depend on processing order
        if (v > grid[o] || (v === grid[o] && (i < bi[o] || (i === bi[o] && j < bj[o])))) { grid[o] = v; bi[o] = i; bj[o] = j; }
      }
    }
    kernel.reset();
    self.postMessage({ ok: true, W, H, grid, bi, bj, pairs: b.count, evaluated: b.evaluated, pruned: b.pruned, backend: kernel.backend,
      ms: { kernel: Math.round(t1 - t0), aggregate: Math.round(performance.now() - t1) } }, [grid.buffer, bi.buffer, bj.buffer]);
  } catch (e) { self.postMessage({ ok: false, error: String(e?.stack ?? e) }); }
};
