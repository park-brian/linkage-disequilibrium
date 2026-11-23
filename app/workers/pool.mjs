// A pool of pyramid workers: the band's first-variant rows are split into shares with equal pair counts, each worker
// gets only its share plus the rows within maxDist after it, and the page merges the pixel grids (max per pixel).
// Rows with uncalled or haploid haplotypes (chrX males, '.') take the kernel's masked path; it splits the same way.
const SIZE = Math.max(1, Math.min(8, (navigator.hardwareConcurrency || 4) - 1));
const workers = [];
const worker = k => (workers[k] ??= new Worker(new URL("./pyramid.mjs", import.meta.url), { type: "module" }));
const call = (k, msg) => new Promise((res, rej) => {
  const w = worker(k);
  w.onmessage = ({ data }) => (data.ok ? res(data) : rej(new Error(data.error))); w.onerror = e => rej(e);
  w.postMessage(msg);
});

/**
 * @param rows [{ pos, haps, missing?, male_haploid? }] sorted by pos
 * @param q { samples, haploidMales, hapIndexes, maxDist, threshold, stat, W, lo, hi, radius }
 * → { W, H, grid, bi, bj, pairs, evaluated, pruned, backend, workers, ms: { kernel (slowest share), aggregate, merge } }
 */
export async function pyramid(rows, q, { size = SIZE } = {}) {
  const n = rows.length, pos = rows.map(r => r.pos);
  // pairs per first row (j > i within maxDist), by two pointers; shares of equal cumulative pairs
  const cnt = new Float64Array(n); for (let i = 0, j = 0; i < n; i++) { if (j < i + 1) j = i + 1; while (j < n && pos[j] - pos[i] <= q.maxDist) j++; cnt[i] = j - i - 1; }
  const total = cnt.reduce((a, b) => a + b, 0), k = Math.max(1, Math.min(size, Math.ceil(total / 2e6))); // ≥ 2M pairs per share
  const cuts = [0]; let acc = 0;
  for (let i = 0; i < n && cuts.length < k; i++) { acc += cnt[i]; if (acc >= (total * cuts.length) / k) cuts.push(i + 1); }
  cuts.push(n);
  const shares = [];
  for (let s = 0; s + 1 < cuts.length; s++) {
    const r0 = cuts[s], r1 = cuts[s + 1]; if (r1 <= r0) continue;
    let end = r1; while (end < n && pos[end] - pos[r1 - 1] <= q.maxDist) end++;
    shares.push({ r0, r1, end });
  }
  const t0 = performance.now();
  const parts = await Promise.all(shares.map(({ r0, r1, end }, s) => call(s, { ...q, rows: rows.slice(r0, end), offset: r0, range: [0, r1 - r0] })));
  const t1 = performance.now(), { W, H } = parts[0], grid = parts[0].grid, bi = parts[0].bi, bj = parts[0].bj;
  for (let s = 1; s < parts.length; s++) { // max per pixel; ties: the smallest (i, j), as within a worker
    const g = parts[s].grid, pi = parts[s].bi, pj = parts[s].bj;
    for (let o = 0; o < g.length; o++) if (g[o] > grid[o] || (g[o] === grid[o] && g[o] >= 0 && (pi[o] < bi[o] || (pi[o] === bi[o] && pj[o] < bj[o])))) { grid[o] = g[o]; bi[o] = pi[o]; bj[o] = pj[o]; }
  }
  const sum = key => parts.reduce((a, p) => a + p[key], 0);
  return { W, H, grid, bi, bj, pairs: sum("pairs"), evaluated: sum("evaluated"), pruned: sum("pruned"), backend: parts[0].backend, workers: parts.length,
    ms: { kernel: Math.max(...parts.map(p => p.ms.kernel)), aggregate: Math.max(...parts.map(p => p.ms.aggregate)), wall: Math.round(t1 - t0), merge: Math.round(performance.now() - t1) } };
}
