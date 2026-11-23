// LD kernel API with two interchangeable backends and identical results (pairs, order, float32 values):
//   WASM  kernel.c → kernel.wasm (scalar i64.popcnt, cache-blocked)      — loadKernel(bytes)
//   JS    kernel_js.js (SWAR popcount, same algorithm)                      — jsKernel(), or loadKernel() fallback
// WASM is 3.5–5.3× faster on large unthresholded bands and 1.7–3× with pruning; at LDlink's module sizes both take
// milliseconds (docs/CORE_REVIEW.md §2).
//
//   const k = await loadKernel(bytes)                       // bytes omitted or no WebAssembly → JS backend
//   const p = k.prepare(rows, { samples, haploidMales, hapIndexes })   // rows: { pos, haps, missing?, male_haploid? }
//   k.band(p, { maxDist, threshold, stat: "r2"|"dprime", prune, strict, rows: [r0, r1), block, maxOut })
//        → { i, j, r2, dprime, count, written, truncated, evaluated, pruned }   order (block, j, i); masked rows: (i, j)
//   k.tile(p, r0, nr, c0, nc)   → { r2, dprime: Uint8Array }  (0..250 = value·250 rounded, 255 = NA; unmasked rows)
//   k.pairCounts(p, I, J)       → Uint32Array |alt_i ∧ alt_j|
//   k.counts(p, I, J)          → [[n, alt_i, alt_j, both]] per pair, over haplotypes valid in both rows
//   k.exact(p, I, J)            → { r2, dprime: Float64Array }  (NaN = NA), exact doubles for 3-decimal output
//   k.reset()                   WASM: frees every prepared set (WASM memory is a bump arena). JS: no-op.
// Thresholds compare the unrounded DOUBLE r²/D′ (LDlink's LD_calcs / LDproxy / LDassoc / SNPclip); stats() is
// bit-identical to Python. strict: true keeps stat > threshold (compares against nextUp(threshold)).
import { nextUp } from "../pyfmt.js";
import * as js from "./kernel_js.js";
import { decodeHaps } from "../haps.js";

/** the WASM binary next to this module (node: readFileSync(KERNEL_WASM); browser: fetch(KERNEL_WASM)) */
export const KERNEL_WASM = new URL("./kernel.wasm", import.meta.url);

const DEFAULT_BLOCK = 128; // measured best for both backends (PARQUET_LAYOUT.md §9)
const OUT_CAP = 1 << 22;   // default output buffer: 4M pairs; grows once to exactly `written` if more pass
const EMPTY = () => ({ count: 0, written: 0, truncated: false, evaluated: 0, pruned: 0, i: new Uint32Array(0), j: new Uint32Array(0), r2: new Float32Array(0), dprime: new Float32Array(0) });


/**
 * Pack the selected haplotypes of every row into a: alt bits, plus a valid mask when some haplotype is uncalled or
 * haploid (chrX males outside the PARs, '.' calls). Rows are w32 = 2·ceil(nh/64) words; a's arrays come zeroed.
 */
function pack(rows, { samples, haploidMales = new Set(), hapIndexes }, a, w32, masked) {
  const nh = hapIndexes.length, hbytes = Math.ceil((2 * samples) / 8), alt8 = new Uint8Array(a.alt.buffer, a.alt.byteOffset, a.alt.byteLength);
  const identity = !masked && hapIndexes.every((k, j) => k === j); // a leading prefix of the panel: copy bytes, clear bits >= nh
  for (let i = 0; i < rows.length; i++) {
    const r = rows[i], bits = decodeHaps(r.haps, hbytes), base = i * w32;
    a.pos[i] = r.pos;
    if (identity) {
      alt8.set(bits.subarray(0, Math.min(hbytes, w32 * 4)), base * 4);
      if (nh % 32) a.alt[base + (nh >> 5)] &= (1 << (nh % 32)) - 1;
      a.alt.fill(0, base + Math.ceil(nh / 32), base + w32);
      continue;
    }
    for (let j = 0; j < nh; j++) {
      const k = hapIndexes[j];
      if (masked) {
        if ((k & 1) && r.male_haploid && haploidMales.has(k >> 1)) continue;
        if (r.missing && (r.missing[k >> 3] >> (k & 7)) & 1) continue;
        a.valid[base + (j >> 5)] |= 1 << (j & 31);
      }
      if ((bits[k >> 3] >> (k & 7)) & 1) a.alt[base + (j >> 5)] |= 1 << (j & 31);
    }
  }
}

/**
 * Masked rows: nv[i] = valid haplotypes in row i; cls[i] = a class id shared exactly by rows with identical valid masks
 * (FNV-1a hash, then a full compare against each earlier class with that hash), so the kernel tests cls[i] == cls[j] alone.
 */
function maskClasses({ valid, nv, cls }, n, w32) {
  const reps = new Map(); let next = 0;
  for (let i = 0; i < n; i++) {
    const b = i * w32; let h = 0x811c9dc5, c = 0;
    for (let k = b; k < b + w32; k++) { const w = valid[k]; h = Math.imul(h ^ w, 0x01000193); c += popcount32(w); }
    nv[i] = c;
    let list = reps.get(h), id = -1;
    if (!list) reps.set(h, (list = []));
    for (const [r, rid] of list) {
      let same = true;
      for (let k = 0; k < w32 && same; k++) same = valid[r * w32 + k] === valid[b + k];
      if (same) { id = rid; break; }
    }
    if (id < 0) list.push([i, (id = next++)]);
    cls[i] = id;
  }
}
const popcount32 = w => { w -= (w >>> 1) & 0x55555555; w = (w & 0x33333333) + ((w >>> 2) & 0x33333333); return (((w + (w >>> 4)) & 0x0f0f0f0f) * 0x01010101) >>> 24; };

/**
 * Shared front end over an engine:
 *   alloc(n, w32, masked) → p with p.views() = { alt, valid, pa, pos } (fresh, zeroed)
 *   rowCounts(p), bandBlocked(p, maxDist, thr, stat, prune, cap, B, r0, r1), bandMasked(p, maxDist, thr, stat, cap, r0, r1),
 *   tile(p, r0, nr, c0, nc), pairCounts(p, I, J), reset()
 */
function makeKernel(engine) {
  function prepare(rows, opts) {
    const nh = opts.hapIndexes.length, w32 = 2 * Math.ceil(nh / 64), masked = rows.some(r => r.missing || r.male_haploid);
    const p = engine.alloc(rows.length, w32, masked);
    Object.assign(p, { n: rows.length, nh, w32, masked });
    pack(rows, opts, p.views(), w32, masked);
    engine.rowCounts(p);
    if (masked) maskClasses(p.views(), rows.length, w32);
    return p;
  }
  /** exact upper bound on band output: pairs i < j within maxDist for rows i in [r0, r1) (two pointers) */
  function pairBound(pos, n, maxDist, r0, r1) {
    let b = 0, j = r0;
    for (let i = r0; i < r1; i++) { if (j < i + 1) j = i + 1; while (j < n && pos[j] - pos[i] <= maxDist) j++; b += j - i - 1; }
    return b;
  }
  function band(p, { maxDist = 500000, threshold = 0, stat = "r2", prune = true, maxOut, block = DEFAULT_BLOCK, rows, strict = false } = {}) {
    if (strict) threshold = nextUp(threshold); // stat > t  ⇔  stat >= nextUp(t) for doubles
    const [r0, r1] = rows ?? [0, p.n];
    if (!Number.isInteger(r0) || !Number.isInteger(r1)) throw new Error("band: rows must be [start, end) integers");
    const st = stat === "r2" ? 0 : 1, lo = Math.max(0, r0), hi = Math.min(p.n, r1);
    if (hi <= lo || p.n < 2) return EMPTY();
    const run = cap => {
      const o = p.masked ? engine.bandMasked(p, maxDist, threshold, st, cap, lo, hi) : engine.bandBlocked(p, maxDist, threshold, st, prune ? 1 : 0, cap, block, lo, hi);
      return { ...o, truncated: o.written > o.count };
    };
    if (maxOut !== undefined) return run(maxOut);
    const first = run(Math.max(1, Math.min(pairBound(p.views().pos, p.n, maxDist, lo, hi), OUT_CAP)));
    return first.truncated ? run(first.written) : first;
  }
  function tile(p, r0, nr, c0, nc) {
    if (p.masked) throw new Error("tile: masked rows are not supported (use band)");
    return engine.tile(p, r0, nr, c0, nc);
  }
  const pairCounts = (p, I, J) => engine.pairCounts(p, I, J);
  /** per pair [n, alt_i, alt_j, alt_i∧alt_j] over the haplotypes valid in both rows (the 2×2 table's margins) */
  function counts(p, I, J) {
    const x = p.masked ? null : pairCounts(p, I, J), v = p.views(); // views AFTER any allocation: memory.grow detaches them
    return Array.from(I, (i, k) => (p.masked ? js.maskedCounts(v.alt, v.valid, v.pa, v.nv, v.cls, i, J[k], p.w32) : [p.nh, v.pa[i], v.pa[J[k]], x[k]]));
  }
  function exact(p, I, J) {
    const m = I.length, r2 = new Float64Array(m), dprime = new Float64Array(m), s = new Float64Array(2), c = counts(p, I, J);
    for (let k = 0; k < m; k++) if (js.stats(...c[k], s)) { r2[k] = s[0]; dprime[k] = s[1]; } else { r2[k] = NaN; dprime[k] = NaN; }
    return { r2, dprime };
  }
  return { backend: engine.name, prepare, band, tile, pairCounts, counts, exact, reset: engine.reset };
}

/** Plain-JS backend (kernel_js.js): ordinary typed arrays, garbage-collected. */
export function jsKernel() {
  return makeKernel({
    name: "js",
    alloc(n, w32, masked) {
      const a = { alt: new Uint32Array(n * w32), valid: masked ? new Uint32Array(n * w32) : null, pa: new Uint32Array(n), pos: new Int32Array(n),
        nv: masked ? new Uint32Array(n) : null, cls: masked ? new Int32Array(n) : null };
      return { views: () => a };
    },
    rowCounts: p => { const v = p.views(); js.rowCounts(v.alt, p.n, p.w32, v.pa); },
    bandBlocked: (p, d, t, s, pr, cap, B, r0, r1) => { const v = p.views(); return js.bandBlocked(v.alt, v.pa, v.pos, p.n, p.w32, p.nh, d, t, s, pr, cap, B, r0, r1); },
    bandMasked: (p, d, t, s, cap, r0, r1) => { const v = p.views(); return js.bandMasked(v.alt, v.valid, v.pa, v.nv, v.cls, v.pos, p.n, p.w32, d, t, s, cap, r0, r1); },
    tile: (p, r0, nr, c0, nc) => { const v = p.views(); return js.tile(v.alt, v.pa, p.w32, p.nh, r0, nr, c0, nc); },
    pairCounts: (p, I, J) => js.pairCounts(p.views().alt, p.w32, I, J),
    reset() {},
  });
}

/** WASM backend over a bump arena in the module's memory; falls back to JS when there are no bytes or no WebAssembly. */
export async function loadKernel(bytes) {
  if (!bytes || typeof WebAssembly === "undefined") return jsKernel();
  const { instance } = await WebAssembly.instantiate(bytes, {});
  const ex = instance.exports, mem = ex.memory, base = (ex.__heap_base.value + 15) & ~15;
  let top = base;
  const alloc = nbytes => {
    const at = top; top = (top + nbytes + 15) & ~15;
    if (top > mem.buffer.byteLength) mem.grow(Math.ceil((top - mem.buffer.byteLength) / 65536));
    return at;
  };
  const scratch = f => { const save = top; try { return f(); } finally { top = save; } }; // per-call buffers
  // memory.grow detaches every view, so arrays are always re-viewed from their offsets
  const u32 = (ptr, n) => new Uint32Array(mem.buffer, ptr, n), i32 = (ptr, n) => new Int32Array(mem.buffer, ptr, n), f32 = (ptr, n) => new Float32Array(mem.buffer, ptr, n);
  const W = p => p.w32 / 2; // the C side counts u64 words
  const bandOut = (cap, call) => scratch(() => {
    const ip = alloc(cap * 4), jp = alloc(cap * 4), rp = alloc(cap * 4), dp = alloc(cap * 4);
    const count = call(ip, jp, rp, dp);
    return { count, written: ex.written(), evaluated: ex.evaluated(), pruned: ex.pruned(),
      i: u32(ip, count).slice(), j: u32(jp, count).slice(), r2: f32(rp, count).slice(), dprime: f32(dp, count).slice() };
  });
  return makeKernel({
    name: "wasm",
    alloc(n, w32, masked) {
      const altPtr = alloc(n * w32 * 4), validPtr = masked ? alloc(n * w32 * 4) : 0, paPtr = alloc(n * 4), posPtr = alloc(n * 4);
      const nvPtr = masked ? alloc(n * 4) : 0, clsPtr = masked ? alloc(n * 4) : 0;
      const p = { altPtr, validPtr, paPtr, posPtr, nvPtr, clsPtr, views: () => ({ alt: u32(altPtr, n * w32), valid: masked ? u32(validPtr, n * w32) : null,
        pa: u32(paPtr, n), pos: i32(posPtr, n), nv: masked ? u32(nvPtr, n) : null, cls: masked ? i32(clsPtr, n) : null }) };
      const v = p.views(); v.alt.fill(0); if (v.valid) v.valid.fill(0); // the arena is reused after reset(): zero it
      return p;
    },
    rowCounts: p => ex.row_counts(p.altPtr, p.n, W(p), p.paPtr),
    bandBlocked: (p, d, t, s, pr, cap, B, r0, r1) => bandOut(cap, (ip, jp, rp, dp) => ex.band_blocked(p.altPtr, p.paPtr, p.posPtr, p.n, W(p), p.nh, d, t, s, pr, ip, jp, rp, dp, cap, B, r0, r1)),
    bandMasked: (p, d, t, s, cap, r0, r1) => bandOut(cap, (ip, jp, rp, dp) => ex.band_masked(p.altPtr, p.validPtr, p.paPtr, p.nvPtr, p.clsPtr, p.posPtr, p.n, W(p), d, t, s, ip, jp, rp, dp, cap, r0, r1)),
    tile: (p, r0, nr, c0, nc) => scratch(() => {
      const a = alloc(nr * nc), b = alloc(nr * nc);
      ex.tile(p.altPtr, p.paPtr, W(p), p.nh, r0, nr, c0, nc, a, b);
      return { r2: new Uint8Array(mem.buffer, a, nr * nc).slice(), dprime: new Uint8Array(mem.buffer, b, nr * nc).slice() };
    }),
    pairCounts: (p, I, J) => scratch(() => {
      const m = I.length, ip = alloc(m * 4), jp = alloc(m * 4), xp = alloc(m * 4);
      u32(ip, m).set(I); u32(jp, m).set(J);
      ex.pair_counts(p.altPtr, W(p), ip, jp, m, xp);
      return u32(xp, m).slice();
    }),
    reset() { top = base; },
  });
}
