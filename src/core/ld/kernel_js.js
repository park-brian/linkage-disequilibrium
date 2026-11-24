// LD kernel, plain-JS backend: the same entry points, algorithm, output order and results as kernel.c (WASM
// backend), over Uint32Array words (2 per u64; w32 = 2·W). Used where WebAssembly is unavailable and as the test oracle;
// about 3.5–5× slower than WASM on large unthresholded bands, the same at LDlink's module sizes (docs/CORE_REVIEW.md §2).
// Popcount is SWAR on 32-bit words (JS has no popcount instruction).

const pop32 = v => {
  v -= (v >>> 1) & 0x55555555;
  v = (v & 0x33333333) + ((v >>> 2) & 0x33333333);
  return Math.imul((v + (v >>> 4)) & 0x0f0f0f0f, 0x01010101) >>> 24;
};
const popAnd = (a, ai, b, bi, w32) => { let c = 0; for (let w = 0; w < w32; w++) c += pop32(a[ai + w] & b[bi + w]); return c; };

/** kernel.c stats(): identical operations on doubles. Writes [r², D′] into s; false when Ms == 0 (NA). */
export function stats(n, pa, pb, x, s) {
  const D = x, B = pa - x, C = pb - x, A = n - pa - pb + x;
  const delta = A * D - B * C, Ms = (A + C) * (B + D) * (A + B) * (C + D);
  if (Ms === 0) return false;
  const dmax = delta < 0 ? ((A + C) * (A + B) < (B + D) * (C + D) ? (A + C) * (A + B) : (B + D) * (C + D))
    : ((A + C) * (C + D) < (A + B) * (B + D) ? (A + C) * (C + D) : (A + B) * (B + D));
  s[0] = delta * delta / Ms; s[1] = (delta < 0 ? -delta : delta) / dmax;
  return true;
}
function r2Bound(pa, pb, n) {
  let m1 = pa / n, m2 = pb / n;
  if (m1 > 0.5) m1 = 1 - m1;
  if (m2 > 0.5) m2 = 1 - m2;
  if (m1 > m2) { const t = m1; m1 = m2; m2 = t; }
  if (m1 <= 0) return 0;
  return m1 * (1 - m2) / (m2 * (1 - m1));
}
function output(maxOut) {
  const m = Math.max(1, maxOut);
  return { written: 0, evaluated: 0, pruned: 0, i: new Uint32Array(m), j: new Uint32Array(m), r2: new Float32Array(m), dprime: new Float32Array(m), s: new Float64Array(2), maxOut };
}
function put(o, i, j, n, pa, pb, x, thr, stat) { // kernel.c put()
  if (!stats(n, pa, pb, x, o.s)) return;
  if ((stat === 0 ? o.s[0] : o.s[1]) < thr) return;
  const k = o.written++;
  if (k < o.maxOut) { o.i[k] = i; o.j[k] = j; o.r2[k] = o.s[0]; o.dprime[k] = o.s[1]; }
}
const done = o => { const count = Math.min(o.written, o.maxOut); return { count, written: o.written, evaluated: o.evaluated, pruned: o.pruned,
  i: o.i.subarray(0, count), j: o.j.subarray(0, count), r2: o.r2.subarray(0, count), dprime: o.dprime.subarray(0, count) }; };

export function rowCounts(alt, n, w32, pa) {
  for (let i = 0; i < n; i++) { let c = 0; for (let w = 0, b = i * w32; w < w32; w++) c += pop32(alt[b + w]); pa[i] = c; }
}

/** kernel.c band_blocked: rows i in [iStart, iEnd), blocks of B, 4-row register tiling, (block, j, i) order */
export function bandBlocked(alt, pa, pos, n, w32, hn, maxDist, thr, stat, prune, maxOut, B, iStart, iEnd) {
  const o = output(maxOut), cut = thr * (1 - 1e-9);
  B = Math.min(1024, Math.max(1, B)); iStart = Math.max(0, iStart); iEnd = Math.min(n, iEnd);
  const cand = new Int32Array(B);
  for (let i0 = iStart; i0 < iEnd; i0 += B) {
    const i1 = Math.min(i0 + B, iEnd);
    for (let j = i0 + 1; j < n && pos[j] - pos[i1 - 1] <= maxDist; j++) {
      let m = 0; const iend = j < i1 ? j : i1, paj = pa[j];
      for (let i = i0; i < iend; i++) {
        if (pos[j] - pos[i] > maxDist) continue;
        if (prune && stat === 0 && r2Bound(pa[i], paj, hn) < cut) { o.pruned++; continue; }
        cand[m++] = i;
      }
      o.evaluated += m;
      const aj = j * w32; let k = 0;
      for (; k + 4 <= m; k += 4) {
        const a0 = cand[k] * w32, a1 = cand[k + 1] * w32, a2 = cand[k + 2] * w32, a3 = cand[k + 3] * w32;
        let c0 = 0, c1 = 0, c2 = 0, c3 = 0;
        for (let w = 0; w < w32; w++) {
          const y = alt[aj + w];
          c0 += pop32(alt[a0 + w] & y); c1 += pop32(alt[a1 + w] & y); c2 += pop32(alt[a2 + w] & y); c3 += pop32(alt[a3 + w] & y);
        }
        put(o, cand[k], j, hn, pa[cand[k]], paj, c0, thr, stat); put(o, cand[k + 1], j, hn, pa[cand[k + 1]], paj, c1, thr, stat);
        put(o, cand[k + 2], j, hn, pa[cand[k + 2]], paj, c2, thr, stat); put(o, cand[k + 3], j, hn, pa[cand[k + 3]], paj, c3, thr, stat);
      }
      for (; k < m; k++) put(o, cand[k], j, hn, pa[cand[k]], paj, popAnd(alt, cand[k] * w32, alt, aj, w32), thr, stat);
    }
  }
  return done(o);
}

/** kernel.c band_masked: per-pair valid intersection (one popcount when cls[i] == cls[j]), (i, j) order; i in [iStart, iEnd) */
export function bandMasked(alt, valid, pa, nv, cls, pos, n, w32, maxDist, thr, stat, maxOut, iStart = 0, iEnd = n) {
  const o = output(maxOut);
  for (let i = iStart; i < iEnd; i++) {
    const ai = i * w32;
    for (let j = i + 1; j < n && pos[j] - pos[i] <= maxDist; j++) {
      const aj = j * w32;
      o.evaluated++;
      if (cls[i] === cls[j]) put(o, i, j, nv[i], pa[i], pa[j], popAnd(alt, ai, alt, aj, w32), thr, stat);
      else put(o, i, j, popAnd(valid, ai, valid, aj, w32), popAnd(alt, ai, valid, aj, w32), popAnd(alt, aj, valid, ai, w32), popAnd(alt, ai, alt, aj, w32), thr, stat);
    }
  }
  return done(o);
}

/** kernel.c tile: rows [r0, r0+nr) × cols [c0, c0+nc), bytes 0..250 (value·250 rounded), 255 = NA */
export function tile(alt, pa, w32, hn, r0, nr, c0, nc) {
  const r2 = new Uint8Array(nr * nc), dprime = new Uint8Array(nr * nc), s = new Float64Array(2);
  for (let a = 0; a < nr; a++) for (let b = 0; b < nc; b++) {
    const k = a * nc + b;
    if (stats(hn, pa[r0 + a], pa[c0 + b], popAnd(alt, (r0 + a) * w32, alt, (c0 + b) * w32, w32), s)) {
      r2[k] = (s[0] * 250 + 0.5) | 0; dprime[k] = (s[1] * 250 + 0.5) | 0;
    } else { r2[k] = 255; dprime[k] = 255; }
  }
  return { r2, dprime };
}

/** masked-path counts for one pair (kernel.c band_masked): [n valid in both, alt_i among valid_j, alt_j among valid_i, both alt] */
export function maskedCounts(alt, valid, pa, nv, cls, i, j, w32) {
  const ai = i * w32, aj = j * w32;
  if (cls[i] === cls[j]) return [nv[i], pa[i], pa[j], popAnd(alt, ai, alt, aj, w32)];
  return [popAnd(valid, ai, valid, aj, w32), popAnd(alt, ai, valid, aj, w32), popAnd(alt, aj, valid, ai, w32), popAnd(alt, ai, alt, aj, w32)];
}

export function pairCounts(alt, w32, I, J) {
  const x = new Uint32Array(I.length);
  for (let k = 0; k < I.length; k++) x[k] = popAnd(alt, I[k] * w32, alt, J[k] * w32, w32);
  return x;
}
