// LD kernel (WASM backend): pairwise haplotype counts by popcount, then r² and D′ with LDlink's formulas
// (server/LDmatrix.py:303-318, LDcommon.LD_calcs; counted by allele index). Built standalone (no libc) by build.sh into
// kernel.wasm (scalar i64.popcnt — V8's SIMD popcount is not faster, docs/CORE_REVIEW.md §2).
// kernel_js.mjs is the plain-JS backend with the same entry points and results; kernel.mjs is the shared API.
//
// Data layout (prepared by kernel.mjs):
//   alt   : n variants × W u64 words (W = ceil(hn / 64)), alt-allele bit per population haplotype (valid haplotypes only)
//   valid : same shape, 1 bit per called haplotype (masked entry point only: chrX males, '.' calls)
//   pa    : alt count per variant (over its valid haplotypes: alt bits are only set where valid)
//   nv,cls: masked only: valid count per variant, and a mask class (equal class <=> identical valid rows; kernel.mjs)
//   pos   : position per variant (ascending)
// Counts are exact integers; stats() is bit-identical to Python (integer products < 2^53, one rounding per op).
#include <stdint.h>

#define EXPORT(name) __attribute__((export_name(#name)))
typedef uint64_t u64;
typedef uint32_t u32;
typedef int32_t i32;

// counters for the last call: written counts every pair that passed, including any beyond maxOut (truncation check)
static u32 g_evaluated, g_pruned, g_written;
EXPORT(evaluated) u32 evaluated(void) { return g_evaluated; }
EXPORT(pruned) u32 pruned(void) { return g_pruned; }
EXPORT(written) u32 written(void) { return g_written; }

static inline u32 popcount_and(const u64* a, const u64* b, i32 W) {
  u32 c = 0;
  for (i32 w = 0; w < W; w++) c += (u32)__builtin_popcountll(a[w] & b[w]);
  return c;
}
// popcount(a_k & b) for 4 rows against one row b, loading each word of b once (register tiling)
static inline void popcount_and4(const u64* a0, const u64* a1, const u64* a2, const u64* a3, const u64* b, i32 W, u32* c) {
  u32 c0 = 0, c1 = 0, c2 = 0, c3 = 0;
  for (i32 w = 0; w < W; w++) {
    u64 y = b[w];
    c0 += (u32)__builtin_popcountll(a0[w] & y); c1 += (u32)__builtin_popcountll(a1[w] & y);
    c2 += (u32)__builtin_popcountll(a2[w] & y); c3 += (u32)__builtin_popcountll(a3[w] & y);
  }
  c[0] = c0; c[1] = c1; c[2] = c2; c[3] = c3;
}
EXPORT(row_counts) void row_counts(const u64* bits, i32 n, i32 W, u32* out) {
  for (i32 i = 0; i < n; i++) { u32 c = 0; for (i32 w = 0; w < W; w++) c += (u32)__builtin_popcountll(bits[(u64)i * W + w]); out[i] = c; }
}

// r² and D′ from counts; returns 0 when Ms == 0 (LDlink "NA")
static inline int stats(double n, double pa, double pb, double x, double* r2, double* dp) {
  double D = x, B = pa - x, C = pb - x, A = n - pa - pb + x;
  double delta = A * D - B * C, Ms = (A + C) * (B + D) * (A + B) * (C + D);
  if (Ms == 0) return 0;
  double dmax = delta < 0 ? ((A + C) * (A + B) < (B + D) * (C + D) ? (A + C) * (A + B) : (B + D) * (C + D))
                          : ((A + C) * (C + D) < (A + B) * (B + D) ? (A + C) * (C + D) : (A + B) * (B + D));
  *r2 = delta * delta / Ms;
  *dp = (delta < 0 ? -delta : delta) / dmax;
  return 1;
}
// upper bound on r² from minor allele frequencies (m1 ≤ m2 ≤ 0.5): r² ≤ m1(1−m2) / (m2(1−m1))
static inline double r2_bound(double pa, double pb, double n) {
  double m1 = pa / n, m2 = pb / n;
  if (m1 > 0.5) m1 = 1 - m1;
  if (m2 > 0.5) m2 = 1 - m2;
  if (m1 > m2) { double t = m1; m1 = m2; m2 = t; }
  if (m1 <= 0) return 0;
  return m1 * (1 - m2) / (m2 * (1 - m1));
}

// keep (i, j) if stat >= thr; write while there is room, count always
static inline void put(i32 i, i32 j, double n, double pa, double pb, u32 x, double thr, i32 stat,
                       u32* outI, u32* outJ, float* outR2, float* outD, u32 maxOut) {
  double r2, dp;
  if (!stats(n, pa, pb, x, &r2, &dp)) return;
  if ((stat == 0 ? r2 : dp) < thr) return;
  if (g_written < maxOut) { outI[g_written] = i; outJ[g_written] = j; outR2[g_written] = (float)r2; outD[g_written] = (float)dp; }
  g_written++;
}

// Banded pairs (j > i, pos[j] - pos[i] <= maxDist) for rows i in [iStart, iEnd), fast path (all rows use the same hn
// haplotypes). Cache-blocked and register-tiled: blocks of B i-rows; each j-row is loaded once per block and ANDed
// against 4 i-rows at a time. Keeps stat >= thr (stat 0 = r², 1 = D′); r²-bound pruning (stat 0) is lossless.
// Output order: (block, j, i). Returns the number written (min(g_written, maxOut)). Empty range → 0; B in [1, 1024].
EXPORT(band_blocked) u32 band_blocked(const u64* alt, const u32* pa, const i32* pos, i32 n, i32 W, u32 hn, i32 maxDist,
                                      double thr, i32 stat, i32 prune, u32* outI, u32* outJ, float* outR2, float* outD,
                                      u32 maxOut, i32 B, i32 iStart, i32 iEnd) {
  g_evaluated = g_pruned = g_written = 0;
  i32 cand[1024];
  if (B > 1024) B = 1024;
  if (B < 1) B = 1;
  if (iStart < 0) iStart = 0;
  if (iEnd > n) iEnd = n;
  if (iEnd <= iStart) return 0;
  const double cut = thr * (1 - 1e-9); // the bound equals r² for perfectly coupled pairs: allow for rounding
  for (i32 i0 = iStart; i0 < iEnd; i0 += B) {
    i32 i1 = i0 + B < iEnd ? i0 + B : iEnd;
    for (i32 j = i0 + 1; j < n && pos[j] - pos[i1 - 1] <= maxDist; j++) { // the block's last row reaches furthest
      i32 m = 0, iend = j < i1 ? j : i1;
      for (i32 i = i0; i < iend; i++) {
        if (pos[j] - pos[i] > maxDist) continue;
        if (prune && stat == 0 && r2_bound(pa[i], pa[j], hn) < cut) { g_pruned++; continue; }
        cand[m++] = i;
      }
      g_evaluated += m;
      const u64* aj = alt + (u64)j * W;
      i32 k = 0;
      for (; k + 4 <= m; k += 4) {
        u32 c[4];
        popcount_and4(alt + (u64)cand[k] * W, alt + (u64)cand[k + 1] * W, alt + (u64)cand[k + 2] * W, alt + (u64)cand[k + 3] * W, aj, W, c);
        for (int q = 0; q < 4; q++) put(cand[k + q], j, hn, pa[cand[k + q]], pa[j], c[q], thr, stat, outI, outJ, outR2, outD, maxOut);
      }
      for (; k < m; k++) put(cand[k], j, hn, pa[cand[k]], pa[j], popcount_and(alt + (u64)cand[k] * W, aj, W), thr, stat, outI, outJ, outR2, outD, maxOut);
    }
  }
  return g_written < maxOut ? g_written : maxOut;
}

// General path: counts over the per-pair valid intersection (chrX across the PAR boundary, '.' calls). Order (i, j).
// Rows with identical masks (same cls) need one popcount: n = nv[i], and the alt counts are pa (alt is within valid).
EXPORT(band_masked) u32 band_masked(const u64* alt, const u64* valid, const u32* pa, const u32* nv, const i32* cls,
                                    const i32* pos, i32 n, i32 W, i32 maxDist, double thr, i32 stat, u32* outI, u32* outJ,
                                    float* outR2, float* outD, u32 maxOut, i32 iStart, i32 iEnd) { // i in [iStart, iEnd), j > i
  g_evaluated = g_pruned = g_written = 0;
  for (i32 i = iStart; i < iEnd; i++) {
    const u64 *ai = alt + (u64)i * W, *vi = valid + (u64)i * W;
    for (i32 j = i + 1; j < n && pos[j] - pos[i] <= maxDist; j++) {
      const u64 *aj = alt + (u64)j * W, *vj = valid + (u64)j * W;
      g_evaluated++;
      if (cls[i] == cls[j]) put(i, j, nv[i], pa[i], pa[j], popcount_and(ai, aj, W), thr, stat, outI, outJ, outR2, outD, maxOut);
      else put(i, j, popcount_and(vi, vj, W), popcount_and(ai, vj, W), popcount_and(aj, vi, W), popcount_and(ai, aj, W),
               thr, stat, outI, outJ, outR2, outD, maxOut);
    }
  }
  return g_written < maxOut ? g_written : maxOut;
}

// Dense tile for matrix/pyramid views: rows [r0, r0+nr) × cols [c0, c0+nc), r² and D′ quantized to bytes
// (0..250 = value*250 rounded; 255 = NA). Fast path only.
EXPORT(tile) void tile(const u64* alt, const u32* pa, i32 W, u32 hn, i32 r0, i32 nr, i32 c0, i32 nc, uint8_t* outR2, uint8_t* outD) {
  for (i32 a = 0; a < nr; a++) {
    const u64* ai = alt + (u64)(r0 + a) * W;
    for (i32 b = 0; b < nc; b++) {
      double r2, dp; u32 k = (u32)a * nc + b;
      if (stats(hn, pa[r0 + a], pa[c0 + b], popcount_and(ai, alt + (u64)(c0 + b) * W, W), &r2, &dp)) {
        outR2[k] = (uint8_t)(r2 * 250 + 0.5); outD[k] = (uint8_t)(dp * 250 + 0.5);
      } else { outR2[k] = 255; outD[k] = 255; }
    }
  }
}

// Exact counts for a list of pairs (exact double stats, LDmatrix 3-decimal output): x[k] = |alt_i ∧ alt_j|
EXPORT(pair_counts) void pair_counts(const u64* alt, i32 W, const u32* I, const u32* J, u32 m, u32* x) {
  for (u32 k = 0; k < m; k++) x[k] = popcount_and(alt + (u64)I[k] * W, alt + (u64)J[k] * W, W);
}
