// LD Score regression (h2, rg) ported from CBIIT/ldsc@ldsc313 (LDlink's dependency), for the configuration LDlink
// uses: one LD Score column, --ref-ld-chr == --w-ld-chr, free intercepts, default --two-step 30, --n-blocks 200,
// M = sum of .l2.M_5_50 (ldsc.py --h2 / --rg; server ldsc_utils.run_herit_command / run_correlation_command).
//   ldscore/sumstats.py   estimate_h2, estimate_rg, merges, allele matching
//   ldscore/regressions.py LD_Score_Regression (two-step), Hsq, Gencov, RG, h2_obs_to_liab, p_z_norm, summaries
//   ldscore/irwls.py      IRWLS (2 reweighting rounds, then weighted fast block jackknife)
//   ldscore/jackknife.py  LstsqJackknifeFast, RatioJackknife, pseudovalues, np.cov(ddof=1)
// Vectors are plain Float64Arrays; the design has p <= 2 columns, so the linear algebra is solved directly.

// ---------- numerics ----------
import { erfc, chi2sf1, normPdf, normIsf, median, VALID_SNPS as VALID, MATCH_ALLELES as MATCH, FLIP } from "../core/stats.js";
const mean = a => { let s = 0; for (const v of a) s += v; return s / a.length; };
export function h2ObsToLiab(h2, P, K) {
  if (Number.isNaN(P) && Number.isNaN(K)) return h2;
  if (K <= 0 || K >= 1) throw new Error("K must be in the range (0,1)");
  if (P <= 0 || P >= 1) throw new Error("P must be in the range (0,1)");
  const thresh = normIsf(K);
  return h2 * (K ** 2 * (1 - K) ** 2 / (P * (1 - P) * normPdf(thresh) ** 2));
}
const pZNorm = (est, se) => { const Z = se === 0 ? Infinity : est / se; return { p: chi2sf1(Z * Z), z: Z }; };

// small dense linear algebra for p <= 2
function solve(A, b) { // A: p x p (array of rows), b: p
  const p = b.length;
  if (p === 1) return [b[0] / A[0][0]];
  // LU with partial pivoting (as LAPACK gesv)
  let [r0, r1] = Math.abs(A[1][0]) > Math.abs(A[0][0]) ? [1, 0] : [0, 1];
  const l = A[r1][0] / A[r0][0], u11 = A[r1][1] - l * A[r0][1], y1 = b[r1] - l * b[r0];
  const x1 = y1 / u11, x0 = (b[r0] - A[r0][1] * x1) / A[r0][0];
  return [x0, x1];
}
// np.floor(np.linspace(0, N, nBlocks + 1)): linspace computes i * (N / nBlocks) in floating point (not exactly) and sets
// the last point to N, so a boundary can land just below an integer (e.g. N = 228: blocks 50, 100, 150 shift by 1)
export const getSeparators = (N, nBlocks) => { const step = N / nBlocks; return Array.from({ length: nBlocks + 1 }, (_, i) => (i === nBlocks ? N : Math.floor(i * step))); };

/** x: array of p Float64Arrays (columns), y: Float64Array */
function lstsqJackknifeFast(x, y, nBlocks, separators) {
  const N = y.length, p = x.length;
  const s = separators ?? getSeparators(N, nBlocks), nb = s.length - 1;
  if (nb > N) throw new Error("More blocks than data points.");
  const xty = [], xtx = [];
  for (let i = 0; i < nb; i++) {
    const a = new Array(p).fill(0), B = Array.from({ length: p }, () => new Array(p).fill(0));
    for (let k = s[i]; k < s[i + 1]; k++) for (let u = 0; u < p; u++) { a[u] += x[u][k] * y[k]; for (let v = 0; v < p; v++) B[u][v] += x[u][k] * x[v][k]; }
    xty.push(a); xtx.push(B);
  }
  const xtyTot = new Array(p).fill(0), xtxTot = Array.from({ length: p }, () => new Array(p).fill(0));
  for (let i = 0; i < nb; i++) for (let u = 0; u < p; u++) { xtyTot[u] += xty[i][u]; for (let v = 0; v < p; v++) xtxTot[u][v] += xtx[i][u][v]; }
  const est = solve(xtxTot, xtyTot);
  const del = [];
  for (let j = 0; j < nb; j++) del.push(solve(xtxTot.map((r, u) => r.map((v, w) => v - xtx[j][u][w])), xtyTot.map((v, u) => v - xty[j][u])));
  return { est, delete_values: del, separators: s, ...jknife(pseudovalues(del, est)) };
}
const pseudovalues = (del, est) => del.map(d => d.map((v, u) => del.length * est[u] - (del.length - 1) * v));
function jknife(pv) {
  const nb = pv.length, p = pv[0].length;
  const m = Array.from({ length: p }, (_, u) => mean(pv.map(r => r[u])));
  const cov = Array.from({ length: p }, () => new Array(p).fill(0));
  for (let u = 0; u < p; u++) for (let v = 0; v < p; v++) { let s = 0; for (const r of pv) s += (r[u] - m[u]) * (r[v] - m[v]); cov[u][v] = s / (nb - 1) / nb; }
  return { jknife_est: m, jknife_cov: cov, jknife_var: cov.map((r, u) => r[u]), jknife_se: cov.map((r, u) => Math.sqrt(r[u])) };
}
function weightCols(x, w) { // IRWLS._weight: w normalized to sum 1
  let sw = 0; for (const v of w) { if (!(v > 0)) throw new Error("Weights must be > 0"); sw += v; }
  return x.map(col => col.map((v, k) => v * (w[k] / sw)));
}
function wls(x, y, w) { // np.linalg.lstsq on weighted data; full-rank p<=2 → normal equations
  const xw = weightCols(x, w), yw = weightCols([y], w)[0], p = x.length;
  const A = Array.from({ length: p }, (_, u) => Array.from({ length: p }, (_, v) => { let s = 0; for (let k = 0; k < y.length; k++) s += xw[u][k] * xw[v][k]; return s; }));
  const b = Array.from({ length: p }, (_, u) => { let s = 0; for (let k = 0; k < y.length; k++) s += xw[u][k] * yw[k]; return s; });
  return solve(A, b);
}
function irwls(x, y, updateFunc, nBlocks, w, separators) {
  w = w.map(Math.sqrt);
  for (let i = 0; i < 2; i++) w = updateFunc(wls(x, y, w)).map(Math.sqrt);
  return lstsqJackknifeFast(weightCols(x, w), weightCols([y], w)[0], nBlocks, separators);
}

// ---------- regression weights ----------
function hsqWeights(ld, wLd, N, M, hsq, intercept = 1) {
  hsq = Math.min(Math.max(hsq, 0), 1);
  return ld.map((l, k) => { const L = Math.max(l, 1), W = Math.max(wLd[k], 1), c = (hsq * N[k]) / M; return (1 / (2 * (intercept + c * L) ** 2)) * (1 / W); });
}
function gencovWeights(ld, wLd, N1, N2, M, h1, h2, rho, igc = 0, ih1 = 1, ih2 = 1) {
  h1 = Math.min(Math.max(h1, 0), 1); h2 = Math.min(Math.max(h2, 0), 1); rho = Math.max(Math.min(rho, 1), -1);
  return ld.map((l, k) => {
    const L = Math.max(l, 1), W = Math.max(wLd[k], 1);
    const a = (N1[k] * (h1 * L)) / M + ih1, b = (N2[k] * (h2 * L)) / M + ih2, c = (Math.sqrt(N1[k] * N2[k]) * (rho * L)) / M + igc;
    return (1 / (a * b + c * c)) * (1 / W);
  });
}

/**
 * LD_Score_Regression with n_annot = 1, free intercept, optional two-step.
 * kind: "hsq" | "gencov"; ctx carries N1/N2/hsq1/hsq2/intercepts for gencov.
 */
function ldScoreRegression(kind, y, x0, w, N, M, nBlocks, step1, ctx = {}) {
  const n = y.length, Mtot = M;
  const nullInt = kind === "hsq" ? 1 : 0;
  const totAgg = (Mtot * (mean(y) - nullInt)) / mean(x0.map((v, k) => v * N[k])); // aggregate()
  const initialW = kind === "hsq" ? hsqWeights(x0, w, N, Mtot, totAgg, 1)
    : gencovWeights(x0, w, ctx.N1, ctx.N2, Mtot, ctx.hsq1, ctx.hsq2, totAgg, 0, ctx.ih1, ctx.ih2);
  const Nbar = mean(N);
  const xs = x0.map((v, k) => (N[k] * v) / Nbar); // x = N * x / Nbar
  const ones = new Float64Array(n).fill(1);
  const upd = (coef, ldCol, wLd, Nsub, intercept, ii) => {
    const est = (Mtot * coef[0]) / Nbar;
    if (kind === "hsq") return hsqWeights(ldCol, wLd, Nsub, Mtot, est, intercept ?? Math.max(coef[1]));
    const N1 = ii ? ctx.N1.filter((_, k) => ii[k]) : ctx.N1, N2 = ii ? ctx.N2.filter((_, k) => ii[k]) : ctx.N2;
    return gencovWeights(ldCol, wLd, N1, N2, Mtot, ctx.hsq1, ctx.hsq2, est, intercept ?? coef[1], ctx.ih1, ctx.ih2);
  };
  let jk, twostepFiltered = null;
  if (step1) {
    const idx = []; step1.forEach((b, k) => { if (b) idx.push(k); });
    twostepFiltered = n - idx.length;
    const sub = a => Float64Array.from(idx, k => a[k]);
    const x1 = [sub(xs), sub(ones)], yp1 = sub(y), w1 = sub(w), N1 = sub(N), iw1 = sub(initialW);
    // step 1: free intercept on SNPs passing the two-step filter (note: ld column = the N-scaled x1, as in ldsc)
    const j1 = irwls(x1, yp1, a => upd(a, x1[0], w1, N1, undefined, step1), nBlocks, iw1);
    const step1Int = j1.est[1];
    const yp = y.map(v => v - step1Int);
    const sep = [0, ...j1.separators.slice(1, -1).map(s => idx[s]), n]; // update_separators
    const j2 = irwls([xs], yp, a => upd(a, x0, w, N, step1Int), nBlocks, initialW, sep);
    let num = 0, den = 0; for (let k = 0; k < n; k++) { num += initialW[k] * xs[k]; den += initialW[k] * xs[k] * xs[k]; }
    const c = num / den;
    const est = [j2.est[0], step1Int];
    const del = j1.delete_values.map((d1, b) => [j2.delete_values[b][0] - c * (d1[1] - step1Int), d1[1]]);
    jk = { est, delete_values: del, ...jknife(pseudovalues(del, est)) };
  } else {
    jk = irwls([xs, ones], y, a => upd(a, x0, w, N, undefined), nBlocks, initialW);
  }
  const coef = jk.est[0] / Nbar, coefCov = jk.jknife_cov[0][0] / Nbar ** 2;
  const tot = M * coef, totSe = Math.sqrt(M * M * coefCov);
  return { tot, tot_se: totSe, intercept: jk.est[1], intercept_se: jk.jknife_se[1], twostep_filtered: twostepFiltered,
    tot_delete_values: jk.delete_values.map(d => (d[0] * M) / Nbar), jk };
}

export function hsq(chisq, ld, wLd, N, M, { nBlocks = 200, twostep = 30 } = {}) {
  const step1 = twostep != null ? Array.from(chisq, v => v < twostep) : null;
  const r = ldScoreRegression("hsq", chisq, ld, wLd, N, M, nBlocks, step1);
  r.mean_chisq = mean(chisq); r.lambda_gc = median(chisq) / 0.4549;
  if (r.mean_chisq > 1) { r.ratio = (r.intercept - 1) / (r.mean_chisq - 1); r.ratio_se = r.intercept_se / (r.mean_chisq - 1); } else r.ratio = r.ratio_se = "NA";
  return r;
}
export function rg(z1, z2, ld, wLd, N1, N2, M, { nBlocks = 200, twostep = 30 } = {}) {
  const h1 = hsq(z1.map(v => v * v), ld, wLd, N1, M, { nBlocks, twostep }), h2 = hsq(z2.map(v => v * v), ld, wLd, N2, M, { nBlocks, twostep });
  const y = z1.map((v, k) => v * z2[k]);
  const step1 = twostep != null ? Array.from(z1, (v, k) => v * v < twostep && z2[k] * z2[k] < twostep) : null;
  const gc = ldScoreRegression("gencov", y, ld, wLd, N1.map((v, k) => Math.sqrt(v * N2[k])), M, nBlocks, step1,
    { N1, N2, hsq1: h1.tot, hsq2: h2.tot, ih1: h1.intercept, ih2: h2.intercept });
  Object.assign(gc, pZNorm(gc.tot, gc.tot_se), { mean_z1z2: mean(y) });
  const out = { hsq1: h1, hsq2: h2, gencov: gc };
  if (h1.tot <= 0 || h2.tot <= 0) return { ...out, rg_ratio: "NA", rg_se: "NA", p: "NA", z: "NA", negative_hsq: true };
  const ratio = gc.tot / Math.sqrt(h1.tot * h2.tot);
  const denom = h1.tot_delete_values.map((d, b) => Math.sqrt(d * h2.tot_delete_values[b])); // NaN if negative
  const nb = denom.length, pv = gc.tot_delete_values.map((num, b) => [nb * ratio - ((nb - 1) * num) / denom[b]]);
  const j = jknife(pv);
  return { ...out, rg_ratio: ratio, rg_se: j.jknife_se[0], rg_jknife: j.jknife_est[0], ...pZNorm(ratio, j.jknife_se[0]) };
}

// ---------- data loading with ldsc's merge semantics ----------
export function parseSumstats(text, { alleles = false, dropna = true } = {}) { // ps.sumstats + drop_duplicates(SNP)
  const lines = text.split(/\r?\n/), h = lines[0].split("\t"), ix = Object.fromEntries(h.map((c, i) => [c, i]));
  const rows = [], seen = new Set();
  for (let i = 1; i < lines.length; i++) {
    if (!lines[i]) continue;
    const t = lines[i].split("\t");
    const num = s => (s === "" || s === undefined || s === "NA" || s === "nan" ? NaN : Number(s));
    const r = { SNP: t[ix.SNP], Z: num(t[ix.Z]), N: num(t[ix.N]) };
    if (alleles) { r.A1 = t[ix.A1] || null; r.A2 = t[ix.A2] || null; }
    if (dropna && (Number.isNaN(r.Z) || Number.isNaN(r.N) || (alleles && (!r.A1 || !r.A2)))) continue;
    if (seen.has(r.SNP)) continue; seen.add(r.SNP); rows.push(r);
  }
  return rows;
}
/** ref panel rows (already sorted by CHR, BP): drop_duplicates(SNP, keep first) → Map SNP → L2, and ordered list */
export function dedupPanel(snps, l2) {
  const order = [], m = new Map();
  for (let i = 0; i < snps.length; i++) if (!m.has(snps[i])) { m.set(snps[i], l2[i]); order.push(snps[i]); }
  return { order, l2: m };
}


/** estimate_h2 (n_annot = 1) */
export function estimateH2(sumstatsText, panel, M, opts = {}) {
  const ss = new Map(parseSumstats(sumstatsText).map(r => [r.SNP, r]));
  const rows = panel.order.filter(s => ss.has(s)).map(s => ss.get(s)); // merge(ref_ld, sumstats) then merge(., w_ld)
  const ld = Float64Array.from(rows, r => panel.l2.get(r.SNP));
  const chisq = Float64Array.from(rows, r => r.Z * r.Z), N = Float64Array.from(rows, r => r.N);
  const nBlocks = Math.min(rows.length, opts.nBlocks ?? 200);
  return { n_snp: rows.length, ...hsq(chisq, ld, ld, N, M, { nBlocks, twostep: 30 }) };
}
/** estimate_rg (two phenotypes, n_annot = 1) */
export function estimateRg(text1, text2, panel, M, opts = {}) {
  const s1 = new Map(parseSumstats(text1, { alleles: true, dropna: true }).map(r => [r.SNP, r]));
  const merged = panel.order.filter(s => s1.has(s)).map(s => s1.get(s));
  const s2 = new Map(parseSumstats(text2, { alleles: true, dropna: false }).map(r => [r.SNP, r]));
  const rows = [];
  for (const r of merged) {
    const q = s2.get(r.SNP); if (!q) continue; // inner merge, left order
    if (Number.isNaN(q.Z) || Number.isNaN(q.N) || !q.A1 || !q.A2) continue; // dropna(how='any')
    const al = r.A1 + r.A2 + q.A1 + q.A2;
    if (!MATCH.has(al)) continue; // _filter_alleles
    rows.push({ ...r, Z2: FLIP.get(al) ? -q.Z : q.Z, N2: q.N });
  }
  const ld = Float64Array.from(rows, r => panel.l2.get(r.SNP));
  const nBlocks = Math.min(opts.nBlocks ?? 200, rows.length);
  // ldsc quirk (sumstats.estimate_rg): args.intercept_h2 is already a list ([None, None]) when it checks
  // `args.intercept_h2 is None`, so the default --two-step 30 is never applied to rg. Reproduce: no two-step.
  return { n_snp: rows.length, ...rg(Float64Array.from(rows, r => r.Z), Float64Array.from(rows, r => r.Z2), ld, ld,
    Float64Array.from(rows, r => r.N), Float64Array.from(rows, r => r.N2), M, { nBlocks, twostep: opts.twostep ?? null }) };
}
