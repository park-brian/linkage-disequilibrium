// Shared numerics and allele sets of the LDscore family (ldsc's regressions.py / sumstats.py / munge_sumstats.py).
/** median of a numeric array (numpy.median) */
export function median(a) { const s = Float64Array.from(a).sort(), n = s.length; return n % 2 ? s[(n - 1) >> 1] : (s[n / 2 - 1] + s[n / 2]) / 2; }
/** erfc with ~1e-16 relative accuracy (W. J. Cody's rational Chebyshev approximations, as in many libms). */
export function erfc(x) {
  if (x === Infinity) return 0;
  if (x === -Infinity) return 2;
  const ax = Math.abs(x);
  let r;
  if (ax < 0.5) {
    const t = x * x;
    const p = [3.16112374387056560e00, 1.13864154151050156e02, 3.77485237685302021e02, 3.20937758913846947e03, 1.85777706184603153e-1];
    const q = [2.36012909523441209e01, 2.44024637934444173e02, 1.28261652607737228e03, 2.84423683343917062e03];
    let num = p[4] * t, den = t;
    for (let i = 0; i < 3; i++) { num = (num + p[i]) * t; den = (den + q[i]) * t; }
    return 1 - x * (num + p[3]) / (den + q[3]);
  } else if (ax < 4) {
    const c = [5.64188496988670089e-1, 8.88314979438837594e00, 6.61191906371416295e01, 2.98635138197400131e02, 8.81952221241769090e02, 1.71204761263407058e03, 2.05107837782607147e03, 1.23033935479799725e03, 2.15311535474403846e-8];
    const d = [1.57449261107098347e01, 1.17693950891312499e02, 5.37181101862009858e02, 1.62138957456669019e03, 3.29079923573345963e03, 4.36261909014324716e03, 3.43936767414372164e03, 1.23033935480374942e03];
    let num = c[8] * ax, den = ax;
    for (let i = 0; i < 7; i++) { num = (num + c[i]) * ax; den = (den + d[i]) * ax; }
    r = (num + c[7]) / (den + d[7]);
  } else {
    const p = [3.05326634961232344e-1, 3.60344899949804439e-1, 1.25781726111229246e-1, 1.60837851487422766e-2, 6.58749161529837803e-4, 1.63153871373020978e-2];
    const q = [2.56852019228982242e00, 1.87295284992346725e00, 5.27905102951428412e-1, 6.05183413124413191e-2, 2.33520497626869185e-3];
    const z = 1 / (ax * ax);
    let num = p[5] * z, den = z;
    for (let i = 0; i < 4; i++) { num = (num + p[i]) * z; den = (den + q[i]) * z; }
    r = z * (num + p[4]) / (den + q[4]);
    r = (0.5641895835477563 - r) / ax;
  }
  // exp(-x^2) computed in two parts to keep precision (as in Cody)
  const xs = Math.trunc(ax * 16) / 16, del = (ax - xs) * (ax + xs);
  r = Math.exp(-xs * xs) * Math.exp(-del) * r;
  return x < 0 ? 2 - r : r;
}
/** chi2.sf(z^2, 1) = erfc(|z| / sqrt(2)) */
export const chi2sf1 = z2 => erfc(Math.sqrt(z2) / Math.SQRT2);
export const normPdf = x => Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);
/** norm.isf(p): Acklam's approximation refined by Newton steps on the accurate erfc */
export function normIsf(p) {
  const q0 = Math.min(p, 1 - p), t = Math.sqrt(-2 * Math.log(q0));
  let x = t - (2.515517 + 0.802853 * t + 0.010328 * t * t) / (1 + 1.432788 * t + 0.189269 * t * t + 0.001308 * t * t * t);
  if (p > 0.5) x = -x;
  for (let i = 0; i < 4; i++) { const f = 0.5 * erfc(x / Math.SQRT2) - p; x += f / normPdf(x); }
  return x;
}

// ---------- alleles (sumstats.py) ----------
const COMP = { A: "T", T: "A", C: "G", G: "C" };
const BASES = ["A", "T", "C", "G"];
const VALID = new Set(BASES.flatMap(a => BASES.map(b => a + b)).filter(x => x[0] !== x[1] && x[0] !== COMP[x[1]]));
const MATCH = new Set(), FLIP = new Map();
for (const a of VALID) for (const b of VALID) {
  const x = a + b;
  const ok = (x[0] === x[2] && x[1] === x[3]) || (x[0] === COMP[x[2]] && x[1] === COMP[x[3]]) || (x[0] === x[3] && x[1] === x[2]) || (x[0] === COMP[x[3]] && x[1] === COMP[x[2]]);
  if (ok) { MATCH.add(x); FLIP.set(x, (x[0] === x[3] && x[1] === x[2]) || (x[0] === COMP[x[3]] && x[1] === COMP[x[2]])); }
}
export { VALID as VALID_SNPS, MATCH as MATCH_ALLELES, FLIP };
