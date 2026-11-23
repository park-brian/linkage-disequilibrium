// munge_sumstats.py port (CBIIT/ldsc@ldsc313), for LDlink's call: --sumstats <file> --merge-alleles w_hm3.snplist
// (all other options at their defaults: --info-min 0.9, --maf-min 0.01, no --N/--signed-sumstats/--a1-inc/--ignore).
// Produces the same .sumstats text (before gzip) and the same log counts.
import { pyRound, pyRoundStr } from "../pyfmt.js";
import { VALID_SNPS, MATCH_ALLELES, erfc, normIsf, normPdf, median } from "./stats.js";

const DEFAULT_CNAMES = {
  SNP: "SNP", MARKERNAME: "SNP", SNPID: "SNP", RS: "SNP", RSID: "SNP", RS_NUMBER: "SNP", RS_NUMBERS: "SNP",
  NSTUDY: "NSTUDY", N_STUDY: "NSTUDY", NSTUDIES: "NSTUDY", N_STUDIES: "NSTUDY",
  P: "P", PVALUE: "P", P_VALUE: "P", PVAL: "P", P_VAL: "P", GC_PVALUE: "P",
  A1: "A1", ALLELE1: "A1", ALLELE_1: "A1", EFFECT_ALLELE: "A1", REFERENCE_ALLELE: "A1", INC_ALLELE: "A1", EA: "A1",
  A2: "A2", ALLELE2: "A2", ALLELE_2: "A2", OTHER_ALLELE: "A2", NON_EFFECT_ALLELE: "A2", DEC_ALLELE: "A2", NEA: "A2",
  N: "N", NCASE: "N_CAS", CASES_N: "N_CAS", N_CASE: "N_CAS", N_CASES: "N_CAS", N_CONTROLS: "N_CON", N_CAS: "N_CAS", N_CON: "N_CON",
  NCONTROL: "N_CON", CONTROLS_N: "N_CON", N_CONTROL: "N_CON", WEIGHT: "N",
  ZSCORE: "Z", "Z-SCORE": "Z", GC_ZSCORE: "Z", Z: "Z", OR: "OR", B: "BETA", BETA: "BETA", LOG_ODDS: "LOG_ODDS", EFFECTS: "BETA", EFFECT: "BETA",
  SIGNED_SUMSTAT: "SIGNED_SUMSTAT", INFO: "INFO", EAF: "FRQ", FRQ: "FRQ", MAF: "FRQ", FRQ_U: "FRQ", F_U: "FRQ",
};
const NULL_VALUES = { LOG_ODDS: 0, BETA: 0, OR: 1, Z: 0 };
// pandas read_csv default NA strings + na_values=['.', 'NA']
const NA = new Set(["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null", "."]);
const cleanHeader = h => h.toUpperCase().replace(/-/g, "_").replace(/\./g, "_").replace(/\n/g, "");

/** chi2.isf(P, 1) = (Q^-1(P/2))^2 with Q = upper normal tail; Newton in log space on the accurate erfc */
const LOG_SQRT_2PI = 0.5 * Math.log(2 * Math.PI);
/** log of the upper normal tail Q(x) = 0.5 erfc(x/√2): erfc while it is representable, the asymptotic series beyond */
function logQ(x) {
  if (x < 30) return Math.log(0.5 * erfc(x / Math.SQRT2));
  const z = 1 / (x * x); // Q(x) ~ φ(x)/x · (1 − z + 3z² − 15z³ + 105z⁴ − 945z⁵ + 10395z⁶)
  const s = 1 + z * (-1 + z * (3 + z * (-15 + z * (105 + z * (-945 + z * 10395)))));
  return -0.5 * x * x - Math.log(x) - LOG_SQRT_2PI + Math.log(s);
}
export function chi2Isf1(P) {
  if (!(P > 0)) return P === 0 ? Infinity : NaN;
  if (P >= 1) return 0;
  const lq = Math.log(P) - Math.LN2; // log(P/2), exact even when P/2 underflows
  let x = lq > -700 ? normIsf(P / 2) : Math.sqrt(-2 * lq - Math.log(-2 * lq) - 2 * LOG_SQRT_2PI);
  for (let i = 0; i < 50; i++) {
    // Newton on f(x) = log Q(x) − lq: f'(x) = −φ(x)/Q(x), so x += (log Q − lq) · Q/φ, with Q/φ = exp(log Q − log φ)
    const lQ = logQ(x), ratio = Math.exp(lQ + 0.5 * x * x + LOG_SQRT_2PI);
    const step = (lQ - lq) * ratio;
    x += step;
    if (Math.abs(step) < 1e-15 * Math.max(1, Math.abs(x))) break;
  }
  return x * x;
}
// pandas Series.quantile (linear interpolation)
// pandas Series.quantile → numpy percentile, method "linear", with numpy's _lerp: for t >= 0.5 it computes
// b - (b - a) * (1 - t) (not a + (b - a) * t), which rounds differently and gives NaN for [x, inf]
function quantile(a, q) {
  const s = Float64Array.from(a).sort(), n = s.length;
  if (!n) return NaN;
  const h = (n - 1) * q, lo = Math.floor(h), t = h - lo, x = s[lo], y = s[Math.min(lo + 1, n - 1)], d = y - x;
  return t >= 0.5 ? y - d * (1 - t) : x + d * t;
}
// Python repr(round(v, 2)) for messages: 0.12, 1.0, -inf, nan
const pyRepr2 = v => (Number.isNaN(v) ? "nan" : !Number.isFinite(v) ? (v > 0 ? "inf" : "-inf") : (r => (Number.isInteger(r) ? r.toFixed(1) : String(r)))(pyRound(v, 2)));
/** numpy's pairwise summation (pairwise_sum in numpy/core/src/umath/loops_utils), used by Series.mean() */
export function npSum(a, lo = 0, n = a.length) {
  if (n < 8) { let s = 0; for (let i = lo; i < lo + n; i++) s += a[i]; return s; } // numpy starts from -0.0; same result here
  if (n <= 128) {
    const r = [0, 1, 2, 3, 4, 5, 6, 7].map(k => a[lo + k]);
    let i = 8;
    for (; i < n - (n % 8); i += 8) for (let k = 0; k < 8; k++) r[k] += a[lo + i + k];
    let res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
    for (; i < n; i++) res += a[lo + i];
    return res;
  }
  let n2 = Math.floor(n / 2); n2 -= n2 % 8;
  return npSum(a, lo, n2) + npSum(a, lo + n2, n - n2);
}
// Python str.split() whitespace (read_header): NOT the BOM (U+FEFF), unlike JS \s / trim()
const PY_WS = new RegExp("[" + [[9, 13], [28, 32], [0x85, 0x85], [0xa0, 0xa0], [0x1680, 0x1680], [0x2000, 0x200a], [0x2028, 0x2029], [0x202f, 0x202f], [0x205f, 0x205f], [0x3000, 0x3000]]
  .map(([a, b]) => String.fromCodePoint(a) + "-" + String.fromCodePoint(b)).join("") + "]+", "u"); // str.isspace() code points
const pySplit = s => s.split(PY_WS).filter(Boolean);
// pandas sep=r'\s+' (C engine: delim_whitespace) splits data rows on spaces and tabs
const pdSplit = s => s.replace(/^[ \t]+|[ \t]+$/g, "").split(/[ \t]+/);
// pandas float parsing: decimal/exponent forms and inf/infinity; hex and other text make the column non-numeric
const pdNumber = v => /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/.test(v) ? Number(v) : /^[+-]?inf(inity)?$/i.test(v) ? (v[0] === "-" ? -Infinity : Infinity) : NaN;

/**
 * @param {string} text   raw summary statistics (whitespace-delimited, header row)
 * @param {{SNP: string[], A1: string[], A2: string[]}} mergeAlleles  w_hm3.snplist columns
 * @returns {{ text: string, log: object }}
 */
export function mungeSumstats(text, mergeAlleles, { infoMin = 0.9, mafMin = 0.01 } = {}) {
  const lines = text.split(/\r?\n/);
  const fileCnames = pySplit(lines[0]); // read_header: f.readline().split()
  const translation = {}; // file column → internal name (keys not cleaned)
  for (const c of fileCnames) { const k = cleanHeader(c); if (k in DEFAULT_CNAMES) translation[c] = DEFAULT_CNAMES[k]; }
  const signCols = Object.keys(translation).filter(c => translation[c] in NULL_VALUES);
  if (signCols.length > 1) throw new Error("Too many signed sumstat columns. Specify which to ignore with the --ignore flag.");
  if (!signCols.length) throw new Error("Could not find a signed summary statistic column.");
  const printed = { ...translation }; // munge_sumstats.py prints cname_translation before the signed column is renamed
  const signName = signCols[0], signNull = NULL_VALUES[translation[signName]];
  translation[signName] = "SIGNED_SUMSTAT";
  for (const c of ["SNP", "P", "SIGNED_SUMSTAT"]) if (!Object.values(translation).includes(c)) throw new Error(`Could not find ${c} column.`);
  // munge_sumstats.py:614-622 sits inside `for field in cname_translation:` and the loop body returns, so only the FIRST
  // mapped field is checked for repeats; the "different ... columns" check sees the count of the LAST mapped value
  const firstField = Object.keys(translation)[0], numk = fileCnames.filter(x => x === firstField).length;
  if (numk > 1) throw new Error(`Found ${numk} columns named ${firstField}`);
  { const vs = Object.values(translation), last = vs.at(-1), numc = vs.filter(v => v === last).length; if (numc > 1) throw new Error(`Found ${numc} different ${last} columns`); }
  const vals = Object.values(translation);
  if (!vals.includes("N") && !(vals.includes("N_CAS") && vals.includes("N_CON"))) throw new Error("Could not determine N.");
  if ((vals.includes("N") || (vals.includes("N_CAS") && vals.includes("N_CON"))) && vals.includes("NSTUDY")) for (const x of Object.keys(translation)) if (translation[x] === "NSTUDY") delete translation[x];
  if (!(Object.values(translation).includes("A1") && Object.values(translation).includes("A2"))) throw new Error("Could not find A1/A2 columns.");
  // usecols, file order; a repeated name is read once: pandas renames later copies to "X.1", which usecols skips
  const used = fileCnames.map((c, i) => [c, i]).filter(([c, i]) => c in translation && fileCnames.indexOf(c) === i);
  const numeric = new Set(["P", "N", "N_CAS", "N_CON", "Z", "OR", "BETA", "LOG_ODDS", "INFO", "FRQ", "SIGNED_SUMSTAT", "NSTUDY"]);
  const merge = new Set(mergeAlleles.SNP);
  const drops = { NA: 0, P: 0, INFO: 0, FRQ: 0, A: 0, MERGE: 0 }, bad = { INFO: 0, FRQ: 0 }; // bad: filter_info / filter_frq warnings
  let dat = [], tot = 0;
  for (let li = 1; li < lines.length; li++) {
    const line = lines[li]; if (!/[^ \t]/.test(line)) continue; // pandas skips blank lines
    tot++;
    const t = pdSplit(line), r = {};
    // extra fields beyond the header are ignored: pandas reads only the usecols positions (checked against the reference)
    let na = false;
    for (const [c, i] of used) {
      const v = t[i], name = translation[c];
      // dropna(subset=columns except one named exactly 'INFO' in the FILE): a lower-case `info` column is not exempt
      if (v === undefined || NA.has(v)) { r[name] = null; if (c !== "INFO") na = true; continue; }
      r[name] = numeric.has(name) ? pdNumber(v) : v;
      if (numeric.has(name) && Number.isNaN(r[name])) {
        // signed columns are read with dtype float64, so pandas fails at read time; others become object → wrong_types
        if (name === "SIGNED_SUMSTAT") throw new Error(`could not convert string to float: '${v}'`);
        throw new Error(`Columns ['${name}'] are expected to be numeric`);
      }
    }
    if (na) { drops.NA++; continue; }
    if (!merge.has(r.SNP)) { drops.MERGE++; continue; }
    if ("INFO" in r) { if (r.INFO > 1.5 || r.INFO < 0) bad.INFO++; if (!(r.INFO >= infoMin)) { drops.INFO++; continue; } }
    if ("FRQ" in r) { const out = r.FRQ < 0 || r.FRQ > 1, maf = Math.min(r.FRQ, 1 - r.FRQ); if (out) bad.FRQ++; if (!(maf > mafMin) || out) { drops.FRQ++; continue; } }
    if (!(r.P > 0 && r.P <= 1)) { drops.P++; continue; }
    r.A1 = r.A1.toUpperCase(); r.A2 = r.A2.toUpperCase();
    if (!VALID_SNPS.has(r.A1 + r.A2)) { drops.A++; continue; }
    dat.push(r);
  }
  if (!dat.length) throw new Error("After applying filters, no SNPs remain.");
  const seen = new Set(), before = dat.length; dat = dat.filter(r => (seen.has(r.SNP) ? false : (seen.add(r.SNP), true)));
  const dupDrops = before - dat.length;
  // process_n
  // dat column order (pandas): usecols in file order with internal names; process_n may append N
  let datCols = used.map(([c]) => translation[c]);
  if (datCols.includes("N_CAS") && datCols.includes("N_CON")) { // munge_sumstats.process_n: overwrites an existing N
    const Ns = dat.map(r => r.N_CAS + r.N_CON), Pc = dat.map((r, k) => r.N_CAS / Ns[k]);
    let Nmax = -Infinity; for (const v of Ns) if (v > Nmax) Nmax = v;
    const atMax = Pc.filter((_, k) => Ns[k] === Nmax), pbar = npSum(atMax) / atMax.length; // P[N == N.max()].mean()
    dat.forEach((r, k) => { r.N = (Ns[k] * Pc[k]) / pbar; delete r.N_CAS; delete r.N_CON; });
    datCols = datCols.filter(c => c !== "N_CAS" && c !== "N_CON");
    if (!datCols.includes("N")) datCols.push("N"); // new column → appended last
  }
  const nMin = quantile(dat.map(r => r.N), 0.9) / 1.5, beforeN = dat.length;
  dat = dat.filter(r => r.N >= nMin);
  // Z from P, signed by SIGNED_SUMSTAT
  const med = median(dat.map(r => r.SIGNED_SUMSTAT));
  if (Math.abs(med - signNull) > 0.1) throw new Error(`WARNING: median value of ${signName} is ${pyRepr2(med)} (should be close to ${signNull}). This column may be mislabeled.`);
  for (const r of dat) { r.Z = Math.sqrt(chi2Isf1(r.P)) * (r.SIGNED_SUMSTAT < signNull ? -1 : 1); }
  // allele_merge: left merge in --merge-alleles order; mismatched alleles → NA
  const bySnp = new Map(dat.map(r => [r.SNP, r]));
  let nMismatch = 0, nWithData = 0;
  const out = mergeAlleles.SNP.map((snp, k) => {
    const r = bySnp.get(snp);
    if (!r) return { SNP: snp };
    nWithData++;
    const ma = (mergeAlleles.A1[k] + mergeAlleles.A2[k]).toUpperCase();
    if (!MATCH_ALLELES.has(r.A1 + r.A2 + ma)) { nMismatch++; return { SNP: snp }; }
    return r;
  });
  if (nMismatch >= nWithData) throw new Error("All SNPs have alleles that do not match --merge-alleles.");
  // columns: SNP first (allele_merge left-merges dat onto the --merge-alleles frame), then dat's columns in order
  // (P renamed to Z in place) among N/Z/A1/A2
  const printCols = ["SNP", ...datCols.map(n => (n === "P" ? "Z" : n)).filter(n => ["N", "Z", "A1", "A2"].includes(n))];
  const fmt = (r, c) => { const v = r[c]; if (v === undefined || v === null) return ""; return c === "Z" || c === "N" ? pyRoundStr(v, 3) : String(v); };
  const body = out.map(r => printCols.map(c => fmt(r, c)).join("\t")).join("\n");
  const z2 = out.filter(r => r.Z !== undefined).map(r => r.Z * r.Z);
  return {
    text: printCols.join("\t") + "\n" + body + "\n",
    log: { read: tot, removed_merge: drops.MERGE, removed_na: drops.NA, removed_info: drops.INFO, removed_maf: drops.FRQ, removed_p: drops.P, removed_alleles: drops.A,
      remain: before, removed_dups: dupDrops, n_min: nMin, removed_n: beforeN - dat.length, median_signed: med, removed_mismatch: nMismatch,
      written: out.length, nonmissing: out.length - out.filter(r => r.Z === undefined).length, mean_chisq: npSum(z2) / z2.length, // pandas mean: numpy pairwise summation
      median_chisq: median(z2), max_chisq: z2.reduce((m, v) => (v > m ? v : m), -Infinity), n_gws: z2.filter(v => v > 29).length, cnames: printed, sign_name: signName, merge_snps: mergeAlleles.SNP.length,
      bad_info: bad.INFO, bad_frq: bad.FRQ, has_info: Object.values(translation).includes("INFO"), has_frq: Object.values(translation).includes("FRQ") },
  };
}
