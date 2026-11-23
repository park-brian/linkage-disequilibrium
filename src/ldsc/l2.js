// LD Score (ldsc --l2) for PLINK bed/bim/fam, ported line-for-line from CBIIT/ldsc@ldsc313 (LDlink's dependency):
//   ldscore/ldscore.py  getBlockLefts, PlinkBEDFile (__read__, __filter_snps_maf__, nextSNPs),
//                       __corSumVarBlocks__ (chunked block algorithm, c = --chunk-size 50), __l2_unbiased__
//   ldsc.py             ldscore(): window coordinates, whole-chromosome check, output (%.3f), M and M_5_50
// LDlink runs: ldsc.py --bfile <x> --l2 --ld-wind-cm|--ld-wind-kb <int> --out <x>  (server/ldscore ldsc_utils.py)
//
// usage (CLI): node ldsc/ldscore_l2.js <bfile prefix> <cm|kb> <window> <out prefix>
import { pyRoundStr } from "../pyfmt.js";
export function readBim(text) {
  return text.split(/\r?\n/).filter(l => l.trim()).map(l => {
    const t = l.trim().split(/\s+/);
    return { CHR: t[0], SNP: t[1], CM: Number(t[2]), BP: Number(t[3]), A1: t[4], A2: t[5] };
  });
}
export const readFamCount = text => text.split(/\r?\n/).filter(l => l.trim()).length;

/** Python's '%.3f' (exact value, ties to even, signed zero): see ../pyfmt.js */
export const pyFixed3 = x => pyRoundStr(x, 3);

/**
 * @param {Uint8Array} bed  .bed bytes (SNP-major)
 * @param {object[]} bim    from readBim
 * @param {number} n        individuals (.fam lines)
 * @param {{ windCm?: number, windKb?: number, mafMin?: number, chunkSize?: number }} opts
 */
export function ldscoreL2(bed, bim, n, { windCm, windKb, mafMin = 0, chunkSize = 50 } = {}) {
  const set = v => v != null && v !== 0 && v !== false; // Python truthiness of the window arguments (NaN is truthy)
  if ([windKb, windCm].filter(set).length !== 1) throw new Error("Must specify exactly one --ld-wind option");
  // ---- PlinkBEDFile.__read__ ----
  if (bed[0] !== 0x6c || bed[1] !== 0x1b) throw new Error("Magic number from Plink .bed file not recognized");
  if (bed[2] !== 0x01) throw new Error("Plink .bed file must be in default SNP-major mode");
  const m0 = bim.length, nru = n % 4 ? n + (4 - (n % 4)) : n, bytesPerSnp = nru / 4;
  if ((bed.length - 3) * 8 !== 2 * m0 * nru) throw new Error(`Plink .bed file has ${(bed.length - 3) * 8} bits, expected ${2 * m0 * nru}`);
  // little-endian bit pairs per individual: (bit0, bit1); ldsc code: 00→0, 01→1, 10→9 (missing), 11→2
  const code = (snp, i) => { const byte = bed[3 + snp * bytesPerSnp + (i >> 2)], s = (i & 3) * 2; const b0 = (byte >> s) & 1, b1 = (byte >> (s + 1)) & 1; return b0 === 0 ? (b1 === 0 ? 0 : 1) : (b1 === 0 ? 9 : 2); };
  // ---- __filter_snps_maf__ (counts over all nru slots; padding slots are 00) ----
  const kept = [], freq = [];
  for (let j = 0; j < m0; j++) {
    let a = 0, b = 0, c = 0;
    for (let i = 0; i < nru; i++) {
      const byte = bed[3 + j * bytesPerSnp + (i >> 2)], s = (i & 3) * 2, b0 = (byte >> s) & 1, b1 = (byte >> (s + 1)) & 1;
      a += b0; b += b1; c += b0 & b1;
    }
    const majorCt = b + c, nNomiss = n - a + c, f = nNomiss > 0 ? majorCt / (2 * nNomiss) : 0, hetMiss = a + b - 2 * c;
    if (Math.min(f, 1 - f) > mafMin && hetMiss < n) { kept.push(j); freq.push(f); }
  }
  const m = kept.length;
  if (!m) throw new Error("After filtering, no SNPs remain");
  const maf = freq.map(f => Math.min(f, 1 - f));
  // ---- nextSNPs: standardized genotypes, one Float64Array of length n per kept SNP ----
  const Y = kept.map(j => {
    const x = new Float64Array(n); let sum = 0, cnt = 0;
    for (let i = 0; i < n; i++) { const g = code(j, i); x[i] = g; if (g !== 9) { sum += g; cnt++; } }
    const avg = sum / cnt; // np.mean over non-missing
    for (let i = 0; i < n; i++) if (x[i] === 9) x[i] = avg;
    let mean = 0; for (let i = 0; i < n; i++) mean += x[i]; mean /= n;
    let v = 0; for (let i = 0; i < n; i++) v += (x[i] - mean) ** 2;
    let denom = Math.sqrt(v / n); if (denom === 0) denom = 1; // np.std (ddof=0)
    for (let i = 0; i < n; i++) x[i] = (x[i] - avg) / denom;
    return x;
  });
  // ---- window: getBlockLefts on cM, or BP with kb*1000 ----
  let coords, maxDist;
  if (set(windCm)) { coords = kept.map(j => bim[j].CM); maxDist = windCm; }
  else if (set(windKb)) { coords = kept.map(j => bim[j].BP); maxDist = windKb * 1000; }
  else throw new Error("Must specify exactly one --ld-wind option");
  const blockLeft = new Int32Array(m);
  for (let i = 0, j = 0; i < m; i++) { while (j < m && Math.abs(coords[j] - coords[i]) > maxDist) j++; blockLeft[i] = j; }
  if (blockLeft[m - 1] === 0)
    throw new Error("Do you really want to compute whole-chomosome LD Score? If so, set the --yes-really flag (warning: it will use a lot of time / memory)");
  // ---- __corSumVarBlocks__ with annot = ones, func = l2_unbiased ----
  const denomN = n > 2 ? n - 2 : n;
  const func = r => { const sq = r * r; return sq - (1 - sq) / denomN; };
  const dotScaled = (a, b) => { let s = 0; for (let i = 0; i < n; i++) s += a[i] * (b[i] / n); return s; }; // A.T @ (B / n)
  const corSum = new Float64Array(m);
  let c = chunkSize;
  const blockSizes = Array.from(blockLeft, (bl, i) => Math.ceil((i - bl) / c) * c);
  let bIdx = 0; while (bIdx < m && !(blockLeft[bIdx] > 0)) bIdx++;
  let b = bIdx < m ? bIdx : m;
  b = Math.ceil(b / c) * c;
  if (b > m) { c = 1; b = m; }
  let cur = 0; const getter = k => { const out = Y.slice(cur, cur + k); cur += k; return out; };
  let lA = 0, A = getter(b);
  for (let lB = 0; lB < b; lB += c) { // chunk inside of block
    const B = A.slice(lB, lB + c);
    for (let ia = 0; ia < A.length; ia++) { let rs = 0; for (const bv of B) rs += func(dotScaled(A[ia], bv)); corSum[lA + ia] += rs; }
  }
  const b0 = b, md = c * Math.floor(m / c), end = md !== m ? md + 1 : md;
  let B = null;
  for (let lB = b0; lB < end; lB += c) { // chunk to right of block
    const oldB = b; b = blockSizes[lB];
    if (lB > b0 && b > 0) { A = [...A.slice(oldB - b + c, oldB), ...B]; lA += oldB - b + c; }
    else if (lB === b0 && b > 0) { A = A.slice(b0 - b, b0); lA = b0 - b; }
    else if (b === 0) { A = []; lA = lB; }
    if (lB === md) c = m - md;
    B = getter(c);
    // rfuncAB (b x c): rows → cor_sum[lA..], columns → cor_sum[lB..]
    for (let ia = 0; ia < A.length; ia++) {
      let rs = 0;
      for (let jb = 0; jb < B.length; jb++) { const v = func(dotScaled(A[ia], B[jb])); rs += v; corSum[lB + jb] += v; }
      corSum[lA + ia] += rs;
    }
    // rfuncBB (c x c)
    for (let ib = 0; ib < B.length; ib++) { let rs = 0; for (const bv of B) rs += func(dotScaled(B[ib], bv)); corSum[lB + ib] += rs; }
  }
  const rows = kept.map((j, k) => ({ CHR: bim[j].CHR, SNP: bim[j].SNP, BP: bim[j].BP, CM: bim[j].CM, MAF: maf[k], L2: corSum[k] }));
  return { rows, M: m, M_5_50: maf.filter(x => x > 0.05).length };
}

/** Text outputs exactly as ldsc writes them (before gzip): .l2.ldscore, .l2.M, .l2.M_5_50 */
export function formatOutputs({ rows, M, M_5_50 }) {
  const ldscore = "CHR\tSNP\tBP\tL2\n" + rows.map(r => `${r.CHR}\t${r.SNP}\t${r.BP}\t${pyFixed3(r.L2)}`).join("\n") + "\n";
  return { ldscore, M: `${M}\n`, M_5_50: `${M_5_50}\n` };
}
