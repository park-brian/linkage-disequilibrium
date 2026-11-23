// VCF genotype parsing shared by the genotype builder and its tests.
/**
 * Allele codes for every haplotype of a VCF record f (tab-split line): codes[2s], codes[2s+1] = allele index, -1 for
 * '.', -2 for an absent second haplotype (haploid GT). GT is the text before the first ':'; '|' or '/' separate.
 */
export function gtCodes(f, nSamples, nAlts) {
  const codes = new Int16Array(2 * nSamples);
  if (nAlts <= 9) { // single-digit indexes: read characters (fast path)
    for (let s = 0; s < nSamples; s++) {
      const g = f[9 + s], c0 = g.charCodeAt(0);
      codes[2 * s] = c0 === 46 ? -1 : c0 - 48;
      if (g.length === 1 || g.charCodeAt(1) === 58) { codes[2 * s + 1] = -2; continue; }
      const c1 = g.charCodeAt(2);
      codes[2 * s + 1] = c1 === 46 ? -1 : c1 - 48;
    }
    return codes;
  }
  for (let s = 0; s < nSamples; s++) {
    const g = f[9 + s], c = g.indexOf(":"), parts = (c < 0 ? g : g.slice(0, c)).split(/[|/]/);
    codes[2 * s] = parts[0] === "." ? -1 : +parts[0];
    codes[2 * s + 1] = parts.length < 2 ? -2 : parts[1] === "." ? -1 : +parts[1];
  }
  return codes;
}
