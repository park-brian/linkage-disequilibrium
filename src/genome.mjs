// Chromosomes, genome keys and the three 1000G builds: the one table every module uses.
// A genome key is chromCode · 1e9 + position (INT64 in the Parquet files, BigInt here).

export const CHROMS = [...Array.from({ length: 22 }, (_, i) => String(i + 1)), "X", "Y"];
export const chromCode = c => (c === "X" ? 23 : c === "Y" ? 24 : c === "MT" ? 25 : +c);
export const chromName = c => (c === 23 ? "X" : c === 24 ? "Y" : c === 25 ? "MT" : String(c));
export const key = (chrom, pos) => BigInt(chromCode(chrom)) * 1_000_000_000n + BigInt(pos);
export const keyPos = k => Number(k) % 1e9; // exact: keys stay below 2^53 (Number or BigInt in)

/**
 * genotypes: genotype directory; positions: dbSNP/annotation position build (High Coverage uses GRCh38 positions);
 * dbsnpKey: column of the rsID index; title / hg: LDlink's labels; vcfPrefix: LDlink's 1000G_chr_prefix
 */
export const BUILDS = {
  grch37: { genotypes: "grch37", positions: "grch37", dbsnpKey: "key37", title: "GRCh37", hg: "hg19", vcfPrefix: "" },
  grch38: { genotypes: "grch38", positions: "grch38", dbsnpKey: "key38", title: "GRCh38", hg: "hg38", vcfPrefix: "" },
  grch38_high_coverage: { genotypes: "grch38_high_coverage", positions: "grch38", dbsnpKey: "key38", title: "GRCh38 High Coverage", hg: "hg38_HC", vcfPrefix: "chr" },
};
export const TITLE = Object.fromEntries(Object.entries(BUILDS).map(([b, v]) => [b, v.title]));
export const tierPath = (build, chrom, tier) => `genotypes/${BUILDS[build].genotypes}/chr${chrom}.${tier}.parquet`;

// ---- LDlink's input parsing (LDcommon.py) ----
// validsnp (LDcommon.py:370): snp_raw.lower().strip().split(), and the modules use snp[0], the first token
export const firstToken = s => String(s ?? "").toLowerCase().trim().split(/\s+/)[0] ?? "";
// LDcommon.py:104-110 — only 1-22, X, Y (and lower-case x, y) are recognized
const LOOKUP = Object.fromEntries([...CHROMS.map(c => [c, c]), ["x", "X"], ["y", "Y"]]);
/**
 * LDlink's chr:pos parsing (validsnp + coord.strip("chr").split(":") + int()): → { chro, pos: BigInt, key } or null.
 * key is null when the position can't be a dbSNP key (≤ 0 or ≥ 1e9).
 */
export function parseCoord(s) {
  const t = firstToken(s).replace(/^[chr]+|[chr]+$/g, "").split(":"); // coord.strip("chr").split(":")
  if (t.length <= 1) return null;
  const chro = LOOKUP[t[0]]; if (!chro) return null;
  if (!/^\s*[+-]?\d+(_\d+)*\s*$/.test(t[1])) return null; // int(temp_coord[1]) (Python allows _ separators)
  const pos = BigInt(t[1].trim().replace(/_/g, ""));
  if (pos <= 0n || pos >= 1_000_000_000n) return { chro, pos, key: null };
  return { chro, pos, key: BigInt(chromCode(chro)) * 1_000_000_000n + pos };
}
