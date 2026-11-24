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
  grch38_high_coverage: { genotypes: "grch38_high_coverage", positions: "grch38", dbsnpKey: "key38", title: "GRCh38 High Coverage",
    hg: "hg38_HC", vcfPrefix: "chr" },
};
export const TITLE = Object.fromEntries(Object.entries(BUILDS).map(([b, v]) => [b, v.title]));
export const tierPath = (build, chrom, tier) => `genotypes/${BUILDS[build].genotypes}/chr${chrom}.${tier}.parquet`;
