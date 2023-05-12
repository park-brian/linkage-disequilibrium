import { RemoteFile } from "generic-filehandle";
import { TabixIndexedFile } from "@gmod/tabix";
import VCF from "@gmod/vcf";
import { sources, populationSamples } from "./sources";

export async function loadTabixIndexedFile(url) {
  const tbiIndexed = new TabixIndexedFile({
    filehandle: new RemoteFile(url),
    tbiFilehandle: new RemoteFile(url + ".tbi"),
  });
  const headerText = await tbiIndexed.getHeader();
  const tbiVCFParser = new VCF({ header: headerText });

  console.log(headerText, tbiIndexed, tbiVCFParser);

  return {
    tbiIndexed,
    tbiVCFParser,
  };
}

/**
 * Retrieves the coordinates of a list of rsids from NCBI.
 * @param {string[]} rsids
 * @returns {Promise<{rsid: string, chromosome: string, grch37Position: number, grch38Position: number}[]}>} snp coordinates
 */
export async function getCoordinates(rsids) {
  const ids = rsids.map((rsid) => rsid.replace(/^rs/, "")).filter((rsid) => rsid.length && !isNaN(rsid));
  const endpoint = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi";
  const params = {
    db: "snp",
    retmode: "json",
    id: ids.join(","),
  };
  const response = await fetch(endpoint, {
    method: "POST",
    headers: {
      "content-type": "application/x-www-form-urlencoded",
    },
    body: new URLSearchParams(params),
  });
  const data = await response.json();
  const results = [];

  for (const uid of data.result.uids) {
    const result = data.result[uid];
    results.push({
      rsid: `rs${uid}`,
      chromosome: result.chr,
      grch37Position: +result.chrpos_prev_assm.split(":")[1],
      grch38Position: +result.chrpos.split(":")[1],
    });
  }

  return results;
}

/**
 * Calculates pairwise linkage disequilibrium between a list of SNPs.
 * @param {string[]} rsids
 * @param {string} genomeBuild
 */
export async function getLinkageDisequilibrium(rsids, populations, genomeBuild) {
  const snpCoordinates = await getCoordinates(rsids);
  validateSnps(snpCoordinates);

  const chromosome = snpCoordinates[0].chromosome;
  const source = sources[genomeBuild];
  const vcfUrl = source.url(chromosome);
  const assembly = source.assembly;
  const variants = await getVariants(snpCoordinates, assembly, vcfUrl, true);
  const samples = [...new Set(populations.map((population) => populationSamples[population]).flat())];

  console.log(variants, samples);
}

/**
 * Validates a list of SNPs.
 * @param {{ chromosome: string; grch37Position: number; grch38Position: number; }[]} snps
 */
export function validateSnps(snps) {
  const chromosomes = new Set(snps.map((snp) => snp.chromosome));

  if (chromosomes.size !== 1) {
    throw new Error("All SNPs must be on the same chromosome");
  }

  if (snps.length < 2) {
    throw new Error("At least two SNPs are required");
  }
}

/**
 * Retrieves a list of variants from a VCF file.
 * @param {{rsid: string, chromosome: string, grch37Position: number, grch38Position: number}[]} snps
 * @param {string} genomeBuild
 * @param {string} vcfUrl
 * @returns {Promise<any[]>} variants
 */
export async function getVariants(snps, assembly, vcfUrl, parallel = false) {
  const variants = [];
  const { tbiIndexed, tbiVCFParser } = await loadTabixIndexedFile(vcfUrl);

  const getVariant = async (snp) => {
    const position = snp[`${assembly}Position`];
    await tbiIndexed.getLines(snp.chromosome, position - 1, position, function (line, fileOffset) {
      const variant = tbiVCFParser.parseLine(line);
      if (variant.ID[0] === snp.rsid) {
        variants.push(variant);
      }
    });
    return variants[0];
  };

  if (parallel) {
    const promises = snps.map(getVariant);
    return Promise.all(promises);
  } else {
    const variants = [];
    for (const snp of snps) {
      variants.push(await getVariant(snp));
    }
    return variants;
  }
}

/**
 * Calculates pairwise linkage disequilibrium between two variants.
 * @param {number} p11 - frequency of haplotype carrying the two variants
 * @param {number} p1 - frequency of one variant at site 1
 * @param {number} p2 -  frequency of one variant at site 2
 * @returns {{dPrime: number, rSquared: number}} LD coefficients
 */
export function calculateLD(p11, p1, p2) {
  // frequencies of other variants at each site
  var q1 = 1 - p1;
  var q2 = 1 - p2;

  // coefficient of Linkage Disequilibrium
  const d = p11 - p1 * p2;
  const rSquared = Math.pow(d, 2) / (p1 * q1 * p2 * q2);

  // maximum possible coefficient of Linkage Disequilibrium
  const dMax = d < 0 ? Math.min(p1 * q2, p2 * q1) : Math.min(p1 * p2, q1 * q2);
  const dPrime = d / dMax;
  
  return { dPrime, rSquared };
}
