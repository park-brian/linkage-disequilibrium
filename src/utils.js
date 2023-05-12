import { RemoteFile } from "generic-filehandle";
import { TabixIndexedFile } from "@gmod/tabix";
import VCF from "@gmod/vcf";
// import chi2gof from "@stdlib/stats/chi2gof";
import { sources, populationSamples } from "./sources";

/**
 * Loads a remote tabix-indexed vcf file
 * @param {string} url
 * @returns {Promise<{tbiIndexed: any; tbiVCFParser: any}>} tabix indexed file and parser
 */
export async function loadTabixIndexedFile(url) {
  const tbiIndexed = new TabixIndexedFile({
    filehandle: new RemoteFile(url),
    tbiFilehandle: new RemoteFile(url + ".tbi"),
  });
  const headerText = await tbiIndexed.getHeader();
  const tbiVCFParser = new VCF({ header: headerText });
  return { tbiIndexed, tbiVCFParser };
}

/**
 * Retrieves the coordinates of a list of rsids from NCBI.
 * @param {string[]} rsids
 * @returns {Promise<{rsid: string, chromosome: string, grch37Position: number, grch38Position: number}[]}>} snp coordinates
 */
export async function getRsidCoordinates(rsids) {
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
  const coords = await getRsidCoordinates(rsids);
  validateSnps(coords);

  const chromosome = coords[0].chromosome;
  const source = sources[genomeBuild];
  const vcfUrl = source.url(chromosome);
  const assembly = source.assembly;
  const variants = await getVariants(coords, assembly, vcfUrl, false);
  const samples = [...new Set(populations.map((population) => populationSamples[population]).flat())];

  // calculate allele frequencies
  for (const variant of variants) {
    const ref = variant.REF;
    const alt = variant.ALT[0];
    variant.haplotypes = [];

    // collect haplotypes for each sample
    for (const sample of samples) {
      const genotype = variant.SAMPLES[sample]?.GT?.at(0);
      if (genotype) {
        // assume phased, bi-allelic variants
        switch (genotype) {
          case "0|0":
            variant.haplotypes.push([ref, ref]);
            break;
          case "0|1":
            variant.haplotypes.push([ref, alt]);
            break;
          case "1|0":
            variant.haplotypes.push([alt, ref]);
            break;
          case "1|1":
            variant.haplotypes.push([alt, alt]);
            break;
          default:
            throw new Error(`Unexpected genotype: ${genotype}`);
        }
      }
    }
  }

  // calculate pairwise linkage disequilibrium
  const ldResults = [];

  for (let i = 0; i < variants.length; i++) {
    for (let j = i; j < variants.length; j++) {
      const variant1 = variants[i];
      const variant2 = variants[j];
      const hap1 = variant1.haplotypes;
      const hap2 = variant2.haplotypes;

      const haplotypes = hap1
        .map((_, i) => [
          [hap1[i][0], hap2[i][0]],
          [hap1[i][1], hap2[i][1]],
        ])
        .flat();

      const haplotypeFrequencies = haplotypes
        .map((a) => a.join(""))
        .reduce(
          (acc, curr) => ({
            ...acc,
            [curr]: (acc[curr] || 0) + 1,
          }),
          {}
        );

      const [p1, p2, q1, q2] = Object.entries(haplotypeFrequencies)
        .sort((a, b) => a[0].localeCompare(b[0]))
        .map((a) => a[1]);

      const { dPrime, rSquared } = calculateLD(p1, p2, q1, q2);

      ldResults.push({
        variant1: variant1,
        variant2: variant2,
        dPrime,
        rSquared,
      });
    }
  }

  return asMatrix(
    ldResults,
    (record) => record.variant1.ID[0],
    (record) => record.dPrime,
    (record) => record.variant2.ID[0],
    (record) => record.rSquared
  );
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
 * Retrieves a list of variants from a VCF file given a list of snp coordinates.
 * @param {{rsid: string, chromosome: string, grch37Position: number, grch38Position: number}[]} snps
 * @param {string} genomeBuild
 * @param {string} vcfUrl
 * @returns {Promise<any[]>} variants
 */
export async function getVariants(snps, assembly, vcfUrl, parallel = false) {
  const { tbiIndexed, tbiVCFParser } = await loadTabixIndexedFile(vcfUrl);

  const getVariant = async (snp) => {
    const matches = [];
    const position = snp[`${assembly}Position`];
    const range = 10;
    await tbiIndexed.getLines(snp.chromosome, position - range, position + range, function (line, fileOffset) {
      const variant = tbiVCFParser.parseLine(line);
      // only add the variant if it matches the rsid (or exact position) and is bi-allelic
      if ((variant.ID?.[0] === snp.rsid || variant.POS === position) && variant.ALT.length === 1) {
        if (!variant.ID) variant.ID = [snp.rsid];
        matches.push(variant);
      }
    });
    return matches[0];
  };

  let variants = [];

  if (parallel) {
    variants = await Promise.all(snps.map(getVariant));
  } else {
    for (const snp of snps) {
      variants.push(await getVariant(snp));
    }
  }

  return variants.sort((a, b) => a.POS - b.POS);
}

/**
 * Calculates dPrime and rSquared between two variants.
 * @param {number} p1 - frequency of one variant at site 1
 * @param {number} p2 - frequency of one variant at site 2
 * @param {number} q1 - frequency of other variant at site 1
 * @param {number} q2 - frequency of other variant at site 2
 * @returns {{dPrime: number, rSquared: number}} LD coefficients
 */
function calculateLD(p1, p2, q1, q2) {
  // Calculate the numerator of the coefficient of linkage disequilibrium
  const d = p1 * q2 - p2 * q1;

  // Calculate the denominator of the coefficient of linkage disequilibrium
  const ms = (p1 + q1) * (p2 + q2) * (p1 + p2) * (q1 + q2);

  // Calculate rSquared, which is the squared coefficient of linkage disequilibrium
  const rSquared = d ** 2 / ms;

  // Calculate the maximum possible coefficient of linkage disequilibrium (dMax)
  const dMax = d < 0 ? Math.min((p1 + q1) * (p1 + p2), (p2 + q2) * (q1 + q2)) : Math.min((p1 + q1) * (q1 + q2), (p1 + p2) * (p2 + q2));

  // Calculate dPrime, which is the absolute value of d divided by dMax
  const dPrime = Math.abs(d / dMax);

  return {
    dPrime: isNaN(dPrime) ? 1 : dPrime,
    rSquared: isNaN(rSquared) ? 1 : rSquared,
  };
}

export function asMatrix(records, xKey, xValue, yKey, yValue) {
  const matrix = {};
  const columns = [];
  const rows = [];
  for (const record of records) {
    const x = xKey(record);
    const y = yKey(record);
    if (!columns.includes(x)) columns.push(x);
    if (!rows.includes(y)) rows.push(y);

    if (!matrix[x]) matrix[x] = {};
    if (!matrix[y]) matrix[y] = {};
    matrix[x][y] = record;
    matrix[y][x] = record;
  }
  return { columns, rows, matrix };
}

export function getTableOptions({ columns, matrix }) {
  const columnDefs = [
    {
      title: "id",
      content: (row) => row.id,
    },
  ].concat(
    columns.map((column) => ({
      title: column,
      content: (row) => `
        <div class="text-end">
          <div class="text-nowrap">D' = ${row[column].dPrime.toFixed(3)}</div>
          <div class="text-nowrap">R² = ${row[column].rSquared.toFixed(3)}</div>
        </div>
      `,
    }))
  );

  const data = Object.entries(matrix).map(([key, value]) => ({ id: key, ...value }));

  // generate data for export
  const header = columnDefs.map((column) => column.title);
  const dPrimeData = [header].concat(data.map((row) => header.map((column) => row[column]?.dPrime?.toFixed(3) || row[column])));
  const rSquaredData = [header].concat(data.map((row) => header.map((column) => row[column]?.rSquared?.toFixed(3) || row[column])));

  return {
    columns: columnDefs,
    data,
    dPrimeData,
    rSquaredData,
  };
}

export function createTable({ columns, data }) {
  const table = document.createElement("table");
  table.className = "table table-striped table-bordered table-hover";

  const thead = document.createElement("thead");
  const tbody = document.createElement("tbody");
  const headerRow = document.createElement("tr");
  const headerCells = columns.map((column) => {
    const cell = document.createElement("th");
    cell.textContent = column.title;
    return cell;
  });
  headerRow.append(...headerCells);
  thead.append(headerRow);
  table.append(thead);

  for (let rowIndex = 0; rowIndex < data.length; rowIndex++) {
    const row = data[rowIndex];
    const tr = document.createElement("tr");
    const cells = columns.map((column, columnIndex) => {
      const cell = document.createElement(columnIndex === 0 ? "th" : "td");
      cell.innerHTML = column.content(row, rowIndex, columnIndex);
      return cell;
    });
    tr.append(...cells);
    tbody.append(tr);
  }
  table.append(tbody);
  return table;
}

export function exportDelimitedTextFile(data, filename, delimiter = "\t") {
  const content = data.map((row) => row.join(delimiter)).join("\r\n");
  const blob = new Blob([content], { type: "text/plain;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}