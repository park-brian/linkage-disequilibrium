// The release's LD score reference panels (ldscore/<pop>.parquet: CHR, SNP, BP, L2 in ldsc's (CHR, BP) order; M_5_50
// per chromosome in the key-value metadata) and HapMap3 SNP list, read as ldsc reads the original files:
//   parse.ldscore → concat chromosomes, sort by (CHR, BP), drop_duplicates(SNP, keep first)
// Only the SNPs the analysis needs are kept (first occurrence, file order = ldsc's merge order), so memory follows the
// sumstats, not the panel (afr has 17.4 M rows). The SNP count ldsc logs comes from ldscore/ldsc_counts.parquet
// (pipeline/release/extras.js): ldsc reads with sep=r'\s+', so a row with an empty SNP field shifts left and its BP
// becomes the SNP; those rows count but never merge (an rsID never equals a bare number), so they are skipped here.
import { parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { required } from "../../core/source.js";

const kvOf = metadata => Object.fromEntries((metadata.key_value_metadata ?? []).map(x => [x.key, x.value]));
const readAll = async (f, columns) => { let r = []; await parquetRead({ file: f.file, metadata: f.metadata, compressors, columns, onComplete: x => { r = x; } }); return r; };

/** ldscore/ldsc_counts.parquet → Map(pop → the SNP count ldsc logs for that panel) */
export async function ldscCounts(file) {
  const f = await required(file)("ldscore/ldsc_counts.parquet");
  return new Map((await readAll(f, ["pop", "rows"])).map(([p, n]) => [p, Number(n)]));
}

/**
 * @param {(rel: string) => Promise<{ file, metadata } | null>} file  a release opener (e.g. planner.file)
 * @param {string} pop  the panel name as LDlink lower-cases it (eur, eas, afr, …)
 * @param {Set<string>} need  the SNPs to keep (the analysis's sumstats)
 * → { order, l2: Map(SNP → L2), M: Σ M_5_50 } or null when there is no such panel
 */
export async function ldscPanel(file, pop, need) {
  const f = await file(`ldscore/${pop}.parquet`);
  if (!f) return null;
  const order = [], l2 = new Map();
  let rowStart = 0;
  for (const rg of f.metadata.row_groups) { // one row group at a time: bounded memory, rows aligned across columns
    const rowEnd = rowStart + Number(rg.num_rows);
    await parquetRead({ file: f.file, metadata: f.metadata, compressors, columns: ["SNP", "L2"], rowStart, rowEnd,
      onComplete: rows => { for (const [snp, v] of rows) if (need.has(snp) && !l2.has(snp)) { l2.set(snp, v); order.push(snp); } } });
    rowStart = rowEnd;
  }
  const M = Object.values(JSON.parse(kvOf(f.metadata).M_5_50)).reduce((a, v) => a + v[0], 0); // np.sum over chromosomes
  return { order, l2, M };
}

/** w_hm3.snplist (ldscore/w_hm3.snplist.parquet) → { SNP, A1, A2 } column arrays, file order */
export async function mergeAlleles(file) {
  const out = { SNP: [], A1: [], A2: [] };
  for (const [s, a, b] of await readAll(await required(file)("ldscore/w_hm3.snplist.parquet"), ["SNP", "A1", "A2"])) { out.SNP.push(s); out.A1.push(a); out.A2.push(b); }
  return out;
}
