// dbSNP reads over the release's index files (LDlink's 2023 mongodump, b155): every doc for an rsID or a position, in
// file order, which is the dump's order (docs/DATA.md "Order"). Choosing among them is the parity layer's job
// (src/parity/resolve.js: LDlink's find_one = the first). Absent index files (a partition or chromosome with no
// docs) read as no docs.
//   <root>/dbsnp/rsid/part-<rsid >> 26>.parquet       rsid · ref_rsid · chrom · key37 · key38
//   <root>/dbsnp/<grch37|grch38>/chr<N>.parquet        key · rsid · is_ref · type · function
//   <root>/dbsnp/duplicated_rsids.parquet              rsIDs with more than one doc
import { compressors } from "hyparquet-compressors";
import { sortedLookup } from "./parquet.js";
import { createSource, once } from "./source.js";
import { groupBy } from "./collect.js";

const partOf = id => Number(BigInt(id) >> 26n);
const partPath = p => `dbsnp/rsid/part-${String(p).padStart(2, "0")}.parquet`;

/**
 * createDbsnp(source) or createDbsnp(root, openFile, { exists }): source = { file(rel) → { file, metadata } | null }
 * (e.g. the planner's, sharing its cache). Without exists(), only ENOENT / 404 count as absent; other errors throw.
 */
export function createDbsnp(root, openFile, { exists } = {}) {
  const { file: open } = typeof root === "object" ? root : createSource(root, openFile, { exists });
  const read = async (path, column, keys, columns) => {
    const f = keys.length ? await open(path) : null;
    return f ? (await sortedLookup({ file: f.file, metadata: f.metadata, column, keys, columns, compressors })).rows : [];
  };
  /** rsIDs (BigInt) → Map(String(rsid) → [docs in dump order]); columns from rsid, ref_rsid, chrom, key37, key38 */
  async function byRsid(ids, columns = ["rsid", "ref_rsid", "chrom", "key37", "key38"]) {
    const parts = await Promise.all([...groupBy([...new Set(ids)], partOf)].map(([p, keys]) => read(partPath(p), "rsid", keys, columns)));
    return groupBy(parts.flat(), r => String(r.rsid));
  }
  /** genome keys (BigInt) on one chromosome of a build's position index (positions = "grch37" | "grch38", chrom a name)
   *  → Map(String(key) → [docs in dump order]); columns from key, rsid, is_ref, type, function */
  async function atKeys(positions, chrom, keys, columns = ["key", "rsid", "is_ref", "function"]) {
    return groupBy(await read(`dbsnp/${positions}/chr${chrom}.parquet`, "key", keys, columns), r => String(r.key));
  }
  /** rsIDs with more than one doc (read once): Set of digit strings */
  const duplicated = once(async () => {
    const f = await open("dbsnp/duplicated_rsids.parquet"); if (!f) return new Set();
    const { parquetRead } = await import("hyparquet"); let rows = [];
    await parquetRead({ file: f.file, metadata: f.metadata, columns: ["rsid"], compressors, onComplete: r => { rows = r; } });
    return new Set(rows.map(r => String(r[0])));
  });
  return { byRsid, atKeys, duplicated };
}
