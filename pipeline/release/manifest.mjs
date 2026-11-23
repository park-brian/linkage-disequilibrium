// Build a manifest for a release directory: one JSON describing every Parquet file so a browser client can
//  - fetch each footer with one exact Range request (no 512 KB guess), or skip it,
//  - pick files / row groups by key range without opening them,
//  - plan point reads vs scans from page-index locations (PLAN_REVIEW.md findings 2 and 7).
// usage: node build_manifest.mjs <releaseRoot>     → <releaseRoot>/manifest.json
import { readdirSync, statSync, writeFileSync } from "node:fs";
import { parquetMetadataAsync, asyncBufferFromFile } from "hyparquet";
import { useDataCwd } from "../../tools/config.mjs";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.mjs)

const [root] = process.argv.slice(2);
const KEY = name => {
  if (name.startsWith("genotypes/")) return "pos";
  if (name.startsWith("dbsnp/rsid/") || /(duplicated_rsids|forge_rsid)\.parquet$/.test(name)) return "rsid";
  if (name.endsWith("gene_names.parquet")) return "name";
  if (name.endsWith("gtex_genes.parquet")) return "gene_id_v26";
  if (/(gtex_tissues|platforms|w_hm3\.snplist)\.parquet$/.test(name)) return null; // small lookup tables, read whole
  if (name.startsWith("ldscore/")) return "CHR"; // sorted by (CHR, BP); a client selects chromosomes by row group
  return "key";
};
const PAYLOAD = name => (name.startsWith("genotypes/") ? "haps" : null);
const walk = d => readdirSync(`${root}/${d}`, { withFileTypes: true }).flatMap(e =>
  e.isDirectory() ? walk(`${d}/${e.name}`) : e.name.endsWith(".parquet") ? [`${d}/${e.name}`] : []);
const num = v => (typeof v === "bigint" ? (v <= BigInt(Number.MAX_SAFE_INTEGER) ? Number(v) : String(v)) : v);

const files = {};
for (const name of readdirSync(root, { withFileTypes: true }).filter(e => e.isDirectory()).flatMap(e => walk(e.name)).sort()) {
  const path = `${root}/${name}`, byteLength = statSync(path).size;
  const m = await parquetMetadataAsync(await asyncBufferFromFile(path));
  const cols = m.row_groups[0]?.columns.map(c => c.meta_data.path_in_schema.join(".")) ?? [];
  const ki = cols.indexOf(KEY(name)), pi = PAYLOAD(name) ? cols.indexOf(PAYLOAD(name)) : -1;
  let first = 0;
  files[name] = {
    byteLength, numRows: Number(m.num_rows),
    footer: { offset: byteLength - 8 - m.metadata_length, length: m.metadata_length + 8 },
    key: KEY(name), columns: cols,
    kv: Object.fromEntries((m.key_value_metadata ?? []).filter(kv => kv.key !== "samples").map(kv => [kv.key, kv.value])),
    rowGroups: m.row_groups.map(rg => {
      const k = ki >= 0 ? rg.columns[ki] : null, st = k?.meta_data.statistics ?? {}, p = pi >= 0 ? rg.columns[pi] : null;
      // compressed bytes per column (+ total) so a planner can cost scans without opening the file
      const colBytes = Object.fromEntries(rg.columns.map(c => [c.meta_data.path_in_schema.join("."), Number(c.meta_data.total_compressed_size)]));
      const o = { firstRow: first, numRows: Number(rg.num_rows), min: num(st.min_value), max: num(st.max_value),
        bytes: Object.values(colBytes).reduce((a, b) => a + b, 0), ...(name.startsWith("genotypes/") ? { colBytes } : {}),
        keyIndex: k?.column_index_length ? [num(k.column_index_offset), k.column_index_length, num(k.offset_index_offset), k.offset_index_length] : null };
      if (p) o.payloadOffsetIndex = p.offset_index_length ? [num(p.offset_index_offset), p.offset_index_length] : null;
      first += Number(rg.num_rows);
      return o;
    }),
  };
}
const manifest = { version: 1, created: new Date().toISOString(), files };
writeFileSync(`${root}/manifest.json`, JSON.stringify(manifest));
console.log(JSON.stringify({ files: Object.keys(files).length, bytes: JSON.stringify(manifest).length }));
