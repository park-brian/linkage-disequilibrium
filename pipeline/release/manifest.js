// Build a manifest for a release directory: one JSON describing every Parquet file so a browser client can
//  - fetch each footer with one exact Range request (no 512 KB guess), or skip it,
//  - pick files / row groups by key range without opening them,
//  - plan point reads vs scans from page-index locations (PLAN_REVIEW.md findings 2 and 7).
// Version 2: row groups stored column-wise, and per-column sizes only for the genotype columns the planner prices
// (src/planner.js PRICED_COLUMNS): 1.9 MB instead of 7.2 MB for release_v4. src/planner.js expandManifest reads both.
// usage: node build_manifest.js <releaseRoot>     → <releaseRoot>/manifest.json
import { readdirSync, statSync, writeFileSync } from "node:fs";
import { parquetMetadataAsync, asyncBufferFromFile } from "hyparquet";
import { useDataCwd } from "../../tools/config.js";
import { PRICED_COLUMNS } from "../../src/planner.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

const [root] = process.argv.slice(2);
const KEY = name => {
  if (name.startsWith("genotypes/")) return "pos";
  if (name.startsWith("dbsnp/rsid/") || /(duplicated_rsids|forge_rsid)\.parquet$/.test(name)) return "rsid";
  if (name.endsWith("gene_names.parquet")) return "name";
  if (name.endsWith("gtex_genes.parquet")) return "gene_id_v26";
  if (/(gtex_tissues|platforms|panels|w_hm3\.snplist)\.parquet$/.test(name)) return null; // small lookup tables, read whole
  if (name.startsWith("ldscore/")) return "CHR"; // sorted by (CHR, BP); a client selects chromosomes by row group
  return "key";
};
const PAYLOAD = name => (name.startsWith("genotypes/") ? "haps" : null);
const walk = d => readdirSync(`${root}/${d}`, { withFileTypes: true }).flatMap(e =>
  e.isDirectory() ? walk(`${d}/${e.name}`) : e.name.endsWith(".parquet") ? [`${d}/${e.name}`] : []);
const num = v => (typeof v === "bigint" ? (v <= BigInt(Number.MAX_SAFE_INTEGER) ? Number(v) : String(v)) : v);

/** [{ min, max, ..., colBytes: { c: n } }] → { min: [...], max: [...], ..., colBytes: { c: [...] } } (keyIndex only if any) */
function columnwise(rgs) {
  const out = { min: rgs.map(r => r.min), max: rgs.map(r => r.max), bytes: rgs.map(r => r.bytes) };
  if (rgs.some(r => r.keyIndex)) out.keyIndex = rgs.map(r => r.keyIndex);
  if (rgs.some(r => "payloadOffsetIndex" in r)) out.payloadOffsetIndex = rgs.map(r => r.payloadOffsetIndex ?? null);
  if (rgs.some(r => r.colBytes)) out.colBytes = Object.fromEntries(Object.keys(rgs[0].colBytes).map(c => [c, rgs.map(r => r.colBytes[c] ?? 0)]));
  return out;
}
const files = {};
for (const name of readdirSync(root, { withFileTypes: true }).filter(e => e.isDirectory()).flatMap(e => walk(e.name)).sort()) {
  const path = `${root}/${name}`, byteLength = statSync(path).size;
  const m = await parquetMetadataAsync(await asyncBufferFromFile(path));
  const cols = m.row_groups[0]?.columns.map(c => c.meta_data.path_in_schema.join(".")) ?? [];
  const ki = cols.indexOf(KEY(name)), pi = PAYLOAD(name) ? cols.indexOf(PAYLOAD(name)) : -1;
  files[name] = {
    byteLength, numRows: Number(m.num_rows),
    footer: { offset: byteLength - 8 - m.metadata_length, length: m.metadata_length + 8 },
    key: KEY(name), columns: cols,
    kv: Object.fromEntries((m.key_value_metadata ?? []).filter(kv => kv.key !== "samples").map(kv => [kv.key, kv.value])),
    rowGroups: columnwise(m.row_groups.map(rg => {
      const k = ki >= 0 ? rg.columns[ki] : null, st = k?.meta_data.statistics ?? {}, p = pi >= 0 ? rg.columns[pi] : null;
      // compressed bytes per column (+ total) so a planner can cost scans without opening the file
      const size = Object.fromEntries(rg.columns.map(c => [c.meta_data.path_in_schema.join("."), Number(c.meta_data.total_compressed_size)]));
      const o = { min: num(st.min_value), max: num(st.max_value), bytes: Object.values(size).reduce((a, b) => a + b, 0),
        keyIndex: k?.column_index_length ? [num(k.column_index_offset), k.column_index_length, num(k.offset_index_offset), k.offset_index_length] : null };
      if (p) o.payloadOffsetIndex = p.offset_index_length ? [num(p.offset_index_offset), p.offset_index_length] : null;
      if (name.startsWith("genotypes/")) o.colBytes = Object.fromEntries(PRICED_COLUMNS.filter(c => c in size).map(c => [c, size[c]]));
      return o;
    })),
  };
}
const manifest = { version: 2, created: new Date().toISOString(), files };
writeFileSync(`${root}/manifest.json`, JSON.stringify(manifest));
console.log(JSON.stringify({ files: Object.keys(files).length, bytes: JSON.stringify(manifest).length }));
