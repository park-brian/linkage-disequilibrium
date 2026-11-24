// Build a manifest for a release directory: the one JSON the client reads first, holding only what it reads:
//  - every Parquet file's size and footer length, so each footer is one exact Range request and no file needs a size
//    probe (src/core/source.js);
//  - the row groups the planner prices its hybrid rsID strategy with (src/parity/planner.js): for the genotype files
//    key range, payload offset index and the label columns' sizes (PRICED_COLUMNS); for the rsID index files key
//    range and key-column index locations. Stored column-wise.
// Version 3 drops what nothing reads (v2 also had row counts, column names, key-value metadata, the dbSNP position,
// annotation and LD score row groups, and per-column sizes of every priced genotype column): 0.91 MB instead of
// 1.92 MB for release_v4 (0.36 / 0.61 MB gzipped). src/parity/planner.js expandManifest reads v1-v3.
// usage: node pipeline/release/manifest.js <releaseRoot>     → <releaseRoot>/manifest.json
import { readdirSync, statSync, writeFileSync } from "node:fs";
import { parquetMetadataAsync, asyncBufferFromFile } from "hyparquet";
import { useDataCwd } from "../../tools/config.js";
import { PRICED_COLUMNS } from "../../src/parity/planner.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

const [root] = process.argv.slice(2);
// the key column of the files that carry row groups (the others are skipped before KEY is asked)
const KEY = name => (name.startsWith("genotypes/") ? "pos" : "rsid");
const PAYLOAD = name => (name.startsWith("genotypes/") ? "haps" : null);
const walk = d => readdirSync(`${root}/${d}`, { withFileTypes: true }).flatMap(e =>
  e.isDirectory() ? walk(`${d}/${e.name}`) : e.name.endsWith(".parquet") ? [`${d}/${e.name}`] : []);
const num = v => (typeof v === "bigint" ? (v <= BigInt(Number.MAX_SAFE_INTEGER) ? Number(v) : String(v)) : v);

/** [{ min, max, ..., colBytes: { c: n } }] → { min: [...], max: [...], ..., colBytes: { c: [...] } } (keyIndex only if any) */
function columnwise(rgs) {
  const out = { min: rgs.map(r => r.min), max: rgs.map(r => r.max) };
  if (rgs.some(r => r.keyIndex)) out.keyIndex = rgs.map(r => r.keyIndex);
  if (rgs.some(r => "payloadOffsetIndex" in r)) out.payloadOffsetIndex = rgs.map(r => r.payloadOffsetIndex ?? null);
  if (rgs.some(r => r.colBytes)) out.colBytes = Object.fromEntries(Object.keys(rgs[0].colBytes).map(c => [c,
    rgs.map(r => r.colBytes[c] ?? 0)]));
  return out;
}
const files = {};
for (const name of readdirSync(root, { withFileTypes: true }).filter(e => e.isDirectory()).flatMap(e => walk(e.name)).sort()) {
  const path = `${root}/${name}`, byteLength = statSync(path).size;
  const m = await parquetMetadataAsync(await asyncBufferFromFile(path));
  const cols = m.row_groups[0]?.columns.map(c => c.meta_data.path_in_schema.join(".")) ?? [];
  const ki = cols.indexOf(KEY(name)), pi = PAYLOAD(name) ? cols.indexOf(PAYLOAD(name)) : -1;
  files[name] = { byteLength, footer: { length: m.metadata_length + 8 } };
  const geno = name.startsWith("genotypes/");
  if (!geno && !name.startsWith("dbsnp/rsid/")) continue; // the planner prices only these
  files[name].rowGroups = columnwise(m.row_groups.map(rg => {
      const k = ki >= 0 ? rg.columns[ki] : null, st = k?.meta_data.statistics ?? {}, p = pi >= 0 ? rg.columns[pi] : null;
      // compressed bytes per column (+ total) so a planner can cost scans without opening the file
      const size = Object.fromEntries(rg.columns.map(c => [c.meta_data.path_in_schema.join("."),
        Number(c.meta_data.total_compressed_size)]));
      const o = { min: num(st.min_value), max: num(st.max_value) };
      if (geno) {
        o.payloadOffsetIndex = p?.offset_index_length ? [num(p.offset_index_offset), p.offset_index_length] : null;
        o.colBytes = Object.fromEntries(PRICED_COLUMNS.filter(c => c in size).map(c => [c, size[c]]));
      } else o.keyIndex = k?.column_index_length ? [num(k.column_index_offset), k.column_index_length, num(k.offset_index_offset),
        k.offset_index_length] : null;
      return o;
    }));
}
const manifest = { version: 3, created: new Date().toISOString(), files };
writeFileSync(`${root}/manifest.json`, JSON.stringify(manifest));
console.log(JSON.stringify({ files: Object.keys(files).length, bytes: JSON.stringify(manifest).length }));
