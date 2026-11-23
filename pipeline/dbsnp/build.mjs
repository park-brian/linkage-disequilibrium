// Write the final dbSNP index files with hyparquet-writer from the DuckDB-sorted intermediates (dbsnp_sort.sql).
// usage: node build_dbsnp_v2.mjs <rsid|pos37|pos38> <outRoot> [partition...]
//   rsid  → <outRoot>/dbsnp/rsid/part-<K>.parquet       rsid · ref_rsid (null unless merged) · chrom · key37 · key38
//   pos37 → <outRoot>/dbsnp/grch37/chr<N>.parquet        key · rsid · is_ref · type · function
//   pos38 → <outRoot>/dbsnp/grch38/chr<N>.parquet
// Row order inside every file is the intermediate's order: (rsid, dump ord) or (key, dump ord) — required for
// LDlink parity (find_one returns the first doc in dump order; find returns all in dump order).
import { mkdirSync, readdirSync, statSync, writeFileSync } from "node:fs";
import { zstdCompressSync, constants } from "node:zlib";
import { parquetMetadataAsync, parquetRead, asyncBufferFromFile } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { ParquetWriter } from "hyparquet-writer";
import { fileWriter } from "hyparquet-writer/src/node.js";
import { useDataCwd } from "../../tools/config.mjs";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.mjs)

const [kind, outRoot, ...only] = process.argv.slice(2);
const ZSTD = input => zstdCompressSync(input, { params: { [constants.ZSTD_c_compressionLevel]: 3 } });
const CHUNK = 1_000_000;
const CHROM_NAME = c => (c === 23 ? "X" : c === 24 ? "Y" : c === 25 ? "MT" : String(c));

const SPEC = {
  rsid: {
    src: "dbsnp_dump/sorted_rsid", partKey: "part", dir: "dbsnp/rsid", name: p => `part-${String(p).padStart(2, "0")}.parquet`,
    pageSize: 8 << 10,
    schema: [
      { name: "root", num_children: 5 },
      { name: "rsid", type: "INT64", repetition_type: "REQUIRED" },
      { name: "ref_rsid", type: "INT64", repetition_type: "OPTIONAL" },
      { name: "chrom", type: "INT32", repetition_type: "OPTIONAL" },
      { name: "key37", type: "INT64", repetition_type: "OPTIONAL" },
      { name: "key38", type: "INT64", repetition_type: "OPTIONAL" },
    ],
    columns: ["rsid", "ref_rsid", "chrom", "key37", "key38"],
    columnData: c => [
      { name: "rsid", data: BigInt64Array.from(c.rsid, BigInt), encoding: "DELTA_BINARY_PACKED", columnIndex: true },
      { name: "ref_rsid", data: c.ref_rsid.map(v => (v == null ? null : BigInt(v))), statistics: false },
      { name: "chrom", data: c.chrom.map(v => (v == null ? null : Number(v))), statistics: false },
      { name: "key37", data: c.key37.map(v => (v == null ? null : BigInt(v))), statistics: false },
      { name: "key38", data: c.key38.map(v => (v == null ? null : BigInt(v))), statistics: false },
    ],
  },
  pos: build => ({
    src: `dbsnp_dump/sorted_pos${build}`, partKey: "chrom", dir: `dbsnp/grch${build}`, name: p => `chr${CHROM_NAME(+p)}.parquet`,
    pageSize: 16 << 10,
    schema: [
      { name: "root", num_children: 5 },
      { name: "key", type: "INT64", repetition_type: "REQUIRED" },
      { name: "rsid", type: "INT64", repetition_type: "REQUIRED" },
      { name: "is_ref", type: "BOOLEAN", repetition_type: "REQUIRED" },
      { name: "type", type: "BYTE_ARRAY", converted_type: "UTF8", repetition_type: "OPTIONAL" },
      { name: "function", type: "BYTE_ARRAY", converted_type: "UTF8", repetition_type: "OPTIONAL" },
    ],
    columns: ["key", "rsid", "is_ref", "type", "function"],
    columnData: c => [
      { name: "key", data: BigInt64Array.from(c.key, BigInt), encoding: "DELTA_BINARY_PACKED", columnIndex: true },
      { name: "rsid", data: BigInt64Array.from(c.rsid, BigInt), statistics: false },
      { name: "is_ref", data: c.is_ref.map(Boolean), statistics: false },
      { name: "type", data: c.type, statistics: false },
      { name: "function", data: c.function, statistics: false },
    ],
  }),
};
const spec = kind === "rsid" ? SPEC.rsid : SPEC.pos(kind.slice(3));
mkdirSync(`${outRoot}/${spec.dir}`, { recursive: true });

async function readChunk(files, metas, cols, start, end) {
  // Read rows [start, end) across the partition's files in order.
  // onChunk may deliver whole pages that extend beyond [rowStart, rowEnd) (and chunks may overlap), so every value
  // is placed by its absolute row number into an exactly-sized buffer; anything outside the window is dropped.
  const len = end - start;
  const out = Object.fromEntries(cols.map(c => [c, new Array(len)]));
  const filled = Object.fromEntries(cols.map(c => [c, new Uint8Array(len)]));
  let base = 0;
  for (let f = 0; f < files.length; f++) {
    const n = Number(metas[f].num_rows), lo = Math.max(start - base, 0), hi = Math.min(end - base, n);
    if (lo < hi) {
      await parquetRead({ file: files[f], metadata: metas[f], columns: cols, compressors, rowStart: lo, rowEnd: hi,
        onChunk: ({ columnName, columnData, rowStart }) => {
          const o = out[columnName], fl = filled[columnName];
          for (let i = 0; i < columnData.length; i++) {
            const r = rowStart + i; // row within this file
            if (r < lo || r >= hi) continue;
            const j = base + r - start; o[j] = columnData[i]; fl[j] = 1;
          }
        } });
    }
    base += n;
  }
  for (const c of cols) if (filled[c].indexOf(0) !== -1) throw new Error(`readChunk: column ${c} missing rows in [${start}, ${end})`);
  return out;
}

const report = {};
const partitions = readdirSync(spec.src).filter(f => f.startsWith(`${spec.partKey}=`) && f.endsWith(".parquet"))
  .map(f => f.slice(spec.partKey.length + 1, -8)).filter(p => !only.length || only.includes(p)).sort((a, b) => +a - +b);
for (const p of partitions) {
  const t0 = performance.now();
  const files = [await asyncBufferFromFile(`${spec.src}/${spec.partKey}=${p}.parquet`)];
  const metas = await Promise.all(files.map(f => parquetMetadataAsync(f)));
  const total = metas.reduce((a, m) => a + Number(m.num_rows), 0);
  if (!total) { console.log(JSON.stringify({ file: spec.name(p), rows: 0, skipped: true })); continue; }
  const outFile = `${outRoot}/${spec.dir}/${spec.name(p)}`;
  const w = new ParquetWriter({ writer: fileWriter(outFile), schema: spec.schema, codec: "ZSTD", compressors: { ZSTD },
    kvMetadata: [{ key: "source", value: "LDlink mongodump 2023-04 dbsnp (b155); order = (key, dump order) for parity" }] });
  // verify (key, dump ord) order — parity depends on it
  let prevKey = null, prevOrd = null, sorted = true;
  for (let s = 0; s < total; s += CHUNK) {
    const c = await readChunk(files, metas, [...spec.columns, "ord"], s, Math.min(s + CHUNK, total));
    const keyName = spec.columns[0];
    c[keyName].forEach((v, i) => {
      const b = BigInt(v), o = BigInt(c.ord[i]);
      if (prevKey !== null && (b < prevKey || (b === prevKey && o <= prevOrd))) sorted = false;
      prevKey = b; prevOrd = o;
    });
    w.write({ columnData: spec.columnData(c), rowGroupSize: CHUNK, pageSize: spec.pageSize });
  }
  w.finish();
  if (!sorted) throw new Error(`${outFile}: input not in (key, dump order) — parity depends on it`); // B4: abort, don't just report
  // the written file must have exactly the intermediate's row count
  const written = Number((await parquetMetadataAsync(await asyncBufferFromFile(outFile))).num_rows);
  if (written !== total) throw new Error(`${outFile}: wrote ${written} rows, expected ${total}`);
  report[spec.name(p)] = { rows: total, bytes: statSync(outFile).size, bytes_per_row: +(statSync(outFile).size / total).toFixed(2), sorted, seconds: +((performance.now() - t0) / 1000).toFixed(1) };
  console.log(JSON.stringify({ file: spec.name(p), ...report[spec.name(p)] }));
}
writeFileSync(`${outRoot}/${spec.dir}/build_report.json`, JSON.stringify(report, null, 2));
