// Build candidate dbSNP index layouts with hyparquet-writer from DuckDB-exported slices
// (out/dbsnp_src/rsid_slice_20Mids.parquet sorted by rsid, out/dbsnp_src/pos37_chr22.parquet sorted by key37).
// usage: node build_dbsnp.mjs <outDir>
// Composite key: chrom code (1-22, X=23, Y=24, MT=25) * 1e9 + position → one INT64 with a total genome order.
import { mkdirSync, statSync, writeFileSync } from "node:fs";
import { zstdCompressSync, constants } from "node:zlib";
import { parquetMetadataAsync, parquetRead, asyncBufferFromFile } from "hyparquet";
import { ParquetWriter } from "hyparquet-writer";
import { fileWriter } from "hyparquet-writer/src/node.js";

const [outDir] = process.argv.slice(2);
mkdirSync(outDir, { recursive: true });
const ZSTD = input => zstdCompressSync(input, { params: { [constants.ZSTD_c_compressionLevel]: 3 } });

async function readColumns(path, columns) {
  const file = await asyncBufferFromFile(path);
  const metadata = await parquetMetadataAsync(file);
  const out = Object.fromEntries(columns.map(c => [c, []]));
  // onChunk may arrive out of row order: keep rowStart and sort before concatenating
  await parquetRead({ file, metadata, columns, onChunk: ({ columnName, columnData, rowStart }) => out[columnName].push({ rowStart, columnData }) });
  return Object.fromEntries(columns.map(c => [c, out[c].sort((a, b) => a.rowStart - b.rowStart).flatMap(x => Array.from(x.columnData))]));
}

const toBig = a => BigInt64Array.from(a, v => (v === null || v === undefined ? 0n : BigInt(v)));
const report = {};
function write(name, schema, columnData, { pageSize, rowGroupSize }) {
  const file = `${outDir}/${name}.parquet`;
  const w = new ParquetWriter({ writer: fileWriter(file), schema, codec: "ZSTD", compressors: { ZSTD } });
  w.write({ columnData, pageSize, rowGroupSize });
  w.finish();
  report[name] = { bytes: statSync(file).size, rows: columnData[0].data.length, bytes_per_row: +(statSync(file).size / columnData[0].data.length).toFixed(2), pageSize, rowGroupSize };
  console.log(name, report[name]);
}

// --- rsid-sorted index: rsid, ref_rsid (null unless merged), key37, key38 ---
{
  const d = await readColumns("out/dbsnp_src/rsid_slice_20Mids.parquet", ["rsid", "ref_rsid", "key37", "key38"]);
  const ref = d.ref_rsid.map((r, i) => (r === d.rsid[i] ? null : BigInt(r)));
  const schema = [
    { name: "root", num_children: 4 },
    { name: "rsid", type: "INT64", repetition_type: "REQUIRED" },
    { name: "ref_rsid", type: "INT64", repetition_type: "OPTIONAL" },
    { name: "key37", type: "INT64", repetition_type: "OPTIONAL" },
    { name: "key38", type: "INT64", repetition_type: "OPTIONAL" },
  ];
  const cols = () => [
    { name: "rsid", data: toBig(d.rsid), encoding: "DELTA_BINARY_PACKED", columnIndex: true },
    { name: "ref_rsid", data: ref, statistics: false },
    { name: "key37", data: d.key37.map(v => (v === null ? null : BigInt(v))), statistics: false },
    { name: "key38", data: d.key38.map(v => (v === null ? null : BigInt(v))), statistics: false },
  ];
  for (const [pageSize, rowGroupSize] of [[8 << 10, 1_000_000], [32 << 10, 1_000_000], [128 << 10, 1_000_000], [32 << 10, 5_000_000]])
    write(`rsid_p${pageSize >> 10}k_rg${rowGroupSize / 1e6}m`, schema, cols(), { pageSize, rowGroupSize });
}

// --- position-sorted index: key37, rsid, type, function ---
{
  const d = await readColumns("out/dbsnp_src/pos37_chr22.parquet", ["key37", "rsid", "type", "function"]);
  const schema = [
    { name: "root", num_children: 4 },
    { name: "key", type: "INT64", repetition_type: "REQUIRED" },
    { name: "rsid", type: "INT64", repetition_type: "REQUIRED" },
    { name: "type", type: "BYTE_ARRAY", converted_type: "UTF8", repetition_type: "OPTIONAL" },
    { name: "function", type: "BYTE_ARRAY", converted_type: "UTF8", repetition_type: "OPTIONAL" },
  ];
  const cols = () => [
    { name: "key", data: toBig(d.key37), encoding: "DELTA_BINARY_PACKED", columnIndex: true },
    { name: "rsid", data: toBig(d.rsid), statistics: false },
    { name: "type", data: d.type, statistics: false },
    { name: "function", data: d.function, statistics: false },
  ];
  for (const [pageSize, rowGroupSize] of [[64 << 10, 1_000_000], [256 << 10, 1_000_000]])
    write(`pos37_chr22_p${pageSize >> 10}k_rg${rowGroupSize / 1e6}m`, schema, cols(), { pageSize, rowGroupSize });
}
writeFileSync(`${outDir}/build_report.json`, JSON.stringify(report, null, 2));
