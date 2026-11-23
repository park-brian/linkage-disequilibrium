// Benchmark v2 tier files: window scans (1/10 Mb) and point reads (scattered positions) through a cached
// AsyncBuffer, plus per-page compressed sizes from the offset index. usage: node bench_v2.mjs <dir>
import { openSync, readSync, fstatSync, readdirSync, writeFileSync } from "node:fs";
import { parquetMetadataAsync, parquetRead, cachedAsyncBuffer, readOffsetIndex } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { sortedLookup } from "../src/parquet.mjs";
import { useDataCwd } from "../tools/config.mjs";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.mjs)

const [dir] = process.argv.slice(2);
const COLUMNS = ["pos", "rsid", "haps", "missing", "male_haploid"];
function countingFile(path) {
  const fd = openSync(path, "r"), byteLength = fstatSync(fd).size, stats = { requests: 0, bytes: 0 };
  return { stats, byteLength, slice(s, e = byteLength) {
    stats.requests++; stats.bytes += e - s; const b = Buffer.allocUnsafe(e - s); readSync(fd, b, 0, e - s, s);
    return Promise.resolve(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength)); } };
}
async function open(path) {
  const raw = countingFile(path), file = cachedAsyncBuffer(raw);
  const metadata = await parquetMetadataAsync(file, { initialFetchSize: 64 << 10 });
  return { raw, file, metadata };
}
async function pageStats(path) {
  const { file, metadata } = await open(path);
  const hapsIdx = metadata.row_groups[0].columns.findIndex(c => c.meta_data.path_in_schema[0] === "haps");
  const sizes = [], rows = [];
  for (const rg of metadata.row_groups) {
    const c = rg.columns[hapsIdx]; if (!c.offset_index_length) continue; const buf = await file.slice(Number(c.offset_index_offset), Number(c.offset_index_offset) + c.offset_index_length);
    const locs = readOffsetIndex({ view: new DataView(buf), offset: 0 }).page_locations;
    locs.forEach((l, i) => { sizes.push(l.compressed_page_size); rows.push((i + 1 < locs.length ? Number(locs[i + 1].first_row_index) : Number(rg.num_rows)) - Number(l.first_row_index)); });
  }
  if (!sizes.length) return { haps_pages: 0, note: "no offset index on haps" };
  sizes.sort((a, b) => a - b); rows.sort((a, b) => a - b);
  return { haps_pages: sizes.length, haps_page_KB_median: +(sizes[sizes.length >> 1] / 1024).toFixed(1), rows_per_page_median: rows[rows.length >> 1] };
}
async function scan(path, lo, hi) {
  const { raw, file, metadata } = await open(path), t0 = performance.now();
  // pos is sorted: find row range from row-group stats, then read with rowStart/rowEnd using the page index via filter
  let n = 0;
  await parquetRead({ file, metadata, compressors, utf8: false, columns: COLUMNS, rowFormat: "object", filter: { pos: { $gte: lo, $lte: hi } }, onComplete: r => { n = r.length; } });
  return { rows: n, requests: raw.stats.requests, KB: Math.round(raw.stats.bytes / 1024), ms: Math.round(performance.now() - t0) };
}
async function points(path, keys) {
  const { raw, file, metadata } = await open(path), t0 = performance.now();
  const r = await sortedLookup({ file, metadata, column: "pos", keys, columns: COLUMNS, compressors, payloadColumn: "haps" });
  return { found: r.rows.length, ranges: r.ranges, requests: raw.stats.requests, KB: Math.round(raw.stats.bytes / 1024), KB_per_point: +(raw.stats.bytes / 1024 / keys.length).toFixed(1), ms: Math.round(performance.now() - t0) };
}
// positions to probe: sample real positions from each file so every probe hits
async function allPositions(path) {
  const { file, metadata } = await open(path); let pos = [];
  await parquetRead({ file, metadata, compressors, utf8: false, columns: ["pos"], onComplete: r => { pos = r.map(x => x[0]); } });
  return pos;
}
let seed = 7; const rnd = () => ((seed = (seed * 1103515245 + 12345) >>> 0) / 2 ** 32);
const sample = (arr, n) => Array.from({ length: n }, () => arr[Math.floor(rnd() * arr.length)]);

const results = [];
const log = o => { results.push(o); console.log(JSON.stringify(o)); };
for (const f of readdirSync(dir).filter(f => f.endsWith(".parquet")).sort()) {
  const path = `${dir}/${f}`;
  log({ file: f, test: "pages", ...(await pageStats(path)) });
  for (const span of [1e6, 1e7]) log({ file: f, test: `window ${span / 1e6}Mb`, ...(await scan(path, 25e6, 25e6 + span)) });
  const pos = await allPositions(path);
  const in10 = pos.filter(p => p >= 25e6 && p <= 35e6);
  for (const n of [30, 300]) log({ file: f, test: `points x${n} across chr22`, ...(await points(path, sample(pos, n))) });
  log({ file: f, test: "points x300 within 10Mb", ...(await points(path, sample(in10, 300))) });
}
writeFileSync(`${dir}/bench_report.json`, JSON.stringify(results, null, 2));
