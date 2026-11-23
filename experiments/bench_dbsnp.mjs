// Benchmark dbSNP index lookups: random rsID batches against the rsid-sorted layouts, and
// window scans against the position-sorted layout. usage: node bench_dbsnp.mjs <dir>
import { openSync, readSync, fstatSync, readdirSync, writeFileSync } from "node:fs";
import { parquetMetadataAsync, parquetReadObjects } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { sortedLookup } from "../src/parquet.mjs";

const [dir] = process.argv.slice(2);
function countingFile(path) {
  const fd = openSync(path, "r"), byteLength = fstatSync(fd).size, stats = { requests: 0, bytes: 0 };
  return { stats, byteLength, slice(s, e = byteLength) {
    stats.requests++; stats.bytes += e - s; const b = Buffer.allocUnsafe(e - s); readSync(fd, b, 0, e - s, s);
    return Promise.resolve(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength)); } };
}
async function lookup(path, column, keys, columns, label) {
  const file = countingFile(path), t0 = performance.now();
  const metadata = await parquetMetadataAsync(file);
  const r = await sortedLookup({ file, metadata, column, keys, columns, compressors });
  const o = { label, file: path.split(/[\/]/).pop(), found: r.rows.length, ranges: r.ranges, rows_decoded: r.pagesRows, requests: file.stats.requests, KB: Math.round(file.stats.bytes / 1024), ms: Math.round(performance.now() - t0) };
  console.log(JSON.stringify(o)); return o;
}
async function run(path, filter, columns, label) {
  const file = countingFile(path), t0 = performance.now();
  const metadata = await parquetMetadataAsync(file);
  const rows = await parquetReadObjects({ file, metadata, compressors, columns, filter });
  const o = { label, file: path.split(/[\\/]/).pop(), rows: rows.length, requests: file.stats.requests, KB: Math.round(file.stats.bytes / 1024), ms: Math.round(performance.now() - t0) };
  console.log(JSON.stringify(o)); return o;
}
const results = [];
const files = readdirSync(dir).filter(f => f.endsWith(".parquet"));
// deterministic pseudo-random rsIDs inside the 20M-id slice [1e9, 1.02e9]
let seed = 42; const rnd = () => ((seed = (seed * 1103515245 + 12345) >>> 0) / 2 ** 32);
for (const n of [1, 16, 300, 3000]) {
  const ids = Array.from({ length: n }, () => BigInt(1_000_000_000 + Math.floor(rnd() * 20_000_000)));
  for (const f of files.filter(f => f.startsWith("rsid_")))
    results.push(await lookup(`${dir}/${f}`, "rsid", ids, ["rsid", "key37", "key38"], `rsid lookup x${n}`));
}
for (const span of [1e6, 1e7]) {
  const lo = 22n * 1_000_000_000n + 25_000_000n, hi = lo + BigInt(span);
  for (const f of files.filter(f => f.startsWith("pos37_")))
    results.push(await run(`${dir}/${f}`, { key: { $gte: lo, $lte: hi } }, ["key", "rsid"], `pos window ${span / 1e6}Mb (key+rsid only)`));
}
writeFileSync(`${dir}/bench_report.json`, JSON.stringify(results, null, 2));
