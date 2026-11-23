// Slice 9 (docs/PLAN.md): can hyparquet 1.33's own page-index pruning replace src/parquet.js sortedLookup? Each case
// looks the same keys up both ways on a cold, byte-counting AsyncBuffer and compares rows (file order), requests,
// bytes and CPU time.
//   ours:      sortedLookup (column index → candidate pages → key column → payload pages via its own offset index)
//   hyparquet: parquetScan({ pruningFilter: { key: { $in } }, usePageIndex: true }) → ranges; readColumn(key) per range,
//              then readColumn(each other column) over the matched rows only
// usage: node bench/pruning.js [release]   (data-relative; default release_v4) → bench/results/pruning.json
import { readFileSync, writeFileSync } from "node:fs";
import { asyncBufferFromFile, parquetMetadataAsync, parquetScan, parquetReadObjects } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { sortedLookup } from "../src/parquet.js";
import { useDataCwd, repo } from "../tools/config.js";
useDataCwd();

const root = process.argv[2] ?? "release_v4";
async function counted(path) {
  const f = await asyncBufferFromFile(`${root}/${path}`), stats = { requests: 0, bytes: 0 };
  return { stats, file: { byteLength: f.byteLength, async slice(s, e) { stats.requests++; const b = await f.slice(s, e); stats.bytes += b.byteLength; return b; } } };
}
const jsonish = rows => JSON.stringify(rows, (k, v) => (typeof v === "bigint" ? `${v}n` : v instanceof Uint8Array ? Array.from(v.subarray(0, 8)) : v));

async function ours(path, column, keys, columns, payloadColumn) {
  const { file, stats } = await counted(path), metadata = await parquetMetadataAsync(file), t0 = performance.now();
  const before = { ...stats };
  const { rows } = await sortedLookup({ file, metadata, column, keys, columns, compressors, payloadColumn });
  return { ms: performance.now() - t0, requests: stats.requests - before.requests, bytes: stats.bytes - before.bytes, rows: rows.map(r => Object.fromEntries(columns.map(c => [c, r[c]]))) };
}
async function theirs(path, column, keys, columns) {
  const { file, stats } = await counted(path), metadata = await parquetMetadataAsync(file), t0 = performance.now();
  const before = { ...stats }, want = new Set(keys.map(String)), out = [];
  const scan = await parquetScan({ file, metadata, compressors, columns, pruningFilter: { [column]: { $in: keys } }, usePageIndex: true });
  for (const { rowStart, rowEnd } of scan.ranges) {
    const k = await scan.readColumn({ column, rowStart, rowEnd }), hits = [];
    for (let i = 0; i < k.length; i++) if (want.has(String(k[i]))) hits.push(rowStart + i);
    for (const row of hits) {
      const r = {};
      for (const c of columns) r[c] = c === column ? k[row - rowStart] : (await scan.readColumn({ column: c, rowStart: row, rowEnd: row + 1 }))[0];
      out.push(r);
    }
  }
  return { ms: performance.now() - t0, requests: stats.requests - before.requests, bytes: stats.bytes - before.bytes, rows: out };
}

// keys drawn from the files themselves (a fixed stride), so every case has hits
async function sampleKeys(path, column, n) {
  const f = await asyncBufferFromFile(`${root}/${path}`), m = await parquetMetadataAsync(f), N = Number(m.num_rows), keys = [];
  for (let i = 0; i < n; i++) { const row = Math.floor(((i + 0.5) / n) * N); keys.push((await parquetReadObjects({ file: f, metadata: m, compressors, columns: [column], rowStart: row, rowEnd: row + 1 }))[0][column]); }
  return [...new Set(keys.map(String))].map(k => (typeof keys[0] === "bigint" ? BigInt(k) : Number(k)));
}

const CASES = [
  { name: "rsID index, 32 keys", path: "dbsnp/rsid/part-00.parquet", column: "rsid", n: 32, columns: ["rsid", "ref_rsid", "chrom", "key37", "key38"] },
  { name: "rsID index, 512 keys", path: "dbsnp/rsid/part-00.parquet", column: "rsid", n: 512, columns: ["rsid", "ref_rsid", "chrom", "key37", "key38"] },
  { name: "dbSNP by position, 200 keys", path: "dbsnp/grch37/chr22.parquet", column: "key", n: 200, columns: ["key", "rsid", "function"] },
  { name: "genotypes, 50 points (haps payload)", path: "genotypes/grch37/chr22.common.parquet", column: "pos", n: 50, columns: ["pos", "rsid", "haps"], payload: "haps" },
];
const report = [];
for (const c of CASES) {
  const keys = await sampleKeys(c.path, c.column, c.n);
  const a = await ours(c.path, c.column, keys, c.columns, c.payload), b = await theirs(c.path, c.column, keys, c.columns);
  const r = { case: c.name, keys: keys.length, sameRows: jsonish(a.rows) === jsonish(b.rows), rows: [a.rows.length, b.rows.length],
    ours: { requests: a.requests, bytes: a.bytes, ms: Math.round(a.ms) }, hyparquet: { requests: b.requests, bytes: b.bytes, ms: Math.round(b.ms) } };
  report.push(r); console.log(JSON.stringify(r));
}
writeFileSync(repo("bench/results/pruning.json"), JSON.stringify({ hyparquet: JSON.parse(readFileSync(repo("node_modules/hyparquet/package.json"), "utf8")).version, root, report }, null, 2));
