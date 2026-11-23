// Same rsID lookups as bench_dbsnp.mjs but through hyparquet's cachedAsyncBuffer (dedupes repeated
// index/footer slices, as a browser client would) and a small initial footer fetch.
import { openSync, readSync, fstatSync } from "node:fs";
import { parquetMetadataAsync, cachedAsyncBuffer } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { sortedLookup } from "../src/parquet.mjs";
function countingFile(path) {
  const fd = openSync(path, "r"), byteLength = fstatSync(fd).size, stats = { requests: 0, bytes: 0 };
  return { stats, byteLength, slice(s, e = byteLength) {
    stats.requests++; stats.bytes += e - s; const b = Buffer.allocUnsafe(e - s); readSync(fd, b, 0, e - s, s);
    return Promise.resolve(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength)); } };
}
let seed = 42; const rnd = () => ((seed = (seed * 1103515245 + 12345) >>> 0) / 2 ** 32);
for (const n of [1, 16, 300, 3000]) {
  const ids = Array.from({ length: n }, () => BigInt(1_000_000_000 + Math.floor(rnd() * 20_000_000)));
  for (const f of ["rsid_p8k_rg1m", "rsid_p32k_rg1m"]) {
    const raw = countingFile(`out/dbsnp/${f}.parquet`), file = cachedAsyncBuffer(raw), t0 = performance.now();
    const metadata = await parquetMetadataAsync(file, { initialFetchSize: 64 << 10 });
    const r = await sortedLookup({ file, metadata, column: "rsid", keys: ids, columns: ["rsid", "key37", "key38"], compressors });
    const o = { label: `rsid lookup x${n} (cached)`, file: f, found: r.rows.length, ranges: r.ranges, requests: raw.stats.requests, KB: Math.round(raw.stats.bytes / 1024), KB_per_key: +(raw.stats.bytes / 1024 / n).toFixed(1), ms: Math.round(performance.now() - t0) };
    console.log(JSON.stringify(o));
  }
}
