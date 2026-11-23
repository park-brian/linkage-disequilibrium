// Step 6: run the planner scenarios over HTTPS against the private S3 bucket, through presigned URLs
// (HTTP Range requests; nothing public). Each scenario is cold (new planner, empty caches).
// usage: node net_bench.js <presigned.json> <localReleaseRoot> <panelDir>
//   presigned.json: { "<relative path>": "<presigned https url>", ... } (see presign.sh)
import { readFileSync, writeFileSync, statSync } from "node:fs";
import { asyncBufferFromFile, parquetMetadataAsync, parquetRead, cachedAsyncBuffer } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { openPlanner } from "../src/planner.js";
import { loadKernel, KERNEL_WASM } from "../src/ld/kernel.js";
import { useDataCwd, repo } from "../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

const [presignedPath, localRoot, panelDir] = process.argv.slice(2);
const urls = JSON.parse(readFileSync(presignedPath, "utf8"));
const NET = "net";

// optional cap on concurrent requests (browser: ~6 per host over HTTP/1.1, which is all S3's REST endpoint speaks;
// HTTP/2 via CloudFront multiplexes ~100)
const MAX_INFLIGHT = +(process.env.MAX_INFLIGHT || Infinity);
let active = 0; const waiters = [];
const acquire = () => (active < MAX_INFLIGHT ? (active++, Promise.resolve()) : new Promise(r => waiters.push(r)));
const release = () => { const w = waiters.shift(); if (w) w(); else active--; };
function counter() {
  const stats = { requests: 0, bytes: 0, maxInFlight: 0, inFlight: 0 };
  const openFile = async path => {
    const rel = path.slice(NET.length + 1), url = urls[rel];
    if (!url) throw new Error(`no URL for ${rel}`);
    const byteLength = statSync(`${localRoot}/${rel}`).size; // known from the manifest in a real client; avoids a HEAD
    return cachedAsyncBuffer({ byteLength, async slice(s, e = byteLength) {
      await acquire();
      stats.requests++; stats.inFlight++; stats.maxInFlight = Math.max(stats.maxInFlight, stats.inFlight);
      try {
        const res = await fetch(url, { headers: { Range: `bytes=${s}-${e - 1}` } });
        if (res.status !== 206 && res.status !== 200) throw new Error(`${rel}: HTTP ${res.status}`);
        const buf = await res.arrayBuffer(); stats.bytes += buf.byteLength;
        return res.status === 200 ? buf.slice(s, e) : buf;
      } finally { stats.inFlight--; release(); }
    } }, { minSize: 0 });
  };
  return { stats, openFile };
}

// baseline round-trip: 10 sequential 1-byte range requests
{
  const url = urls["manifest.json"], t = [];
  for (let i = 0; i < 10; i++) { const t0 = performance.now(); await (await fetch(url, { headers: { Range: "bytes=0-0" } })).arrayBuffer(); t.push(performance.now() - t0); }
  t.sort((a, b) => a - b);
  console.log(JSON.stringify({ test: "rtt", median_ms: Math.round(t[5]), min_ms: Math.round(t[0]), max_ms: Math.round(t[9]) }));
}

// inputs drawn from the local copy (same procedure as verify_planner.js)
async function labels(path, lo, hi) {
  const file = await asyncBufferFromFile(`${localRoot}/${path}`), metadata = await parquetMetadataAsync(file);
  let rows = []; await parquetRead({ file, metadata, compressors, rowFormat: "object", columns: ["pos", "rsid"], filter: { pos: { $gte: lo, $lte: hi } }, onComplete: r => { rows = r; } });
  return rows.filter(r => r.rsid != null);
}
let seed = 5; const rnd = () => ((seed = (seed * 1103515245 + 12345) >>> 0) / 2 ** 32);
const pick = (arr, n) => Array.from({ length: n }, () => arr[Math.floor(rnd() * arr.length)]);
const g = "genotypes/grch37";
const clustered = [...await labels(`${g}/chr22.common.parquet`, 30e6, 32e6), ...await labels(`${g}/chr22.low.parquet`, 30e6, 32e6), ...await labels(`${g}/chr22.rare.parquet`, 30e6, 32e6)];
const chr22all = await labels(`${g}/chr22.common.parquet`, 0, 1e9), chrXall = await labels(`${g}/chrX.common.parquet`, 0, 1e9);

const kernel = await loadKernel(readFileSync(KERNEL_WASM));
const eur = new Set(readFileSync(`${panelDir}/EUR.txt`, "utf8").split(/\s+/).filter(Boolean));

const results = [];
async function run(label, q, ld = false) {
  const c = counter();
  const t0 = performance.now();
  const pl = await openPlanner(NET, c.openFile);
  const tManifest = performance.now() - t0;
  const r = await pl.query(q);
  const tQuery = performance.now() - t0;
  let ldInfo = null;
  if (ld) {
    // LD on the fetched rows: EUR, 500 kb band, r² ≥ 0.8, pruned, blocked kernel, single thread
    const meta = (await (await c.openFile(`${NET}/${g}/chr22.common.parquet`))); // samples come from file metadata
    const md = await parquetMetadataAsync(meta);
    const samples = md.key_value_metadata.find(x => x.key === "samples").value.split(",");
    const idx = samples.flatMap((s, i) => (eur.has(s) ? [2 * i, 2 * i + 1] : []));
    const t1 = performance.now();
    const p = kernel.prepare(r.rows, { samples: samples.length, hapIndexes: idx });
    const b = kernel.band(p, { maxDist: 500000, threshold: 0.8, prune: true, block: 128 });
    ldInfo = { ld_ms: Math.round(performance.now() - t1), pairs_evaluated: b.evaluated, pairs_pruned: b.pruned, kept: b.count };
  }
  const o = { label, rows: r.rows.length, found: r.variants.filter(v => v?.found).length, rsid: r.plan.rsid?.strategy ?? null,
    fetch: Object.fromEntries(Object.entries(r.plan.fetch.reduce((a, f) => ((a[f.kind] = (a[f.kind] || 0) + 1), a), {}))),
    requests: c.stats.requests, max_in_flight: c.stats.maxInFlight, MB: +(c.stats.bytes / 1e6).toFixed(2), predicted_MB: +(r.plan.predictedBytes / 1e6).toFixed(2),
    manifest_ms: Math.round(tManifest), query_ms: Math.round(tQuery), MB_per_s: +(c.stats.bytes / 1e3 / tQuery).toFixed(1), ...(ldInfo ?? {}) };
  results.push(o); console.log(JSON.stringify(o));
}
await run("window 1 Mb, all tiers", { build: "grch37", regions: [{ chrom: "22", start: 25e6, end: 26e6 }] });
await run("window 10 Mb, common (+ LD EUR r2>=0.8)", { build: "grch37", regions: [{ chrom: "22", start: 25e6, end: 35e6 }], tiers: ["common"] }, true);
await run("window 10 Mb, all tiers", { build: "grch37", regions: [{ chrom: "22", start: 25e6, end: 35e6 }] });
await run("300 rsIDs clustered in 2 Mb", { build: "grch37", inputs: pick(clustered, 300).map(r => `rs${r.rsid}`) });
await run("3000 rsIDs clustered in 2 Mb", { build: "grch37", inputs: pick(clustered, 3000).map(r => `rs${r.rsid}`) });
await run("300 rsIDs scattered chr22 + chrX", { build: "grch37", inputs: [...pick(chr22all, 150), ...pick(chrXall, 150)].map(r => `rs${r.rsid}`) });
await run("30 coords scattered on chr22", { build: "grch37", inputs: pick(chr22all, 30).map(r => `chr22:${r.pos}`) });
await run("300 coords in 2 Mb", { build: "grch37", inputs: pick(clustered, 300).map(r => `chr22:${r.pos}`) });
writeFileSync(repo(`bench/results/v2_net_bench${process.env.MAX_INFLIGHT ? "_inflight" + process.env.MAX_INFLIGHT : ""}.json`), JSON.stringify(results, null, 2));
