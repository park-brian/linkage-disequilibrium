// Benchmark reading genotype Parquet layouts with hyparquet and computing LD with bitsets.
// usage: node bench_genotypes.mjs <dir with *.parquet> <panelDir>
// Reports per query: HTTP-like requests, bytes read, wall time (local disk, so time ≈ CPU decode),
// and LD throughput. Network time ≈ requests in sequential phases × RTT + bytes / bandwidth.
import { openSync, readSync, fstatSync, readFileSync, readdirSync, writeFileSync } from "node:fs";
import { parquetMetadataAsync, parquetReadObjects } from "hyparquet";
import { compressors } from "hyparquet-compressors";

const [dir, panelDir] = process.argv.slice(2);

function countingFile(path) {
  const fd = openSync(path, "r"), byteLength = fstatSync(fd).size;
  const stats = { requests: 0, bytes: 0 };
  return {
    stats, byteLength,
    slice(start, end = byteLength) {
      stats.requests++; stats.bytes += end - start;
      const buf = Buffer.allocUnsafe(end - start); readSync(fd, buf, 0, end - start, start);
      return Promise.resolve(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength));
    },
  };
}

const POPCNT = new Uint8Array(1 << 16); for (let i = 1; i < POPCNT.length; i++) POPCNT[i] = (i & 1) + POPCNT[i >> 1];
const pc32 = v => POPCNT[v & 0xffff] + POPCNT[v >>> 16];

// Gather the selected population's haplotype bits into a dense Uint32 vector (+ valid mask).
function packPopulation(rows, hapIdx) {
  const W = Math.ceil(hapIdx.length / 32);
  return rows.map(r => {
    const alt = new Uint32Array(W), valid = new Uint32Array(W), h = r.haps, m = r.missing;
    for (let j = 0; j < hapIdx.length; j++) {
      const k = hapIdx[j];
      if (m && (m[k >> 3] >> (k & 7)) & 1) continue;
      valid[j >>> 5] |= 1 << (j & 31);
      if ((h[k >> 3] >> (k & 7)) & 1) alt[j >>> 5] |= 1 << (j & 31);
    }
    return { pos: r.pos, alt, valid, full: !m };
  });
}

function ldPair(a, b) {
  let n = 0, pa = 0, pb = 0, x = 0;
  const W = a.alt.length;
  if (a.full && b.full) {
    for (let w = 0; w < W; w++) { const aa = a.alt[w], bb = b.alt[w]; pa += pc32(aa); pb += pc32(bb); x += pc32(aa & bb); n += pc32(a.valid[w]); }
  } else {
    for (let w = 0; w < W; w++) { const v = a.valid[w] & b.valid[w], aa = a.alt[w] & v, bb = b.alt[w] & v; n += pc32(v); pa += pc32(aa); pb += pc32(bb); x += pc32(aa & bb); }
  }
  // haplotype counts with A = ref/ref, B = ref/alt, C = alt/ref, D = alt/alt (index order, not string order)
  const D = x, B = pa - x, C = pb - x, A = n - pa - pb + x;
  const delta = A * D - B * C, Ms = (A + C) * (B + D) * (A + B) * (C + D);
  if (!Ms) return [NaN, NaN];
  const dmax = delta < 0 ? Math.min((A + C) * (A + B), (B + D) * (C + D)) : Math.min((A + C) * (C + D), (A + B) * (B + D));
  return [(delta * delta) / Ms, Math.abs(delta / dmax)];
}

async function query(path, filter, label, extra = {}) {
  const file = countingFile(path);
  const t0 = performance.now();
  const metadata = await parquetMetadataAsync(file);
  const metaReq = { ...file.stats };
  const rows = await parquetReadObjects({ file, metadata, compressors, columns: ["pos", "haps", "missing"], filter, ...extra });
  const t1 = performance.now();
  return { label, file: path.split(/[\\/]/).pop(), rows: rows.length, requests: file.stats.requests, meta_requests: metaReq.requests,
    MB: +(file.stats.bytes / 1e6).toFixed(2), read_ms: Math.round(t1 - t0), data: rows };
}

const files = readdirSync(dir).filter(f => f.endsWith(".parquet")).map(f => `${dir}/${f}`);
const meta0 = await parquetMetadataAsync(countingFile(files[0]));
const samples = meta0.key_value_metadata.find(kv => kv.key === "samples").value.split(",");
const pop = name => new Set(readFileSync(`${panelDir}/${name}.txt`, "utf8").split(/\s+/).filter(Boolean));
const hapIdxOf = set => samples.flatMap((s, i) => (set.has(s) ? [2 * i, 2 * i + 1] : []));
const POPS = { CEU: hapIdxOf(pop("CEU")), EUR: hapIdxOf(pop("EUR")), ALL: samples.flatMap((_, i) => [2 * i, 2 * i + 1]) };

const results = [];
const START = 25_000_000;
function report(r, extra = {}) { const { data, ...rest } = r; const o = { ...rest, ...extra }; results.push(o); console.log(JSON.stringify(o)); }

// 1) Region scans: 1 Mb and 10 Mb, every layout
for (const span of [1_000_000, 10_000_000]) {
  for (const f of files) {
    const r = await query(f, { pos: { $gte: START, $lte: START + span } }, `region ${span / 1e6}Mb`);
    report(r);
  }
}

// 2) Scattered lists: 300 / 3000 common SNPs sampled from the 10 Mb window (positions from the common tier)
const common = files.find(f => f.includes("tier_common"));
const all10 = (await query(common, { pos: { $gte: START, $lte: START + 10_000_000 } }, "seed")).data.map(r => r.pos);
for (const n of [300, 3000]) {
  const step = Math.floor(all10.length / n), sel = all10.filter((_, i) => i % step === 0).slice(0, n);
  for (const f of files.filter(f => /tier_common|bits_zstd3_p(64k|256k)|bits_zstd3_p16k/.test(f))) {
    report(await query(f, { pos: { $in: sel } }, `list ${n} SNPs over 10Mb`));
  }
}

// 3) LD compute on the 10 Mb common tier: banded (≤500 kb) pyramid pairs, per population
const region = await query(common, { pos: { $gte: START, $lte: START + 10_000_000 } }, "ld input");
for (const [name, idx] of Object.entries(POPS)) {
  const t0 = performance.now();
  const packed = packPopulation(region.data, idx);
  const t1 = performance.now();
  let pairs = 0, strong = 0;
  for (let i = 0; i < packed.length; i++) {
    for (let j = i + 1; j < packed.length && packed[j].pos - packed[i].pos <= 500_000; j++) {
      const [r2] = ldPair(packed[i], packed[j]); pairs++; if (r2 >= 0.8) strong++;
    }
  }
  const t2 = performance.now();
  const o = { label: `LD banded 500kb, 10Mb common tier, ${name}`, variants: packed.length, haplotypes: idx.length, pairs,
    pairs_r2_ge_0_8: strong, pack_ms: Math.round(t1 - t0), ld_ms: Math.round(t2 - t1), Mpairs_per_s: +(pairs / (t2 - t1) / 1e3).toFixed(1) };
  results.push(o); console.log(JSON.stringify(o));
}
writeFileSync(`${dir}/bench_report.json`, JSON.stringify(results, null, 2));
