// Build candidate genotype Parquet layouts from one 1000G VCF.gz in a single pass.
// usage: node build_genotypes.mjs <vcf.gz> <outDir> <panelDir> [maxRecords]
// Also writes MAF-tiered files keyed on maf_max = max minor-allele frequency over the 5 super-populations
// (so population-specific variants stay in the common tier).
// Each output file has: pos INT32 (sorted), id, ref, alt, ac, haps (one bit per haplotype), missing (bitset or null).
// Multi-allelic records are split into one row per ALT (like the GRCh38 releases). Symbolic SVs are skipped.
import { createReadStream, mkdirSync, statSync, writeFileSync } from "node:fs";
import { createGunzip, zstdCompressSync, gzipSync, constants } from "node:zlib";
import { createInterface } from "node:readline";
import { ParquetWriter } from "hyparquet-writer";
import { fileWriter } from "hyparquet-writer/src/node.js";

const [vcfPath, outDir, panelDir, maxArg] = process.argv.slice(2);
const SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"];
const { readFileSync } = await import("node:fs");
const panels = Object.fromEntries(SUPERPOPS.map(p => [p, new Set(readFileSync(`${panelDir}/${p}.txt`, "utf8").split(/\s+/).filter(Boolean))]));
const MAX = maxArg ? +maxArg : Infinity;
mkdirSync(outDir, { recursive: true });

const zstd = level => input => zstdCompressSync(input, { params: { [constants.ZSTD_c_compressionLevel]: level } });
const gzip = input => gzipSync(input, { level: 6 });

// Candidate layouts: [name, encoding, codec, pageSize (uncompressed bytes), rowGroupRows, tier filter on maf_max]
const T = { common: m => m >= 0.05, low: m => m >= 0.01 && m < 0.05, rare: m => m < 0.01 };
const VARIANTS = [
  ["tier_common_zstd3_p256k", "bits", "ZSTD", 256 << 10, 50000, T.common],
  ["tier_low_zstd3_p256k", "bits", "ZSTD", 256 << 10, 50000, T.low],
  ["tier_rare_zstd3_p256k", "bits", "ZSTD", 256 << 10, 50000, T.rare],
  ["tier_rare_sparse_zstd3_p256k", "sparse", "ZSTD", 256 << 10, 50000, T.rare],
  ["bits_snappy_p64k", "bits", "SNAPPY", 64 << 10, 50000],
  ["bits_gzip_p64k", "bits", "GZIP", 64 << 10, 50000],
  ["bits_zstd3_p16k", "bits", "ZSTD", 16 << 10, 50000],
  ["bits_zstd3_p64k", "bits", "ZSTD", 64 << 10, 50000],
  ["bits_zstd3_p256k", "bits", "ZSTD", 256 << 10, 50000],
  ["bits_zstd3_p1m", "bits", "ZSTD", 1 << 20, 50000],
  ["sparse_zstd3_p64k", "sparse", "ZSTD", 64 << 10, 50000],
];

let H = 0, nSamples = 0, samples = [], popIdx = [];
const writers = [];
function openWriters() {
  const bytes = H >> 3;
  for (const [name, enc, codec, pageSize, rg, tier] of VARIANTS) {
    const schema = [
      { name: "root", num_children: 8 },
      { name: "pos", type: "INT32", repetition_type: "REQUIRED" },
      { name: "id", type: "BYTE_ARRAY", converted_type: "UTF8", repetition_type: "OPTIONAL" },
      { name: "ref", type: "BYTE_ARRAY", converted_type: "UTF8", repetition_type: "REQUIRED" },
      { name: "alt", type: "BYTE_ARRAY", converted_type: "UTF8", repetition_type: "REQUIRED" },
      { name: "ac", type: "INT32", repetition_type: "REQUIRED" },
      enc === "bits"
        ? { name: "haps", type: "FIXED_LEN_BYTE_ARRAY", type_length: bytes, repetition_type: "REQUIRED" }
        : { name: "haps", type: "BYTE_ARRAY", repetition_type: "REQUIRED" },
      { name: "missing", type: "BYTE_ARRAY", repetition_type: "OPTIONAL" },
      { name: "maf_max", type: "FLOAT", repetition_type: "REQUIRED" },
    ];
    const file = `${outDir}/${name}.parquet`;
    const w = new ParquetWriter({
      writer: fileWriter(file), schema, codec: codec === "ZSTD9" ? "ZSTD" : codec,
      compressors: { ZSTD: codec === "ZSTD9" ? zstd(9) : zstd(3), GZIP: gzip },
      kvMetadata: [{ key: "samples", value: samples.join(",") }, { key: "haps_encoding", value: enc }],
    });
    writers.push({ name, enc, pageSize, rg, w, file, tier, rows: 0 });
  }
}

// bits -> sparse: [flag(minor is ALT?1:0... stored as which allele is listed), uint16 deltas of listed haplotype indexes]
function toSparse(bits, ac) {
  const listRef = ac > H - ac; // list the minor allele's haplotypes
  const idx = [];
  for (let k = 0; k < H; k++) {
    const b = (bits[k >> 3] >> (k & 7)) & 1;
    if (b === (listRef ? 0 : 1)) idx.push(k);
  }
  const out = new Uint8Array(1 + 3 * idx.length); // varint deltas, worst case 3 bytes (H < 2^21)
  out[0] = listRef ? 1 : 0;
  let o = 1, prev = 0;
  for (const k of idx) { let d = k - prev; prev = k; while (d >= 0x80) { out[o++] = (d & 0x7f) | 0x80; d >>= 7; } out[o++] = d; }
  return out.subarray(0, o);
}

let chunk = null;
function newChunk() { return { pos: [], id: [], ref: [], alt: [], ac: [], haps: [], sparse: [], missing: [], maf: [] }; }
function flush() {
  if (!chunk || !chunk.pos.length) return;
  for (const wr of writers) {
    const { enc, pageSize, rg, w, tier } = wr;
    const keep = tier ? chunk.maf.map((m, i) => (tier(m) ? i : -1)).filter(i => i >= 0) : null;
    const sel = arr => (keep ? keep.map(i => arr[i]) : arr);
    if (keep && !keep.length) continue;
    wr.rows += keep ? keep.length : chunk.pos.length;
    w.write({
      rowGroupSize: rg, pageSize,
      columnData: [
        { name: "pos", data: Int32Array.from(sel(chunk.pos)), columnIndex: true, encoding: "DELTA_BINARY_PACKED" },
        { name: "id", data: sel(chunk.id), statistics: false },
        { name: "ref", data: sel(chunk.ref), statistics: false },
        { name: "alt", data: sel(chunk.alt), statistics: false },
        { name: "ac", data: Int32Array.from(sel(chunk.ac)), statistics: false },
        { name: "haps", data: sel(enc === "bits" ? chunk.haps : chunk.sparse), statistics: false },
        { name: "missing", data: sel(chunk.missing), statistics: false },
        { name: "maf_max", data: Float32Array.from(sel(chunk.maf)), statistics: false },
      ],
    });
  }
  chunk = newChunk();
}

const t0 = performance.now();
let rows = 0, skipped = 0, records = 0;
const rl = createInterface({ input: createReadStream(vcfPath).pipe(createGunzip()), crlfDelay: Infinity });
for await (const line of rl) {
  if (line.charCodeAt(0) === 35) {
    if (line.startsWith("#CHROM")) {
      samples = line.split("\t").slice(9); nSamples = samples.length; H = 2 * nSamples;
      H = Math.ceil(H / 8) * 8; // pad to whole bytes
      popIdx = SUPERPOPS.map(p => samples.map((x, i) => (panels[p].has(x) ? i : -1)).filter(i => i >= 0));
      openWriters(); chunk = newChunk();
    }
    continue;
  }
  if (++records > MAX) break;
  const f = line.split("\t");
  const alts = f[4].split(",");
  if (alts.some(a => a.startsWith("<"))) { skipped++; continue; }
  for (let a = 0; a < alts.length; a++) {
    const want = 49 + a; // char code of allele index '1' + a (indexes > 9 not present in 1000G)
    const bits = new Uint8Array(H >> 3);
    let miss = null, ac = 0;
    for (let s = 0; s < nSamples; s++) {
      const g = f[9 + s];
      const c0 = g.charCodeAt(0), haploid = g.length === 1 || g.charCodeAt(1) === 58; // "0" or "0:..."
      const c1 = haploid ? -1 : g.charCodeAt(2);
      const k = 2 * s;
      if (c0 === want) { bits[k >> 3] |= 1 << (k & 7); ac++; }
      else if (c0 === 46) { (miss ??= new Uint8Array(H >> 3))[k >> 3] |= 1 << (k & 7); }
      if (haploid || c1 === 46) { (miss ??= new Uint8Array(H >> 3))[(k + 1) >> 3] |= 1 << ((k + 1) & 7); }
      else if (c1 === want) { bits[(k + 1) >> 3] |= 1 << ((k + 1) & 7); ac++; }
    }
    chunk.pos.push(+f[1]); chunk.id.push(f[2] === "." ? null : f[2]); chunk.ref.push(f[3]); chunk.alt.push(alts[a]);
    let mafMax = 0;
    for (const idx of popIdx) {
      let alt = 0, n = 0;
      for (const s of idx) for (const k of [2 * s, 2 * s + 1]) {
        if (miss && (miss[k >> 3] >> (k & 7)) & 1) continue;
        n++; alt += (bits[k >> 3] >> (k & 7)) & 1;
      }
      if (n) mafMax = Math.max(mafMax, Math.min(alt, n - alt) / n);
    }
    chunk.maf.push(mafMax);
    chunk.ac.push(ac); chunk.haps.push(bits); chunk.sparse.push(toSparse(bits, ac)); chunk.missing.push(miss);
    rows++;
  }
  if (chunk.pos.length >= 50000) flush();
}
flush();
for (const { w } of writers) w.finish();
const report = {
  vcf: vcfPath.split(/[\\/]/).pop(), records: Math.min(records, MAX), rows, skipped_sv: skipped, samples: nSamples,
  seconds: +((performance.now() - t0) / 1000).toFixed(1),
  files: Object.fromEntries(writers.map(({ name, file, rows }) => [name, { bytes: statSync(file).size, rows, bytes_per_row: +(statSync(file).size / rows).toFixed(1) }])),
};
writeFileSync(`${outDir}/build_report.json`, JSON.stringify(report, null, 2));
console.log(JSON.stringify(report, null, 2));
