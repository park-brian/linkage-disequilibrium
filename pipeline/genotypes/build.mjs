// v2 genotype builder (PLAN_REVIEW.md step 1). One pass over a 1000G VCF.gz → MAF-tiered Parquet files.
//
// usage: node --max-old-space-size=16000 build_genotypes_v2.mjs <vcf.gz> <outDir> <panelDir> <dbsnpPos.parquet> <chromCode>
//   dbsnpPos.parquet: DuckDB export sorted by (key, rsid) with columns key37|key38, rsid, ref_rsid, type, function
//                     (e.g. dbsnp_src/pos37_chr22.parquet); key = chromCode * 1e9 + pos
//
// Changes from v1 (see docs/PLAN_REVIEW.md):
//  - tiers use maf_max over the 26 sub-populations (not the 5 super-populations)
//  - rsid INT64 labelled from the 2023 dbSNP dump: VCF rsID (merged → current) else position + type, lowest rsID wins;
//    rsid_candidates = number of type-compatible current rsIDs at the position; rsid_source records which rule fired
//  - dbSNP function joined by position (via the chosen rsID)
//  - multiallelic_source flag on rows split from multi-allelic records; VCF order kept within a position (parity)
//  - vcf_order (source VCF row index) for exact parity order across tiers; ALLELES=<NCBI slice> adds allele_rsids
//  - vcf_id: the raw VCF ID when it contains "rs" and is not "rs<rsid>" (merged or unlisted rsIDs; LDlink prints it in
//    LDproxy's RS_Number and its "does not match RS number" warnings); null otherwise
//  - male_haploid row flag + male sample indexes in key-value metadata instead of per-row missing bitsets;
//    `missing` holds only true '.' calls
import { createReadStream, mkdirSync, readFileSync, statSync, writeFileSync, readdirSync } from "node:fs";
import { createGunzip, zstdCompressSync, constants } from "node:zlib";
import { createInterface } from "node:readline";
import { parquetMetadataAsync, parquetRead, asyncBufferFromFile } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { ParquetWriter } from "hyparquet-writer";
import { fileWriter } from "hyparquet-writer/src/node.js";
import { TYPES, minimalForm, makeLabeler } from "./rsid_label.mjs";
import { gtCodes } from "./vcf_gt.mjs";
import { useDataCwd } from "../../tools/config.mjs";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.mjs)

const [vcfPath, outDir, panelDir, dbsnpPath, chromArg] = process.argv.slice(2);
const CHROM = BigInt(chromArg), KEYBASE = CHROM * 1_000_000_000n;
mkdirSync(outDir, { recursive: true });
const ZSTD = input => zstdCompressSync(input, { params: { [constants.ZSTD_c_compressionLevel]: 3 } });

// ---------- populations ----------
const SUPER = ["AFR", "AMR", "EAS", "EUR", "SAS"];
const SUBPOPS = readdirSync(panelDir).map(f => f.replace(/\.txt$/, "")).filter(p => !SUPER.includes(p) && p !== "ALL");
const panel = Object.fromEntries(SUBPOPS.map(p => [p, new Set(readFileSync(`${panelDir}/${p}.txt`, "utf8").split(/\s+/).filter(Boolean))]));

// ---------- dbSNP positions (sorted by key, rsid) ----------
// Load the dbSNP slice straight into typed arrays (chrX has ~42M rows). Values are placed by absolute row number:
// onChunk may deliver chunks in any order and may overlap (see PARQUET_LAYOUT.md §7).
const dbFile = await asyncBufferFromFile(dbsnpPath), dbMeta = await parquetMetadataAsync(dbFile);
const keyCol = dbMeta.schema.find(s => /^key/.test(s.name)).name;
const N = Number(dbMeta.num_rows);
const dbPos = new Int32Array(N), dbRsid = new Float64Array(N), dbRef = new Float64Array(N), dbType = new Uint8Array(N).fill(255);
const dbFunc = new Array(N).fill(null);
const seen = { [keyCol]: new Uint8Array(N), rsid: new Uint8Array(N), ref_rsid: new Uint8Array(N), type: new Uint8Array(N), function: new Uint8Array(N) };
await parquetRead({ file: dbFile, metadata: dbMeta, columns: [keyCol, "rsid", "ref_rsid", "type", "function"], compressors,
  onChunk: ({ columnName, columnData, rowStart }) => {
    const s = seen[columnName];
    for (let i = 0; i < columnData.length; i++) {
      const r = rowStart + i; if (r >= N) break;
      const v = columnData[i]; s[r] = 1;
      if (columnName === keyCol) dbPos[r] = Number(BigInt(v) - KEYBASE);
      else if (columnName === "rsid") dbRsid[r] = Number(v);
      else if (columnName === "ref_rsid") dbRef[r] = Number(v);
      else if (columnName === "type") dbType[r] = TYPES.indexOf(v);
      else dbFunc[r] = v ?? null;
    }
  } });
for (const [c, s] of Object.entries(seen)) if (s.indexOf(0) !== -1) throw new Error(`dbSNP slice: column ${c} has unfilled rows`);
for (let i = 1; i < N; i++) if (dbPos[i] < dbPos[i - 1]) throw new Error(`dbSNP slice not sorted at row ${i}`);
console.error(`dbSNP rows loaded: ${dbPos.length}`);
// ---------- optional NCBI alleles (step 8 allele-aware mode): ALLELES=<slice.parquet> ----------
// slice: one chromosome of an NCBI dbSNP VCF (columns rsid, pos, ref, alt), sorted by pos (VCF order). Each 1000G row
// gets allele_rsids = the rsIDs with an ALT whose minimal form (trim common suffix, then common prefix; POS moves
// right) equals the row's minimal (pos, ref, alt). Stored as typed arrays + one latin1 buffer (chrX: 42.6M records).
let alleleIndex = null;
if (process.env.ALLELES) {
  const f = await asyncBufferFromFile(process.env.ALLELES), m = await parquetMetadataAsync(f);
  const n = Number(m.num_rows), apos = new Int32Array(n), aend = new Int32Array(n), arsid = new Float64Array(n), off = new Uint32Array(n + 1);
  let buf = new Uint8Array(64 << 20), used = 0, row = 0, maxRef = 1;
  for (const rg of m.row_groups) { // row group by row group, in order, so the byte buffer is filled sequentially
    const rowEnd = row + Number(rg.num_rows);
    await parquetRead({ file: f, metadata: m, compressors, columns: ["rsid", "pos", "ref", "alt"], rowStart: row, rowEnd, rowFormat: "array",
      onComplete: rows => {
        for (const [rs, p, ref, alt] of rows) {
          const s = `${ref}\t${alt}`;
          if (used + s.length > buf.length) { const nb = new Uint8Array(buf.length * 2); nb.set(buf.subarray(0, used)); buf = nb; }
          for (let k = 0; k < s.length; k++) buf[used++] = s.charCodeAt(k);
          apos[row] = p; aend[row] = p + ref.length - 1; arsid[row] = Number(rs); off[++row] = used; if (ref.length > maxRef) maxRef = ref.length;
        }
      } });
    if (row !== rowEnd) throw new Error(`alleles slice: expected ${rowEnd} rows, got ${row}`);
  }
  for (let i = 1; i < n; i++) if (apos[i] < apos[i - 1]) throw new Error(`alleles slice not sorted at row ${i}`);
  const dec = new TextDecoder("latin1");
  const text = i => dec.decode(buf.subarray(off[i], off[i + 1]));
  const first = p => { let a = 0, b = n; while (a < b) { const h = (a + b) >> 1; if (apos[h] < p) a = h + 1; else b = h; } return a; };
  alleleIndex = {
    records: n, maxRef,
    /** rsIDs (ascending, unique) whose ALT has the same minimal form as (pos, ref, alt) */
    match(pos, ref, alt) {
      const [mp, mr, ma] = minimalForm(pos, ref, alt), out = new Set();
      // a record at Q has minimal position in [Q, Q + len(ref) - 1], so only Q in [mp - maxRef + 1, mp] can match
      for (let i = first(mp - maxRef + 1); i < n && apos[i] <= mp; i++) {
        if (aend[i] < mp) continue; // REF span ends before mp: cannot reach it (skips the decode)
        const [r, alts] = text(i).split("\t");
        for (const a of alts.split(",")) { const [qp, qr, qa] = minimalForm(apos[i], r, a); if (qp === mp && qr === mr && qa === ma) { out.add(arsid[i]); break; } }
      }
      return [...out].sort((x, y) => x - y);
    },
  };
  console.error(`NCBI allele records loaded: ${n} (max REF ${maxRef})`);
}

// rsID labels (rsid_label.mjs): VCF rsID, else NCBI allele match, else position + type, searching the dump at the
// position(s) LDlink stores the variant at (POS+1 for deletions/delins); rsid_pos_delta = dump pos − POS
const labelRsid = makeLabeler({ dbPos, dbRsid, dbRef, dbType, dbFunc }, alleleIndex);

// ---------- writers ----------
// FINAL=1 writes only the chosen layouts (PARQUET_LAYOUT.md §6); otherwise the page-size sweep
const LAYOUTS = process.env.FINAL ? [["common", "bits", 64 << 10], ["low", "bits", 64 << 10], ["rare", "sparse", 16 << 10]] : [
  ["common", "bits", 16 << 10], ["common", "bits", 64 << 10], ["common", "bits", 256 << 10],
  ["low", "bits", 64 << 10], ["low", "bits", 256 << 10],
  ["rare", "sparse", 16 << 10], ["rare", "sparse", 64 << 10], ["rare", "sparse", 256 << 10],
];
const tierOf = m => (m >= 0.05 ? "common" : m >= 0.01 ? "low" : "rare");
let H = 0, nSamples = 0, samples = [], subIdx = [];
const writers = [];
function openWriters() {
  for (const [tier, enc, pageSize] of LAYOUTS) {
    const schema = [
      { name: "root", num_children: 12 },
      { name: "pos", type: "INT32", repetition_type: "REQUIRED" },
      { name: "rsid", type: "INT64", repetition_type: "OPTIONAL" },
      { name: "rsid_pos_delta", type: "INT32", repetition_type: "OPTIONAL" },
      { name: "rsid_candidates", type: "INT32", repetition_type: "REQUIRED" },
      { name: "rsid_source", type: "INT32", repetition_type: "REQUIRED" },
      { name: "ref", type: "BYTE_ARRAY", converted_type: "UTF8", repetition_type: "REQUIRED" },
      { name: "alt", type: "BYTE_ARRAY", converted_type: "UTF8", repetition_type: "REQUIRED" },
      { name: "multiallelic_source", type: "BOOLEAN", repetition_type: "REQUIRED" },
      { name: "male_haploid", type: "BOOLEAN", repetition_type: "REQUIRED" },
      { name: "ac", type: "INT32", repetition_type: "REQUIRED" },
      enc === "bits"
        ? { name: "haps", type: "FIXED_LEN_BYTE_ARRAY", type_length: H >> 3, repetition_type: "REQUIRED" }
        : { name: "haps", type: "BYTE_ARRAY", repetition_type: "REQUIRED" },
      { name: "missing", type: "BYTE_ARRAY", repetition_type: "OPTIONAL" },
      { name: "maf_max", type: "FLOAT", repetition_type: "REQUIRED" },
      { name: "function", type: "BYTE_ARRAY", converted_type: "UTF8", repetition_type: "OPTIONAL" },
      { name: "vcf_order", type: "INT32", repetition_type: "REQUIRED" },
      { name: "vcf_id", type: "BYTE_ARRAY", converted_type: "UTF8", repetition_type: "OPTIONAL" },
      ...(alleleIndex ? [{ name: "allele_is_label", type: "BOOLEAN", repetition_type: "REQUIRED" },
        { name: "allele_rsids", type: "BYTE_ARRAY", converted_type: "UTF8", repetition_type: "OPTIONAL" }] : []),
    ];
    schema[0].num_children = schema.length - 1;
    const name = `${tier}_${enc}_p${pageSize >> 10}k`, file = `${outDir}/${name}.parquet`;
    const w = new ParquetWriter({ writer: fileWriter(file), schema, codec: "ZSTD", compressors: { ZSTD } });
    writers.push({ name, tier, enc, pageSize, w, file, rows: 0 });
  }
}
function toSparse(bits, ac) {
  const listRef = ac > H - ac, idx = [];
  for (let k = 0; k < H; k++) if (((bits[k >> 3] >> (k & 7)) & 1) === (listRef ? 0 : 1)) idx.push(k);
  const out = new Uint8Array(1 + 3 * idx.length); out[0] = listRef ? 1 : 0;
  let o = 1, prev = 0;
  for (const k of idx) { let d = k - prev; prev = k; while (d >= 0x80) { out[o++] = (d & 0x7f) | 0x80; d >>= 7; } out[o++] = d; }
  return out.subarray(0, o);
}
const COLS = ["pos", "rsid", "rsid_pos_delta", "rsid_candidates", "rsid_source", "ref", "alt", "multiallelic_source", "male_haploid", "ac", "bits", "sparse", "missing", "maf_max", "function", "vcf_order", "vcf_id", "allele_is_label", "allele_rsids"];
const newChunk = () => Object.fromEntries(COLS.map(c => [c, []]));
let chunk;
function flush() {
  if (!chunk || !chunk.pos.length) return;
  const byTier = { common: [], low: [], rare: [] };
  chunk.maf_max.forEach((m, i) => byTier[tierOf(m)].push(i));
  for (const wr of writers) {
    const keep = byTier[wr.tier]; if (!keep.length) continue;
    const sel = a => keep.map(i => a[i]);
    wr.rows += keep.length;
    wr.w.write({
      rowGroupSize: 50000, pageSize: wr.pageSize,
      columnData: [
        { name: "pos", data: Int32Array.from(sel(chunk.pos)), columnIndex: true, encoding: "DELTA_BINARY_PACKED" },
        { name: "rsid", data: sel(chunk.rsid).map(v => (v === null ? null : BigInt(v))), statistics: false },
        { name: "rsid_pos_delta", data: sel(chunk.rsid_pos_delta), statistics: false },
        { name: "rsid_candidates", data: Int32Array.from(sel(chunk.rsid_candidates)), statistics: false },
        { name: "rsid_source", data: Int32Array.from(sel(chunk.rsid_source)), statistics: false },
        { name: "ref", data: sel(chunk.ref), statistics: false },
        { name: "alt", data: sel(chunk.alt), statistics: false },
        { name: "multiallelic_source", data: sel(chunk.multiallelic_source), statistics: false },
        { name: "male_haploid", data: sel(chunk.male_haploid), statistics: false },
        { name: "ac", data: Int32Array.from(sel(chunk.ac)), statistics: false },
        { name: "haps", data: sel(wr.enc === "bits" ? chunk.bits : chunk.sparse), statistics: false },
        { name: "missing", data: sel(chunk.missing), statistics: false },
        { name: "maf_max", data: Float32Array.from(sel(chunk.maf_max)), statistics: false },
        { name: "function", data: sel(chunk.function), statistics: false },
        { name: "vcf_order", data: Int32Array.from(sel(chunk.vcf_order)), statistics: false, encoding: "DELTA_BINARY_PACKED" },
        { name: "vcf_id", data: sel(chunk.vcf_id), statistics: false },
        ...(alleleIndex ? [{ name: "allele_is_label", data: sel(chunk.allele_is_label), statistics: false },
          { name: "allele_rsids", data: sel(chunk.allele_rsids), statistics: false }] : []),
      ],
    });
  }
  chunk = newChunk();
}

// ---------- VCF pass ----------
const stats = { records: 0, rows: 0, skipped_sv: 0, rsid_source: [0, 0, 0, 0, 0], ambiguous_pos_type: 0,
  vcf_rs_rows: 0, vcf_rs_found_at_pos: 0, vcf_id_rows: 0, pos_type_agrees_with_vcf: 0, multiallelic_rows: 0, male_haploid_rows: 0, max_ref_len: 1,
  ...(alleleIndex ? { allele_rsids: [0, 0, 0], allele_label_agree: 0, allele_label_disagree: 0 } : {}) }; // allele_rsids: rows with 0 / 1 / 2+ matches
const haploidSamples = new Set();
const t0 = performance.now();
const rl = createInterface({ input: createReadStream(vcfPath).pipe(createGunzip()), crlfDelay: Infinity });
for await (const line of rl) {
  if (line.charCodeAt(0) === 35) {
    if (line.startsWith("#CHROM")) {
      samples = line.split("\t").slice(9); nSamples = samples.length; H = Math.ceil((2 * nSamples) / 8) * 8;
      subIdx = SUBPOPS.map(p => samples.map((x, i) => (panel[p].has(x) ? i : -1)).filter(i => i >= 0));
      openWriters(); chunk = newChunk();
    }
    continue;
  }
  stats.records++;
  const f = line.split("\t"), alts = f[4].split(",");
  if (alts.some(a => a.startsWith("<"))) { stats.skipped_sv++; continue; }
  const pos = +f[1];
  // GT per haplotype, parsed once per record: allele index, -1 for '.', -2 for an absent 2nd haplotype (haploid call).
  // Single-digit indexes (<= 9 ALTs) by character; otherwise split the GT (1000G has records with 12 ALTs).
  const codes = gtCodes(f, nSamples, alts.length);
  let anyHaploid = false;
  for (let s = 0; s < nSamples; s++) if (codes[2 * s + 1] === -2) { anyHaploid = true; haploidSamples.add(s); }
  for (let a = 0; a < alts.length; a++) {
    const want = a + 1, bits = new Uint8Array(H >> 3);
    let miss = null, ac = 0;
    for (let k = 0; k < 2 * nSamples; k++) {
      const c = codes[k];
      if (c === want) { bits[k >> 3] |= 1 << (k & 7); ac++; } else if (c === -1) (miss ??= new Uint8Array(H >> 3))[k >> 3] |= 1 << (k & 7);
    }
    let mafMax = 0;
    for (const idx of subIdx) {
      let alt = 0, n = 0;
      for (const s of idx) for (let k = 2 * s; k <= 2 * s + 1; k++) {
        const c = codes[k]; if (c < 0) continue; // missing or absent (haploid)
        n++; alt += c === want ? 1 : 0;
      }
      if (n) mafMax = Math.max(mafMax, Math.min(alt, n - alt) / n);
    }
    const lab = labelRsid(pos, f[2] === "." ? null : f[2], f[3], alts[a]);
    stats.rsid_source[lab.source]++;
    if (lab.candidates > 1) stats.ambiguous_pos_type++;
    if (lab.vcfRs !== null) {
      stats.vcf_rs_rows++;
      if (lab.source === 1 || lab.source === 2) { stats.vcf_rs_found_at_pos++; if (lab.posType === lab.rsid) stats.pos_type_agrees_with_vcf++; }
    }
    if (alts.length > 1) stats.multiallelic_rows++;
    if (anyHaploid) stats.male_haploid_rows++;
    chunk.pos.push(pos); chunk.rsid.push(lab.rsid); chunk.rsid_pos_delta.push(lab.delta); chunk.rsid_candidates.push(lab.candidates); chunk.rsid_source.push(lab.source);
    chunk.ref.push(f[3]); chunk.alt.push(alts[a]); chunk.multiallelic_source.push(alts.length > 1); chunk.male_haploid.push(anyHaploid);
    chunk.ac.push(ac); chunk.bits.push(bits); chunk.sparse.push(toSparse(bits, ac)); chunk.missing.push(miss);
    chunk.maf_max.push(mafMax); chunk.function.push(lab.func ?? null);
    if (f[3].length > stats.max_ref_len) stats.max_ref_len = f[3].length;
    chunk.vcf_order.push(stats.rows); // row index in the source VCF (split ALTs in ALT order): exact VCF order across tiers
    const vid = f[2].includes("rs") && f[2] !== `rs${lab.rsid}` ? f[2] : null;
    chunk.vcf_id.push(vid); if (vid !== null) stats.vcf_id_rows++;
    if (alleleIndex) {
      const am = lab.alleles; // computed once by the labeller
      // matches == [rsid label] (the common case) is stored as a flag; the list only when it says something else
      const isLabel = am.length === 1 && am[0] === lab.rsid;
      chunk.allele_is_label.push(isLabel); chunk.allele_rsids.push(isLabel || !am.length ? null : am.join(","));
      stats.allele_rsids[Math.min(am.length, 2)]++;
      if (lab.rsid !== null && am.length) (am.includes(lab.rsid) ? stats.allele_label_agree++ : stats.allele_label_disagree++);
    }
    stats.rows++;
  }
  if (chunk.pos.length >= 50000) flush();
}
flush();
const kv = [
  { key: "samples", value: samples.join(",") },
  { key: "male_haploid_samples", value: [...haploidSamples].sort((a, b) => a - b).join(",") },
  { key: "tier_rule", value: "maf_max over 26 sub-populations: common >= 0.05, low >= 0.01, rare < 0.01" },
  { key: "rsid_rule", value: "dbSNP 2023 dump: VCF rsID (merged->current) else position+type, lowest rsID; source 0 none,1 vcf,2 vcf-merged,3 pos+type,4 NCBI allele match; searched at POS and POS+1 (dump stores deletions/delins at POS+1); rsid_pos_delta = dump position - pos" },
  { key: "max_ref_len", value: String(stats.max_ref_len) }, // readers: rows overlapping p start in [p - max_ref_len + 1, p]
  { key: "vcf_id_rule", value: "vcf_id = the source VCF ID when it contains 'rs' and differs from 'rs' + rsid, else null" },
  { key: "vcf_order_rule", value: "row index in the source VCF (multi-allelic records split in ALT order); parity picks the last overlapping row by vcf_order" },
  ...(alleleIndex ? [{ key: "allele_rule", value: `allele_rsids: current rsIDs from ${process.env.ALLELES.split(/[\\/]/).slice(-2).join("/")} (NCBI dbSNP VCF) with an ALT whose minimal (pos, ref, alt) equals the row's; allele_is_label = the matches are exactly [rsid]; otherwise allele_rsids = comma list, ascending (null = no match)` }] : []),
];
for (const { w } of writers) { w.kvMetadata = kv; w.finish(); }
const report = { vcf: vcfPath.split(/[\\/]/).pop(), ...stats, seconds: +((performance.now() - t0) / 1000).toFixed(1),
  files: Object.fromEntries(writers.map(({ name, file, rows }) => [name, { bytes: statSync(file).size, rows, bytes_per_row: +(statSync(file).size / rows).toFixed(1) }])) };
writeFileSync(`${outDir}/build_report.json`, JSON.stringify(report, null, 2));
console.log(JSON.stringify(report, null, 2));
