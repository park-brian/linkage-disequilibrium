// Phase B integrity check over a genotype tree (<root>/<build>/chr<N>/{common_bits_p64k,low_bits_p64k,rare_sparse_p16k}.parquet
// + build_report.json, as written by build_all.sh): per chromosome
//   - tier row counts sum to the build report's rows; the expected columns and key-value metadata are present
//   - pos non-decreasing and vcf_order strictly increasing within each file; vcf_order across tiers is a permutation of
//     0..rows-1 (every source row exactly once)
//   - every row's maf_max lies in its tier (common >= 0.05 > low >= 0.01 > rare)
//   - rsid_pos_delta is set exactly when rsid is set; allele_is_label implies an rsid
// usage: node pipeline/release/check.js <root> [out.json]
import { readdirSync, readFileSync, existsSync, writeFileSync } from "node:fs";
import { asyncBufferFromFile, parquetMetadataAsync, parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { useDataCwd } from "../../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

const [root, outPath] = process.argv.slice(2);
const TIERS = { common: ["common_bits_p64k", m => m >= 0.05], low: ["low_bits_p64k", m => m >= 0.01 && m < 0.05], rare: ["rare_sparse_p16k", m => m < 0.01] };
const COLS = ["pos", "rsid", "rsid_pos_delta", "rsid_candidates", "rsid_source", "ref", "alt", "multiallelic_source", "male_haploid", "ac", "haps", "missing", "maf_max", "function", "vcf_order", "allele_is_label", "allele_rsids"];
const KV = ["samples", "male_haploid_samples", "tier_rule", "rsid_rule", "max_ref_len", "vcf_order_rule", "allele_rule"];
const report = { root, chromosomes: [], failures: 0 };
for (const build of readdirSync(root, { withFileTypes: true }).filter(d => d.isDirectory()).map(d => d.name).sort()) {
  for (const chr of readdirSync(`${root}/${build}`, { withFileTypes: true }).filter(d => d.isDirectory()).map(d => d.name)) {
    const dir = `${root}/${build}/${chr}`, r = { build, chr, problems: [] };
    if (!existsSync(`${dir}/build_report.json`)) { r.problems.push("no build_report.json (build incomplete)"); report.chromosomes.push(r); report.failures++; continue; }
    const br = JSON.parse(readFileSync(`${dir}/build_report.json`, "utf8"));
    r.rows = br.rows; r.bytes = 0; r.tiers = {};
    const seen = new Uint8Array(br.rows);
    for (const [tier, [name, inTier]] of Object.entries(TIERS)) {
      const path = `${dir}/${name}.parquet`;
      if (!existsSync(path)) { r.problems.push(`missing ${name}`); continue; }
      const file = await asyncBufferFromFile(path), m = await parquetMetadataAsync(file);
      r.bytes += file.byteLength;
      const cols = m.schema.slice(1).map(s => s.name);
      for (const c of COLS) if (!cols.includes(c)) r.problems.push(`${tier}: no column ${c}`);
      const kv = Object.fromEntries((m.key_value_metadata ?? []).map(x => [x.key, x.value]));
      for (const k of KV) if (!(k in kv)) r.problems.push(`${tier}: no metadata ${k}`);
      let n = 0, prevPos = -1, prevOrd = -1, bad = { order: 0, ord: 0, tier: 0, delta: 0, label: 0, dupOrd: 0, range: 0 };
      await parquetRead({ file, metadata: m, compressors, utf8: false, columns: ["pos", "vcf_order", "maf_max", "rsid", "rsid_pos_delta", "allele_is_label"], rowFormat: "array",
        onComplete: rows => {
          for (const [pos, ord, maf, rsid, delta, isLabel] of rows) {
            if (pos < prevPos) bad.order++;
            if (ord <= prevOrd) bad.ord++;
            if (!inTier(maf)) bad.tier++;
            if ((rsid == null) !== (delta == null)) bad.delta++;
            if (isLabel && rsid == null) bad.label++;
            if (ord < 0 || ord >= br.rows) bad.range++; else if (seen[ord]++) bad.dupOrd++;
            prevPos = pos; prevOrd = ord; n++;
          }
        } });
      r.tiers[tier] = n;
      for (const [k, v] of Object.entries(bad)) if (v) r.problems.push(`${tier}: ${v} rows fail ${k}`);
      if (n !== Number(m.num_rows)) r.problems.push(`${tier}: read ${n} of ${m.num_rows}`);
    }
    const total = Object.values(r.tiers).reduce((a, b) => a + b, 0);
    if (total !== br.rows) r.problems.push(`tiers hold ${total} rows, build report says ${br.rows}`);
    const missing = seen.reduce((a, v) => a + (v === 0 ? 1 : 0), 0);
    if (missing) r.problems.push(`${missing} vcf_order values never appear`);
    if (r.problems.length) report.failures++;
    report.chromosomes.push(r);
    console.log(`${build} ${chr}: rows ${r.rows}, ${(r.bytes / 1e6).toFixed(1)} MB${r.problems.length ? "  PROBLEMS: " + r.problems.join("; ") : ""}`);
  }
}
report.totals = Object.fromEntries([...new Set(report.chromosomes.map(c => c.build))].map(b => {
  const cs = report.chromosomes.filter(c => c.build === b);
  return [b, { chromosomes: cs.length, rows: cs.reduce((a, c) => a + (c.rows ?? 0), 0), gb: +(cs.reduce((a, c) => a + (c.bytes ?? 0), 0) / 1e9).toFixed(2) }];
}));
console.log(JSON.stringify({ totals: report.totals, failures: report.failures }));
if (outPath) writeFileSync(outPath, JSON.stringify(report, null, 2));
process.exitCode = report.failures ? 1 : 0;
