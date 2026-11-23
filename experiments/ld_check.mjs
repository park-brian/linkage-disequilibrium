// End-to-end check: LD (r², D′) computed from the v2 Parquet bitsets must equal the reference LDlink math
// (scripts/profile/ldmatrix_reference.mjs) run directly on the VCF, for the same SNPs and population.
// usage: node ld_check.mjs <v2 dir> <vcf.gz> <tbi> <vcfChromName> <panelFile> <windowStart> [n]
import { readFileSync, writeFileSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { asyncBufferFromFile, parquetMetadataAsync, parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { sortedLookup } from "../src/parquet.mjs";

const [dir, vcf, tbi, vcfChrom, panelFile, startArg, nArg] = process.argv.slice(2);
const START = +startArg, N = +(nArg || 25);
const file = await asyncBufferFromFile(`${dir}/common_bits_p64k.parquet`), metadata = await parquetMetadataAsync(file);
const kv = Object.fromEntries(metadata.key_value_metadata.map(x => [x.key, x.value]));
const samples = kv.samples.split(","), haploidMales = new Set((kv.male_haploid_samples || "").split(",").filter(Boolean).map(Number));
const panel = new Set(readFileSync(panelFile, "utf8").split(/\s+/).filter(Boolean));

// candidate SNVs: unique position, biallelic, in [START, START+2Mb), evenly spaced
let rows = [];
await parquetRead({ file, metadata, compressors, utf8: false, columns: ["pos", "ref", "alt", "multiallelic_source"], rowFormat: "object", onComplete: r => { rows = r; } });
const count = new Map(); for (const r of rows) count.set(r.pos, (count.get(r.pos) || 0) + 1);
const cand = rows.filter(r => r.pos >= START && r.pos < START + 2e6 && count.get(r.pos) === 1 && !r.multiallelic_source && r.ref.length === 1 && r.alt.length === 1);
const step = Math.max(1, Math.floor(cand.length / N)), picks = cand.filter((_, i) => i % step === 0).slice(0, N).map(r => r.pos);

// --- Parquet side: bitsets → counts by allele index ---
const got = (await sortedLookup({ file, metadata, column: "pos", keys: picks, columns: ["pos", "haps", "missing", "male_haploid"], compressors, payloadColumn: "haps" })).rows;
got.sort((a, b) => a.pos - b.pos);
const popSamples = samples.map((s, i) => (panel.has(s) ? i : -1)).filter(i => i >= 0);
const bit = (b, k) => (b[k >> 3] >> (k & 7)) & 1;
function hapList(row) { // [haplotype index, allele] for valid haplotypes of the population
  const out = [];
  for (const s of popSamples) for (const k of [2 * s, 2 * s + 1]) {
    if (k === 2 * s + 1 && row.male_haploid && haploidMales.has(s)) continue; // haploid male: no 2nd haplotype
    if (row.missing && bit(row.missing, k)) continue;
    out.push([k, bit(row.haps, k)]);
  }
  return new Map(out);
}
const H = got.map(hapList);
const round3 = x => Math.round(x * 1000) / 1000;
function ld(a, b) {
  let A = 0, B = 0, C = 0, D = 0;
  for (const [k, x] of a) { if (!b.has(k)) continue; const y = b.get(k); if (!x && !y) A++; else if (!x && y) B++; else if (x && !y) C++; else D++; }
  const delta = A * D - B * C, Ms = (A + C) * (B + D) * (A + B) * (C + D);
  if (!Ms) return ["NA", "NA"];
  const dmax = delta < 0 ? Math.min((A + C) * (A + B), (B + D) * (C + D)) : Math.min((A + C) * (C + D), (A + B) * (B + D));
  return [round3(delta * delta / Ms), round3(Math.abs(delta / dmax))];
}
const mine = { r2: [], dprime: [] };
for (let i = 0; i < H.length; i++) { mine.r2.push([]); mine.dprime.push([]); for (let j = 0; j < H.length; j++) { const [r, d] = i === j ? [1, 1] : ld(H[i], H[j]); mine.r2[i].push(r); mine.dprime[i].push(d); } }

// --- reference side: LDlink math on the VCF ---
const q = got.map(r => `${vcfChrom}:${r.pos}:p${r.pos}`);
const ref = JSON.parse(execFileSync("node", ["../scripts/profile/ldmatrix_reference.mjs", "../../convert-vcf-json/node_modules", vcf, panelFile, "--", ...q],
  { env: { ...process.env, TBI: tbi }, maxBuffer: 1 << 28 }).toString());
// diagonal: reference computes self-pairs (1 unless monomorphic); compare off-diagonal only
let cmp = 0, bad = [];
for (let i = 0; i < got.length; i++) for (let j = 0; j < got.length; j++) {
  if (i === j) continue; cmp++;
  for (const k of ["r2", "dprime"]) if (String(mine[k][i][j]) !== String(ref[k][i][j])) bad.push({ i, j, k, mine: mine[k][i][j], ref: ref[k][i][j], pos: [got[i].pos, got[j].pos] });
}
const out = { dir, population: panelFile.split(/[\\/]/).pop(), snps: got.length, ref_snps: ref.snps.length, ref_samples: ref.n_samples, pop_samples: popSamples.length, pairs_compared: cmp, mismatches: bad.length, examples: bad.slice(0, 5) };
writeFileSync(`${dir}/ld_check_${out.population.replace(".txt", "")}.json`, JSON.stringify(out, null, 2));
console.log(JSON.stringify(out));
