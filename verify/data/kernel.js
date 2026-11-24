// V1 (ROADMAP §2.2): the WASM LD kernel itself vs LDlink's LD math run on the VCF (verify/reference/ldmatrix_reference.js).
// Every kernel path is compared on the same SNPs and populations:
//   WASM + JS backends × band (blocked; pruned / unpruned; whole and split into worker ranges), tile (bytes),
//   band_masked (chunks with male haploid calls or '.' — chrX), and populations whose haplotypes are a leading
//   prefix of the sample order (the shape that hid R1).
// SNPs: unique-position, non-multiallelic rows from all three tiers (bitset and sparse decode), SNVs and indels.
// Indels: LDlink builds its 2×2 table from sorted() allele strings (LDmatrix.py:270-307); where that order differs
// from allele-index order the reference table is permuted (ROADMAP K4) — such pairs are reported separately.
// usage: node verify/data/kernel.js <genoPrefix e.g. release_v4/genotypes/grch37/chr22 (data-relative)> <vcf.gz> <tbi> <vcfChrom> <start> <out.json>
import { readFileSync, writeFileSync, mkdirSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { asyncBufferFromFile, parquetMetadataAsync, parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { loadKernel, jsKernel, KERNEL_WASM } from "../../src/ld/kernel.js";
import { sortedLookup } from "../../src/parquet.js";
import * as CONFIG from "../../tools/config.js";
import { useDataCwd, repo } from "../../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

const [prefix, vcf, tbi, vcfChrom, startArg, outPath] = process.argv.slice(2);
const START = +startArg, SPAN = 1_000_000, PANELS = CONFIG.PANELS;
const tmp = "kernel_check"; mkdirSync(tmp, { recursive: true });

// ---------- rows: picks from each tier ----------
const tiers = await Promise.all(["common", "low", "rare"].map(async t => {
  const file = await asyncBufferFromFile(`${prefix}.${t}.parquet`); return { t, file, metadata: await parquetMetadataAsync(file) };
}));
const kv = Object.fromEntries(tiers[0].metadata.key_value_metadata.map(x => [x.key, x.value]));
const samples = kv.samples.split(","), haploidMales = new Set((kv.male_haploid_samples || "").split(",").filter(Boolean).map(Number));
const all = [];
for (const { t, file, metadata } of tiers) await parquetRead({ file, metadata, compressors, columns: ["pos", "ref", "alt", "multiallelic_source"], rowFormat: "object",
  onComplete: r => { for (const x of r) all.push({ ...x, tier: t }); } });
const nAt = new Map(); for (const r of all) nAt.set(r.pos, (nAt.get(r.pos) ?? 0) + 1);
const ok = r => r.pos >= START && r.pos < START + SPAN && nAt.get(r.pos) === 1 && !r.multiallelic_source;
const want = { common: 22, low: 12, rare: 8 };
const picks = [];
for (const t of Object.keys(want)) {
  const c = all.filter(r => r.tier === t && ok(r)), snv = c.filter(r => r.ref.length === 1 && r.alt.length === 1), ind = c.filter(r => !(r.ref.length === 1 && r.alt.length === 1));
  const take = (a, k) => { const s = Math.max(1, Math.floor(a.length / Math.max(1, k))); return a.filter((_, i) => i % s === 0).slice(0, k); };
  const ni = Math.min(ind.length, Math.round(want[t] / 3));
  picks.push(...take(snv, want[t] - ni).map(r => ({ ...r })), ...take(ind, ni).map(r => ({ ...r })));
}
const rows = [];
for (const { t, file, metadata } of tiers) {
  const keys = picks.filter(p => p.tier === t).map(p => p.pos);
  const r = await sortedLookup({ file, metadata, column: "pos", keys, columns: ["pos", "ref", "alt", "haps", "missing", "male_haploid"], compressors });
  rows.push(...r.rows.map(x => ({ ...x, tier: t })));
}
rows.sort((a, b) => a.pos - b.pos);
if (rows.length !== picks.length) throw new Error(`looked up ${rows.length} of ${picks.length} picks`);

// ---------- populations ----------
const pops = [];
for (const p of ["CEU", "EUR", "ALL"]) pops.push({ name: p, file: `${PANELS}/${p}.txt` });
for (const k of [20, 32, 100]) { // leading prefixes of the VCF sample order: nh = 40, 64 (multiple of 32), 200
  const f = `${tmp}/prefix${k}.txt`; writeFileSync(f, samples.slice(0, k).join("\n") + "\n"); pops.push({ name: `prefix${k}`, file: f });
}

// ---------- kernel builds ----------
const kernels = { wasm: await loadKernel(readFileSync(KERNEL_WASM)), js: jsKernel() }; // both backends, one API
const round3 = x => Math.round(x * 1000) / 1000;
// float32 output: a value within 1e-6 of a 3-decimal rounding boundary may round either way vs the double reference
const nearBoundary = x => Math.abs(((x * 1000) % 1) - 0.5) < 1e-3;
function setAlleles(a1, a2) { // LDmatrix.py allele strings (as in the reference)
  if (a1.length === 1 && a2.length === 1) return [a1, a2];
  if (a1.length === 1) return ["-", a2.slice(1)];
  if (a2.length === 1) return [a1.slice(1), "-"];
  return [a1.slice(1), a2.slice(1)];
}
// does LDlink's sorted-string 2×2 table equal the allele-index table for this pair?
function indexConsistent(i, j) {
  const [a1, a2] = setAlleles(rows[i].ref, rows[i].alt), [b1, b2] = setAlleles(rows[j].ref, rows[j].alt);
  const byIndex = [a1 + b1, a1 + b2, a2 + b1, a2 + b2];
  return new Set(byIndex).size === 4 && JSON.stringify([...byIndex].sort()) === JSON.stringify(byIndex);
}

const report = { genotypes: prefix, vcf: vcf.split(/[\\/]/).pop(), start: START, snps: rows.length,
  by_tier: Object.fromEntries(Object.keys(want).map(t => [t, rows.filter(r => r.tier === t).length])), indels: rows.filter(r => !(r.ref.length === 1 && r.alt.length === 1)).length,
  masked_rows: rows.filter(r => r.missing || r.male_haploid).length, populations: [] };
const q = rows.map(r => `${vcfChrom}:${r.pos}:p${r.pos}`);
for (const pop of pops) {
  const ref = JSON.parse(execFileSync("node", [repo("verify/reference/ldmatrix_reference.js"), CONFIG.REF_NODE_MODULES, vcf, pop.file, "--", ...q],
    { env: { ...process.env, TBI: tbi }, maxBuffer: 1 << 28 }).toString());
  if (ref.snps.length !== rows.length || ref.snps.some((s, i) => s.pos !== rows[i].pos)) throw new Error(`${pop.name}: reference returned ${ref.snps.length} SNPs (${ref.warnings.join("; ")})`);
  const panel = new Set(readFileSync(pop.file, "utf8").split(/\s+/).filter(Boolean));
  const popSamples = samples.map((s, i) => (panel.has(s) ? i : -1)).filter(i => i >= 0);
  const hapIndexes = popSamples.flatMap(s => [2 * s, 2 * s + 1]);
  const res = { population: pop.name, samples: popSamples.length, nh: hapIndexes.length, paths: {} };
  const n = rows.length, expect = (i, j, k) => ref[k][i][j];
  function compare(name, getPair) { // getPair(i, j) → { r2, dprime } (floats) or null (NA / not emitted)
    const c = { pairs: 0, match: 0, k4_divergent: 0, float32_boundary: 0, mismatch: 0, examples: [] };
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      c.pairs++;
      const g = getPair(i, j), er = expect(i, j, "r2"), ed = expect(i, j, "dprime");
      const same = g === null ? er === "NA" && ed === "NA" : round3(g.r2) === er && round3(g.dprime) === ed;
      if (same) { c.match++; continue; }
      if (!indexConsistent(i, j)) { c.k4_divergent++; continue; }
      if (g && (nearBoundary(g.r2) || nearBoundary(g.dprime)) && Math.abs(g.r2 - er) <= 0.0011 && Math.abs(g.dprime - ed) <= 0.0011) { c.float32_boundary++; continue; }
      c.mismatch++; if (c.examples.length < 5) c.examples.push({ i, j, pos: [rows[i].pos, rows[j].pos], alleles: [`${rows[i].ref}>${rows[i].alt}`, `${rows[j].ref}>${rows[j].alt}`], got: g, ref: [er, ed] });
    }
    res.paths[name] = c;
  }
  const asMap = b => { const m = new Map(); for (let k = 0; k < b.count; k++) m.set(`${b.i[k]},${b.j[k]}`, { r2: b.r2[k], dprime: b.dprime[k] }); return m; };
  for (const [kname, k] of Object.entries(kernels)) {
    const p = k.prepare(rows, { samples: samples.length, haploidMales, hapIndexes });
    res.masked = p.masked;
    const big = { maxDist: 1e9, maxOut: 1 << 14 };
    const full = asMap(k.band(p, { ...big, prune: true }));
    compare(`${kname}.band`, (i, j) => full.get(`${i},${j}`) ?? null);
    const unpruned = asMap(k.band(p, { ...big, prune: false }));
    compare(`${kname}.band_unpruned`, (i, j) => unpruned.get(`${i},${j}`) ?? null);
    if (!p.masked) {
      const blocked = asMap(k.band(p, { ...big, block: 4 }));
      compare(`${kname}.band_block4`, (i, j) => blocked.get(`${i},${j}`) ?? null);
      // worker split: rows [0, m) and [m, n) must together give the same pairs, each once
      const m = Math.floor(n / 2), a = k.band(p, { ...big, block: 4, rows: [0, m] }), b = k.band(p, { ...big, block: 4, rows: [m, n] });
      const split = new Map(); let dup = 0;
      for (const x of [a, b]) for (let t = 0; t < x.count; t++) { const key = `${x.i[t]},${x.j[t]}`; if (split.has(key)) dup++; split.set(key, { r2: x.r2[t], dprime: x.dprime[t] }); }
      compare(`${kname}.band_split`, (i, j) => split.get(`${i},${j}`) ?? null);
      res.paths[`${kname}.band_split`].duplicate_pairs = dup;
      const t = k.tile(p, 0, n, 0, n);
      const tc = { pairs: 0, within_quantization: 0, na_agree: 0, bad: 0, examples: [] };
      for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
        if (!indexConsistent(i, j)) continue;
        tc.pairs++;
        const vr = t.r2[i * n + j], vd = t.dprime[i * n + j], er = expect(i, j, "r2"), ed = expect(i, j, "dprime");
        if (vr === 255 || er === "NA") { if (vr === 255 && er === "NA") tc.na_agree++; else { tc.bad++; if (tc.examples.length < 5) tc.examples.push({ i, j, vr, er }); } continue; }
        if (Math.abs(vr / 250 - er) <= 0.0025 && Math.abs(vd / 250 - ed) <= 0.0025) tc.within_quantization++; else { tc.bad++; if (tc.examples.length < 5) tc.examples.push({ i, j, vr, vd, er, ed }); }
      }
      res.paths[`${kname}.tile`] = tc;
    }
  }
  { const a = kernels.wasm.band(kernels.wasm.prepare(rows, { samples: samples.length, haploidMales, hapIndexes }), { maxDist: 1e9 });
    const b = kernels.js.band(kernels.js.prepare(rows, { samples: samples.length, haploidMales, hapIndexes }), { maxDist: 1e9 });
    let same = a.count === b.count; for (let t = 0; same && t < a.count; t++) same = a.i[t] === b.i[t] && a.j[t] === b.j[t] && a.r2[t] === b.r2[t] && a.dprime[t] === b.dprime[t];
    res.backends_identical = same; if (!same) res.paths.backends = { bad: 1 }; kernels.wasm.reset(); }
  const worst = Object.values(res.paths).reduce((a, c) => a + (c.mismatch ?? 0) + (c.bad ?? 0) + (c.duplicate_pairs ?? 0), 0);
  res.failures = worst;
  report.populations.push(res);
  console.log(JSON.stringify({ population: pop.name, nh: res.nh, masked: res.masked, failures: worst,
    paths: Object.fromEntries(Object.entries(res.paths).map(([k, c]) => [k, c.mismatch ?? c.bad])), k4: res.paths["wasm.band"].k4_divergent, f32: res.paths["wasm.band"].float32_boundary }));
}
report.failures = report.populations.reduce((a, p) => a + p.failures, 0);
writeFileSync(outPath, JSON.stringify(report, null, 2));
if (report.failures > 0) process.exitCode = 1;
console.log(JSON.stringify({ snps: report.snps, by_tier: report.by_tier, indels: report.indels, masked_rows: report.masked_rows, failures: report.failures }));
