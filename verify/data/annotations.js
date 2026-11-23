// Verify the annotation Parquet files reproduce LDlink's query answers (expected values from parity_annotations.sql,
// computed on the raw dump tables). usage: node verify_annotations.js <releaseRoot>
import { readFileSync, writeFileSync } from "node:fs";
import { asyncBufferFromFile, parquetMetadataAsync, parquetRead, cachedAsyncBuffer } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { sortedLookup } from "../../src/parquet.js";
import { useDataCwd } from "../../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

const [root] = process.argv.slice(2);
const load = f => JSON.parse(readFileSync(`annot_src/${f}`, "utf8"));
const files = new Map();
async function open(rel) {
  if (!files.has(rel)) files.set(rel, (async () => { const file = cachedAsyncBuffer(await asyncBufferFromFile(`${root}/${rel}`)); return { file, metadata: await parquetMetadataAsync(file) }; })());
  return files.get(rel);
}
const code = c => { const s = String(c).replace(/^chr/i, "").toUpperCase(); return s === "X" ? 23 : s === "Y" ? 24 : s === "MT" ? 25 : +s; };
const key = (c, p) => BigInt(code(c)) * 1_000_000_000n + BigInt(p);
async function range(rel, lo, hi, columns) {
  const { file, metadata } = await open(rel); let rows = [];
  await parquetRead({ file, metadata, compressors, rowFormat: "object", columns, filter: { key: { $gte: lo, $lte: hi } }, onComplete: r => { rows = r; } });
  return rows; // file order = (key, dump order)
}
async function points(rel, column, keys, columns) {
  const { file, metadata } = await open(rel);
  return (await sortedLookup({ file, metadata, column, keys, columns, compressors })).rows; // file order
}
const report = {};
const tally = (name, cases, bad) => { report[name] = { cases, mismatches: bad.length, examples: bad.slice(0, 3) }; console.log(JSON.stringify({ check: name, cases, mismatches: bad.length, example: bad[0] ?? null })); };

for (const b of ["37", "38"]) {
  const dir = `annotations/grch${b}`;
  { // RegulomeDB find_one
    const cases = load(`parity_regulome${b}.json`), bad = [];
    const rows = await points(`${dir}/regulome.parquet`, "key", cases.map(c => key(c.c, c.p)), ["key", "score"]);
    const first = new Map(); for (const r of rows) { const k = String(r.key); if (!first.has(k)) first.set(k, r.score); }
    for (const c of cases) { const got = first.get(String(key(c.c, c.p))) ?? null; if (got !== (c.expected ?? null)) bad.push({ ...c, got }); }
    tally(`regulome${b}`, cases.length, bad);
  }
  { // recombination ranges
    const cases = load(`parity_recomb${b}.json`), bad = [];
    for (const c of cases) {
      const rows = await range(`${dir}/recomb.parquet`, key(c.c, Math.max(0, c.lo)), key(c.c, c.hi), ["key", "rate"]);
      const pos = rows.map(r => Number(BigInt(r.key) % 1_000_000_000n)), sum = rows.reduce((a, r) => a + r.rate, 0);
      const exp = c.positions ?? [];
      if (rows.length !== c.n || pos.some((p, i) => p !== exp[i]) || Math.abs(sum - (c.rate_sum ?? 0)) > 1e-6 * Math.max(1, Math.abs(c.rate_sum ?? 0))) bad.push({ c: c.c, lo: c.lo, n: c.n, got: rows.length });
    }
    tally(`recomb${b}`, cases.length, bad);
  }
  { // refGene overlap (set comparison; LDlink's query has no sort)
    const cases = load(`parity_genes${b}.json`), bad = [];
    for (const c of cases) {
      // overlap: txStart <= hi and txEnd >= lo; key = chrom*1e9 + txStart, so scan keys up to hi and filter txEnd
      const rows = (await range(`${dir}/genes.parquet`, key(c.chrom, 0), key(c.chrom, c.hi), ["key", "name", "txStart", "txEnd"])).filter(r => r.txEnd >= c.lo);
      const got = rows.map(r => `${r.name}@${r.txStart}-${r.txEnd}`).sort(), exp = [...(c.transcripts ?? [])].sort();
      if (got.length !== exp.length || got.some((g, i) => g !== exp[i])) bad.push({ chrom: c.chrom, lo: c.lo, expected: exp.length, got: got.length });
    }
    tally(`genes${b}`, cases.length, bad);
  }
  { // GWAS window (index order = (position, dump order))
    const cases = load(`parity_gwas${b}.json`), bad = [];
    for (const c of cases) {
      const rows = await range(`${dir}/gwas.parquet`, key(c.c, c.lo), key(c.c, c.hi), ["key", "pubmedid", "snp_id_current", "trait"]);
      const got = rows.map(r => `${r.pubmedid ?? ""}|${r.snp_id_current ?? ""}|${r.trait ?? ""}`), exp = c.records ?? [];
      if (got.length !== exp.length || got.some((g, i) => g !== exp[i])) bad.push({ c: c.c, lo: c.lo, expected: exp.length, got: got.length });
    }
    tally(`gwas${b}`, cases.length, bad);
  }
  { // GTEx: match tissue set + position, inner join tissues and genes
    const cases = load(`parity_gtex${b}.json`), bad = [];
    const t = await open("annotations/gtex_tissues.parquet"); let trows = []; await parquetRead({ ...t, compressors, rowFormat: "object", onComplete: r => { trows = r; } });
    const tissueName = new Map(trows.map(r => [r.tissueSiteDetailId, r.tissueSiteDetail]));
    const g = await open("annotations/gtex_genes.parquet"); let grows = []; await parquetRead({ ...g, compressors, rowFormat: "object", columns: ["gene_id_v26", "gene_name_v26"], onComplete: r => { grows = r; } });
    const geneName = new Map(grows.map(r => [r.gene_id_v26, r.gene_name_v26]));
    const rows = await points(`${dir}/gtex_eqtl.parquet`, "key", cases.map(c => key(c.c, c.p)), ["key", "gene_id", "tissue", "slope", "pval_nominal"]);
    const at = new Map(); for (const r of rows) { const k = String(r.key); if (!at.has(k)) at.set(k, []); at.get(k).push(r); }
    for (const c of cases) {
      const T = new Set(c.tissues);
      const got = (at.get(String(key(c.c, c.p))) ?? []).filter(r => T.has(r.tissue) && tissueName.has(r.tissue) && geneName.has(r.gene_id))
        .map(r => [r.gene_id, r.tissue, r.slope, r.pval_nominal, geneName.get(r.gene_id), tissueName.get(r.tissue)]);
      const exp = (c.rows ?? []).map(s => { const p = s.split("|"); return [p[0], p[1], Number(p[2]), Number(p[3]), p[4], p[5]]; });
      if (got.length !== exp.length || got.some((x, i) => x.some((v, j) => v !== exp[i][j]))) bad.push({ c: c.c, p: c.p, expected: exp.length, got: got.length, first: [got[0], exp[0]] });
    }
    tally(`gtex${b}`, cases.length, bad);
  }
  { // SNPchip: every doc at the position, platform lists in dump order
    const cases = load(`parity_arrays${b}.json`), bad = [];
    const rows = await points(`${dir}/arrays.parquet`, "key", cases.map(c => key(c.c, c.p)), ["key", "platforms"]);
    const at = new Map(); for (const r of rows) { const k = String(r.key); if (!at.has(k)) at.set(k, []); at.get(k).push(r.platforms ?? ""); }
    for (const c of cases) { const got = at.get(String(key(c.c, c.p))) ?? [], exp = c.platforms ?? []; if (got.length !== exp.length || got.some((x, i) => x !== exp[i])) bad.push({ c: c.c, p: c.p, expected: exp, got }); }
    tally(`arrays${b}`, cases.length, bad);
  }
}
{ // FORGEdb find_one by rsID (rsid-sorted parity file)
  const cases = load("parity_forge.json"), bad = [];
  const ids = cases.map(c => c.snp_id.replace(/^rs/, "")).filter(x => /^\d+$/.test(x)).map(BigInt);
  const rows = await points("annotations/forge_rsid.parquet", "rsid", ids, ["rsid", "score"]);
  const first = new Map(); for (const r of rows) { const k = String(r.rsid); if (!first.has(k)) first.set(k, r.score); }
  for (const c of cases) { const got = first.get(c.snp_id.replace(/^rs/, "")); const exp = c.expected == null ? undefined : Number(c.expected); if (got !== exp) bad.push({ ...c, got }); }
  tally("forge_rsid", cases.length, bad);
}
{ // the region view's range readers (src/annotations.js forgeRange / regulomeRange) against hyparquet's own filtered read
  const { createAnnotations } = await import("../../src/annotations.js");
  const ann = createAnnotations(open), bad = [];
  const WINDOWS = [["grch37", "22", 30_000_000, 31_000_000], ["grch37", "X", 2_600_000, 3_200_000], ["grch38", "8", 127_700_000, 128_200_000], ["grch38", "1", 0, 600_000]];
  let cases = 0;
  for (const [b, c, lo, hi] of WINDOWS) {
    const base = `annotations/${b}`, klo = key(c, lo), khi = key(c, hi), pos = k => Number(BigInt(k) % 1_000_000_000n);
    const fExp = (await range(`${base}/forge.parquet`, klo, khi, ["key", "rsid", "score"])).map(r => ({ pos: pos(r.key), rsid: r.rsid, score: r.score }));
    const fGot = await ann.forgeRange(b, c, lo, hi);
    const seen = new Set(), rExp = [];
    for (const r of await range(`${base}/regulome.parquet`, klo, khi, ["key", "score"])) if (!seen.has(pos(r.key))) { seen.add(pos(r.key)); rExp.push({ pos: pos(r.key), score: r.score }); }
    const rGot = await ann.regulomeRange(b, c, lo, hi);
    cases += fExp.length + rExp.length;
    if (JSON.stringify(fGot, (k, v) => (typeof v === "bigint" ? String(v) : v)) !== JSON.stringify(fExp, (k, v) => (typeof v === "bigint" ? String(v) : v))) bad.push({ track: "forge", b, c, lo, hi, got: fGot.length, expected: fExp.length });
    if (JSON.stringify(rGot) !== JSON.stringify(rExp)) bad.push({ track: "regulome", b, c, lo, hi, got: rGot.length, expected: rExp.length });
    if (!fExp.length || !rExp.length) bad.push({ empty: true, b, c, lo, hi, forge: fExp.length, regulome: rExp.length });
  }
  tally("range_tracks", cases, bad);
}
writeFileSync("verify/annotations_report.json", JSON.stringify(report, null, 2)); // $LD_DATA/verify, not the release
if (Object.values(report).some(r => r.mismatches > 0)) process.exitCode = 1;
