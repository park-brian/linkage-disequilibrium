// Validate the planner: (1) hybrid rsID resolution gives exactly the direct (LDlink-parity) answer, with its bytes
// and its predicted bytes; (2) select (parity / allele) picks the same rows as on a full-chromosome read; (3)
// non-canonical inputs resolve like the direct resolver. Every scenario runs cold (new planner, empty caches).
// usage: node verify/data/planner.js <releaseRoot>
import { openSync, readSync, fstatSync, writeFileSync, existsSync } from "node:fs";
import { asyncBufferFromFile, parquetMetadataAsync, parquetRead, cachedAsyncBuffer } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { openPlanner } from "../../src/parity/planner.js";
import { pickAllele } from "../../src/core/genotypes.js";
import { pickParity } from "../../src/parity/select.js";
import { useDataCwd } from "../../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

const [root] = process.argv.slice(2);
const group = path => (/genotypes\//.test(path) ? "genotypes" : /dbsnp\/rsid|duplicated/.test(path) ? "rsids" : /dbsnp\//.test(path) ? "coords" : "other");
function counter() {
  const stats = { requests: 0, bytes: 0, by: { genotypes: 0, rsids: 0, coords: 0, other: 0 } };
  const openFile = async path => {
    const fd = openSync(path, "r"), byteLength = fstatSync(fd).size, gname = group(path);
    return cachedAsyncBuffer({ byteLength, slice(s, e = byteLength) {
      stats.requests++; stats.bytes += e - s; stats.by[gname] += e - s; const b = Buffer.allocUnsafe(e - s); readSync(fd, b, 0, e - s, s);
      return Promise.resolve(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength)); } }, { minSize: 0 });
  };
  return { stats, openFile };
}
async function run(label, q) {
  const c = counter(), pl = await openPlanner(root, c.openFile), base = { ...c.stats, by: { ...c.stats.by } };
    // manifest read counted separately
  const t0 = performance.now(), r = await pl.query(q), ms = performance.now() - t0;
  const measured = c.stats.bytes - base.bytes;
  const o = { label, inputs: q.inputs?.length ?? 0, found: r.variants.filter(v => v?.found).length, rsid: r.plan.rsid,
    measured_MB: +(measured / 1e6).toFixed(2), requests: c.stats.requests - base.requests, ms: Math.round(ms),
    by: Object.fromEntries(["genotypes", "rsids", "coords"].map(k => [k, +((c.stats.by[k] - (base.by?.[k] ?? 0)) / 1e6).toFixed(2)])) };
  console.log(JSON.stringify(o));
  return { o, r };
}

// ---- inputs drawn from the data ----
async function labels(path, lo, hi) {
  const file = await asyncBufferFromFile(`${root}/${path}`), metadata = await parquetMetadataAsync(file);
  let rows = []; await parquetRead({ file, metadata, compressors, rowFormat: "object", columns: ["pos", "rsid", "rsid_source",
    "rsid_pos_delta"], filter: { pos: { $gte: lo, $lte: hi } }, onComplete: r => { rows = r; } });
  return rows.filter(r => r.rsid != null);
}
let seed = 5; const rnd = () => ((seed = (seed * 1103515245 + 12345) >>> 0) / 2 ** 32);
const pick = (arr, n) => Array.from({ length: n }, () => arr[Math.floor(rnd() * arr.length)]);
const g = "genotypes/grch37";
const clustered = [...await labels(`${g}/chr22.common.parquet`, 30e6, 32e6), ...await labels(`${g}/chr22.low.parquet`, 30e6, 32e6),
  ...await labels(`${g}/chr22.rare.parquet`, 30e6, 32e6)];
const chr22all = await labels(`${g}/chr22.common.parquet`, 0, 1e9), chrXall = await labels(`${g}/chrX.common.parquet`, 0, 1e9);
const dupFile = await asyncBufferFromFile(`${root}/dbsnp/duplicated_rsids.parquet`);
const dupIds = (await (async () => { let r = []; await parquetRead({ file: dupFile, metadata: await parquetMetadataAsync(dupFile),
  compressors, rowFormat: "object", columns: ["rsid"], rowStart: 0, rowEnd: 2000, onComplete: x => { r = x; } }); return r;
  })()).map(x => `rs${x.rsid}`);
const extras = [...pick(dupIds, 20), "rs52815562", "rs2179", "rs9999999999", "rs1", "rs999"]; // duplicates, merged, absent

const results = [];
// every field a caller sees (not `via`, which says how it was resolved)
const key = v => JSON.stringify(v == null ? null : { found: !!v.found, rsid: v.rsid ?? null, ref_rsid: v.ref_rsid ?? null,
  chrom: v.chrom ?? null, pos: v.pos ?? null, warning: v.warning ?? null });
// 1) hybrid vs direct must agree exactly (found, chrom, pos) for every input
for (const [label, inputs] of [
  ["clustered 300 rsIDs in 2 Mb (+extras)", [...pick(clustered, 300).map(r => `rs${r.rsid}`), ...extras]],
  ["clustered 3000 rsIDs in 2 Mb (+extras)", [...pick(clustered, 3000).map(r => `rs${r.rsid}`), ...extras]],
  ["scattered 300 rsIDs chr22 + chrX", [...pick(chr22all, 150), ...pick(chrXall, 150)].map(r => `rs${r.rsid}`)],
]) {
  const auto = await run(`${label} [auto]`, { build: "grch37", inputs, tiers: ["common", "low", "rare"] });
  const dir = await run(`${label} [force direct]`, { build: "grch37", inputs, tiers: ["common", "low", "rare"], forceRsid: "direct" });
  const diff = inputs.map((x, i) => [x, auto.r.variants[i], dir.r.variants[i]]).filter(([, a, b]) => key(a) !== key(b)).map(([x, a,
    b]) => ({ input: x, auto: key(a), direct: key(b) }));
  results.push({ check: "hybrid == direct", label, inputs: inputs.length, strategy: auto.o.rsid?.strategy, differences: diff.length,
    examples: diff.slice(0, 5), by: auto.o.by,
    auto_MB: auto.o.measured_MB, direct_MB: dir.o.measured_MB });
  console.log(JSON.stringify(results.at(-1)));
}
// 2) select: the planner's per-variant choice == select_rows on the WHOLE chromosome (all tiers), independent of
//    the planner's overlap windows, tier files and merged-ID handling
const SELCOLS = ["pos", "rsid", "ref", "alt", "vcf_order", "multiallelic_source", "allele_is_label", "allele_rsids"];
// the reference reads each variant's OWN chromosome in full (inputs include e.g. rs52815562 on chr11), lazily, once
const wholeByChrom = new Map();
async function wholeChrom(chrom) {
  if (wholeByChrom.has(chrom)) return wholeByChrom.get(chrom);
  const rows = [];
  for (const t of ["common", "low", "rare"]) {
    const path = `${root}/${g}/chr${chrom}.${t}.parquet`; if (!existsSync(path)) continue;
    const file = await asyncBufferFromFile(path);
    await parquetRead({ file, metadata: await parquetMetadataAsync(file), compressors, utf8: false, rowFormat: "object", columns: SELCOLS,
      onComplete: r => { for (const x of r) rows.push(x); } });
  }
  rows.sort((a, b) => a.pos - b.pos); wholeByChrom.set(chrom, rows); return rows;
}
const reach = 1000, windowOf = (whole, p) => { let a = 0, b = whole.length; while (a < b) { const h = (a + b) >> 1;
  if (whole[h].pos < p - reach) a = h + 1; else b = h; } const o = []; for (let j = a; j < whole.length && whole[j].pos <= p;
  j++) o.push(whole[j]); return o; };
const rowKey = r => (r ? `${r.pos}:${r.ref}>${r.alt}#${r.vcf_order}` : null);
const indelLabels = clustered.filter(r => r.rsid_pos_delta === 1);
const selInputs = [...pick(clustered, 400).map(r => `rs${r.rsid}`), ...pick(indelLabels, 100).map(r => `rs${r.rsid}`), "rs52815562",
  ...pick(clustered, 100).map(r => `chr22:${r.pos}`)];
for (const mode of ["parity", "allele"]) {
  for (const forceRsid of [undefined, "direct"]) {
    const r = await run(`select ${mode} [${forceRsid ?? "auto"}]`, { build: "grch37", inputs: selInputs, select: mode, forceRsid });
    const bad = [];
    for (const [i, v] of r.r.variants.entries()) {
      if (v?.pos == null) { if (v?.selected) bad.push({ input: selInputs[i], why: "selected without position" }); continue; }
      const w = windowOf(await wholeChrom(v.chrom), v.pos), e = mode === "parity" ? pickParity(w, v.pos) : pickAllele(w, v.pos,
        v.ref_rsid ?? v.rsid ?? "");
      const got = v.selected;
      if (!got || got.status !== e.status || rowKey(got.row) !== rowKey(e.row)) bad.push({ input: selInputs[i], chrom: v.chrom, pos: v.pos,
        got: got && [got.status, rowKey(got.row)], want: [e.status, rowKey(e.row)] });
    }
    const statuses = {};
      for (const v of r.r.variants) if (v?.selected) statuses[v.selected.status] = (statuses[v.selected.status] ?? 0) + 1;
    results.push({ check: "select == select_rows(whole chromosome)", label: `${mode} [${forceRsid ?? "auto"}]`, inputs: selInputs.length,
      differences: bad.length, statuses, examples: bad.slice(0, 5) });
    console.log(JSON.stringify(results.at(-1)));
  }
}
// 3) non-canonical inputs: planner (auto) == direct resolver for each, on every field
const someRs = pick(clustered, 3).map(r => String(r.rsid)), somePos = pick(clustered, 3).map(r => r.pos);
const odd = [
  ...someRs.flatMap(id => [`RS${id}`, ` rs${id} `, `rs${id}rs`, `rsrs${id}`, `${id}`]),
  ...somePos.flatMap(p => [`CHR22:${p}`, `chr22:${p} `, `22:${p}`, `chr22:${String(p).replace(/(\d{3})$/, "_$1")}`, `chr22:+${p}`,
    `chrchr22:${p}`]),
  "chr22:abc", "chr22:1.5", "chr22:0", "chr22:-5", "chr23:5", "chrx:2700000", "x:2700000", "rs52815562", "rs9999999999", "rs", "", "rs1 2",
];
for (const forceRsid of [undefined, "direct"]) {
  const a = await run(`non-canonical inputs [${forceRsid ?? "auto"}]`, { build: "grch37", inputs: odd, forceRsid });
  const c = counter(), pl = await openPlanner(root, c.openFile);
  // the resolver alone (no planner): LDlink's whole chain, resolve.js resolveInputs
  const { createResolver } = await import("../../src/parity/resolve.js");
  const want = await createResolver(root, c.openFile).resolveInputs(odd, "grch37");
  const bad = odd.map((x, i) => [x, key(a.r.variants[i]), key(want[i])]).filter(([, g, w]) => g !== w).map(([input, got, want]) => ({ input,
    got, want }));
  results.push({ check: "non-canonical inputs == resolver", label: forceRsid ?? "auto", inputs: odd.length, differences: bad.length,
    examples: bad.slice(0, 8) });
  console.log(JSON.stringify(results.at(-1)));
}
writeFileSync("verify/planner_report.json", JSON.stringify(results, null, 2)); // $LD_DATA/verify, not the release
if (results.some(r => r.differences > 0)) process.exitCode = 1;
