// SNPchip: local LDlink (REST POST /LDlinkRest/snpchip → details text) vs modules/snpchip.mjs, line by line.
// Platform requests: all codes, Illumina only, Affymetrix only, one code, and "" (which crashes LDlink).
// usage (from experiments/): ROOT=out/release_v4 node compare_snpchip.mjs <ldmatrix_cases.json> <out.json> [ldlinkBase]
import { readFileSync, writeFileSync } from "node:fs";
import { asyncBufferFromFile, cachedAsyncBuffer, parquetMetadataAsync, parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { loadKernel, KERNEL_WASM } from "../src/ld/kernel.mjs";
import { createContext } from "../src/ldlink/context.mjs";
import { snpchip } from "../src/ldlink/snpchip.mjs";

const [casesPath, outPath, base = "http://127.0.0.1:8090"] = process.argv.slice(2);
const ROOT = process.env.ROOT ?? "out/release_v4", PANELS = "D:/Projects/LDlink/db/data/1000G_population_samples";
const ctx = await createContext({ root: ROOT, open: async p => cachedAsyncBuffer(await asyncBufferFromFile(p)),
  readPanel: async p => readFileSync(`${PANELS}/${p}.txt`, "utf8").split(/\s+/).filter(Boolean), kernel: await loadKernel(readFileSync(KERNEL_WASM)) });
const plats = await ctx.platforms(), codes = plats.map(p => p.code);
const PLATS = { all: codes.join("+"), illumina: plats.filter(p => /^Illumina/.test(p.platform)).map(p => p.code).join("+"),
  affymetrix: plats.filter(p => /^Affymetrix/.test(p.platform)).map(p => p.code).join("+"), one: codes[0] };

const cases = [];
for (const c of JSON.parse(readFileSync(casesPath, "utf8"))) for (const [k, q] of Object.entries(PLATS)) cases.push({ name: `${c.name} [${k}]`, build: c.build, snps: c.inputs.join("\n"), platforms: q });
{ const ids = []; const f = await asyncBufferFromFile(`${ROOT}/genotypes/grch37/chr22.common.parquet`);
  await parquetRead({ file: f, metadata: await parquetMetadataAsync(f), compressors, columns: ["pos", "rsid"], rowFormat: "object", filter: { pos: { $gte: 30e6, $lte: 33e6 } }, onComplete: r => { for (const x of r) if (x.rsid != null) ids.push(`rs${x.rsid}`); } });
  for (const b of ["grch37", "grch38"]) cases.push({ name: `chr22 dense 2000 [all] ${b}`, build: b, snps: ids.slice(0, 2000).join("\n"), platforms: PLATS.all }); }
const edge = (name, snps, platforms = PLATS.all, build = "grch37") => cases.push({ name, build, snps: snps.join("\n"), platforms });
edge("several chromosomes, unsorted", ["rs12465506", "rs3", "rs4", "rs4892901", "chr13:32446842", "foo", "rs9999999999"]);
edge("empty platform request", ["rs3", "rs4"], "");
edge("nothing valid", ["foo", "rs9999999999"]);
edge("unknown platform code", ["rs3", "rs4"], "NOPE");

const report = { base, root: ROOT, cases: [], totals: { cases: 0, identical: 0, differ: 0 } };
for (const c of cases) {
  const res = await fetch(`${base}/LDlinkRest/snpchip`, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ snps: c.snps, platforms: c.platforms, genome_build: c.build }) });
  const body = await res.text(), l = body.trimStart().startsWith("{") ? { error: JSON.parse(body).error } : { text: body };
  const o = await snpchip(ctx, c), same = JSON.stringify(l) === JSON.stringify(o);
  let detail = {};
  if (!same && l.text && o.text) {
    const a = l.text.split("\n"), b = o.text.split("\n"), d = [];
    for (let i = 0; i < Math.max(a.length, b.length); i++) if (a[i] !== b[i]) d.push({ line: i, ldlink: a[i], ours: b[i] });
    detail = { lines: [a.length, b.length], differing_lines: d.length, examples: d.slice(0, 5) };
  } else if (!same) detail = { ldlink: l.error ?? l.text?.slice(0, 300), ours: o.error ?? o.text?.slice(0, 300) };
  report.totals.cases++; report.totals[same ? "identical" : "differ"]++;
  report.cases.push({ name: c.name, verdict: same ? "identical" : "differ", lines: l.text ? l.text.split("\n").length : 0, ...detail });
  console.log(JSON.stringify(report.cases.at(-1)).slice(0, 900));
}
writeFileSync(outPath, JSON.stringify(report, null, 2));
console.log(JSON.stringify(report.totals));
