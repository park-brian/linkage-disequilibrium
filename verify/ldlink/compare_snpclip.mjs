// SNPclip: local LDlink (REST POST /LDlinkRest/snpclip → details text) vs modules/snpclip.mjs, line by line.
// Cases: every LDmatrix case at three threshold settings, a dense 1,500-rsID chr22 list, and edge cases.
// usage: node verify/ldlink/compare_snpclip.mjs <ldmatrix_cases.json> <out.json> [ldlinkBase]
import { readFileSync, writeFileSync } from "node:fs";
import { asyncBufferFromFile, cachedAsyncBuffer, parquetMetadataAsync, parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { loadKernel, KERNEL_WASM } from "../../src/ld/kernel.mjs";
import { createContext } from "../../src/ldlink/context.mjs";
import { snpclip } from "../../src/ldlink/snpclip.mjs";
import * as CONFIG from "../../tools/config.mjs";
import { useDataCwd, repo } from "../../tools/config.mjs";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.mjs)

const [casesPath, outPath, base = CONFIG.LDLINK_URL] = process.argv.slice(2);
const ROOT = process.env.ROOT ?? CONFIG.RELEASE, PANELS = CONFIG.PANELS;
const ctx = await createContext({ root: ROOT, open: async p => cachedAsyncBuffer(await asyncBufferFromFile(p)),
  readPanel: async p => readFileSync(`${PANELS}/${p}.txt`, "utf8").split(/\s+/).filter(Boolean), kernel: await loadKernel(readFileSync(KERNEL_WASM)) });

const cases = [];
for (const c of JSON.parse(readFileSync(casesPath, "utf8"))) for (const [r2, maf] of [[0.1, 0.01], [0.5, 0.05], [0.8, 0]])
  cases.push({ name: `${c.name} r2≥${r2} maf≥${maf}`, build: c.build, pop: c.pop, snps: c.inputs.join("\n"), r2_threshold: r2, maf_threshold: maf });
{ // dense list: labelled chr22 rows in 30–31 Mb, all tiers
  const ids = [];
  for (const t of ["common", "low"]) {
    const f = await asyncBufferFromFile(`${ROOT}/genotypes/grch37/chr22.${t}.parquet`);
    await parquetRead({ file: f, metadata: await parquetMetadataAsync(f), compressors, columns: ["pos", "rsid"], rowFormat: "object", filter: { pos: { $gte: 30e6, $lte: 31e6 } },
      onComplete: r => { for (const x of r) if (x.rsid != null) ids.push(`rs${x.rsid}`); } });
  }
  cases.push({ name: "chr22 1 Mb dense (1500 rsIDs) EUR", build: "grch37", pop: "EUR", snps: ids.slice(0, 1500).join("\n"), r2_threshold: 0.2, maf_threshold: 0.01 });
}
const edge = (name, snps, extra = {}) => cases.push({ name, build: "grch37", pop: "CEU", snps: snps.join("\n"), r2_threshold: 0.1, maf_threshold: 0.01, ...extra });
edge("different chromosomes", ["rs3", "rs4", "rs12465506"]);
edge("junk lines and duplicates", ["  rs3", "RS4", "foo", "", "rs3", "chr13:32446842", "rs", "rs9999999999", "chrx"]);
edge("chrY GRCh38", ["rs2032597", "rs2032599", "rs3"], { build: "grch38", pop: "EUR" });

const report = { base, root: ROOT, cases: [], totals: { cases: 0, identical: 0, differ: 0 } };
for (const c of cases) {
  const res = await fetch(`${base}/LDlinkRest/snpclip`, { method: "POST", headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ snps: c.snps, pop: c.pop, r2_threshold: String(c.r2_threshold), maf_threshold: String(c.maf_threshold), genome_build: c.build }) });
  const body = await res.text(), l = body.trimStart().startsWith("{") ? { error: JSON.parse(body).error } : { text: body };
  const o = await snpclip(ctx, c), same = JSON.stringify(l) === JSON.stringify(o);
  let detail = {};
  if (!same && l.text && o.text) {
    const a = l.text.split("\n"), b = o.text.split("\n"), d = [];
    for (let i = 0; i < Math.max(a.length, b.length); i++) if (a[i] !== b[i]) d.push({ line: i, ldlink: a[i], ours: b[i] });
    detail = { lines: [a.length, b.length], differing_lines: d.length, examples: d.slice(0, 5) };
  } else if (!same) detail = { ldlink: l.error ?? l.text?.slice(0, 300), ours: o.error ?? o.text?.slice(0, 300) };
  const kept = l.text ? (l.text.match(/Variant kept\./g) ?? []).length : 0;
  report.totals.cases++; report.totals[same ? "identical" : "differ"]++;
  report.cases.push({ name: c.name, verdict: same ? "identical" : "differ", lines: l.text ? l.text.split("\n").length - 2 : 0, kept, ...detail });
  console.log(JSON.stringify(report.cases.at(-1)).slice(0, 1000));
}
writeFileSync(outPath, JSON.stringify(report, null, 2));
console.log(JSON.stringify(report.totals));
