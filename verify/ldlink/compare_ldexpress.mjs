// LDexpress: local LDlink (REST POST /LDlinkRest/ldexpress → text) vs modules/ldexpress.mjs, line by line.
// Cases: chunks of ≤ 5 inputs from every LDmatrix case, five tissues (r² ≥ 0.1, p < 0.1), variants and edges.
// usage: node verify/ldlink/compare_ldexpress.mjs <ldmatrix_cases.json> <out.json> [ldlinkBase]
import { readFileSync, writeFileSync } from "node:fs";
import { asyncBufferFromFile, cachedAsyncBuffer } from "hyparquet";
import { loadKernel, KERNEL_WASM } from "../../src/ld/kernel.mjs";
import { createContext } from "../../src/ldlink/context.mjs";
import { ldexpress } from "../../src/ldlink/ldexpress.mjs";
import * as CONFIG from "../../tools/config.mjs";
import { useDataCwd, repo } from "../../tools/config.mjs";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.mjs)

const [casesPath, outPath, base = CONFIG.LDLINK_URL] = process.argv.slice(2);
const ROOT = process.env.ROOT ?? CONFIG.RELEASE, PANELS = CONFIG.PANELS;
const ctx = await createContext({ root: ROOT, open: async p => cachedAsyncBuffer(await asyncBufferFromFile(p)),
  readPanel: async p => readFileSync(`${PANELS}/${p}.txt`, "utf8").split(/\s+/).filter(Boolean), kernel: await loadKernel(readFileSync(KERNEL_WASM)) });
const cases = [];
for (const c of JSON.parse(readFileSync(casesPath, "utf8"))) {
  const ins = c.inputs.map(s => s.replace(/\r$/, ""));
  for (let i = 0; i < ins.length; i += 10) cases.push({ name: `${c.name} [${i}..]`, build: c.build, pop: c.pop, snps: ins.slice(i, i + 10).join("\n") });
}
const edge = (name, snps, extra = {}) => cases.push({ name, build: "grch37", pop: "CEU", snps: snps.join("\n"), ...extra });
const known = ["rs12465506", "rs3", "rs7837688", "rs4242382", "rs6983267"];
edge("GWAS hits r2>=0.1", known);
edge("GWAS hits r2>=0.8", known, { r2_d_threshold: 0.8 });
edge("GWAS hits p<1e-5, all tissues", known, { p_threshold: 1e-5, tissues: "ALL" });
edge("blank line", ["rs3", "", "rs4"]); edge("too many", Array.from({ length: 11 }, (_, i) => `rs${i + 3}`));
edge("GWAS hits D'>=0.9", known, { r2_d: "d", r2_d_threshold: 0.9 });
edge("GWAS hits window 50 kb", known, { window: 50000 });
edge("GWAS hits EUR+EAS", known, { pop: "EUR+EAS" });
edge("GWAS hits GRCh38", known, { build: "grch38" });
edge("GWAS hits High Coverage", known, { build: "grch38_high_coverage", pop: "YRI" });
edge("junk", ["foo", "rs9999999999", "chr13:32446842", "rs3"]);
edge("nothing found", ["rs553970249"]);
edge("bad window", ["rs3"], { window: 2000000 });
const TISSUES = "Adipose_Subcutaneous+Whole_Blood+Testis+Lung+Brain_Cerebellum";

const report = { base, root: ROOT, cases: [], totals: { cases: 0, identical: 0, differ: 0 } };
for (const c of cases) {
  const args = { r2_d: "r2", r2_d_threshold: 0.1, p_threshold: 0.1, window: 500000, tissues: TISSUES, ...c };
  if (args.tissues === "ALL") args.tissues = (await ctx.gtexTissueIds()).join("+");
  const res = await fetch(`${base}/LDlinkRest/ldexpress`, { method: "POST", headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ snps: c.snps, pop: c.pop, tissues: args.tissues, r2_d: args.r2_d, r2_d_threshold: String(args.r2_d_threshold), p_threshold: String(args.p_threshold), window: String(args.window), genome_build: c.build }) });
  const body = await res.text(), l = body.trimStart().startsWith("{") ? { error: JSON.parse(body).error } : { text: body };
  const o = await ldexpress(ctx, args), same = JSON.stringify(l) === JSON.stringify(o);
  let detail = {};
  if (!same && l.text && o.text) {
    const a = l.text.split("\n"), b = o.text.split("\n"), d = [];
    for (let i = 0; i < Math.max(a.length, b.length); i++) if (a[i] !== b[i]) d.push({ line: i, ldlink: a[i], ours: b[i] });
    detail = { lines: [a.length, b.length], differing_lines: d.length, examples: d.slice(0, 5) };
  } else if (!same) detail = { ldlink: l.error ?? l.text?.slice(0, 300), ours: o.error ?? o.text?.slice(0, 300) };
  report.totals.cases++; report.totals[same ? "identical" : "differ"]++;
  report.cases.push({ name: c.name, verdict: same ? "identical" : "differ", lines: l.text ? l.text.split("\n").length : 0, ...detail });
  console.log(JSON.stringify(report.cases.at(-1)).slice(0, 1500));
}
writeFileSync(outPath, JSON.stringify(report, null, 2));
console.log(JSON.stringify(report.totals));
