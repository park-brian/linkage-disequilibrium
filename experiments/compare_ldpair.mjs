// LDpair: local LDlink (REST GET /LDlinkRest/ldpair, json_out) vs modules/ldpair.mjs, the whole JSON per pair.
// Pairs: consecutive inputs of every LDmatrix case (indels, allele-mismatch rsIDs, coordinates, chrX, merged, 12-ALT
// site, three builds) plus error and edge cases. LDlink's tie order among equal haplotype counts is hash-randomized
// (see ldpair.mjs) and reported separately as "tie_order".
// LDpop (module=ldpop): the API's text table, verbatim; every 4th pair also with pop=ALL (all 32 populations).
// usage (from experiments/): node compare_ldpair.mjs <ldmatrix_cases.json> <out.json> [ldpair|ldpop] [ldlinkBase]
import { readFileSync, writeFileSync } from "node:fs";
import { asyncBufferFromFile, cachedAsyncBuffer } from "hyparquet";
import { loadKernel, KERNEL_WASM } from "../src/ld/kernel.mjs";
import { createContext } from "../src/ldlink/context.mjs";
import { ldpair } from "../src/ldlink/ldpair.mjs";
import { ldpop } from "../src/ldlink/ldpop.mjs";

const [casesPath, outPath, module = "ldpair", base = "http://127.0.0.1:8090"] = process.argv.slice(2);
const PANELS = "D:/Projects/LDlink/db/data/1000G_population_samples";
const ctx = await createContext({ root: process.env.ROOT ?? "out/release_v4", open: async p => cachedAsyncBuffer(await asyncBufferFromFile(p)),
  readPanel: async p => readFileSync(`${PANELS}/${p}.txt`, "utf8").split(/\s+/).filter(Boolean), kernel: await loadKernel(readFileSync(KERNEL_WASM)) });

const pairs = [];
for (const c of JSON.parse(readFileSync(casesPath, "utf8")))
  for (let i = 0; i + 1 < c.inputs.length; i++) pairs.push({ case: c.name, build: c.build, pop: c.pop, var1: c.inputs[i], var2: c.inputs[i + 1] });
const edge = (case_, var1, var2, pop = "CEU", build = "grch37") => pairs.push({ case: case_, build, pop, var1, var2 });
edge("different chromosomes", "rs3", "rs12465506");
edge("coordinate + upper case", "chr13:32446842", "RS4");
edge("invalid second", "rs3", "xx1"); edge("invalid first", "foo", "rs4"); edge("too short", "rs", "rs4");
edge("not in dbSNP", "rs9999999999", "rs4"); edge("not in dbSNP (second)", "rs3", "rs9999999999");
edge("bad population", "rs3", "rs4", "FOO"); edge("two populations", "rs3", "rs4", "CEU+YRI"); edge("ALL", "rs3", "rs4", "ALL");
edge("coordinate not in dbSNP", "chr13:1", "rs4"); edge("inner whitespace", "rs3 rs5", "rs4");
edge("chrY on GRCh38", "rs2032597", "rs2032599", "EUR", "grch38"); edge("chrY GRCh37 (haploid)", "rs2032597", "rs2032599", "EUR", "grch37");
edge("rare pair (monomorphic in CEU)", "rs553970249", "rs536672276");
if (module === "ldpop") pairs.push(...pairs.filter((_, i) => i % 4 === 0).map(q => ({ ...q, case: `${q.case} [ALL]`, pop: "ALL" })));

const canon = o => (Array.isArray(o) ? o.map(canon) : o && typeof o === "object" ? Object.fromEntries(Object.keys(o).sort().map(k => [k, canon(o[k])])) : o); // json.dumps(sort_keys=True)
const drop = o => (Array.isArray(o) ? o.map(x => { const { request, ...r } = x; return r; }) : o);
const tieOnly = (a, b) => { // same haplotype multiset and the same count sequence: only tied entries permuted
  const h = x => Object.values(x?.[0]?.haplotypes ?? {}), ha = h(a), hb = h(b), k = x => x.map(y => JSON.stringify(y)).sort().join();
  if (!ha.length || k(ha) !== k(hb) || ha.map(y => y.count).join() !== hb.map(y => y.count).join()) return false;
  const strip = x => JSON.stringify(x.map(({ haplotypes, ...r }) => r));
  return strip(a) === strip(b);
};
const report = { base, pairs: [], totals: { pairs: 0, identical: 0, tie_order: 0, differ: 0 } };
for (const q of pairs) {
  const url = `${base}/LDlinkRest/${module}?var1=${encodeURIComponent(q.var1)}&var2=${encodeURIComponent(q.var2)}&pop=${encodeURIComponent(q.pop)}&genome_build=${q.build}&json_out=true&r2_d=r2`;
  const body = await (await fetch(url)).text(), json = body.trimStart().startsWith("{") || body.trimStart().startsWith("[");
  const l = json ? canon(drop(JSON.parse(body))) : { text: body }, o = canon(await (module === "ldpop" ? ldpop : ldpair)(ctx, q));
  const same = JSON.stringify(l) === JSON.stringify(o), verdict = same ? "identical" : tieOnly(l, o) ? "tie_order" : "differ";
  report.totals.pairs++; report.totals[verdict]++;
  report.pairs.push({ ...q, verdict, ...(verdict === "differ" ? { ldlink: l, ours: o } : {}) });
  if (verdict === "differ") console.log(JSON.stringify({ ...q, ldlink: l, ours: o }).slice(0, 1500));
}
writeFileSync(outPath, JSON.stringify(report, null, 2));
console.log(JSON.stringify(report.totals));
