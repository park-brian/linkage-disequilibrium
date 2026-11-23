// LDmatrix: local LDlink (REST, ../ldlink-local) vs our pipeline, cell by cell.
// Ours = src/ldlink/ldmatrix.mjs (LDlink's input chain via the planner → parity rows → kernel exact → round(x, 3)).
// usage: node verify/ldlink/compare_ldmatrix.mjs <cases.json> <out.json> [ldlinkBase=$LDLINK_URL]
//   cases.json: [{ name, build, pop, inputs: [...] }]
import { readFileSync, writeFileSync } from "node:fs";
import { asyncBufferFromFile, cachedAsyncBuffer } from "hyparquet";
import { loadKernel, KERNEL_WASM } from "../../src/ld/kernel.mjs";
import { createContext } from "../../src/ldlink/context.mjs";
import { ldmatrix } from "../../src/ldlink/ldmatrix.mjs";
import * as CONFIG from "../../tools/config.mjs";
import { useDataCwd, repo } from "../../tools/config.mjs";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.mjs)

const [casesPath, outPath, base = CONFIG.LDLINK_URL] = process.argv.slice(2);
const ROOT = process.env.ROOT ?? CONFIG.RELEASE, PANELS = CONFIG.PANELS;
const ctx = await createContext({ root: ROOT, open: async p => cachedAsyncBuffer(await asyncBufferFromFile(p)),
  readPanel: async p => readFileSync(`${PANELS}/${p}.txt`, "utf8").split(/\s+/).filter(Boolean), kernel: await loadKernel(readFileSync(KERNEL_WASM)) });

async function ldlink(c, stat) {
  const q = `snps=${encodeURIComponent(c.inputs.join("\n"))}&pop=${encodeURIComponent(c.pop)}&r2_d=${stat}&genome_build=${c.build}`;
  const t = await (await fetch(`${base}/LDlinkRest/ldmatrix?${q}`)).text();
  if (t.trim().startsWith("{")) return { error: JSON.parse(t).error ?? t };
  const rows = t.trim().split(/\r?\n/).map(l => l.split("\t"));
  return { ids: rows[0].slice(1), cells: rows.slice(1).map(r => r.slice(1)) };
}
const cases = JSON.parse(readFileSync(casesPath, "utf8")), report = { base, cases: [] };
for (const c of cases) {
  const o = await ldmatrix(ctx, c), res = { name: c.name, build: c.build, pop: c.pop, inputs: c.inputs.length, ours_snps: o.ids.length };
  for (const stat of ["r2", "d"]) {
    const l = await ldlink(c, stat);
    if (l.error) { res[stat] = { ldlink_error: l.error, ours_snps: o.ids.length }; continue; }
    const idsSame = JSON.stringify(l.ids) === JSON.stringify(o.ids);
    let cells = 0, diff = 0; const ex = [];
    if (idsSame) for (let i = 0; i < l.ids.length; i++) for (let j = 0; j < l.ids.length; j++) {
      cells++; const a = l.cells[i][j], b = o[stat][i][j];
      if (a !== b && !(a !== "NA" && b !== "NA" && Number(a) === Number(b))) { diff++; if (ex.length < 5) ex.push({ i: l.ids[i], j: l.ids[j], ldlink: a, ours: b }); }
    }
    res[stat] = { ldlink_snps: l.ids.length, ids_identical: idsSame, cells, mismatches: diff, examples: ex,
      ...(idsSame ? {} : { only_ldlink: l.ids.filter(x => !o.ids.includes(x)), only_ours: o.ids.filter(x => !l.ids.includes(x)), order_differs: l.ids.length === o.ids.length }) };
  }
  report.cases.push(res);
  console.log(JSON.stringify(res, (k, v) => (k === "examples" && !v.length ? undefined : v)));
}
report.totals = { cases: report.cases.length, matrices: report.cases.length * 2,
  identical: report.cases.reduce((a, c) => a + ["r2", "d"].filter(s => c[s]?.ids_identical && c[s].mismatches === 0).length, 0) };
writeFileSync(outPath, JSON.stringify(report, null, 2));
console.log(JSON.stringify(report.totals));
