// LDhap: LDlink's calculate_hap called inside the local backend container (tools/ldlink-local/drive.py: full JSON with
// warnings, plus the API's text) vs modules/ldhap.mjs. Cases: chunks of ≤ 30 inputs from every LDmatrix case, plus edges.
// usage: node verify/ldlink/compare_ldhap.mjs <ldmatrix_cases.json> <out.json>
import { readFileSync, writeFileSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { asyncBufferFromFile, cachedAsyncBuffer } from "hyparquet";
import { loadKernel, KERNEL_WASM } from "../../src/ld/kernel.mjs";
import { createContext } from "../../src/ldlink/context.mjs";
import { ldhap, ldhapText } from "../../src/ldlink/ldhap.mjs";
import * as CONFIG from "../../tools/config.mjs";
import { useDataCwd, repo } from "../../tools/config.mjs";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.mjs)

const [casesPath, outPath] = process.argv.slice(2);
const LOCAL = CONFIG.LDLINK_HOME, PANELS = CONFIG.PANELS;
const ctx = await createContext({ root: process.env.ROOT ?? CONFIG.RELEASE, open: async p => cachedAsyncBuffer(await asyncBufferFromFile(p)),
  readPanel: async p => readFileSync(`${PANELS}/${p}.txt`, "utf8").split(/\s+/).filter(Boolean), kernel: await loadKernel(readFileSync(KERNEL_WASM)) });

const cases = [];
for (const c of JSON.parse(readFileSync(casesPath, "utf8"))) {
  const ins = c.inputs.map(s => s.replace(/\r$/, ""));
  for (let i = 0; i < ins.length; i += 10) cases.push({ name: `${c.name} [${i}..]`, build: c.build, pop: c.pop, snps: ins.slice(i, i + 10) });
  cases.push({ name: `${c.name} [first 30]`, build: c.build, pop: c.pop, snps: ins.slice(0, 30) });
}
const edge = (name, snps, pop = "CEU", build = "grch37") => cases.push({ name, build, pop, snps });
edge("over 30 inputs", Array.from({ length: 31 }, (_, i) => `rs${i + 3}`));
edge("different chromosomes", ["rs3", "rs4", "rs12465506"]);
edge("coordinate and its rsID", ["chr13:32446842", "rs3", "rs4"]);
edge("duplicates and blanks", ["rs3", "", "rs3", "RS4", " rs4 ", "rs3 extra"]);
edge("nothing valid", ["foo", "rs9999999999"]);
edge("single multi-allelic", ["rs62224611"]);
edge("> 1 Mb range", ["rs3", "rs9999", "rs4"]);
edge("bad population", ["rs3", "rs4"], "FOO");
edge("chrY on GRCh38", ["rs2032597", "rs2032599", "rs3"], "EUR", "grch38");
edge("chrY GRCh37 (haploid)", ["rs2032597", "rs2032599"], "EUR");

writeFileSync(`${CONFIG.LDLINK_TMP}/ldhap_cases.json`, JSON.stringify(cases.map(c => ({ module: "ldhap", ...c }))));
const lines = execFileSync("docker", ["compose", "exec", "-T", "backend", "python3.13", "-", "/data/tmp/ldhap_cases.json"],
  { cwd: LOCAL, input: readFileSync(repo("tools/ldlink-local/drive.py")), env: { ...process.env, MSYS_NO_PATHCONV: "1" }, maxBuffer: 1 << 28 }).toString().trim().split("\n");
const canon = o => (Array.isArray(o) ? o.map(canon) : o && typeof o === "object" ? Object.fromEntries(Object.keys(o).sort().map(k => [k, canon(o[k])])) : o);
const report = { cases: [], totals: { cases: 0, identical: 0, differ: 0, ldlink_exception: 0 } };
for (const [i, c] of cases.entries()) {
  const l = JSON.parse(lines[i]), o = await ldhap(ctx, c);
  let verdict;
  if (l.exception) verdict = "ldlink_exception";
  else {
    const jsonSame = JSON.stringify(canon(l.json)) === JSON.stringify(canon(o));
    const textSame = l.text === null ? o.error !== undefined : !o.error && ldhapText(o, c.build) === l.text;
    verdict = jsonSame && textSame ? "identical" : "differ";
  }
  report.totals.cases++; report.totals[verdict]++;
  report.cases.push({ name: c.name, build: c.build, pop: c.pop, inputs: c.snps.length, verdict, ...(verdict !== "identical" ? { ldlink: l, ours: o } : {}) });
  if (verdict !== "identical") console.log(JSON.stringify({ name: c.name, ldlink: l.exception ?? l.json, ours: o }).slice(0, 1500));
}
writeFileSync(outPath, JSON.stringify(report, null, 2));
console.log(JSON.stringify(report.totals));
