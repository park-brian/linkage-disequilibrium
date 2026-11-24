// Faithful port of LDlink's LDmatrix LD math (server/LDmatrix.py:191-372, LDcommon.py:486-499)
// used ONLY to produce reference test vectors from local 1000G VCF + tabix index. Read-only.
// usage: node verify/reference/ldmatrix_reference.js <node_modules dir with @gmod/tabix> <vcf.gz> <pop file>... -- chr:pos:rsid ...
import { createRequire } from "node:module";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { useDataCwd } from "../../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)
const [modDir, vcfPath, ...rest] = process.argv.slice(2);
const sep = rest.indexOf("--");
const popFiles = rest.slice(0, sep), queries = rest.slice(sep + 1).map(q => q.split(":"));
const require = createRequire(pathToFileURL(resolve(modDir) + "/").href);
const { TabixIndexedFile } = require("@gmod/tabix");
const tbi = new TabixIndexedFile({ path: vcfPath, tbiPath: process.env.TBI || vcfPath + ".tbi" });
const header = (await tbi.getHeader()).trim().split("\n").pop().split("\t");
const popIds = new Set(popFiles.flatMap(f => readFileSync(f, "utf8").split(/\s+/).filter(Boolean)));
const index = header.map((h, i) => (i >= 9 && popIds.has(h) ? i : -1)).filter(i => i >= 0);

function setAlleles(a1, a2) {
  if (a1.length === 1 && a2.length === 1) return [a1, a2];
  if (a1.length === 1) return ["-", a2.slice(1)];
  if (a2.length === 1) return [a1.slice(1), "-"];
  return [a1.slice(1), a2.slice(1)];
}
const hap1 = index.map(() => []), hap2 = index.map(() => []);
const snps = [], alleles = [], warnings = [];
queries.sort((a, b) => +a[1] - +b[1]);
for (const [chr, pos, rs] of queries) {
  const lines = [];
  await tbi.getLines(chr, +pos - 1, +pos, l => lines.push(l));
  // parse_vcf: last record wins, unless one has POS == query pos (LDcommon.py:592-605)
  let match = null;
  for (const l of lines) { match = l; if (l.split("\t", 2)[1] === pos) match = l; }
  if (!match) { warnings.push(`${rs} missing from 1000G`); continue; }
  const g = match.split("\t");
  if (g[3].includes(",") || g[4].includes(",")) { warnings.push(`${rs} variant is not a biallelic`); continue; }
  const [a1, a2] = setAlleles(g[3], g[4]);
  index.forEach((col, i) => {
    const gt = g[col];
    const m = { "0|0": [a1, a1], "0|1": [a1, a2], "1|0": [a2, a1], "1|1": [a2, a2], "0": [a1, "."], "1": [a2, "."] }[gt] || [".", "."];
    hap1[i].push(m[0]); hap2[i].push(m[1]);
  });
  snps.push({ rs, chr, pos: +pos, vcfId: g[2], ref: g[3], alt: g[4] }); alleles.push([a1, a2]);
}
const allHaps = hap1.concat(hap2), n = snps.length;
const round3 = x => Math.round(x * 1000) / 1000; // note: python round() is banker's rounding on exact .5 ties
const dprime = [...Array(n)].map(() => Array(n)), r2 = [...Array(n)].map(() => Array(n)), counts = {};
for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
  const hap = {};
  for (const h of allHaps) { const k = h[i] + h[j]; hap[k] = (hap[k] || 0) + 1; }
  for (const k of Object.keys(hap)) if (k.includes(".")) delete hap[k];
  if (Object.keys(hap).length !== 4)
    for (const h of [alleles[i][0] + alleles[j][0], alleles[i][0] + alleles[j][1], alleles[i][1] + alleles[j][0], alleles[i][1] + alleles[j][1]])
      if (!(h in hap)) hap[h] = 0;
  const ks = Object.keys(hap).sort(); // python sorted(): codepoint order, same as JS default sort for ASCII
  const [A, B, C, D] = ks.map(k => hap[k]);
  const delta = A * D - B * C, Ms = (A + C) * (B + D) * (A + B) * (C + D);
  let dp = "NA", rr = "NA";
  if (Ms !== 0) {
    dp = round3(Math.abs(delta / (delta < 0 ? Math.min((A + C) * (A + B), (B + D) * (C + D)) : Math.min((A + C) * (C + D), (A + B) * (B + D)))));
    rr = round3(delta ** 2 / Ms);
  }
  dprime[i][j] = dprime[j][i] = dp; r2[i][j] = r2[j][i] = rr;
  if (i !== j) counts[`${snps[i].rs}|${snps[j].rs}`] = Object.fromEntries(ks.map(k => [k, hap[k]]));
}
console.log(JSON.stringify({ vcf: vcfPath.split(/[\/]/).pop(), pops: popFiles.map(f => f.split(/[\/]/).pop()), n_samples: index.length,
  snps: snps.map((s, i) => ({ ...s, alleles: alleles[i].join("/") })), warnings, dprime, r2, counts }));
