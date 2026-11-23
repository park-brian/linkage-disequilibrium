// Read-only: compare per-sample GTs for a few SNPs between GRCh37 phase 3 and GRCh38 High Coverage. VCF paths are relative to LD_VCF; @gmod/tabix comes from LD_REF_NODE_MODULES.
import { createRequire } from "node:module"; import { pathToFileURL } from "node:url"; import { resolve } from "node:path";
import { VCF, REF_NODE_MODULES } from "../../tools/config.js";
const require = createRequire(pathToFileURL(REF_NODE_MODULES + "/").href);
const { TabixIndexedFile } = require("@gmod/tabix");
async function get(path, tbi, chr, pos) {
  const t = new TabixIndexedFile({ path, tbiPath: tbi }); const h = (await t.getHeader()).trim().split("\n").pop().split("\t");
  let line; await t.getLines(chr, pos - 1, pos, l => { if (l.split("\t")[1] == pos) line = l; });
  const g = line.split("\t"); return Object.fromEntries(h.slice(9).map((s, i) => [s, g[9 + i]]));
}
const pairs = [[217894234, 217029511], [217936015, 217071292], [217926328, 217061605]];
for (const [p37, p38] of pairs) {
  const a = await get(VCF + "/GRCh37/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz", VCF + "/GRCh37/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz.tbi", "2", p37);
  const f = "CCDG_14151_B01_GRM_WGS_2020-08-05_chr2.filtered.shapeit2-duohmm-phased.vcf.gz";
  const c = await get(VCF + "/GRCh38_High_Coverage/" + f, VCF + "/GRCh38_High_Coverage/" + f + ".tbi", "chr2", p38);
  let same = 0, flip = 0, diff = 0;
  for (const s of Object.keys(a)) { if (!(s in c)) continue; if (a[s] === c[s]) same++; else if (a[s].split("|").reverse().join("|") === c[s]) flip++; else diff++; }
  console.log(p37, { same, phase_flipped: flip, different: diff });
}
