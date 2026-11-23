// LDexpress: the API text (POST /ldexpress) over five tissues (all 54 in one case), chunks of every LDmatrix list,
// known hits under variants, and input edges.
import { ldexpress } from "../../../src/ldlink/ldexpress.mjs";
import { matrixCases, clean, post, asText } from "../lib.mjs";

const TISSUES = "Adipose_Subcutaneous+Whole_Blood+Testis+Lung+Brain_Cerebellum";
export async function cases(ctx) {
  const cs = [];
  for (const c of matrixCases()) {
    const ins = c.inputs.map(clean);
    for (let i = 0; i < ins.length; i += 10) cs.push({ name: `${c.name} [${i}..]`, build: c.build, pop: c.pop, snps: ins.slice(i, i + 10).join("\n") });
  }
  const edge = (name, snps, extra = {}) => cs.push({ name, build: "grch37", pop: "CEU", snps: snps.join("\n"), ...extra });
  const known = ["rs12465506", "rs3", "rs7837688", "rs4242382", "rs6983267"];
  edge("GWAS hits r2>=0.1", known);
  edge("GWAS hits r2>=0.8", known, { r2_d_threshold: 0.8 });
  edge("GWAS hits p<1e-5, all tissues", known, { p_threshold: 1e-5, tissues: (await ctx.gtexTissueIds()).join("+") });
  edge("blank line", ["rs3", "", "rs4"]); edge("too many", Array.from({ length: 11 }, (_, i) => `rs${i + 3}`));
  edge("GWAS hits D'>=0.9", known, { r2_d: "d", r2_d_threshold: 0.9 });
  edge("GWAS hits window 50 kb", known, { window: 50000 });
  edge("GWAS hits EUR+EAS", known, { pop: "EUR+EAS" });
  edge("GWAS hits GRCh38", known, { build: "grch38" });
  edge("GWAS hits High Coverage", known, { build: "grch38_high_coverage", pop: "YRI" });
  edge("junk", ["foo", "rs9999999999", "chr13:32446842", "rs3"]);
  edge("nothing found", ["rs553970249"]);
  edge("bad window", ["rs3"], { window: 2000000 });
  return cs.map(c => ({ r2_d: "r2", r2_d_threshold: 0.1, p_threshold: 0.1, window: 500000, tissues: TISSUES, ...c }));
}
export const ldlink = cs => Promise.all(cs.map(async c => asText(await post("ldexpress", { snps: c.snps, pop: c.pop, tissues: c.tissues, r2_d: c.r2_d,
  r2_d_threshold: String(c.r2_d_threshold), p_threshold: String(c.p_threshold), window: String(c.window), genome_build: c.build }))));
export const ours = (ctx, c) => ldexpress(ctx, c);
