// SNPchip: the API's details text (POST /snpchip) for four platform sets, dense lists on two builds, and edges.
import { snpchip } from "../../../src/ldlink/snpchip.mjs";
import { matrixCases, post, asText } from "../lib.mjs";
import { denseIds } from "./snpclip.mjs";

export async function cases(ctx) {
  const plats = await ctx.platforms(), codes = plats.map(p => p.code);
  const P = { all: codes.join("+"), illumina: plats.filter(p => /^Illumina/.test(p.platform)).map(p => p.code).join("+"),
    affymetrix: plats.filter(p => /^Affymetrix/.test(p.platform)).map(p => p.code).join("+"), one: codes[0] };
  const cs = [];
  for (const c of matrixCases()) for (const [k, q] of Object.entries(P)) cs.push({ name: `${c.name} [${k}]`, build: c.build, snps: c.inputs.join("\n"), platforms: q });
  const ids = (await denseIds(["common"], 30e6, 33e6)).slice(0, 2000).join("\n");
  for (const b of ["grch37", "grch38"]) cs.push({ name: `chr22 dense 2000 [all] ${b}`, build: b, snps: ids, platforms: P.all });
  const edge = (name, snps, platforms = P.all, build = "grch37") => cs.push({ name, build, snps: snps.join("\n"), platforms });
  edge("several chromosomes, unsorted", ["rs12465506", "rs3", "rs4", "rs4892901", "chr13:32446842", "foo", "rs9999999999"]);
  edge("empty platform request", ["rs3", "rs4"], "");
  edge("nothing valid", ["foo", "rs9999999999"]);
  edge("unknown platform code", ["rs3", "rs4"], "NOPE");
  return cs;
}
export const ldlink = cs => Promise.all(cs.map(async c => asText(await post("snpchip", { snps: c.snps, platforms: c.platforms, genome_build: c.build }))));
export const ours = (ctx, c) => snpchip(ctx, c);
