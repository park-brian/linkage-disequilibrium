// SNPclip: the API's details text (POST /snpclip): every LDmatrix list at three threshold settings, a dense list, edges.
import { parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { snpclip } from "../../../src/parity/ldlink/snpclip.js";
import { matrixCases, post, asText } from "../lib.js";

/** labelled rsIDs of GRCh37 chr22 rows in [lo, hi], in file order (tiers in turn), read through the context's files
 *  (local or over HTTP) */
export async function denseIds(ctx, tiers, lo, hi) {
  const ids = [];
  for (const t of tiers) {
    const { file, metadata } = await ctx.planner.file(`genotypes/grch37/chr22.${t}.parquet`);
    await parquetRead({ file, metadata, compressors, columns: ["pos", "rsid"], rowFormat: "object",
      filter: { pos: { $gte: lo, $lte: hi } },
      onComplete: r => { for (const x of r) if (x.rsid != null) ids.push(`rs${x.rsid}`); } });
  }
  return ids;
}
export async function cases(ctx) {
  const cs = [];
  for (const c of matrixCases()) for (const [r2, maf] of [[0.1, 0.01], [0.5, 0.05], [0.8, 0]])
    cs.push({ name: `${c.name} r2≥${r2} maf≥${maf}`, build: c.build, pop: c.pop, snps: c.inputs.join("\n"), r2_threshold: r2,
      maf_threshold: maf });
  cs.push({ name: "chr22 1 Mb dense (1500 rsIDs) EUR", build: "grch37", pop: "EUR", snps: (await denseIds(ctx, ["common", "low"], 30e6,
    31e6)).slice(0, 1500).join("\n"), r2_threshold: 0.2, maf_threshold: 0.01 });
  const edge = (name, snps, extra = {}) => cs.push({ name, build: "grch37", pop: "CEU", snps: snps.join("\n"), r2_threshold: 0.1,
    maf_threshold: 0.01, ...extra });
  edge("different chromosomes", ["rs3", "rs4", "rs12465506"]);
  edge("junk lines and duplicates", ["  rs3", "RS4", "foo", "", "rs3", "chr13:32446842", "rs", "rs9999999999", "chrx"]);
  edge("chrY GRCh38", ["rs2032597", "rs2032599", "rs3"], { build: "grch38", pop: "EUR" });
  return cs;
}
export const ldlink = cs => Promise.all(cs.map(async c => asText(await post("snpclip", { snps: c.snps, pop: c.pop,
  r2_threshold: String(c.r2_threshold), maf_threshold: String(c.maf_threshold), genome_build: c.build }))));
export const ours = (ctx, c) => snpclip(ctx, c);
