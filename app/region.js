// The region view (pyramid): the window's data read on the page, LD + pixels in the worker pool.
import { pyramid } from "./workers/pool.js";
import { join } from "./render.js";

const MAX_SPAN = 5_000_000, MAX_VARIANTS = 60000;
const num = s => +s.replace(/[,_]/g, "");

/** a = the form: region ("chr:start-end" or a variant), width, pop, build, r2_d, tiers, maxDist, threshold, workers, tissues */
export async function region(ctx, a, width) {
  let chrom, lo, hi;
  const m = /^\s*(?:chr)?([0-9]{1,2}|X|Y)\s*:\s*([\d,_]+)\s*-\s*([\d,_]+)\s*$/i.exec(a.region);
  if (m) { chrom = m[1].toUpperCase(); lo = num(m[2]); hi = num(m[3]); }
  else {
    const v = (await ctx.planner.query({ build: a.build, inputs: [a.region.trim()] })).variants[0];
    if (!v?.found || v.pos == null) return { error: `${a.region} is not a region (chr:start-end) or a variant in dbSNP 155.` };
    chrom = v.chrom; lo = Math.max(0, v.pos - Math.round(a.width / 2)); hi = v.pos + Math.round(a.width / 2);
  }
  if (!(hi > lo)) return { error: "The region's end must be after its start." };
  if (hi - lo > MAX_SPAN) return { error: "Regions up to 5 Mb, please." };
  if (!a.tiers.length) return { error: "Choose at least one MAF tier." };
  const t0 = performance.now();
  const [rows, genes, recomb, gwas, eqtl, smp, ids, forge, regulome] = await Promise.all([
    ctx.windowRows(a.build, chrom, lo, hi, a.tiers), ctx.genes(a.build, chrom, lo, hi), ctx.recomb(a.build, chrom, lo, hi),
    ctx.gwas(a.build, chrom, lo, hi), ctx.eqtlRange(a.build, chrom, lo, hi, a.tissues), ctx.samples(a.build, chrom),
      ctx.popIds(join(a.pop)),
    a.scores === "no" ? null : ctx.forgeRange(a.build, chrom, lo, hi), a.scores === "no" ? null : ctx.regulomeRange(a.build, chrom, lo,
      hi)]);
  const use = rows.filter(r => r.pos >= lo && r.pos <= hi).sort((x, y) => x.pos - y.pos || x.vcf_order - y.vcf_order);
  if (use.length < 2) return { error: "Fewer than two 1000G variants in this region and tier." };
  if (use.length > MAX_VARIANTS)
    return { error: `${use.length.toLocaleString()} variants: narrow the region or the tiers (up to 60,000).` };
  const tData = Math.round(performance.now() - t0), W = Math.min(1100, Math.max(600, width - 40));
  const res = await pyramid(use.map(r => ({ pos: r.pos, haps: r.haps, missing: r.missing, male_haploid: r.male_haploid })), {
    samples: smp.ids.length, haploidMales: [...smp.haploidMales], hapIndexes: smp.ids.flatMap((s, i) => (ids.has(s) ? [2 * i,
      2 * i + 1] : [])),
    maxDist: Math.min(a.maxDist, hi - lo), threshold: a.threshold, stat: a.r2_d, W, lo, hi, radius: Math.max(0, Math.min(6,
      Math.floor(W / use.length / 2))),
  }, a.workers > 0 ? { size: a.workers } : {});
  return { ...res, chrom, lo, hi, rows: use, labels: use.map(r => (r.rsid != null ? `rs${r.rsid}` : `${chrom}:${r.pos}`)), genes, recomb,
    gwas: gwas.map(r => ({ ...r })), forge, regulome,
    eqtl: { rows: eqtl, tissues: a.tissues.length ? a.tissues.join(", ") : "none chosen" }, tData };
}
