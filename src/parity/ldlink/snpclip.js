// SNPclip on the core: LDlink's API text (details<ref>.txt; server/SNPclip.py calculate_clip + the LDlink.py route)
// or { error }. LD pruning in input order: each kept variant removes every later one with r² ≥ threshold.
// LDlink's mechanics, reproduced:
//  - the route keeps only lines starting rs/chr (after lstrip), lower-cased; validsnp then dedupes token lists
//  - MAF from genotype strings (haploid counts one allele); kept when maf ≥ maf_threshold (unrounded)
//  - r² pairs genotype strings per individual: "0|1-1|0" etc., and haploid with haploid; an individual diploid at one
//    variant and haploid at the other (chrX PAR vs non-PAR males) is skipped for that pair — a second prepared set
//    without the haploid males serves those pairs
//  - a multi-allelic record is listed with its raw REF and ALT list ("A=NA, G,T=NA")
import { pyRound, pyStr } from "../pyfmt.js";
import { decodeHaps } from "../../core/haps.js";
import { TITLE } from "../../genome.js";
import { dedupeLists, vcfAlleles, addWarning, check, failFirst, routeLines, screen, WARNED, sameChromError } from "./common.js";
import { stats } from "../../core/ld/kernel_js.js";

const INTERNAL = { error: "Internal server error. Please contact LDlink admin." };

export async function snpclip(ctx, { snps: text, pop, r2_threshold = 0.1, maf_threshold = 0.01, build = "grch37" }) {
  const bad = failFirst(check.build(build));
  if (bad) return bad;
  // route: lstrip, keep rs/chr lines, lower-case; validsnp: limit on lines, token lists, dedupe
  const lines = routeLines(text);
  const bad2 = failFirst(check.list(lines.length, 5000));
  if (bad2) return bad2;
  const lists = dedupeLists(lines);
  const bad3 = failFirst(check.pop(pop));
  if (bad3) return bad3;
  const o = {}, details = new Map(), set = (k, v) => details.set(k, v); // OrderedDict: re-setting keeps the first position
  const toks = lists.filter(t => t.length).map(t => t[0]);
  const { items, kept, warned } = await screen(ctx, build, toks, o), labels = items.map(x => x.label);
  const REMOVED = { label: "Not a RS number, query removed.", dbsnp: `Variant not found in dbSNP155 (${TITLE[build]}), variant removed.`,
    chrY: "Chromosome Y variants are unavailable for GRCh38, only available for GRCh37." };
  for (const { label, bad: why } of items) if (why) set(label, ["NA", "NA", REMOVED[why]]);
  if (warned.length) addWarning(o, WARNED(warned));
  const W = () => o.warning ?? "";
  if (!kept.length) return { error: `Input SNP list does not contain any valid RS numbers or coordinates. ${W()}` };
  // check_same_chromosome sets the error; parse_vcf (tabix on the first SNP's chromosome) may replace it
  const diff = kept.findIndex(k => k.v.chrom !== kept[0].v.chrom), chrom = kept[0].v.chrom;
  const found = kept.filter(k => k.v.chrom === chrom && k.v.selected?.row);
  if (found.length === 0) return { error: `Input variant list does not contain any valid RS numbers or coordinates. ${W()}` };
  if (diff > 0) return { error: `${sameChromError(kept, diff)} ${W()}` };
  const entries = [...new Map(found.map(k => [`${k.v.pos}_${k.label}`, k])).values()]; // parse_vcf, unsorted

  const { ids: samples, haploidMales } = await ctx.samples(build, chrom), hbytes = Math.ceil((2 * samples.length) / 8);
  const panel = await ctx.popIds(pop), idx = samples.flatMap((s, i) => (panel.has(s) ? [i] : []));
  const rows = [], order = [];
  for (const { label, v } of entries) {
    const row = { ...(await ctx.row(build, chrom, v.selected.row)), chrom };
    if (v.selected.status === "not_biallelic") {
      const alts = [row.alt];
      for (let k = row.vcf_order - 1; ; k--) { // the record's other ALTs: the split rows just before, same pos and REF
        const w = (await ctx.windowRows(build, chrom, row.pos, row.pos)).find(x => x.vcf_order === k && x.pos === row.pos && x.ref === row.ref && x.multiallelic_source);
        if (!w) break; alts.unshift(w.alt);
      }
      set(label, [`chr${chrom}:${v.pos}`, `${row.ref}=NA, ${alts.join(",")}=NA`, "Variant is not biallelic, variant removed."]);
      continue;
    }
    const bits = decodeHaps(row.haps, hbytes), bit = k => (bits[k >> 3] >> (k & 7)) & 1;
    let z = 0, n1 = 0;
    for (const s of idx) {
      if (row.missing && (row.missing[(2 * s) >> 3] >> ((2 * s) & 7)) & 1) continue;
      const ks = row.male_haploid && haploidMales.has(s) ? [2 * s] : [2 * s, 2 * s + 1];
      for (const k of ks) bit(k) ? n1++ : z++;
    }
    if (z + n1 === 0) return INTERNAL; // ZeroDivisionError in calc_maf
    const f0 = z / (z + n1), f1 = n1 / (z + n1), maf = Math.min(f0, f1), [a0, a1] = vcfAlleles(row.ref, row.alt);
    const d = [`chr${chrom}:${v.pos}`, `${a0}=${pyStr(pyRound(f0, 3))}, ${a1}=${pyStr(pyRound(f1, 3))}`];
    set(label, d);
    if (maf_threshold <= maf) { rows.push(row); order.push(label); } else d.push(`Variant MAF is ${pyStr(pyRound(maf, 4))}, variant removed.`);
  }
  for (const k of kept) if (!order.includes(k.label) && !details.has(k.label)) set(k.label, [`chr${chrom}:${k.v.pos}-${k.v.pos}`, "NA", "Variant not in 1000G data, variant removed."]);
  // greedy pruning; pairs whose rows differ in haploid-ness use the set without haploid males
  const alive = rows.map(() => true), s2 = new Float64Array(2);
  const failed = await ctx.withPrepared(build, rows, pop, p => { // synchronous: both prepared sets live until the reset
    const noMales = [...haploidMales].some(i => panel.has(samples[i])) && rows.some(r => r.male_haploid) && rows.some(r => !r.male_haploid)
      ? ctx.kernel.prepare(rows, { samples: samples.length, haploidMales, hapIndexes: idx.filter(i => !haploidMales.has(i)).flatMap(i => [2 * i, 2 * i + 1]) }) : null;
    for (let i = 0; i < rows.length; i++) {
      if (!alive[i]) continue;
      details.get(order[i]).push("Variant kept.");
      const J = []; for (let j = i + 1; j < rows.length; j++) if (alive[j]) J.push(j);
      const mixed = j => noMales && rows[i].male_haploid !== rows[j].male_haploid;
      const cA = ctx.kernel.counts(p, J.map(() => i), J), cB = noMales ? ctx.kernel.counts(noMales, J.map(() => i), J) : null;
      for (const [k, j] of J.entries()) {
        if (!stats(...(mixed(j) ? cB[k] : cA[k]), s2)) return true; // r2 None: TypeError in LDlink
        if (r2_threshold <= s2[0]) { alive[j] = false; details.get(order[j]).push(`Variant in LD with ${order[i]} (R2=${pyStr(pyRound(s2[0], 4))}), variant removed.`); }
      }
    }
    return false;
  });
  if (failed) return INTERNAL;
  const out = ["RS Number\tPosition\tAlleles\tDetails"];
  for (const t of labels) { const d = details.get(t); if (!d) return INTERNAL; out.push(`${t}\t${d.join("\t")}`); }
  return { text: out.join("\n") + "\n" };
}
