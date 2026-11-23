// SNPclip on the core: LDlink's API text (details<ref>.txt; server/SNPclip.py calculate_clip + the LDlink.py route)
// or { error }. LD pruning in input order: each kept variant removes every later one with r² ≥ threshold.
// LDlink's mechanics, reproduced:
//  - the route keeps only lines starting rs/chr (after lstrip), lower-cased; validsnp then dedupes token lists
//  - MAF from genotype strings (haploid counts one allele); kept when maf ≥ maf_threshold (unrounded)
//  - r² pairs genotype strings per individual: "0|1-1|0" etc., and haploid with haploid; an individual diploid at one
//    variant and haploid at the other (chrX PAR vs non-PAR males) is skipped for that pair — a second prepared set
//    without the haploid males serves those pairs
//  - a multi-allelic record is listed with its raw REF and ALT list ("A=NA, G,T=NA")
import { pyRound, pyStr } from "../pyfmt.mjs";
import { decodeHaps } from "../haps.mjs";
import { TITLE } from "../genome.mjs";
import { POPS, POP_ERROR } from "./common.mjs";
import { vcfAlleles, addCoordWarning } from "./common.mjs";
import { stats } from "../ld/kernel_js.mjs";

const gate = t => t.length > 2 && (t.startsWith("rs") || t.startsWith("chr")) && /\d$/.test(t);
const INTERNAL = { error: "Internal server error. Please contact LDlink admin." };

export async function snpclip(ctx, { snps: text, pop, r2_threshold = 0.1, maf_threshold = 0.01, build = "grch37" }) {
  if (!TITLE[build]) return { error: `Invalid genome build. Please specify either ${Object.keys(TITLE).join(", ")}.` };
  // route: lstrip, keep rs/chr lines, lower-case; validsnp: limit on lines, token lists, dedupe
  const lines = text.split(/\r\n|\n|\r/).map(s => s.replace(/^\s+/, "")).filter(s => /^(rs|chr)/i.test(s)).map(s => s.toLowerCase());
  if (lines.length > 5000) return { error: `Maximum variant list is 5000  RS numbers or coordinates. Your list contains ${lines.length} entries.` };
  const lists = [];
  for (const l of lines) { const t = l.trim().split(/\s+/).filter(Boolean); if (!lists.some(x => x.join(" ") === t.join(" "))) lists.push(t); }
  for (const p of String(pop).split("+")) if (!POPS.includes(p)) return { error: POP_ERROR(p) };
  const o = {}, details = new Map(), set = (k, v) => details.set(k, v); // OrderedDict: re-setting keeps the first position
  const toks = lists.filter(t => t.length).map(t => t[0]);
  const vs = (await ctx.planner.query({ build, inputs: toks, select: "parity" })).variants;
  vs.forEach(v => addCoordWarning(o, v.warning));
  const labels = toks.map((t, i) => (t.startsWith("rs") ? t : (vs[i].rsid ?? t))); // replace_coords_rsid_list
  const kept = [], warn = [], join = m => { o.warning = o.warning === undefined ? m : `${o.warning}. ${m}`; };
  labels.forEach((label, i) => {
    const v = vs[i];
    if (!gate(label)) { warn.push(label); return set(label, ["NA", "NA", "Not a RS number, query removed."]); }
    if (!v.found || v.pos == null) { warn.push(label); return set(label, ["NA", "NA", `Variant not found in dbSNP155 (${TITLE[build]}), variant removed.`]); }
    if (v.chrom === "Y" && build !== "grch37") {
      join(`Input variants on chromosome Y are unavailable for GRCh38, only available for GRCh37 (${v.rsid} = chr${v.chrom}:${v.pos})`); warn.push(label);
      return set(label, ["NA", "NA", "Chromosome Y variants are unavailable for GRCh38, only available for GRCh37."]);
    }
    kept.push({ label, v });
  });
  if (warn.length) join(`The following RS number(s) or coordinate(s) inputs have warnings: ${warn.join(", ")}`);
  const W = () => o.warning ?? "";
  if (!kept.length) return { error: `Input SNP list does not contain any valid RS numbers or coordinates. ${W()}` };
  // check_same_chromosome sets the error; parse_vcf (tabix on the first SNP's chromosome) may replace it
  const diff = kept.findIndex(k => k.v.chrom !== kept[0].v.chrom), chrom = kept[0].v.chrom;
  const found = kept.filter(k => k.v.chrom === chrom && k.v.selected?.row);
  if (found.length === 0) return { error: `Input variant list does not contain any valid RS numbers or coordinates. ${W()}` };
  if (diff > 0) return { error: `Not all input variants are on the same chromosome: ${kept[diff - 1].label}=chr${kept[diff - 1].v.chrom}:${kept[diff - 1].v.pos}, ${kept[diff].label}=chr${kept[diff].v.chrom}:${kept[diff].v.pos}. ${W()}` };
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
  const p = await ctx.prepare(build, rows, pop);
  const noMales = [...haploidMales].some(i => panel.has(samples[i])) && rows.some(r => r.male_haploid) && rows.some(r => !r.male_haploid)
    ? ctx.kernel.prepare(rows, { samples: samples.length, haploidMales, hapIndexes: idx.filter(i => !haploidMales.has(i)).flatMap(i => [2 * i, 2 * i + 1]) }) : null;
  const alive = rows.map(() => true), s2 = new Float64Array(2);
  for (let i = 0; i < rows.length; i++) {
    if (!alive[i]) continue;
    details.get(order[i]).push("Variant kept.");
    const J = []; for (let j = i + 1; j < rows.length; j++) if (alive[j]) J.push(j);
    const mixed = j => noMales && rows[i].male_haploid !== rows[j].male_haploid;
    const cA = ctx.kernel.counts(p, J.map(() => i), J), cB = noMales ? ctx.kernel.counts(noMales, J.map(() => i), J) : null;
    for (const [k, j] of J.entries()) {
      if (!stats(...(mixed(j) ? cB[k] : cA[k]), s2)) { ctx.kernel.reset(); return INTERNAL; } // r2 None: TypeError in LDlink
      if (r2_threshold <= s2[0]) { alive[j] = false; details.get(order[j]).push(`Variant in LD with ${order[i]} (R2=${pyStr(pyRound(s2[0], 4))}), variant removed.`); }
    }
  }
  ctx.kernel.reset();
  const out = ["RS Number\tPosition\tAlleles\tDetails"];
  for (const t of labels) { const d = details.get(t); if (!d) return INTERNAL; out.push(`${t}\t${d.join("\t")}`); }
  return { text: out.join("\n") + "\n" };
}
