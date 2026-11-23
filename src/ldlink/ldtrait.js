// LDtrait on the core: LDlink's API text (trait_variants_annotated<ref>.txt; server/LDtrait.py, LDtrait_ld_sub.py,
// the LDlink.py route) or { error }. For each query, the GWAS Catalog records within ±window, kept when the pair's
// r² (or D′) ≥ threshold. LDlink's mechanics, reproduced:
//  - records in index order: position, then insertion order (= our file's (key, dump order))
//  - pair LD from one tabix row per position (the last overlapping row); NA (never matched) when that row is
//    multi-allelic or its POS differs from the queried position (indels listed at POS+1), or Ms = 0
//  - R2 and D′ printed unrounded (str(float)); Alleles are the GWAS variant's, "A=0.692, G=0.308" ("NA" for a
//    frequency of 1); columns: RAF ("Risk Allele"), float(OR or BETA) ("Effect Size"), the CI with brackets, letters and
//    spaces removed ("Beta or OR"), P-VALUE verbatim
import { pyRound, pyStr } from "../pyfmt.js";
import { TITLE } from "../genome.js";
import { POPS, POP_ERROR, chrYMessage, dedupeLists, vcfAlleles, addCoordWarning, gate, table2x2 } from "./common.js";

const PYFLOAT = /^\s*[+-]?((\d(_?\d)*(\.(\d(_?\d)*)?)?|\.\d(_?\d)*)([eE][+-]?\d(_?\d)*)?|inf(inity)?|nan)\s*$/i;
/** castFloat: str(float(v)) when Python's float() accepts v, else v */
const castFloat = v => (PYFLOAT.test(v) ? pyStr(/inf/i.test(v) ? (v.includes("-") ? -Infinity : Infinity) : /nan/i.test(v) ? NaN : Number(v.trim().replace(/_/g, ""))) : v);
const rangeString = v => v.replace(/\[*\]*[a-zA-Z]*\s*/g, "") || "NA"; // findRangeString
const nz = v => (v != null && v.length > 0 ? v : "NA");

export async function ldtrait(ctx, { snps: text, pop, r2_d = "r2", r2_d_threshold = 0.1, window = 500000, build = "grch37" }) {
  if (!["r2", "d"].includes(r2_d)) return { error: "Invalid r2_d parameter." };
  if (!Number.isInteger(+window) || window < 0 || window > 1000000) return { error: "Invalid window parameter." };
  if (!TITLE[build]) return { error: `Invalid genome build. Please specify either ${Object.keys(TITLE).join(", ")}.` };
  const lines = text.split(/\r\n|\n|\r/).map(s => s.replace(/^\s+/, "")).filter(s => /^(rs|chr)/i.test(s)).map(s => s.toLowerCase());
  if (lines.length > 50) return { error: `Maximum variant list is 50  RS numbers or coordinates. Your list contains ${lines.length} entries.` };
  const lists = dedupeLists(lines);
  for (const p of String(pop).split("+")) if (!POPS.includes(p)) return { error: POP_ERROR(p) };
  const o = {}, join = m => { o.warning = o.warning === undefined ? m : `${o.warning}. ${m}`; };
  const toks = lists.map(t => t[0]);
  const vs = (await ctx.planner.query({ build, inputs: toks, select: "parity" })).variants;
  vs.forEach(v => addCoordWarning(o, v.warning));
  const queries = [], warn = [];
  toks.forEach((t, i) => {
    const v = vs[i], label = t.startsWith("rs") ? t : (v.rsid ?? t);
    if (!gate(label) || !v.found || v.pos == null) return warn.push(label);
    if (v.chrom === "Y" && build !== "grch37") { join(chrYMessage(v)); return warn.push(label); }
    queries.push({ label, v });
  });
  if (warn.length) join(`The following RS number(s) or coordinate(s) inputs have warnings: ${warn.join(", ")}`);
  if (!queries.length) return { error: `Input SNP list does not contain any valid RS numbers or coordinates. ${o.warning ?? ""}` };

  const out = [];
  const usable = (row, pos) => row && !row.multiallelic_source && row.pos === pos;
  for (const { label, v } of queries) {
    const lo = Math.max(0, v.pos - window), hi = v.pos + window, recs = await ctx.gwas(build, v.chrom, lo, hi);
    // one read of the window; rowAt(p) = the last 1000G row overlapping p (tabix p-p), as each LDtrait_ld_sub lookup
    const win = await ctx.windowRows(build, v.chrom, Math.min(lo, v.pos), Math.max(hi, v.pos));
    const rowAt = p => { let r = null; for (const x of win) if (x.pos <= p && x.pos + x.ref.length - 1 >= p) r = x; return r && { ...r, chrom: v.chrom }; };
    const q = rowAt(v.pos), ld = new Map(), pairs = []; // ldInfo[query]["rs" + SNP_ID_CURRENT]
    for (const rec of recs) {
      const key = `rs${rec.snp_id_current ?? ""}`;
      if (ld.has(key)) continue;
      const g = rowAt(rec.pos);
      ld.set(key, null);
      if (usable(q, v.pos) && usable(g, rec.pos)) pairs.push({ key, g });
    }
    if (pairs.length) { // every pair against the query in one kernel call
      const p = await ctx.prepare(build, [q, ...pairs.map(x => x.g)], pop), counts = ctx.kernel.counts(p, pairs.map(() => 0), pairs.map((_, i) => i + 1)); ctx.kernel.reset();
      const [a1, b1] = vcfAlleles(q.ref, q.alt);
      pairs.forEach(({ key, g }, i) => {
        const [n, p1, p2, x] = counts[i], [a2, b2] = vcfAlleles(g.ref, g.alt);
        const hap = { [`${b1}_${b2}`]: x, [`${b1}_${a2}`]: p1 - x, [`${a1}_${b2}`]: p2 - x, [`${a1}_${a2}`]: n - p1 - p2 + x };
        const s = Object.keys(hap).sort(), [A, B, C, D] = s.map(k => hap[k]), N = A + B + C + D;
        const t = table2x2(A, B, C, D);
        if (!t) return;
        const f = c => (N > c ? pyStr(pyRound(c / N, 3)) : "NA");
        ld.set(key, { r2: t.r2, dprime: t.dprime, alleles: `${s[0].split("_")[1]}=${f(A + C)}, ${s[1].split("_")[1]}=${f(B + D)}` });
      });
    }
    const rows = [];
    for (const rec of recs) {
      const l = ld.get(`rs${rec.snp_id_current ?? ""}`);
      if (!l || !((r2_d === "r2" ? l.r2 : l.dprime) >= r2_d_threshold)) continue;
      rows.push([label, rec.trait ?? "", rec.pubmedid ?? "", `rs${rec.snp_id_current ?? ""}`, `chr${rec.chromosome}:${rec.pos}`, l.alleles, pyStr(l.r2), pyStr(l.dprime),
        nz(rec.risk_allele_frequency), rec.or_beta != null && rec.or_beta.length ? castFloat(rec.or_beta) : "NA",
        rec.ci_95_text != null && rec.ci_95_text.length ? rangeString(rec.ci_95_text) : "NA", nz(rec.p_value)].join("\t"));
    }
    if (rows.length) out.push(...rows);
  }
  if (!out.length) return { error: "No entries in the GWAS Catalog are identified using the LDtrait search criteria." };
  const head = `Query\tGWAS Trait\tPMID\tRS Number\tPosition (${TITLE[build]})\tAlleles\tR2\tD'\tRisk Allele\tEffect Size (95% CI)\tBeta or OR\tP-value`;
  return { text: [head, ...out].join("\n") + "\n" + (o.warning !== undefined ? `Warning(s):\n${o.warning}` : "") };
}
