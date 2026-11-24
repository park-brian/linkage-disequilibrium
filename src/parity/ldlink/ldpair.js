// LDpair on the core, with LDlink's exact JSON (server/LDpair.py calculate_pair, one GET pair; "request" omitted).
// Returns { error } for a request-level or single-pair error (as LDlink's API does), else [pairResult].
import { pyRoundRepr as r } from "../pyfmt.js";
import { TITLE } from "../../core/genome.js";
import { chrYMessage, vcfAlleles, ldStats, queryVariant, check, failFirst, resolveOne, labelOf } from "./common.js";

const VALID = /^(rs\d+|chr\d+:\d+|chr[X|Y]:\d+)/i; // LDpair.py: re.match (a prefix match)

/**
 * Haplotype table in LDlink's layout from the 2×2 margins [n, alt1, alt2, both] (kernel.counts). A..D follow the
 * sorted "allele1_allele2" keys (Python sorted() on ASCII == JS sort). LDlink orders hap1..hap4 by count, descending;
 * ties keep dict insertion order: first appearance while iterating a Python set of sample IDs (hash-randomized per
 * process), then the zero-filled keys. Ours: observed ties in key order, zero-filled keys last in LDlink's fill order.
 */
function pairTable(label1, label2, row1, row2, [n, p1, p2, x]) {
  const [a1, b1] = vcfAlleles(row1.ref, row1.alt), [a2, b2] = vcfAlleles(row2.ref, row2.alt);
  const hap = { [`${b1}_${b2}`]: x, [`${b1}_${a2}`]: p1 - x, [`${a1}_${b2}`]: p2 - x, [`${a1}_${a2}`]: n - p1 - p2 + x };
  const fill = [`${a1}_${a2}`, `${a1}_${b2}`, `${b1}_${a2}`, `${b1}_${b2}`]; // LDpair.py zero-fill order
  const s = Object.keys(hap).sort(), [A, B, C, D] = s.map(k => hap[k]), N = A + B + C + D, f = v => r(v / N, 3);
  const order = [...s.filter(k => hap[k] > 0), ...fill.filter(k => hap[k] === 0)].sort((u, v) => hap[v] - hap[u]); // stable
  const st = ldStats(A, B, C, D), half = k => s[k].split("_");
  const corr = st && st.r2 > 0.1
    ? ((A * D) / Math.max(B * C, 0.01) > 1 ? [0, 3] : [1, 2]).map(k => `${label1}(${half(k)[0]}) allele is correlated with ${label2}(${half(k)[1]}) allele`)
    : [`${label1} and ${label2} are in linkage equilibrium`];
  return {
    corr_alleles: corr,
    haplotypes: Object.fromEntries(order.map((k, i) => [`hap${i + 1}`, { alleles: k, count: String(hap[k]), frequency: f(hap[k]) }])),
    snp1: { allele_1: { allele: half(0)[0], count: String(A + B), frequency: f(A + B) }, allele_2: { allele: half(2)[0], count: String(C + D), frequency: f(C + D) }, rsnum: label1 },
    snp2: { allele_1: { allele: half(0)[1], count: String(A + C), frequency: f(A + C) }, allele_2: { allele: half(1)[1], count: String(B + D), frequency: f(B + D) }, rsnum: label2 },
    statistics: st ? { chisq: r(st.chisq, 4), d_prime: r(st.dprime, 4), p: st.p >= 0.0001 ? r(st.p, 4) : "<0.0001", r2: r(st.r2, 4) } : { chisq: "NA", d_prime: "NA", p: "NA", r2: "NA" },
    two_by_two: { cells: { c11: String(A), c12: String(B), c21: String(C), c22: String(D) }, total: String(N) },
  };
}

export async function ldpair(ctx, { var1 = "", var2 = "", pop, build = "grch37" }) {
  const bad = failFirst(check.build(build), check.pop(pop));
  if (bad) return bad;
  const o = {}, W = () => o.warning ?? "", fail = error => ({ error }); // the API returns only the error of a failed pair
  if (var1.length < 3 || var2.length < 3) return fail(`Missing or additional SNPs in pair. ${W()}`);
  const snps = [];
  for (const raw of [var1, var2]) { // each SNP is validated and resolved before the next is looked at
    const s = raw.toLowerCase().trim();
    if (!VALID.test(s)) return fail(`${s} is not a valid SNP. ${W()}`);
    const v = await resolveOne(ctx, build, s, o), label = labelOf(s, v);
    if (!v.found || v.pos == null) return fail(`${label} is not in dbSNP build 155 (${TITLE[build]}). ${W()}`);
    snps.push({ label, v });
  }
  const [s1, s2] = snps;
  if (s1.v.chrom !== s2.v.chrom) o.warning = `${W()}${s1.label} and ${s2.label} are on different chromosomes. `;
  for (const { v } of snps) if (v.chrom === "Y" && build !== "grch37")
    return fail(`${chrYMessage(v, "-")}. ${W()}`);
  const rows = [];
  for (const { label, v } of snps) rows.push(await queryVariant(ctx, build, label, v, o));
  if (o.error !== undefined) return fail(o.error);
  snps.forEach(({ v }, i) => { if (rows[i].pos !== v.pos) o.warning = `${W()}VCF File does not match variant coordinates for SNP${i + 1}. `; });
  const c = await ctx.withPrepared(build, rows, pop, p => ctx.kernel.counts(p, [0], [1])[0]);
  const t = pairTable(s1.label, s2.label, rows[0], rows[1], c);
  t.snp1.coord = `chr${s1.v.chrom}:${s1.v.pos}`; t.snp2.coord = `chr${s2.v.chrom}:${s2.v.pos}`;
  return [{ ...t, pair: [var1, var2], ...(o.warning !== undefined ? { warning: o.warning } : {}) }];
}
