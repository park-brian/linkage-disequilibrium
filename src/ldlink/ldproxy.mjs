// LDproxy on the core: LDlink's API text (server/LDproxy.py + LDproxy_sub.py; proxy<ref>.txt) and its warning, or
// { error }. LDlink's own mechanics, reproduced:
//  - the window is cut into NUM_SUBPROCESSES chunks (LDcommon.chunkWindow), each read with tabix (overlap
//    semantics), so a variant spanning a chunk boundary is listed twice; the last chunk ends a few bp past pos + window
//  - rows: every biallelic 1000G row (VCF order per chunk) with Ms ≠ 0; sort |distance| ascending, then the statistic
//    descending (stable); the first row is printed as the query; proxies need R² > 0.01 and RS ≠ the query label
//  - RS_Number: the VCF ID when it starts with "rs", else the first dbSNP doc at the VCF position; Function from
//    db.dbsnp.find_one by that ID / position; FORGEdb by RS (empty → "NA"); RegulomeDB at the VCF position
//  - MAF is the proxy's minor allele frequency over haplotypes called in both variants
import { pyRoundRepr as r } from "../pyfmt.mjs";
import { decodeHaps } from "../haps.mjs";
import { TITLE } from "../genome.mjs";
import { POPS, POP_ERROR, chrYMessage, vcfAlleles, vcfId, addCoordWarning, rsAnnotations, table2x2, chunkedWindow, countsAgainst } from "./common.mjs";

export async function ldproxy(ctx, { var: input = "", pop, r2_d = "r2", window = 500000, build = "grch37", chunks = 4 }) {
  const snp0 = input.toLowerCase().trim(), o = {};
  // the REST layer's checks (LDlink.py _validate_common_parameters) come first
  if (!["r2", "d"].includes(r2_d)) return { error: "Invalid r2_d parameter." };
  if (!Number.isInteger(+window) || window < 0 || window > 1000000) return { error: "Invalid window parameter." };
  const v = /\s/.test(snp0) || !snp0 ? { found: false } : (await ctx.planner.query({ build, inputs: [snp0], select: "parity" })).variants[0];
  addCoordWarning(o, v.warning);
  let snp = snp0.startsWith("rs") ? snp0 : (v.rsid ?? snp0);
  if (!v.found || v.pos == null) return { error: `${snp} is not in dbSNP 155 (${TITLE[build]}).` };
  if (v.chrom === "Y" && build !== "grch37") return { error: chrYMessage(v) };
  for (const p of String(pop).split("+")) if (!POPS.includes(p)) return { error: POP_ERROR(p) };
  // get_query_variant_c, then LDproxy's pass over its messages
  const sel = v.selected;
  if (!sel?.row) return { error: `${snp} Variant is not in 1000G reference panel.` };
  const q = { ...(await ctx.row(build, v.chrom, sel.row)), chrom: v.chrom }, qid = vcfId(q);
  const { ids: samples, haploidMales } = await ctx.samples(build, v.chrom), hbytes = Math.ceil((2 * samples.length) / 8), panel = await ctx.popIds(pop);
  const mismatch = qid !== snp && qid.includes("rs") ? `Genomic position does not match RS number at 1000G position (chr${v.chrom}:${q.pos} = ${qid}).` : null;
  if (mismatch) o.warning = mismatch + (o.warning ?? "");
  if (sel.status === "not_biallelic") return { error: `${snp} variant is not a biallelic.${o.warning ?? ""}Variant is not a biallelic.` };
  { // monoallelic: the '0' and '1' characters of the population's GT strings
    const bits = decodeHaps(q.haps, hbytes); let n0 = 0, n1 = 0;
    samples.forEach((s, i) => {
      if (!panel.has(s) || (q.missing && (q.missing[(2 * i) >> 3] >> ((2 * i) & 7)) & 1)) return;
      for (const k of q.male_haploid && haploidMales.has(i) ? [2 * i] : [2 * i, 2 * i + 1]) (bits[k >> 3] >> (k & 7)) & 1 ? n1++ : n0++;
    });
    if (!n0 || !n1) return { error: "Variant is monoallelic in the chosen population(s)." };
  }
  if (mismatch) { o.warning += mismatch; snp = qid; } // the loop appends the warning again and the query label becomes the VCF ID
  // the window, chunk by chunk (tabix per chunk), biallelic rows only
  const { rows } = await chunkedWindow(ctx, build, v.chrom, v.pos, window, chunks), cnt = await countsAgainst(ctx, build, v.chrom, q, rows, pop);
  const [qa0, qa1] = vcfAlleles(q.ref, q.alt), out = [];
  for (const x of rows) {
    const [A, B, C, D] = cnt.get(x), N = A + B + C + D, t = table2x2(A, B, C, D); if (!t) continue;
    const { dprime, r2 } = t, [pa0, pa1] = vcfAlleles(x.ref, x.alt);
    const match = r2 > 0.1 ? ((A * D) / Math.max(B * C, 0.01) > 1 ? `${qa0}=${pa0},${qa1}=${pa1}` : `${qa0}=${pa1},${qa1}=${pa0}`) : "NA";
    out.push({ x, dist: x.pos - q.pos, dprime, r2, match, maf: Math.min((A + C) / N, (B + D) / N), id: vcfId(x) });
  }
  const byStat = r2_d === "d" ? "dprime" : "r2";
  out.sort((a, b) => Math.abs(a.dist) - Math.abs(b.dist)).sort((a, b) => b[byStat] - a[byStat]);
  // annotations: RS by VCF ID or by position (the first dbSNP doc there, which also gives Function), FORGEdb, RegulomeDB
  const noRs = [...new Set(out.filter(e => !e.id.startsWith("rs")).map(e => e.x.pos))], at = await ctx.dbsnpAt(build, v.chrom, noRs);
  for (const e of out) if (!e.id.startsWith("rs")) { const d = at.get(e.x.pos); e.rs = d ? `rs${d.rsid}` : "."; e.func = d ? d.function : "."; } else e.rs = e.id;
  // where each rsID's doc is: a VCF-labelled row's label at POS + rsid_pos_delta; a position-lookup rsID at the VCF POS
  const docPos = e => (!e.id.startsWith("rs") ? e.x.pos : e.x.rsid_source === 1 && e.id === `rs${e.x.rsid}` && e.x.rsid_pos_delta != null ? e.x.pos + e.x.rsid_pos_delta : null);
  // unknown doc position: try the VCF POS (a single-doc rsID found there is exact; a miss falls back to find_one)
  const items = out.filter(e => /^rs\d+$/.test(e.rs)).map(e => ({ id: String(BigInt(e.rs.slice(2))), pos: docPos(e) ?? e.x.pos }));
  const ann = await rsAnnotations(ctx, build, v.chrom, items), rg = await ctx.regulome(build, v.chrom, [...new Set(out.map(e => e.x.pos))]);
  for (const e of out) if (e.id.startsWith("rs")) { const a = /^rs\d+$/.test(e.id) ? ann.get(String(BigInt(e.id.slice(2)))) : null; e.func = a && a.func !== undefined ? a.func : "."; }
  const forgeOf = e => (/^rs\d+$/.test(e.rs) ? ann.get(String(BigInt(e.rs.slice(2))))?.forge : undefined);
  // the first line takes RS, distance and statistics from the top row but Coord and Alleles from the query (query_snp)
  const line = (e, top) => [e.rs, `chr${v.chrom}:${(top ? q : e.x).pos}`, `(${vcfAlleles((top ? q : e.x).ref, (top ? q : e.x).alt).join("/")})`, r(e.maf, 4), String(e.dist), r(e.dprime, 4), r(e.r2, 4), e.match,
    forgeOf(e)?.toString() ?? "NA", rg.get(e.x.pos) ?? ".", e.func ?? "NA"].join("\t");
  const text = ["RS_Number\tCoord\tAlleles\tMAF\tDistance\tDprime\tR2\tCorrelated_Alleles\tFORGEdb\tRegulomeDB\tFunction", line(out[0], true),
    ...out.slice(1).filter(e => e.r2 > 0.01 && e.rs !== snp).map(e => line(e, false))].join("\n") + "\n";
  return { text, ...(o.warning !== undefined ? { warning: o.warning } : {}) };
}
