// LDexpress on the core: LDlink's API text (express_variants_annotated<ref>.txt; server/LDexpress.py,
// LDexpress_ld_sub.py, LDexpress_tissues_sub.py, the LDlink.py route) or { error }.
// For each query: the window rows in LD (LDproxy's chunked tabix window, r² or D′ ≥ threshold), then every GTEx eQTL at
// each row's position in the chosen tissues with pval_nominal < p_threshold. LDlink's mechanics, reproduced:
//  - inputs: "+".join(line.strip().lower()), at most 10 entries (blank lines count); blank entries are then dropped,
//    and a 1–2 character entry ends the request ("Input list of RS numbers is empty")
//  - window rows in chunk order, VCF order within a chunk (boundary-spanning variants twice); RS ID is the VCF rsID or
//    the first dbSNP doc at the VCF position; R2/D′ printed unrounded
//  - eQTLs in index order, joined to gtex_tissues and gtex_genes as inner joins; Effect Size / P-value are GTEx's
//    %.6g strings; allele columns are the window variant's ("A=0.616"; "NA" for a frequency of 1)
//  - a query not in 1000G deletes the request's population file (get_query_variant_c cleanup), so every later query
//    returns nothing; a multi-allelic query is computed on its GT digits ("0|2": only 0/1 haplotypes count)
import { pyRound, pyStr, cFormatG } from "../pyfmt.js";
import { decodeHaps } from "../../core/haps.js";
import { dedupeLists, vcfAlleles, vcfId, addWarning, check, failFirst, screen, WARNED, table2x2, chunkedWindow,
  countsAgainst } from "./common.js";

/**
 * A multi-allelic query record as LDexpress_ld_sub sees it: GT digits per haplotype, counted only when "0" or "1".
 * → a row whose ALT bits are the record's first ALT and whose haplotypes carrying ALT ≥ 2 are uncalled.
 */
function multiQuery(q, all, nSamples) {
  const recs = all.filter(x => x.pos === q.pos && x.ref === q.ref && x.multiallelic_source && x.vcf_order <= q.vcf_order).sort((a,
    b) => a.vcf_order - b.vcf_order);
  const mine = []; for (let k = recs.length - 1, o = q.vcf_order; k >= 0 && recs[k].vcf_order === o; k--, o--) mine.unshift(recs[k]);
    // contiguous split rows
  const hbytes = Math.ceil((2 * nSamples) / 8), missing = new Uint8Array(hbytes);
  if (q.missing) missing.set(q.missing.subarray(0, hbytes));
  for (const r of mine.slice(1)) { const b = decodeHaps(r.haps, hbytes); for (let i = 0; i < hbytes; i++) missing[i] |= b[i]; }
  return { ...mine[0], chrom: q.chrom, haps: decodeHaps(mine[0].haps, hbytes).slice(), missing, multiallelic_source: false };
}

export async function ldexpress(ctx, { snps: text, pop, tissues, r2_d = "r2", r2_d_threshold = 0.1, p_threshold = 0.1, window = 500000,
  build = "grch37", chunks = 4 }) {
  const bad = failFirst(check.r2d(r2_d), check.window(window), check.build(build));
  if (bad) return bad;
  const raw = text.split(/\r\n|\n|\r/).map(s => s.trim().toLowerCase()).join("+").split("+");
  const bad2 = failFirst(check.list(raw.length, 10));
  if (bad2) return bad2;
  const lists = dedupeLists(raw);
  const bad3 = failFirst(check.pop(pop));
  if (bad3) return bad3;
  const o = {};
  const toks = lists.filter(t => t.length).map(t => t[0]); // replace_coords_rsid_list drops empty entries
  // LDexpress stops at the first entry whose label is too short (an empty one included): the whole list fails
  if (toks.some(t => !t)) return { error: "Input list of RS numbers is empty" };
  const { items, kept: queries, warned } = await screen(ctx, build, toks, o);
  if (items.some(x => x.label.length <= 2)) return { error: "Input list of RS numbers is empty" };
  if (warned.length) addWarning(o, WARNED(warned));
  if (!queries.length) return { error: `Input SNP list does not contain any valid RS numbers or coordinates. ${o.warning}` };

  const tissueIds = String(tissues).split("+"), out = [], stat = r2_d === "d" ? "dprime" : "r2";
  let popsGone = false;
  for (const { label, v } of queries) {
    // get_query_variant_c(..., is_output=True) deletes pops_<ref>.txt when a query is not in 1000G: every later
    // query's subprocesses then fail to read it and return nothing
    if (!v.selected?.row) { popsGone = true; continue; }
    if (popsGone) continue;
    let q = { ...(await ctx.row(build, v.chrom, v.selected.row)), chrom: v.chrom };
    const { all, rows } = await chunkedWindow(ctx, build, v.chrom, v.pos, window, chunks);
    if (q.multiallelic_source) q = multiQuery(q, all, (await ctx.samples(build, v.chrom)).ids.length);
    const cnt = await countsAgainst(ctx, build, v.chrom, q, rows, pop), ld = [];
    for (const x of rows) {
      const [A, B, C, D] = cnt.get(x), N = A + B + C + D, t = table2x2(A, B, C, D); if (!t) continue;
      const { r2, dprime } = t;
      if (!((stat === "r2" ? r2 : dprime) >= r2_d_threshold)) continue;
      const [pa0, pa1] = vcfAlleles(x.ref, x.alt), f = c => (N > c ? pyStr(pyRound(c / N, 3)) : "NA");
      ld.push({ x, r2, dprime, al1: `${pa0}=${f(A + C)}`, al2: `${pa1}=${f(B + D)}`, id: vcfId(x) });
    }
    const noRs = [...new Set(ld.filter(e => !e.id.startsWith("rs")).map(e => e.x.pos))], at = await ctx.dbsnpAt(build, v.chrom, noRs);
    const eq = await ctx.eqtl(build, v.chrom, ld.map(e => e.x.pos), tissueIds);
    for (const e of ld) {
      const rs = e.id.startsWith("rs") ? e.id : at.get(e.x.pos) ? `rs${at.get(e.x.pos).rsid}` : ".";
      for (const t of eq.get(e.x.pos) ?? []) if (t.pval_nominal < p_threshold)
        out.push([label, rs, `chr${v.chrom}:${e.x.pos}`, pyStr(e.r2), pyStr(e.dprime), t.gene_name_v26, t.gene_id, t.tissueSiteDetail,
          e.al1, e.al2, cFormatG(t.slope), cFormatG(t.pval_nominal)].join("\t"));
    }
  }
  if (!out.length) return { error: "No entries in GTEx are identified using the LDexpress search criteria." };
  const head = "Query\tRS ID\tPosition\tR2\tD'\tGene Symbol\tGencode ID\tTissue\tNon-effect Allele Freq\tEffect Allele Freq\tEffect Size\tP-value";
  return { text: [head, ...out].join("\n") + "\n" + (o.warning !== undefined ? `Warning(s):\n${o.warning}` : "") };
}
