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
import { pyRound, pyStr, cFormatG } from "../pyfmt.mjs";
import { decodeHaps } from "../haps.mjs";
import { TITLE } from "../genome.mjs";
import { POPS, POP_ERROR } from "./common.mjs";
import { vcfAlleles, vcfId, addCoordWarning } from "./common.mjs";
import { chunkWindow } from "./ldproxy.mjs";

const gate = t => t.length > 2 && (t.startsWith("rs") || t.startsWith("chr")) && /\d$/.test(t);

/**
 * A multi-allelic query record as LDexpress_ld_sub sees it: GT digits per haplotype, counted only when "0" or "1".
 * → a row whose ALT bits are the record's first ALT and whose haplotypes carrying ALT ≥ 2 are uncalled.
 */
function multiQuery(q, all, nSamples) {
  const recs = all.filter(x => x.pos === q.pos && x.ref === q.ref && x.multiallelic_source && x.vcf_order <= q.vcf_order).sort((a, b) => a.vcf_order - b.vcf_order);
  const mine = []; for (let k = recs.length - 1, o = q.vcf_order; k >= 0 && recs[k].vcf_order === o; k--, o--) mine.unshift(recs[k]); // contiguous split rows
  const hbytes = Math.ceil((2 * nSamples) / 8), missing = new Uint8Array(hbytes);
  if (q.missing) missing.set(q.missing.subarray(0, hbytes));
  for (const r of mine.slice(1)) { const b = decodeHaps(r.haps, hbytes); for (let i = 0; i < hbytes; i++) missing[i] |= b[i]; }
  return { ...mine[0], chrom: q.chrom, haps: decodeHaps(mine[0].haps, hbytes).slice(), missing, multiallelic_source: false };
}

export async function ldexpress(ctx, { snps: text, pop, tissues, r2_d = "r2", r2_d_threshold = 0.1, p_threshold = 0.1, window = 500000, build = "grch37", chunks = 4 }) {
  if (!["r2", "d"].includes(r2_d)) return { error: "Invalid r2_d parameter." };
  if (!Number.isInteger(+window) || window < 0 || window > 1000000) return { error: "Invalid window parameter." };
  if (!TITLE[build]) return { error: `Invalid genome build. Please specify either ${Object.keys(TITLE).join(", ")}.` };
  const raw = text.split(/\r\n|\n|\r/).map(s => s.trim().toLowerCase()).join("+").split("+");
  if (raw.length > 10) return { error: `Maximum variant list is 10  RS numbers or coordinates. Your list contains ${raw.length} entries.` };
  const lists = [];
  for (const l of raw) { const t = l.trim().split(/\s+/).filter(Boolean); if (!lists.some(x => x.join(" ") === t.join(" "))) lists.push(t); }
  for (const p of String(pop).split("+")) if (!POPS.includes(p)) return { error: POP_ERROR(p) };
  const o = {}, join = m => { o.warning = o.warning === undefined ? m : `${o.warning}. ${m}`; };
  const toks = lists.filter(t => t.length).map(t => t[0]); // replace_coords_rsid_list drops empty entries
  const vs = (await ctx.planner.query({ build, inputs: toks.filter(Boolean), select: "parity" })).variants;
  let vi = 0; const byTok = toks.map(t => (t ? vs[vi++] : null));
  byTok.forEach(v => v && addCoordWarning(o, v.warning)); // replace_coords_rsid_list (all entries first)
  const queries = [], warn = [];
  for (const [i, t] of toks.entries()) {
    const v = byTok[i];
    if (!t || (t.startsWith("rs") ? t : (v.rsid ?? t)).length <= 2) return { error: "Input list of RS numbers is empty" };
    const label = t.startsWith("rs") ? t : (v.rsid ?? t);
    if (!gate(label) || !v.found || v.pos == null) { warn.push(label); continue; }
    if (v.chrom === "Y" && build !== "grch37") { join(`Input variants on chromosome Y are unavailable for GRCh38, only available for GRCh37 (${v.rsid} = chr${v.chrom}:${v.pos})`); warn.push(label); continue; }
    queries.push({ label, v });
  }
  if (warn.length) join(`The following RS number(s) or coordinate(s) inputs have warnings: ${warn.join(", ")}`);
  if (!queries.length) return { error: `Input SNP list does not contain any valid RS numbers or coordinates. ${o.warning}` };

  const tissueIds = String(tissues).split("+"), out = [], stat = r2_d === "d" ? "dprime" : "r2";
  let popsGone = false;
  for (const { label, v } of queries) {
    // get_query_variant_c(..., is_output=True) deletes pops_<ref>.txt when a query is not in 1000G: every later
    // query's subprocesses then fail to read it and return nothing
    if (!v.selected?.row) { popsGone = true; continue; }
    if (popsGone) continue;
    let q = { ...(await ctx.row(build, v.chrom, v.selected.row)), chrom: v.chrom };
    const parts = chunkWindow(v.pos, window, chunks), all = await ctx.windowRows(build, v.chrom, parts[0][0], parts.at(-1)[1]);
    if (q.multiallelic_source) q = multiQuery(q, all, (await ctx.samples(build, v.chrom)).ids.length);
    const rows = parts.flatMap(([lo, hi]) => all.filter(x => x.pos <= hi && x.pos + x.ref.length - 1 >= lo && !x.multiallelic_source));
    const uniq = [...new Set(rows)], p = await ctx.prepare(build, [q, ...uniq.map(x => ({ ...x, chrom: v.chrom }))], pop);
    const counts = ctx.kernel.counts(p, uniq.map(() => 0), uniq.map((_, i) => i + 1)); ctx.kernel.reset();
    const cnt = new Map(uniq.map((x, i) => [x, counts[i]])), ld = [];
    for (const x of rows) {
      const [n, a1, a2, both] = cnt.get(x), D = both, C = a1 - both, B = a2 - both, A = n - a1 - a2 + both, N = n;
      const Ms = (A + C) * (B + D) * (A + B) * (C + D); if (Ms === 0) continue;
      const delta = A * D - B * C, dprime = Math.abs(delta / (delta < 0 ? Math.min((A + C) * (A + B), (B + D) * (C + D)) : Math.min((A + C) * (C + D), (A + B) * (B + D))));
      const r2 = (delta * delta) / Ms;
      if (!((stat === "r2" ? r2 : dprime) >= r2_d_threshold)) continue;
      const [pa0, pa1] = vcfAlleles(x.ref, x.alt), f = c => (N > c ? pyStr(pyRound(c / N, 3)) : "NA");
      ld.push({ x, r2, dprime, al1: `${pa0}=${f(A + C)}`, al2: `${pa1}=${f(B + D)}`, id: vcfId(x) });
    }
    const noRs = [...new Set(ld.filter(e => !e.id.startsWith("rs")).map(e => e.x.pos))], at = await ctx.dbsnpAt(build, v.chrom, noRs);
    const eq = await ctx.eqtl(build, v.chrom, ld.map(e => e.x.pos), tissueIds);
    for (const e of ld) {
      const rs = e.id.startsWith("rs") ? e.id : at.get(e.x.pos) ? `rs${at.get(e.x.pos).rsid}` : ".";
      for (const t of eq.get(e.x.pos) ?? []) if (t.pval_nominal < p_threshold)
        out.push([label, rs, `chr${v.chrom}:${e.x.pos}`, pyStr(e.r2), pyStr(e.dprime), t.gene_name_v26, t.gene_id, t.tissueSiteDetail, e.al1, e.al2, cFormatG(t.slope), cFormatG(t.pval_nominal)].join("\t"));
    }
  }
  if (!out.length) return { error: "No entries in GTEx are identified using the LDexpress search criteria." };
  const head = "Query\tRS ID\tPosition\tR2\tD'\tGene Symbol\tGencode ID\tTissue\tNon-effect Allele Freq\tEffect Allele Freq\tEffect Size\tP-value";
  return { text: [head, ...out].join("\n") + "\n" + (o.warning !== undefined ? `Warning(s):\n${o.warning}` : "") };
}
