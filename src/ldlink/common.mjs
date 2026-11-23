// Pieces several LDlink modules share (server/LDcommon.py and the per-module copies of the same code).
import { erfc } from "../ldsc/stats.mjs";

// LDcommon.get_population: the populations LDlink accepts
export const POPS = ["ALL", "AFR", "AMR", "EAS", "EUR", "SAS", "ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"];
export const POP_ERROR = p => `${p} is not an ancestral population. Choose one of the following ancestral populations: AFR, AMR, EAS, EUR, or SAS; or one of the following sub-populations: ACB, ASW, BEB, CDX, CEU, CHB, CHS, CLM, ESN, FIN, GBR, GIH, GWD, IBS, ITU, JPT, KHV, LWK, MSL, MXL, PEL, PJL, PUR, STU, TSI, or YRI.`;

/** LDcommon.check_allele / set_alleles: indels as "-" and the allele without its anchor base */
export const vcfAlleles = (ref, alt) => (ref.length === 1 && alt.length === 1 ? [ref, alt] : ref.length === 1 ? ["-", alt.slice(1)] : alt.length === 1 ? [ref.slice(1), "-"] : [ref.slice(1), alt.slice(1)]);

/**
 * The 2×2 statistics LDpair / LDpop / LDhap print (A..D = hap[sorted(hap)[0..3]]): D′, r², chi-square (Python's exact
 * integer numerator, then one division) and p = 2·(1 − Φ(√χ²)) written as LDlink does; null when Ms == 0 ("NA").
 */
export function ldStats(A, B, C, D) {
  const N = A + B + C + D, delta = A * D - B * C, Ms = (A + C) * (B + D) * (A + B) * (C + D);
  if (Ms === 0) return null;
  const dprime = Math.abs(delta / (delta < 0 ? Math.min((A + C) * (A + B), (B + D) * (C + D)) : Math.min((A + C) * (C + D), (A + B) * (B + D))));
  const r2 = (delta * delta) / Ms, chisq = Number(BigInt(N) * BigInt(delta) ** 2n) / Ms;
  const p = 2 * (1 - 0.5 * (1 + (1 - erfc(Math.sqrt(chisq) / Math.SQRT2))));
  return { dprime, r2, chisq, p };
}

/**
 * The VCF ID LDlink sees on the row (its "does not match RS number" warnings, LDproxy's RS_Number): vcf_id when the VCF
 * ID contains "rs" and is not the label, else "rs<rsid>" for a VCF-labelled row (rsid_source 1), else none (".").
 * Files built before vcf_id existed cannot tell a merged rsID (source 2) or an unlisted one (docs/CORE_REVIEW.md §7).
 */
export const vcfId = row => row.vcf_id ?? (row.rsid_source === 1 ? `rs${row.rsid}` : row.vcf_id === undefined && row.rsid_source === 2 ? "rs?(merged)" : ".");

/**
 * LDcommon.get_query_variant_c for one resolved SNP: the last overlapping 1000G row (planner select "parity"), with
 * LDlink's error/warning side effects on `o`. → the full genotype row ({ chrom } added), or null when not in 1000G.
 */
export async function queryVariant(ctx, build, label, v, o) {
  const sel = v.selected;
  if (!sel?.row) { o.error = `${label} Variant is not in 1000G reference panel.${o.error ?? ""}`; return null; }
  const row = { ...(await ctx.row(build, v.chrom, sel.row)), chrom: v.chrom }, id = vcfId(row);
  if (id !== label && id.includes("rs")) o.warning = `Genomic position does not match RS number at 1000G position (chr${v.chrom}:${row.pos} = ${id}).${o.warning ?? ""}`;
  if (sel.status === "not_biallelic") o.error = `${label} variant is not a biallelic.${o.warning ?? ""}`; // sic: LDlink appends the warning
  return row;
}

/** replace_coord_rsid's warning, joined as LDlink joins them */
export const addCoordWarning = (o, w) => { if (w) o.warning = o.warning === undefined ? w : `${o.warning}. ${w}`; };

/**
 * dbSNP function (db.dbsnp.find_one({id})) and FORGEdb score (db.forge_score.find_one({snp_id})) for rsIDs whose doc
 * position is known: items [{ id: digits, pos }] → Map(id → { func: string | null | undefined, forge: number | undefined }).
 * Fast path, by position: for an rsID with a single dbSNP doc (not in duplicated_rsids) that doc is the one at `pos`,
 * and the per-build FORGEdb file holds each rsID once at its first doc's position — so a positional match is exact.
 * Anything else (duplicated rsIDs, not found where expected) takes the find_one path (scattered rsID-index reads).
 * func: undefined = no dbSNP doc ("."); null = the dump's "NA".
 */
export async function rsAnnotations(ctx, build, chrom, items) {
  const dup = await ctx.duplicatedRsids(), out = new Map(), slow = [];
  const fast = items.filter(i => i.pos != null && !dup.has(i.id)), positions = [...new Set(fast.map(i => i.pos))];
  const [docs, fg] = await Promise.all([ctx.dbsnpAllAt(build, chrom, positions), ctx.forgeAt(build, chrom, positions)]);
  for (const i of fast) {
    const d = docs.get(i.pos)?.find(r => String(r.rsid) === i.id);
    if (!d) { slow.push(i); continue; }
    const f = fg.get(i.pos)?.find(r => String(r.rsid) === i.id);
    out.set(i.id, { func: d.function, forge: f?.score });
  }
  const rest = [...new Set([...items.filter(i => i.pos == null || dup.has(i.id)), ...slow].map(i => i.id))].filter(id => !out.has(id));
  if (rest.length) {
    const [byId, forge] = await Promise.all([ctx.dbsnpById(rest), ctx.forge(rest)]);
    for (const id of rest) out.set(id, { func: byId.get(id) ? byId.get(id).function : undefined, forge: forge.get(id) });
  }
  return out;
}
