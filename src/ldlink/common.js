// Pieces several LDlink modules share (server/LDcommon.py and the per-module copies of the same code).
import { erfc } from "../ldsc/stats.js";

// LDcommon.get_population: the populations LDlink accepts
export const POPS = ["ALL", "AFR", "AMR", "EAS", "EUR", "SAS", "ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"];
export const POP_ERROR = p => `${p} is not an ancestral population. Choose one of the following ancestral populations: AFR, AMR, EAS, EUR, or SAS; or one of the following sub-populations: ACB, ASW, BEB, CDX, CEU, CHB, CHS, CLM, ESN, FIN, GBR, GIH, GWD, IBS, ITU, JPT, KHV, LWK, MSL, MXL, PEL, PJL, PUR, STU, TSI, or YRI.`;

/** LDcommon.check_allele / set_alleles: indels as "-" and the allele without its anchor base */
export const vcfAlleles = (ref, alt) => (ref.length === 1 && alt.length === 1 ? [ref, alt] : ref.length === 1 ? ["-", alt.slice(1)] : alt.length === 1 ? [ref.slice(1), "-"] : [ref.slice(1), alt.slice(1)]);

/** The chrY-on-GRCh38 message; LDpair and LDpop write "rsid - chrY:pos", the rest "rsid = chrY:pos" */
export const chrYMessage = (v, sep = "=") =>
  `Input variants on chromosome Y are unavailable for GRCh38, only available for GRCh37 (${v.rsid} ${sep} chr${v.chrom}:${v.pos})`;

/** The input gate every module applies to a token (e.g. LDmatrix.py:71-72): longer than 2, starts rs/chr, ends in a digit. */
export const gate = t => t.length > 2 && (t.startsWith("rs") || t.startsWith("chr")) && /\d$/.test(t);

/** validsnp: one whitespace token list per line, duplicate lists removed (first kept; empty lists stay, once) */
export function dedupeLists(lines) {
  const seen = new Set(), lists = [];
  for (const l of lines) { const t = l.trim().split(/\s+/).filter(Boolean), k = t.join(" "); if (!seen.has(k)) { seen.add(k); lists.push(t); } }
  return lists;
}

/** kernel counts [n, a1, a2, both] (alt counts) → LDlink's 2×2 [A, B, C, D] (ref/ref … alt/alt) */
export const cells = ([n, a1, a2, both]) => [n - a1 - a2 + both, a2 - both, a1 - both, both];

/**
 * LD_calcs' core on a 2×2 (A..D = hap[sorted(hap)[0..3]]): D′ with the sign-dependent Dmax, r² = Δ²/Ms;
 * null when Ms == 0 (LDlink skips the pair or prints "NA").
 */
export function table2x2(A, B, C, D) {
  const delta = A * D - B * C, Ms = (A + C) * (B + D) * (A + B) * (C + D);
  if (Ms === 0) return null;
  const dmax = delta < 0 ? Math.min((A + C) * (A + B), (B + D) * (C + D)) : Math.min((A + C) * (C + D), (A + B) * (B + D));
  return { delta, Ms, dprime: Math.abs(delta / dmax), r2: (delta * delta) / Ms };
}

/**
 * The 2×2 statistics LDpair / LDpop / LDhap print: table2x2 plus chi-square (Python's exact integer numerator, then
 * one division) and p = 2·(1 − Φ(√χ²)) written as LDlink does; null when Ms == 0 ("NA").
 */
export function ldStats(A, B, C, D) {
  const t = table2x2(A, B, C, D);
  if (!t) return null;
  const chisq = Number(BigInt(A + B + C + D) * BigInt(t.delta) ** 2n) / t.Ms;
  const p = 2 * (1 - 0.5 * (1 + (1 - erfc(Math.sqrt(chisq) / Math.SQRT2))));
  return { dprime: t.dprime, r2: t.r2, chisq, p };
}

/**
 * LDlink's chunked tabix window (LDproxy / LDexpress): the window cut into `chunks` pieces, each read with overlap
 * semantics, so a variant spanning a chunk boundary is listed twice; biallelic rows only, chunk order.
 * → { all: every row read, rows }
 */
export async function chunkedWindow(ctx, build, chrom, pos, window, chunks) {
  const parts = chunkWindow(pos, window, chunks), all = await ctx.windowRows(build, chrom, parts[0][0], parts.at(-1)[1]);
  return { all, rows: parts.flatMap(([lo, hi]) => all.filter(x => x.pos <= hi && x.pos + x.ref.length - 1 >= lo && !x.multiallelic_source)) };
}

/** LDcommon.chunkWindow → [[start, stop]] (Python floats, math.ceil; the last chunk ends a few bp past pos + window) */
export function chunkWindow(pos, window, n) {
  const min = pos - window <= 0 ? 0 : pos - window, range = pos + window - min, out = [];
  let lo = min;
  for (let i = 0; i < n; i++) { const hi = lo + range / n; out.push([Math.ceil(lo), Math.ceil(hi)]); lo = hi + 1; }
  return out;
}

/** the query row against each row (duplicates counted once) in one kernel call → Map(row → [A, B, C, D]) */
export async function countsAgainst(ctx, build, chrom, q, rows, pop) {
  const uniq = [...new Set(rows)], p = await ctx.prepare(build, [q, ...uniq.map(x => ({ ...x, chrom }))], pop);
  const counts = ctx.kernel.counts(p, uniq.map(() => 0), uniq.map((_, i) => i + 1)); ctx.kernel.reset();
  return new Map(uniq.map((x, i) => [x, cells(counts[i])]));
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
