// Genotype-row selection for a resolved variant (rsID → chrom/pos), two modes:
//   parity — LDlink: tabix chr:p-p returns every VCF row whose REF overlaps p; LDlink keeps the LAST one
//            (LDcommon.py:590-603, the POS check there is dead code). With MAF tiers split into files, "last" is the
//            highest vcf_order. A source record that was multi-allelic → LDlink's "Variant is not a biallelic." error.
//   allele — opt-in: the row whose allele_rsids (NCBI dbSNP alleles, build time) contains the query's current rsID.
//            Reports when no 1000G row carries the rsID's alleles instead of returning a different variant.
// Rows must have pos, ref, vcf_order (+ rsid, allele_is_label, allele_rsids for allele mode), and cover [p - max_ref_len + 1, p] (file
// key-value metadata `max_ref_len`), across all tiers, for a complete candidate set.
import { rangeScan, pool } from "./parquet.js";
import { tierPath } from "./genome.js";
import { kv, required } from "./source.js";
export { decodeHaps } from "./haps.js";
import { parquetMetadataAsync } from "hyparquet";
import { compressors } from "hyparquet-compressors";

export const overlaps = (r, p) => r.pos <= p && r.pos + r.ref.length - 1 >= p;
/** rows → (lo, hi) → the rows with lo ≤ pos ≤ hi in the input order (a binary search over a pos-sorted index, not a
 *  scan of every row per query: uploads ask about hundreds of thousands of positions) */
export function posRange(rows) {
  const ix = rows.map((_, i) => i).sort((i, j) => rows[i].pos - rows[j].pos || i - j);
  return (lo, hi) => {
    let a = 0, b = ix.length;
    while (a < b) { const m = (a + b) >> 1; if (rows[ix[m]].pos < lo) a = m + 1; else b = m; }
    const hit = [];
    for (let k = a; k < ix.length && rows[ix[k]].pos <= hi; k++) hit.push(ix[k]);
    return hit.sort((i, j) => i - j).map(i => rows[i]);
  };
}
/** rows → p → the rows overlapping p (overlaps(), tabix p-p) in the input order */
export function overlapping(rows) {
  const range = posRange(rows), reach = rows.reduce((m, r) => Math.max(m, r.ref.length), 1) - 1;
  return p => range(p - reach, p).filter(x => overlaps(x, p));
}

/** LDlink parity: last overlapping row in VCF order, or null ("Variant is not in 1000G reference panel.") */
export function pickParity(rows, p) {
  let best = null;
  for (const r of rows) if (overlaps(r, p) && (best === null || r.vcf_order > best.vcf_order)) best = r;
  if (!best) return { status: "not_in_1000g", row: null };
  return { status: best.multiallelic_source ? "not_biallelic" : "ok", row: best };
}

// allele_is_label: the row's allele matches are exactly [rsid]; else allele_rsids lists them (null = none)
const hasRsid = (r, id) => (r.allele_is_label ? String(r.rsid) === id : r.allele_rsids != null && r.allele_rsids.split(",").includes(id));

/**
 * Allele-aware: rows in the window whose allele_rsids contain `rsid` (current rsID digits, e.g. "4672393"; for a merged
 * input use the dbSNP doc's ref_id). status: "ok" | "ok_split_multiallelic" (one ALT of a multi-allelic record) |
 * "multiallelic_in_1000g" (the matches are the split ALTs of one multi-allelic record; LDlink: "not biallelic") |
 * "ambiguous" (several distinct rows; first by VCF order returned, all in `candidates`) |
 * "alleles_not_in_1000g" (rows overlap p, none carry the rsID's alleles; `parity` is what LDlink would have used) |
 * "not_in_1000g".
 */
export function pickAllele(rows, p, rsid) {
  const id = String(rsid).replace(/^rs/i, "");
  const m = rows.filter(r => hasRsid(r, id)).sort((a, b) => a.vcf_order - b.vcf_order);
  const parity = pickParity(rows, p);
  if (m.length === 1) return { status: m[0].multiallelic_source ? "ok_split_multiallelic" : "ok", row: m[0], parity: parity.row };
  // several matches that are the split ALTs of one multi-allelic 1000G record (same pos/ref): LDlink's "not biallelic"
  if (m.length > 1 && m.every(r => r.multiallelic_source && r.pos === m[0].pos && r.ref === m[0].ref))
    return { status: "multiallelic_in_1000g", row: m[0], candidates: m, parity: parity.row };
  if (m.length > 1) return { status: "ambiguous", row: m[0], candidates: m, parity: parity.row };
  return { status: parity.row ? "alleles_not_in_1000g" : "not_in_1000g", row: null, parity: parity.row };
}

/**
 * Fetch candidate rows for positions on one chromosome from the tier files (page-level range reads), merging
 * windows closer than `gap` into one read.
 * @param files [{ tier, file, metadata }]  (metadata may be omitted; it is read)
 */
export async function candidateRows(files, positions, { columns, gap = 10_000 } = {}) {
  const fs = await Promise.all(files.map(async f => ({ ...f, metadata: f.metadata ?? await parquetMetadataAsync(f.file) })));
  const reach = Math.max(1, ...fs.map(f => +((f.metadata.key_value_metadata ?? []).find(k => k.key === "max_ref_len")?.value ?? 1)));
  const ps = [...new Set(positions)].sort((a, b) => a - b), wins = [];
  for (const p of ps) { const lo = p - reach + 1, w = wins.at(-1); if (w && lo - w.hi <= gap) w.hi = p; else wins.push({ lo, hi: p }); }
  const out = [];
  await pool(fs.flatMap(f => wins.map(w => ({ f, w }))), 8, async ({ f, w }) => {
    const r = await rangeScan({ file: f.file, metadata: f.metadata, column: "pos", lo: w.lo, hi: w.hi, columns, compressors });
    for (const row of r.rows ?? r) out.push({ ...row, tier: f.tier });
  });
  return { rows: out, reach };
}


/** the columns a module needs from a genotype row (vcf_id: files built with it, GRCh37 from release v4) */
const ROW_COLS = ["pos", "vcf_order", "rsid", "rsid_source", "rsid_pos_delta", "ref", "alt", "multiallelic_source", "haps", "missing", "male_haploid"];
const TIERS = ["common", "low", "rare"];

/** Genotype rows and samples over a source's file(rel) → { file, metadata } (e.g. planner.file). */
export function createGenotypes(fileOrNull) {
  const file = required(fileOrNull); // every genotype file a query names must exist
  const colsOf = metadata => ROW_COLS.concat(metadata.schema.some(c => c.name === "vcf_id") ? ["vcf_id"] : []);
  /** the full row for a selected row reference { tier, pos, vcf_order } (planner select) */
  async function row(build, chrom, ref) {
    const { file: f, metadata } = await file(tierPath(build, chrom, ref.tier));
    const got = await rangeScan({ file: f, metadata, column: "pos", lo: ref.pos, hi: ref.pos, columns: colsOf(metadata), compressors });
    return got.find(x => x.vcf_order === ref.vcf_order);
  }
  /** every row overlapping [lo, hi] (tabix semantics: pos ≤ hi and pos + len(ref) − 1 ≥ lo), in the given MAF tiers (all by default), VCF order */
  async function windowRows(build, chrom, lo, hi, tiers = TIERS) {
    const out = [];
    for (const tier of tiers) {
      const { file: f, metadata } = await file(tierPath(build, chrom, tier)), reach = +kv(metadata).max_ref_len - 1;
      for (const r of await rangeScan({ file: f, metadata, column: "pos", lo: lo - reach, hi, columns: colsOf(metadata), compressors }))
        if (r.pos + r.ref.length - 1 >= lo) out.push(r);
    }
    return out.sort((a, b) => a.vcf_order - b.vcf_order);
  }
  /** sample IDs and haploid-male sample indexes of a chromosome (key-value metadata of its common tier) */
  async function samples(build, chrom) {
    const m = kv((await file(tierPath(build, chrom, "common"))).metadata);
    return { ids: m.samples.split(","), haploidMales: new Set((m.male_haploid_samples || "").split(",").filter(Boolean).map(Number)) };
  }
  return { row, windowRows, samples };
}
