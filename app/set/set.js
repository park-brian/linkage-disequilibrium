// The variant set: what the page's one input means. A region ("chr22:30,000,000-31,000,000") is every common 1000G
// variant in it; anything else is a list of rsIDs / chr:pos, resolved as LDlink resolves them (the planner with LDlink's
// row choice, as LDmatrix does). Each usable row carries its genotype row, so a lens can compute on it directly, and
// its allele frequency in the chosen populations (the kernel's counts of a row against itself).
import { vcfAlleles } from "../../src/parity/ldlink/common.js";

export const REGION = /^(?:chr)?([0-9]{1,2}|X|Y)\s*:\s*([\d,_]+)\s*-\s*([\d,_]+)$/i;
const TOKEN = /^(rs\d+|chr([0-9]{1,2}|x|y):\d+)$/;
export const MAX_SPAN = 5_000_000, MAX_REGION_ROWS = 60_000;
const num = s => +s.replace(/[,_]/g, "");

/** → { kind: "list" | "region", rows: [{ label, chrom?, pos?, alleles?, maf?, g?, status? }], ok: rows with genotypes,
 *      chrom?, lo?, hi?, error? } */
export async function resolveSet(ctx, { q, build, pop }) {
  const text = q.trim(), m = REGION.exec(text);
  if (!text) return { kind: "list", rows: [], ok: [] };
  const set = m ? await regionSet(ctx, build, m[1].toUpperCase(), num(m[2]), num(m[3])) : await listSet(ctx, build, text);
  set.ok = set.rows.filter(r => r.g);
  if (set.ok.length) {
    const I = set.ok.map((_, i) => i), c = await ctx.withPrepared(build, set.ok.map(r => r.g), pop, p => ctx.kernel.counts(p, I, I));
    set.ok.forEach((r, i) => { const f = c[i][0] ? c[i][1] / c[i][0] : NaN; r.af = f; r.maf = Math.min(f, 1 - f); });
  }
  return set;
}

async function regionSet(ctx, build, chrom, lo, hi) {
  const bad = error => ({ kind: "region", rows: [], chrom, lo, hi, error });
  if (!(hi > lo)) return bad("The region's end must come after its start.");
  if (hi - lo > MAX_SPAN) return bad("Regions up to 5 Mb for now.");
  if (chrom === "Y" && build !== "grch37") return bad("Chromosome Y is only in the GRCh37 panel.");
  const got = (await ctx.windowRows(build, chrom, lo, hi, ["common"])).filter(r => r.pos >= lo && r.pos <= hi)
    .sort((a, b) => a.pos - b.pos || a.vcf_order - b.vcf_order);
  if (got.length > MAX_REGION_ROWS) return bad(`${got.length.toLocaleString()} common variants: narrow the region.`);
  return { kind: "region", chrom, lo, hi, rows: got.map(r => ({ label: r.rsid != null ? `rs${r.rsid}` : `chr${chrom}:${r.pos}`, chrom,
    pos: r.pos, alleles: vcfAlleles(r.ref, r.alt).join("/"), g: { ...r, chrom } })) };
}

async function listSet(ctx, build, text) {
  const toks = [...new Set(text.toLowerCase().split(/[\s,;]+/).filter(Boolean))], good = toks.filter(t => TOKEN.test(t));
  const res = good.length ? (await ctx.planner.query({ build, inputs: good, select: "parity" })).variants : [], seen = new Set();
  const byTok = new Map(good.map((t, i) => [t, res[i]]));
  const rows = await Promise.all(toks.map(async t => {
    const v = byTok.get(t);
    if (!v) return { label: t, status: "not an rsID or chr:pos" };
    const label = v.rsid ?? t, at = { label, chrom: v.chrom, pos: v.pos };
    if (!v.found || v.pos == null) return { label, status: "not in dbSNP 155 for this build" };
    if (v.chrom === "Y" && build !== "grch37") return { ...at, status: "chromosome Y is GRCh37 only" };
    if (!v.selected?.row) return { ...at, status: "not in 1000 Genomes" };
    if (v.selected.status === "not_biallelic") return { ...at, status: "not biallelic" };
    const g = await ctx.row(build, v.chrom, v.selected.row);
    return { ...at, alleles: vcfAlleles(g.ref, g.alt).join("/"), g: { ...g, chrom: v.chrom } };
  }));
  return { kind: "list", rows: rows.filter(r => !seen.has(r.label) && seen.add(r.label)) };
}

/** the genomic window a set is drawn in: the region itself, or a list's span with a margin (one variant: ± 250 kb) */
export function windowOf(set) {
  if (set.kind === "region") return { chrom: set.chrom, lo: set.lo, hi: set.hi };
  const chroms = [...new Set(set.ok.map(r => r.chrom))];
  if (!chroms.length) return { error: "No variant in this set is in 1000 Genomes." };
  if (chroms.length > 1) return { error: `This set spans chromosomes ${chroms.join(", ")}. `
    + `A locus is one chromosome: select a variant to see its own.` };
  const pos = set.ok.map(r => r.pos), a = Math.min(...pos), b = Math.max(...pos), pad = b > a ? Math.max(20_000, Math.round((b - a) * 0.15))
    : 250_000;
  if (b - a + 2 * pad > MAX_SPAN) return { error: "These variants span more than 5 Mb. Select a closer group to see their locus." };
  return { chrom: chroms[0], lo: Math.max(1, a - pad), hi: b + pad };
}
