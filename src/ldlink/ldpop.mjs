// LDpop on the core: LDlink's API text table (server/LDpop.py calculate_pop; the API returns LDpop_<ref>.txt) or
// { error }. One row per population: the requested ones plus, for a super-population, its sub-populations.
// LDpop's own quirks, reproduced: allele frequencies count a haploid call ("0" → "0|.") as one REF *and* one ALT allele
// over 2 alleles, and SNP2's frequencies divide by SNP1's allele total.
import { pyRoundRepr as r } from "../pyfmt.mjs";
import { decodeHaps } from "../haps.mjs";
import { TITLE } from "../genome.mjs";
import { POPS, POP_ERROR } from "./common.mjs";
import { ldStats, queryVariant, addCoordWarning } from "./common.mjs";

export const GROUPS = { AFR: ["YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"], AMR: ["MXL", "PUR", "CLM", "PEL"], EAS: ["CHB", "JPT", "CHS", "CDX", "KHV"], EUR: ["CEU", "TSI", "FIN", "GBR", "IBS"], SAS: ["GIH", "PJL", "BEB", "STU", "ITU"] };
// LDpop.py location_data, in pop_order (the table's row order): [location, super-population, latitude, longitude]
export const LOCATIONS = {
  ALL: ["All Populations"], AFR: ["African"],
  YRI: ["Yoruba in Ibadan, Nigeria", "AFR", 7.40026, 3.910742], LWK: ["Luhya in Webuye, Kenya", "AFR", 0.59738, 34.777227],
  GWD: ["Gambian in Western Divisions in the Gambia", "AFR", 13.474133, -16.394272], MSL: ["Mende in Sierra Leone", "AFR", 8.176076, -11.040253],
  ESN: ["Esan in Nigeria", "AFR", 6.687988, 6.212868], ASW: ["Americans of African Ancestry in SW USA", "AFR", 35.310647, -107.975885],
  ACB: ["African Caribbeans in Barbados", "AFR", 13.172483, -59.552779], AMR: ["Ad Mixed American"],
  MXL: ["Mexican Ancestry from Los Angeles USA", "AMR", 34.113837, -118.440427], PUR: ["Puerto Ricans from Puerto Rico", "AMR", 18.234429, -66.418775],
  CLM: ["Colombians from Medellin, Colombia", "AMR", 6.252089, -75.594652], PEL: ["Peruvians from Lima, Peru", "AMR", -12.046543, -77.046155],
  EAS: ["East Asian"], CHB: ["Han Chinese in Beijing, China", "EAS", 39.906802, 116.407323], JPT: ["Japanese in Tokyo, Japan", "EAS", 35.709444, 139.731815],
  CHS: ["Southern Han Chinese", "EAS", 24.719998, 113.043464], CDX: ["Chinese Dai in Xishuangbanna, China", "EAS", 22.008264, 100.796045],
  KHV: ["Kinh in Ho Chi Minh City, Vietnam", "EAS", 10.812236, 106.633978], EUR: ["European"],
  CEU: ["Utah Residents (CEPH) with Northern and Western European Ancestry", "EUR", 39.250493, -111.631295], TSI: ["Toscani in Italia", "EUR", 43.444187, 11.117199],
  FIN: ["Finnish in Finland", "EUR", 63.112, 26.770837], GBR: ["British in England and Scotland", "EUR", 54.55902, -2.143222],
  IBS: ["Iberian Population in Spain", "EUR", 40.482057, -4.088383], SAS: ["South Asian"],
  GIH: ["Gujarati Indian from Houston, Texas", "SAS", 29.760619, -95.361356], PJL: ["Punjabi from Lahore, Pakistan", "SAS", 31.515188, 74.357703],
  BEB: ["Bengali from Bangladesh", "SAS", 24.013458, 90.233561], STU: ["Sri Lankan Tamil from the UK", "SAS", 7.595905, 80.843382],
  ITU: ["Indian Telugu from the UK", "SAS", 15.489823, 78.487081],
};

/** per individual of `ids` (sample indexes): GT class of one row — "d" diploid, "h" haploid, null uncalled */
function calls(row, ids, haploidMales, hbytes) {
  const bits = decodeHaps(row.haps, hbytes), bit = k => (bits[k >> 3] >> (k & 7)) & 1, miss = k => row.missing && (row.missing[k >> 3] >> (k & 7)) & 1;
  return ids.map(i => (miss(2 * i) ? null : row.male_haploid && haploidMales.has(i) ? { h: true, a: bit(2 * i) } : { h: false, a: bit(2 * i) + bit(2 * i + 1) }));
}

export async function ldpop(ctx, { var1 = "", var2 = "", pop, build = "grch37" }) {
  if (!TITLE[build]) return { error: `Invalid genome build. Please specify either ${Object.keys(TITLE).join(", ")}.` };
  const in1 = var1.toLowerCase().trim(), in2 = var2.toLowerCase().trim(), o = {};
  const snps = [];
  for (const s of [in1, in2]) { // both replace_coord_rsid first, then both get_coords
    const v = /\s/.test(s) || !s ? { found: false } : (await ctx.planner.query({ build, inputs: [s], select: "parity" })).variants[0];
    addCoordWarning(o, v.warning);
    snps.push({ label: s.startsWith("rs") ? s : (v.rsid ?? s), v });
  }
  for (const { label, v } of snps) if (!v.found || v.pos == null) return { error: `${label} is not in dbSNP build 155 (${TITLE[build]}).` };
  const [s1, s2] = snps;
  if (s1.v.chrom !== s2.v.chrom) o.warning = `${o.warning ?? ""}${s1.label} and ${s2.label} are on different chromosomes. `;
  for (const { v } of snps) if (v.chrom === "Y" && build !== "grch37")
    return { error: `Input variants on chromosome Y are unavailable for GRCh38, only available for GRCh37 (${v.rsid} - chr${v.chrom}:${v.pos})` };
  let pops = String(pop).split("+");
  pops = pops.includes("ALL") ? [...pops, ...Object.values(GROUPS).flat(), "ALL", ...Object.keys(GROUPS)] : [...pops, ...pops.flatMap(p => GROUPS[p] ?? [])];
  pops = [...new Set(pops)];
  for (const p of pops) if (!POPS.includes(p)) return { error: POP_ERROR(p) };
  const rows = [];
  for (const { label, v } of snps) rows.push(await queryVariant(ctx, build, label, v, o));
  if (o.error !== undefined) return { error: o.error };

  const { ids: samples, haploidMales } = await ctx.samples(build, s1.v.chrom), hbytes = Math.ceil((2 * samples.length) / 8);
  const males = new Set([...haploidMales, ...(await ctx.samples(build, s2.v.chrom)).haploidMales]);
  const out = [];
  for (const P of Object.keys(LOCATIONS).filter(p => pops.includes(p))) {
    const panel = await ctx.popIds(P), idx = samples.flatMap((s, i) => (panel.has(s) ? [i] : []));
    const [g1, g2] = rows.map(row => calls(row, idx, males, hbytes));
    // allele counts as LDpop.py counts genotype strings; total from SNP1 only
    const tot = 2 * g1.filter(Boolean).length;
    if (!tot) return { error: `Insufficient haplotype data for ${s1.label} and ${s2.label} in 1000G reference panel.` };
    const freq = g => { let ref = 0, alt = 0; for (const c of g) if (c) { if (c.h) { ref++; alt++; } else { ref += 2 - c.a; alt += c.a; } } return [r((ref / tot) * 100, 2), r((alt / tot) * 100, 2)]; };
    const [f1, f2] = [freq(g1), freq(g2)];
    const p = await ctx.prepare(build, rows, P), [n, a1, a2, x] = ctx.kernel.counts(p, [0], [1])[0]; ctx.kernel.reset();
    const st = ldStats(n - a1 - a2 + x, a2 - x, a1 - x, x), num = v => (st ? r(v, 4) : "NA");
    out.push([LOCATIONS[P][0], P, panel.size, `${rows[0].ref}: ${f1[0]}%, ${rows[0].alt}: ${f1[1]}%`, `${rows[1].ref}: ${f2[0]}%, ${rows[1].alt}: ${f2[1]}%`,
      num(st?.r2), num(st?.dprime), num(st?.chisq), num(st?.p)].join("\t"));
  }
  return { text: [["Population", "Abbrev", "N", `${in1} Allele Freq`, `${in2} Allele Freq`, "R2", "D'", "Chisq", "P"].join("\t"), ...out].join("\n") + "\n" + (o.warning !== undefined ? `\n${o.warning}` : "") };
}
