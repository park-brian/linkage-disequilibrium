// LDhap on the core: LDlink's JSON (server/LDhap.py calculate_hap) and the API's text (snps + haplotypes files).
// Haplotypes are whole phased chromosomes over the kept SNPs; one containing a haploid or uncalled allele is dropped
// (chrX males lose their second haplotype). Order: count descending, ties by haplotype string — deterministic.
import { pyRound, pyStr } from "../pyfmt.js";
import { decodeHaps } from "../haps.js";
import { TITLE, BUILDS } from "../genome.js";
import { POPS, POP_ERROR, chrYMessage, dedupeLists, vcfAlleles, addCoordWarning, gate } from "./common.js";

const pad = (i, n) => String(i).padStart(String(n).length, "0");

export async function ldhap(ctx, { snps: lines, pop, build = "grch37" }) {
  if (!TITLE[build]) return { error: `Invalid genome build. Please specify either ${Object.keys(TITLE).join(", ")}.` };
  // validsnp: one token list per line, duplicate lists removed (before coordinates are replaced)
  const lists = dedupeLists(lines.map(l => l.toLowerCase()));
  if (lists.length > 30) return { error: `Maximum variant list is 30  RS numbers or coordinates. Your list contains ${lists.length} entries.` };
  for (const p of String(pop).split("+")) if (!POPS.includes(p)) return { error: POP_ERROR(p) };
  const o = {}, W = () => o.warning ?? "", warnAdd = s => { o.warning = W() + s; };
  const fail = error => ({ error, ...(o.warning !== undefined ? { warning: o.warning } : {}) }); // calculate_hap returns the whole output
  const toks = lists.filter(t => t.length).map(t => t[0]);
  const vs = (await ctx.planner.query({ build, inputs: toks, select: "parity" })).variants;
  vs.forEach(v => addCoordWarning(o, v.warning)); // replace_coords_rsid_list, all inputs first
  const kept = [], warn = [];
  toks.forEach((t, i) => {
    const v = vs[i], label = t.startsWith("rs") ? t : (v.rsid ?? t);
    if (!gate(label) || !v.found || v.pos == null) return warn.push(label);
    if (v.chrom === "Y" && build !== "grch37") { warnAdd(`${chrYMessage(v)}. `); return warn.push(label); }
    kept.push({ label, v });
  });
  if (warn.length) warnAdd(`The following RS number(s) or coordinate(s) inputs have warnings: ${warn.join(", ")}. `);
  if (!kept.length) return fail(`Input variant list does not contain any valid RS numbers or coordinates. ${W()}`);
  const diff = kept.findIndex(k => k.v.chrom !== kept[0].v.chrom); // check_same_chromosome: reported after parse_vcf
  if (diff > 0) o.error = `Not all input variants are on the same chromosome: ${kept[diff - 1].label}=chr${kept[diff - 1].v.chrom}:${kept[diff - 1].v.pos}, ${kept[diff].label}=chr${kept[diff].v.chrom}:${kept[diff].v.pos}. ${W()}`;
  const ps = kept.map(k => k.v.pos), range = Math.max(...ps) - Math.min(...ps);
  if (range > 1e6) warnAdd(`Switch rate errors become more common as distance between query variants increases (Query range = ${range} bp). `);
  // parse_vcf: tabix on the first SNP's chromosome; entries keyed chr:pos_label (dups collapse), sorted by position
  const chrom = kept[0].v.chrom, found = kept.filter(k => k.v.chrom === chrom && k.v.selected?.row);
  const missing = kept.filter(k => !found.includes(k)).map(k => k.label);
  if (missing.length === kept.length) return fail(`Input variant list does not contain any valid RS numbers or coordinates. ${W()}`);
  if (missing.length) o.warning = `Query variant ${missing.join(" ")} is missing from 1000G (${TITLE[build]}) data. ${W()}`;
  if (o.error !== undefined) return fail(o.error);
  const entries = [...new Map(found.map(k => [`${k.v.pos}_${k.label}`, k])).values()].sort((a, b) => a.v.pos - b.v.pos);

  const { ids: samples, haploidMales } = await ctx.samples(build, chrom), hbytes = Math.ceil((2 * samples.length) / 8);
  const panel = await ctx.popIds(pop), idx = samples.flatMap((s, i) => (panel.has(s) ? [i] : []));
  const h1 = idx.map(() => []), h2 = idx.map(() => []), snps = [];
  for (const { label, v } of entries) {
    const row = await ctx.row(build, chrom, v.selected.row);
    if (row.pos !== v.pos) warnAdd(`Genomic position (${row.pos}) in 1000G data does not match dbSNP155 (${TITLE[build]}) search coordinates for query variant ${label}. `);
    if (v.selected.status === "not_biallelic") {
      if (entries.length > 1) { warnAdd(`${label} variant is not a biallelic`); continue; }
      return fail(`${o.error ?? ""}${label} variant is not a biallelic`);
    }
    const [a1, a2] = vcfAlleles(row.ref, row.alt), bits = decodeHaps(row.haps, hbytes), bit = k => (bits[k >> 3] >> (k & 7)) & 1;
    let c0 = 0, c1 = 0;
    idx.forEach((s, j) => {
      if (row.missing && (row.missing[(2 * s) >> 3] >> ((2 * s) & 7)) & 1) { h1[j].push("."); h2[j].push("."); return; }
      const x = bit(2 * s), hap = row.male_haploid && haploidMales.has(s), y = hap ? null : bit(2 * s + 1);
      h1[j].push(x ? a2 : a1); h2[j].push(y === null ? "." : y ? a2 : a1);
      for (const b of y === null ? [x] : [x, y]) b ? c1++ : c0++;
    });
    const f0 = pyRound(c0 / (c0 + c1), 4), f1 = pyRound(c1 / (c0 + c1), 4), s3 = f => pyStr(pyRound(f, 3));
    snps.push({ Alleles: f0 >= f1 ? `${a1}=${s3(f0)}, ${a2}=${s3(f1)}` : `${a2}=${s3(f1)}, ${a1}=${s3(f0)}`, Coord: `chr${chrom}:${v.pos}`, RS: label });
  }
  const counts = new Map();
  for (const h of [h1, h2].flatMap(x => x)) { const k = h.join("_"); counts.set(k, (counts.get(k) ?? 0) + 1); }
  for (const k of [...counts.keys()]) if (k.includes(".")) counts.delete(k);
  const total = [...counts.values()].reduce((a, b) => a + b, 0);
  const haps = [...counts].sort((a, b) => (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0)).sort((a, b) => b[1] - a[1]);
  const out = {
    haplotypes: Object.fromEntries(haps.map(([k, n], i) => [`haplotype_${pad(i + 1, haps.length)}`, { Count: n, Frequency: pyRound(n / total, 4), Haplotype: k }])),
    snps: Object.fromEntries(snps.map((s, i) => [`snp_${pad(i + 1, snps.length)}`, s])),
    ...(o.warning !== undefined ? { warning: o.warning } : {}),
  };
  return out;
}

/** the API's text: snps_<ref>.txt + separator + haplotypes_<ref>.txt (LDhap.py file writers, LDlink.py ldhap route) */
export function ldhapText(out, build) {
  const s = ["RS_Number\tPosition (" + BUILDS[build].hg + ")\tAllele Frequency", ...Object.keys(out.snps).sort().map(k => {
    const x = out.snps[k], [a, b] = x.Alleles.trim().split(",");
    return [x.RS, x.Coord, `${a}${"0".repeat(Math.max(0, 7 - a.length))},${b}${"0".repeat(Math.max(0, 8 - b.length))}`].join("\t");
  })].join("\n") + "\n";
  const h = ["Haplotype\tCount\tFrequency", ...Object.keys(out.haplotypes).sort().map(k => { const x = out.haplotypes[k]; return [x.Haplotype, x.Count, pyStr(x.Frequency)].join("\t"); })].join("\n") + "\n";
  return `${s}\n${"#".repeat(85)}\n\n${h}`;
}
