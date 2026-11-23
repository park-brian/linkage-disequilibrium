// LDassoc on the core: LDlink's table (assoc<ref>.txt) and the error/warning of assoc<ref>.json (server/LDassoc.py
// calculate_assoc + LDassoc_sub.py; the web route is the only entry point, there is no API). For a region of an
// association file: LD of every 1000G row at the file's positions with an origin variant, sorted by p-value.
// LDlink's mechanics, reproduced:
//  - region: variant (origin ± window), gene (genes_name_coords ± window), or region (start..end, window 0);
//    max 3 Mb; "region" compares the coordinate *strings* (coord_s >= coord_e)
//  - no origin: the lowest-p position whose tabix lookup is in 1000G, biallelic (first row) and polymorphic; on High
//    Coverage that lookup asks tabix for "chrchr8:…" and always fails (LDlink bug, reproduced)
//  - the origin genotype is the FIRST row overlapping its position (snp_no_dups line 0)
//  - one tabix region per association position, in file order: a row is listed under every position it overlaps,
//    and kept when its own POS is one of the file's positions (with that position's last p-value)
//  - sort |distance| then p ascending (stable); first line: the top row with the origin's Coord/Alleles and MAF
//  - Correlated_Alleles is written proxy-first, with the allele indexes crossed (allele_n[q] = allele[p], sic)
//  - a coordinate origin raises KeyError in LDlink (var_coord has 'position', not 'position_grchNN')
import { pyRound, pyStr } from "../pyfmt.mjs";
import { decodeHaps } from "../haps.mjs";
import { TITLE, BUILDS, CHROMS } from "../genome.mjs";
import { POPS, POP_ERROR } from "./common.mjs";
import { vcfAlleles, vcfId, rsAnnotations } from "./common.mjs";

const INTERNAL = { error: "Internal server error. Please contact LDlink admin." };
const stripChr = s => s.replace(/^[chr]+|[chr]+$/g, ""); // str.strip("chr")
const PYFLOAT = /^\s*[+-]?((\d(_?\d)*(\.(\d(_?\d)*)?)?|\.\d(_?\d)*)([eE][+-]?\d(_?\d)*)?|inf(inity)?|nan)\s*$/i;
const pyFloat = s => (/inf/i.test(s) ? (s.includes("-") ? -Infinity : Infinity) : /nan/i.test(s) ? NaN : Number(s.trim().replace(/_/g, "")));
const PYINT = /^\s*[+-]?\d(_?\d)*\s*$/;

/**
 * @param {{ text: string, columns: [chr, bp, pval], region: "variant"|"gene"|"region", origin?, window?, name?, start?, end?, pop, build }} a
 * @returns {{ text?, error?, warning? }}
 */
export async function ldassoc(ctx, a) {
  const { build, pop } = a, o = {}, W = m => { o.warning = o.warning === undefined ? m : `${o.warning}. ${m}`; };
  const fail = error => ({ error, ...(o.warning !== undefined ? { warning: o.warning } : {}) });
  if (a.region === "variant" && a.origin == null) return fail("--origin required when --variant is specified.");
  let chromosome, org, snp = a.origin;
  if (a.origin != null) {
    if (a.origin.slice(0, 2) === "rs") {
      const v = (await ctx.planner.query({ build, inputs: [a.origin] })).variants[0];
      if (!v.found) return fail(`${a.origin} is not in dbSNP 155 (${TITLE[build]}).`);
      if (v.chrom === "Y" && build !== "grch37") return fail(`Input variants on chromosome Y are unavailable for GRCh38, only available for GRCh37 (${v.rsid} = chr${v.chrom}:${v.pos})`);
      if (v.pos == null) return INTERNAL; // int("NA")
      chromosome = v.chrom; org = v.pos;
    } else if (CHROMS.includes(stripChr(a.origin.split(":")[0])) && a.origin.split(":").length === 2) return INTERNAL;
    else return fail(`--origin (${a.origin}) is not a RS number (ex: rs12345) or chromosomal position (ex: chr22:25855459).`);
  }
  const lines = a.text.split("\n"), header = (lines[0] ?? "").trim().split(/\s+/).filter(Boolean), first = (lines[1] ?? "").trim().split(/\s+/).filter(Boolean);
  if (header.length !== first.length) return fail(`Header has ${header.length} elements and first line has ${first.length} elements.`);
  if (!a.columns.every(c => header.includes(c))) return fail("Variables mapping is not listed in the association file header.");
  const [ci, pi, vi] = a.columns.map(c => header.indexOf(c));
  const window = a.window ?? (a.region === "variant" ? 500000 : a.region === "gene" ? 100000 : 0);
  let c1, c2;
  if (a.region === "variant") { c1 = Math.max(0, org - window); c2 = org + window; }
  else if (a.region === "gene") {
    if (a.name == null) return fail("Gene name (--name) is needed when --gene option is used.");
    const g = await ctx.gene(a.name.toUpperCase(), build);
    if (!g || g.begin == null) return fail(`Gene name ${a.name} is not in RefSeq database.`);
    c1 = Math.max(0, g.begin - window); c2 = g.end + window;
    if (a.origin != null) {
      if (g.chrom !== chromosome) return fail(`Origin variant ${a.origin} is not on the same chromosome as ${a.name} (chr${chromosome} is not equal to chr${g.chrom}).`);
      if (c1 > org || org > c2) return fail(`Origin variant ${a.origin} (chr${chromosome}:${org}) is not in the coordinate range chr${g.chrom}:${c1}-${c2}.`);
    } else chromosome = g.chrom;
  } else if (a.region === "region") {
    if (a.start == null) return fail("Start coordinate is needed when --region option is used.");
    if (a.end == null) return fail("End coordinate is needed when --region option is used.");
    if (a.start.split(":").length !== 2) return fail("Start coordinate is not in correct format (ex: chr22:25855459).");
    if (a.end.split(":").length !== 2) return fail("End coordinate is not in correct format (ex: chr22:25855459).");
    const up = c => (c === "x" || c === "y" ? c.toUpperCase() : c);
    const cs = up(stripChr(a.start).split(":")[0]), ps = a.start.split(":")[1], ce = up(stripChr(a.end).split(":")[0]), pe = a.end.split(":")[1];
    if (!CHROMS.includes(cs)) return fail(`Start chromosome (chr${cs}) is not an autosome (chr1-chr22) or sex chromosome (chrX or chrY).`);
    if (!CHROMS.includes(ce)) return fail(`End chromosome (chr${ce}) is not an autosome (chr1-chr22) or sex chromosome (chrX or chrY).`);
    if (cs !== ce) return fail(`Start and end chromosome must be the same (chr${cs} is not equal to chr${ce}).`);
    if (ps >= pe) return fail(`End coordinate (${a.end}) must be greater than start coordinate(${a.start}).`);
    if (!PYINT.test(ps) || !PYINT.test(pe)) return INTERNAL;
    c1 = Math.max(0, +ps.replace(/_/g, "") - window); c2 = +pe.replace(/_/g, "") + window;
    if (a.origin != null) {
      if (cs !== chromosome) return fail(`Origin variant ${a.origin} is not on the same chromosome as start and stop coordinates (chr${chromosome} is not equal to chr${ce}).`);
      if (c1 > org || org > c2) return fail(`Origin variant ${a.origin} is not in the coordinate range ${a.start} to ${a.end} -/+ a ${window} bp window.`);
    } else chromosome = cs;
  }
  if (c2 - c1 > 3000000) return fail(`Queried region is ${c2 - c1} base pairs. Max size is 3000000 base pairs.`);
  // the association file, line by line (line 1 is the header)
  const assoc = [], pval = new Map(); // assoc: [pos, p] in file order; pval: pos → last p string
  lines.forEach((line, n) => {
    if (n === lines.length - 1 && line === "") return; // no line after the final newline
    const col = line.trim().split(/\s+/).filter(Boolean);
    if (col.length !== header.length) { o.warning = `Line ${n + 1} of association data file has a different number of elements than the header`; return; }
    if (stripChr(col[ci]) !== chromosome || !PYINT.test(col[pi])) return;
    const p = +col[pi].replace(/_/g, "");
    if (p < c1 || p > c2 || !PYFLOAT.test(col[vi])) return;
    assoc.push([p, pyFloat(col[vi])]); pval.set(p, col[vi]);
  });
  if (!assoc.length) return fail("There are no variants in the association file with genomic coordinates inside the plotting window.");
  for (const p of String(pop).split("+")) if (!POPS.includes(p)) return fail(POP_ERROR(p));

  const all = await ctx.windowRows(build, chromosome, c1, c2), at = p => all.filter(x => x.pos <= p && x.pos + x.ref.length - 1 >= p);
  const { ids: samples, haploidMales } = await ctx.samples(build, chromosome), panel = await ctx.popIds(pop), hbytes = Math.ceil((2 * samples.length) / 8);
  const polymorphic = row => { // get_query_variant_c / LDassoc: '0' and '1' characters of the population's GT strings
    const bits = decodeHaps(row.haps, hbytes); let n0 = 0, n1 = 0;
    samples.forEach((s, i) => {
      if (!panel.has(s) || (row.missing && (row.missing[(2 * i) >> 3] >> ((2 * i) & 7)) & 1)) return;
      for (const k of row.male_haploid && haploidMales.has(i) ? [2 * i] : [2 * i, 2 * i + 1]) (bits[k >> 3] >> (k & 7)) & 1 ? n1++ : n0++;
    });
    return n0 > 0 && n1 > 0;
  };
  let q;
  if (org === undefined) { // the lowest p-value position usable as origin
    let missing = 0, found = false;
    for (const [p] of [...assoc].sort((x, y) => x[1] - y[1])) {
      snp = `${BUILDS[build].vcfPrefix}${chromosome}:${p}`;
      const rows = BUILDS[build].vcfPrefix ? [] : at(p); // sic: this lookup prepends the prefix twice ("chrchr8:…"), so High Coverage never finds one
      if (!rows.length) { missing++; W(`Lowest P-value variant (${snp}) is not in 1000G reference panel, using next lowest P-value variant`); continue; }
      if (rows.length > 1) W(`Multiple variants map to lowest P-value variant (${snp}), using first variant in VCF file`);
      const r = rows[0];
      if (r.multiallelic_source) { W(`Lowest P-value variant (${snp}) is not a biallelic variant, using next lowest P-value variant`); continue; }
      if (!polymorphic(r)) { W(`Lowest P-value variant (${snp}) is monoallelic in the ${pop} population, using next lowest P-value variant`); continue; }
      org = p; q = r; found = true; break;
    }
    if (missing === assoc.length) return fail(`Association file Lowest P-value variant (${snp}) is not in 1000G reference panel`);
    if (!found) return INTERNAL; // org_coord unset: NameError
  } else {
    if (!pval.has(org)) return fail(`Association file is missing a p-value for origin variant ${snp}.`);
    const rows = at(org);
    if (!rows.length) return fail("No LDassoc variants could be processed from the association file and 1000G reference panel.");
    const last = rows.at(-1), id = vcfId(last); // get_query_variant_c checks the LAST row; LDassoc_sub uses the FIRST
    if (id !== snp && id.includes("rs")) o.warning = `Genomic position does not match RS number at 1000G position (chr${chromosome}:${last.pos} = ${id}).${o.warning ?? ""}`;
    if (last.multiallelic_source) o.error = `${snp} variant is not a biallelic.${o.warning ?? ""}`;
    q = rows[0];
  }
  // LD of the origin with every row listed under each association position (file order), Ms ≠ 0, POS in the file
  const listed = assoc.flatMap(([p]) => at(p)).filter(x => !x.multiallelic_source);
  const uniq = [...new Set([q, ...listed])], idx = new Map(uniq.map((x, i) => [x, i]));
  const prep = await ctx.prepare(build, uniq.map(x => ({ ...x, chrom: chromosome })), pop);
  const counts = ctx.kernel.counts(prep, uniq.map(() => 0), uniq.map((_, i) => i)); ctx.kernel.reset();
  const [qa0, qa1] = vcfAlleles(q.ref, q.alt), qal = { 0: qa0, 1: qa1 }, out = [];
  for (const x of listed) {
    const [n, a1, a2, both] = counts[idx.get(x)], D = both, C = a1 - both, B = a2 - both, A = n - a1 - a2 + both, N = n;
    if (!N) return INTERNAL; // ZeroDivisionError in LD_calcs (maf before the Ms test)
    const Ms = (A + C) * (B + D) * (A + B) * (C + D); if (Ms === 0) continue;
    if (!pval.has(x.pos)) continue;
    const delta = A * D - B * C, dprime = Math.abs(delta / (delta < 0 ? Math.min((A + C) * (A + B), (B + D) * (C + D)) : Math.min((A + C) * (C + D), (A + B) * (B + D))));
    const r2 = (delta * delta) / Ms, [pa0, pa1] = vcfAlleles(x.ref, x.alt), pal = { 0: pa0, 1: pa1 };
    const m = (k1, k2) => `${pal[k1[0]]}=${qal[k1[1]]},${pal[k2[0]]}=${qal[k2[1]]}`;
    const match = r2 > 0.1 ? ((A * D) / Math.max(B * C, 0.01) > 1 ? m("00", "11") : m("01", "10")) : "NA";
    out.push({ x, dist: x.pos - q.pos, dprime, r2, match, mafq: Math.min((A + B) / N, (C + D) / N), mafp: Math.min((A + C) / N, (B + D) / N), id: vcfId(x), p: pyFloat(pval.get(x.pos)) });
  }
  if (!out.length) return fail("No LDassoc variants could be processed from the association file and 1000G reference panel.");
  out.sort((u, v) => Math.abs(u.dist) - Math.abs(v.dist)).sort((u, v) => u.p - v.p);
  // annotations, as LDproxy (get_coords for rsIDs: int-normalized)
  const noRs = [...new Set(out.filter(e => !e.id.startsWith("rs")).map(e => e.x.pos))], dAt = await ctx.dbsnpAt(build, chromosome, noRs);
  const digits = e => { const t = e.id.replace(/^[rs]+|[rs]+$/g, ""); return PYINT.test(t) ? String(BigInt(t.trim().replace(/_/g, ""))) : null; };
  for (const e of out) {
    if (e.id.startsWith("rs")) e.rs = e.id;
    else { const d = dAt.get(e.x.pos); e.rs = d ? `rs${d.rsid}` : "."; e.func = d ? d.function : "."; }
  }
  // function (find_one by the int-normalized ID) and FORGEdb (find_one by RS), positionally where exact (common.rsAnnotations)
  const docPos = e => (!e.id.startsWith("rs") ? e.x.pos : e.x.rsid_source === 1 && e.id === `rs${e.x.rsid}` && e.x.rsid_pos_delta != null ? e.x.pos + e.x.rsid_pos_delta : null);
  const items = [...out.filter(e => e.id.startsWith("rs") && digits(e)).map(e => ({ id: digits(e), pos: docPos(e) })),
    ...out.filter(e => /^rs\d+$/.test(e.rs)).map(e => ({ id: String(BigInt(e.rs.slice(2))), pos: e.rs === e.id ? docPos(e) : e.x.pos }))];
  const ann = await rsAnnotations(ctx, build, chromosome, items), rg = await ctx.regulome(build, chromosome, [...new Set(out.map(e => e.x.pos))]);
  for (const e of out) if (e.id.startsWith("rs")) { const a = digits(e) && ann.get(digits(e)); e.func = a && a.func !== undefined ? a.func : "."; }
  const forgeOf = e => (/^rs\d+$/.test(e.rs) ? ann.get(String(BigInt(e.rs.slice(2))))?.forge : undefined);
  const line = (e, top) => [e.rs, `chr${chromosome}:${(top ? q : e.x).pos}`, `(${vcfAlleles((top ? q : e.x).ref, (top ? q : e.x).alt).join("/")})`,
    pyStr(pyRound(top ? e.mafq : e.mafp, 4)), String(e.dist), pyStr(pyRound(e.dprime, 4)), pyStr(pyRound(e.r2, 4)), e.match, pyStr(e.p),
    forgeOf(e)?.toString() ?? "NA", rg.get(e.x.pos) ?? ".", e.func ?? "NA"].join("\t");
  if (o.error !== undefined) return fail(o.error);
  const text = ["RS_Number\tCoord\tAlleles\tMAF\tDistance\tDprime\tR2\tCorrelated_Alleles\tP-value\tFORGEdb\tRegulomeDB\tFunction", line(out[0], true),
    ...out.slice(1).filter(e => e.rs !== snp).map(e => line(e, false))].join("\n") + "\n";
  return { text, ...(o.warning !== undefined ? { warning: o.warning } : {}) };
}
