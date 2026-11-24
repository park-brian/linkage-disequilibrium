// The core's own LD API (REVIEW.md §4): numbers and typed errors, without LDlink's conventions. The parity ports
// (src/parity/ldlink) reproduce LDlink's outputs byte for byte; this is for callers who want the statistics.
//   const ld = createLD({ root, open, kernel, manifest? })          open(path, byteLength?) → AsyncBuffer
//   await ld.variant("rs3", "grch37")                   → { input, rsid, chrom, pos, ref, alt } (the row it uses)
//   await ld.pair("rs3", "rs4", { pop: "CEU" })         → { a, b, r2, dprime, n, haplotypes }
//   await ld.matrix(["rs3", "rs4", …], { pop: "CEU" })  → { variants, r2: [[…]], dprime: [[…]] }
//   await ld.proxy("rs3", { pop: "CEU", window: 500000, minR2: 0.1 }) → { query, proxies: [{ variant, r2, dprime,
//                                                        distance }] } (r² ≥ minR2, r² descending, then nearest)
// Options: pop "CEU" or "CEU+YRI" (union of panels), build "grch37" (default) | "grch38" | "grch38_high_coverage".
// r2 / dprime are doubles (NaN when undefined: a monomorphic variant in the population); haplotypes are counts keyed
// "<allele a><allele b>" over the population's haplotypes valid at both variants (n of them).
// Choices, where LDlink's are fixed by its code:
//   - an rsID's position: its dbSNP doc; docs at several positions → "ambiguous" (LDlink takes the first in dump order)
//   - its 1000G row: the row carrying the rsID's alleles (NCBI b155, core/genotypes.js pickAllele), not LDlink's last
//     overlapping row; no such row → "not_in_1000g"; the split ALTs of one multi-allelic record → "multiallelic"
//   - a coordinate (chr7:12345): the row that starts there; several → "ambiguous"
// Errors are LDError { code, message, input? }: bad_input | not_found | ambiguous | not_in_1000g | multiallelic |
// unknown_population | chromosomes_differ.
import { parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { BUILDS, chromName, keyPos, tierPath } from "./genome.js";
import { createSource, once, required } from "./source.js";
import { createDbsnp } from "./dbsnp.js";
import { createGenotypes, candidateRows, pickAllele } from "./genotypes.js";
import { groupBy } from "./collect.js";

export class LDError extends Error {
  constructor(code, message, extra = {}) { super(message); this.name = "LDError"; this.code = code; Object.assign(this, extra); }
}
const TIERS = ["common", "low", "rare"];
const COLUMNS = ["pos", "vcf_order", "rsid", "ref", "alt", "multiallelic_source", "allele_is_label", "allele_rsids",
  "haps", "missing", "male_haploid"];

export function createLD({ root, open, kernel, manifest }) {
  const source = createSource(root, open, manifest ? { exists: rel => !!manifest.files[rel], size: rel => manifest.files[rel] } : {});
  const need = required(source.file), db = createDbsnp(source), geno = createGenotypes(source.file);
  const panels = once(async () => {
    const { file, metadata } = await need("annotations/panels.parquet"); let rows = [];
    await parquetRead({ file, metadata, compressors, columns: ["pop", "sample"], rowFormat: "object", onComplete: r => { rows = r; } });
    return new Map([...groupBy(rows, r => r.pop)].map(([p, xs]) => [p, new Set(xs.map(x => x.sample))]));
  });
  const buildOf = b => { if (!BUILDS[b]) throw new LDError("bad_input", `unknown build ${b}`); return b; };

  /** an input → { input, rsid, chrom, pos } (no row yet) */
  async function locate(input, build) {
    const s = String(input).trim(), rs = /^rs(\d+)$/i.exec(s), co = /^(?:chr)?([1-9]|1\d|2[0-2]|X|Y):(\d+)$/i.exec(s);
    if (co) return { input: s, rsid: null, chrom: co[1].toUpperCase(), pos: Number(co[2]) };
    if (!rs) throw new LDError("bad_input", `${s}: neither an rsID (rs123) nor a coordinate (chr7:12345)`, { input: s });
    const docs = (await db.byRsid([BigInt(rs[1])])).get(rs[1]) ?? [], col = BUILDS[build].dbsnpKey;
    const at = [...new Set(docs.filter(d => d[col] != null).map(d => String(d[col])))];
    if (!at.length) throw new LDError("not_found", `${s} has no ${build} position in dbSNP`, { input: s });
    if (at.length > 1) throw new LDError("ambiguous", `${s} has ${at.length} positions in dbSNP`, { input: s, keys: at });
    const d = docs.find(x => String(x[col]) === at[0]);
    // a merged rsID's alleles are listed under the rsID it was merged into
    return { input: s, rsid: String(d.ref_rsid ?? d.rsid), chrom: chromName(Number(d.chrom)), pos: keyPos(at[0]) };
  }
  /** located variants (one chromosome) → their 1000G rows, in input order */
  async function rows(build, chrom, vs) {
    const files = await Promise.all(TIERS.map(async tier => ({ tier, ...(await need(tierPath(build, chrom, tier))) })));
    const { rows: cand } = await candidateRows(files, vs.map(v => v.pos), { columns: COLUMNS });
    return vs.map(v => {
      if (v.rsid === null) { // a coordinate: the row that starts there
        const at = cand.filter(r => r.pos === v.pos);
        if (!at.length) throw new LDError("not_in_1000g", `${v.input}: no 1000G variant starts there`, { input: v.input });
        if (at.length > 1) throw new LDError("ambiguous", `${v.input}: ${at.length} 1000G variants start there`, { input: v.input });
        return at[0];
      }
      const s = pickAllele(cand, v.pos, v.rsid);
      if (s.status === "ok" || s.status === "ok_split_multiallelic") return s.row;
      const at = { input: v.input };
      if (s.status === "multiallelic_in_1000g") throw new LDError("multiallelic", `${v.input} is a multi-allelic 1000G record`, at);
      if (s.status === "ambiguous") throw new LDError("ambiguous", `${v.input} matches ${s.candidates.length} 1000G rows`, at);
      throw new LDError("not_in_1000g", `${v.input}: no 1000G row carries its alleles`, at);
    });
  }
  const describe = (v, r) => ({ input: v.input, rsid: r.rsid == null ? v.rsid && `rs${v.rsid}` : `rs${r.rsid}`, chrom: v.chrom,
    pos: r.pos, ref: r.ref, alt: r.alt });
  /** inputs on one chromosome → { vs (described), rows } */
  async function resolveAll(inputs, build) {
    const vs = await Promise.all(inputs.map(x => locate(x, build)));
    const chroms = [...new Set(vs.map(v => v.chrom))];
    if (chroms.length > 1) throw new LDError("chromosomes_differ", `the variants are on chromosomes ${chroms.join(", ")}`);
    const rs = await rows(build, vs[0].chrom, vs);
    return { vs: vs.map((v, i) => describe(v, rs[i])), rows: rs, chrom: vs[0].chrom };
  }
  /** the kernel options for a chromosome and population (CEU or CEU+YRI) */
  async function population(build, chrom, pop) {
    const m = await panels(), ids = new Set();
    for (const p of String(pop).split("+")) {
      if (!m.has(p)) throw new LDError("unknown_population", `unknown population ${p}`);
      for (const s of m.get(p)) ids.add(s);
    }
    const ss = await geno.samples(build, chrom);
    const hapIndexes = ss.ids.flatMap((s, i) => (ids.has(s) ? [2 * i, 2 * i + 1] : []));
    return { samples: ss.ids.length, haploidMales: ss.haploidMales, hapIndexes };
  }
  /** rows → kernel counts and exact stats for pairs (I[k], J[k]); prepare … reset without yielding (shared arena) */
  function measure(rowsIn, opts, I, J) {
    const p = kernel.prepare(rowsIn, opts);
    try { return { counts: kernel.counts(p, I, J), ...kernel.exact(p, I, J) }; } finally { kernel.reset(); }
  }

  return {
    async variant(input, build = "grch37") {
      const { vs } = await resolveAll([input], buildOf(build));
      return vs[0];
    },
    async pair(a, b, { pop, build = "grch37" } = {}) {
      const { vs, rows: rs, chrom } = await resolveAll([a, b], buildOf(build));
      const m = measure(rs, await population(build, chrom, pop), [0], [1]), [n, a1, a2, both] = m.counts[0];
      const [A, B] = vs;
      return { a: A, b: B, r2: m.r2[0], dprime: m.dprime[0], n,
        haplotypes: { [A.ref + B.ref]: n - a1 - a2 + both, [A.ref + B.alt]: a2 - both, [A.alt + B.ref]: a1 - both,
          [A.alt + B.alt]: both } };
    },
    async matrix(inputs, { pop, build = "grch37" } = {}) {
      const { vs, rows: rs, chrom } = await resolveAll(inputs, buildOf(build)), n = rs.length, I = [], J = [];
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { I.push(i); J.push(j); }
      const m = measure(rs, await population(build, chrom, pop), I, J);
      const grid = a => Array.from({ length: n }, (_, i) => Array.from(a.subarray(i * n, (i + 1) * n)));
      return { variants: vs, r2: grid(m.r2), dprime: grid(m.dprime) };
    },
    async proxy(input, { pop, build = "grch37", window = 500_000, minR2 = 0.1 } = {}) {
      const { vs, rows: [q], chrom } = await resolveAll([input], buildOf(build));
      const self = r => r.pos === q.pos && r.vcf_order === q.vcf_order && r.ref === q.ref && r.alt === q.alt;
      const near = (await geno.windowRows(build, chrom, Math.max(1, q.pos - window), q.pos + window)).filter(r => !self(r));
      if (!near.length) return { query: vs[0], proxies: [] };
      const m = measure([q, ...near], await population(build, chrom, pop), near.map(() => 0), near.map((_, i) => i + 1));
      const variant = r => ({ rsid: r.rsid == null ? null : `rs${r.rsid}`, chrom, pos: r.pos, ref: r.ref, alt: r.alt });
      const proxies = near.map((r, i) => ({ variant: variant(r), r2: m.r2[i], dprime: m.dprime[i], distance: r.pos - q.pos }))
        .filter(x => x.r2 >= minR2).sort((x, y) => y.r2 - x.r2 || Math.abs(x.distance) - Math.abs(y.distance));
      return { query: vs[0], proxies };
    },
  };
}
