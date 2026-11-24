// LDlink's variant resolution (server/LDcommon.py) over the core's dbSNP reads (src/core/dbsnp.js, every doc in dump
// order):
//   get_coords(db, rsid)           → db.dbsnp.find_one({"id": rsid})          → the first doc in dump order
//   get_rsnum + replace_coord_rsid → db.dbsnp.find({chromosome, position: str(pos)}) → all docs in dump order, then the
//                                    first `id == ref_id` doc (a warning when several), else the first doc
import { BUILDS, chromName, key as gkey, keyPos } from "../core/genome.js";
import { createDbsnp } from "../core/dbsnp.js";
import { groupBy } from "../core/collect.js";
import { firstToken, parseCoord } from "./input.js";

/** createResolver(source) or createResolver(root, openFile, { exists }) (createDbsnp's arguments) */
export function createResolver(root, openFile, opts) {
  const db = createDbsnp(root, openFile, opts);

  /** rsIDs → [{ input, rsid, chrom, pos, ref_rsid, found }] in input order (one entry per input, parity: first doc) */
  async function resolveRsids(inputs, build) {
    const ids = inputs.map(s => { const t = firstToken(s).replace(/^[rs]+|[rs]+$/g, "").trim();
      return /^[+]?\d+(_\d+)*$/.test(t) ? BigInt(t.replace(/[+_]/g, "")) : null; }); // rsid.strip("rs"); int()
    const docs = await db.byRsid(ids.filter(id => id !== null));
    return inputs.map((input, i) => {
      const row = ids[i] === null ? null : docs.get(String(ids[i]))?.[0]; // find_one: the first doc
      if (!row) return { input, found: false };
      const key = row[BUILDS[build].dbsnpKey];
      return { input, found: true, rsid: `rs${row.rsid}`, ref_rsid: row.ref_rsid == null ? null : `rs${row.ref_rsid}`,
        chrom: row.chrom == null ? null : chromName(Number(row.chrom)),
        pos: key == null ? null : keyPos(key) }; // null ⇔ position "NA" for this build
    });
  }

  /** "chr7:12345" → { input, rsid | null (input kept), warning? }  (replace_coord_rsid semantics) */
  async function resolveCoords(inputs, build) {
    const parsed = inputs.map(parseCoord), docs = new Map();
    await Promise.all([...groupBy(parsed.filter(p => p?.key != null), p => p.chro)].map(async ([chro, ps]) => {
      for (const [k, list] of await db.atKeys(BUILDS[build].positions, chro, ps.map(p => p.key), ["key", "rsid", "is_ref"])) docs.set(k,
        list);
    }));
    return inputs.map((input, i) => {
      const p = parsed[i];
      const list = p?.key == null ? [] : docs.get(String(p.key)) ?? [];
      if (!list.length) return { input, rsid: null };
      if (list.length === 1) return { input, rsid: `rs${list[0].rsid}` };
      const refs = list.filter(d => d.is_ref).map(d => `rs${d.rsid}`);
      if (refs.length > 1) return { input, rsid: refs[0],
        warning: `Multiple rsIDs (${refs.join(", ")}) map to genomic coordinates ${input}` };
      if (refs.length === 0) return { input, rsid: `rs${list[0].rsid}`, warning: `Multiple rsIDs () map to genomic coordinates ${input}` };
        // sic: LDlink joins the empty ref list
      return { input, rsid: refs[0] };
    });
  }
  /**
   * LDlink's whole input chain for one list (LDmatrix and friends): first token of each input; coordinates →
   * replace_coord_rsid (rsID + "Multiple rsIDs" warning, or the token unchanged); then get_coords on the rsID —
   * so a coordinate's position is its rsID's find_one position, as in LDlink. → [{ input, token, found, rsid,
   * ref_rsid, chrom, pos, warning? }]
   */
  async function resolveInputs(inputs, build) {
    const toks = inputs.map(firstToken);
    const coIdx = toks.map((t, i) => (t && t.slice(0, 2) !== "rs" ? i : -1)).filter(i => i >= 0);
    const co = await resolveCoords(coIdx.map(i => toks[i]), build);
    const asRs = [...toks], warn = new Map();
    coIdx.forEach((i, k) => { if (co[k].rsid) asRs[i] = co[k].rsid; if (co[k].warning) warn.set(i, co[k].warning); });
    const rsIdx = asRs.map((t, i) => (t && t.slice(0, 2) === "rs" ? i : -1)).filter(i => i >= 0);
    const rs = await resolveRsids(rsIdx.map(i => asRs[i]), build);
    const byIdx = new Map(rsIdx.map((i, k) => [i, rs[k]]));
    return inputs.map((input, i) => ({ found: false, rsid: null, ref_rsid: null, chrom: null, pos: null, ...byIdx.get(i), input,
      token: toks[i], ...(warn.has(i) ? { warning: warn.get(i) } : {}) }));
  }
  /** db.dbsnp.find_one({chromosome, position}) per position → Map(pos → { rsid, function } | null) */
  async function findAt(build, chrom, positions) {
    const m = await db.atKeys(BUILDS[build].positions, chrom, positions.map(p => gkey(chrom, p)), ["key", "rsid", "function"]);
    return new Map(positions.map(p => [p, m.get(String(gkey(chrom, p)))?.[0] ?? null]));
  }
  /** db.dbsnp.find_one({id}) per rsID digits → Map(id → { function } | null): the first doc, then its row at its position */
  async function findById(ids) {
    const all = await db.byRsid(ids.map(id => BigInt(id)), ["rsid", "chrom", "key37", "key38"]);
    const docs = new Map([...all].map(([k, list]) => [k, list[0]]));
    const out = new Map(ids.map(id => [String(id), null]));
    for (const [b, col] of [["grch37", "key37"], ["grch38", "key38"]]) { // the doc's function, from whichever build has its position
      const byChrom = groupBy([...docs].filter(([k, d]) => out.get(k) === null && d[col] != null).map(([k, d]) => [k, BigInt(d[col]),
        Number(d.chrom)]), x => x[2]);
      for (const [c, list] of byChrom) {
        const at = await db.atKeys(b, chromName(c), list.map(x => x[1]), ["key", "rsid", "function"]);
        for (const [k, kk] of list) { const r = (at.get(String(kk)) ?? []).find(x => String(x.rsid) === k); if (r) out.set(k,
          { function: r.function }); }
      }
    }
    return out;
  }
  /** db.dbsnp.find({chromosome, position}) per position → Map(pos → [{ rsid, function }] in dump order) */
  async function findAllAt(build, chrom, positions) {
    const m = await db.atKeys(BUILDS[build].positions, chrom, positions.map(p => gkey(chrom, p)), ["key", "rsid", "function"]);
    return new Map(positions.map(p => [p, m.get(String(gkey(chrom, p))) ?? []]));
  }
  return { resolveRsids, resolveCoords, resolveInputs, findAt, findById, findAllAt, duplicated: db.duplicated };
}
