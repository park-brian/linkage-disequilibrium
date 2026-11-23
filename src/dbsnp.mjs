// LDlink-parity variant resolution over the Parquet dbSNP indexes (built from LDlink's 2023 mongodump).
// Mirrors server/LDcommon.py:
//   get_coords(db, rsid)           → db.dbsnp.find_one({"id": rsid})          → first doc in dump order
//   get_rsnum + replace_coord_rsid → db.dbsnp.find({chromosome, position: str(pos)}) → all docs in dump order, then
//                                    first `id == ref_id` doc (warning if several), else first doc
// Index files (see build_dbsnp_v2.mjs): <root>/dbsnp/rsid/part-<rsid>>26>.parquet, <root>/dbsnp/grch<37|38>/chr<N>.parquet
import { compressors } from "hyparquet-compressors";
import { sortedLookup } from "./parquet.mjs";
import { BUILDS, chromName, firstToken, parseCoord, key as gkey } from "./genome.mjs";
import { createSource } from "./source.mjs";

/**
 * @param {(path: string) => Promise<AsyncBuffer>} openFile
 * @param {{ exists?: (relPath: string) => boolean }} [opts] exists: e.g. manifest membership. An absent file is "not
 *   found" without a request; with no exists(), only ENOENT/404 errors count as absent and every other error throws.
 */
export function createResolver(root, openFile, { exists } = {}) {
  const { file: open } = createSource(root, openFile, { exists }); // → { file, metadata } or null when absent

  /** rsIDs → [{ input, rsid, chrom, pos, ref_rsid, found }] in input order (one entry per input, parity: first doc) */
  async function resolveRsids(inputs, build) {
    const ids = inputs.map(s => { const t = firstToken(s).replace(/^[rs]+|[rs]+$/g, "").trim(); return /^[+]?\d+(_\d+)*$/.test(t) ? BigInt(t.replace(/[+_]/g, "")) : null; }); // rsid.strip("rs"); int()
    const byPart = new Map();
    ids.forEach(id => { if (id !== null) { const p = Number(id >> 26n); if (!byPart.has(p)) byPart.set(p, []); byPart.get(p).push(id); } });
    const first = new Map();
    await Promise.all([...byPart].map(async ([p, keys]) => {
      const f = await open(`dbsnp/rsid/part-${String(p).padStart(2, "0")}.parquet`);
      if (!f) return; // no partition for these IDs: not in dbSNP
      const { file, metadata } = f;
      const r = await sortedLookup({ file, metadata, column: "rsid", keys, columns: ["rsid", "ref_rsid", "chrom", "key37", "key38"], compressors });
      for (const row of r.rows) { const k = String(row.rsid); if (!first.has(k)) first.set(k, row); } // rows are in file order
    }));
    return inputs.map((input, i) => {
      const row = ids[i] === null ? null : first.get(String(ids[i]));
      if (!row) return { input, found: false };
      const key = row[BUILDS[build].dbsnpKey];
      return { input, found: true, rsid: `rs${row.rsid}`, ref_rsid: row.ref_rsid == null ? null : `rs${row.ref_rsid}`,
        chrom: row.chrom == null ? null : chromName(Number(row.chrom)),
        pos: key == null ? null : Number(BigInt(key) % 1_000_000_000n) }; // null ⇔ position "NA" for this build
    });
  }

  /** "chr7:12345" → { input, rsid | null (input kept), warning? }  (replace_coord_rsid semantics) */
  async function resolveCoords(inputs, build) {
    const parsed = inputs.map(parseCoord);
    const byChrom = new Map();
    parsed.forEach(p => { if (p?.key != null) { if (!byChrom.has(p.chro)) byChrom.set(p.chro, []); byChrom.get(p.chro).push(p.key); } });
    const docs = new Map();
    await Promise.all([...byChrom].map(async ([chro, keys]) => {
      const f = await open(`dbsnp/${BUILDS[build].positions}/chr${chro}.parquet`);
      if (!f) return; // no position index for this chromosome: no rsIDs there
      const { file, metadata } = f;
      const r = await sortedLookup({ file, metadata, column: "key", keys, columns: ["key", "rsid", "is_ref"], compressors });
      for (const row of r.rows) { const k = String(row.key); if (!docs.has(k)) docs.set(k, []); docs.get(k).push(row); } // file (dump) order
    }));
    return inputs.map((input, i) => {
      const p = parsed[i];
      const list = p?.key == null ? [] : docs.get(String(p.key)) ?? [];
      if (!list.length) return { input, rsid: null };
      if (list.length === 1) return { input, rsid: `rs${list[0].rsid}` };
      const refs = list.filter(d => d.is_ref).map(d => `rs${d.rsid}`);
      if (refs.length > 1) return { input, rsid: refs[0], warning: `Multiple rsIDs (${refs.join(", ")}) map to genomic coordinates ${input}` };
      if (refs.length === 0) return { input, rsid: `rs${list[0].rsid}`, warning: `Multiple rsIDs () map to genomic coordinates ${input}` }; // sic: LDlink joins the empty ref list
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
    return inputs.map((input, i) => ({ found: false, rsid: null, ref_rsid: null, chrom: null, pos: null, ...byIdx.get(i), input, token: toks[i], ...(warn.has(i) ? { warning: warn.get(i) } : {}) }));
  }
  /** first row per key in file (= dump) order: Map(String(key) → row) */
  async function firstBy(path, column, keys, columns) {
    const m = new Map(), f = keys.length ? await open(path) : null; if (!f) return m;
    for (const r of (await sortedLookup({ file: f.file, metadata: f.metadata, column, keys, columns, compressors })).rows) { const k = String(r[column]); if (!m.has(k)) m.set(k, r); }
    return m;
  }
  /** db.dbsnp.find_one({chromosome, position}) per position → Map(pos → { rsid, function } | null) */
  async function findAt(build, chrom, positions) {
    const m = await firstBy(`dbsnp/${BUILDS[build].positions}/chr${chrom}.parquet`, "key", positions.map(p => gkey(chrom, p)), ["key", "rsid", "function"]);
    return new Map(positions.map(p => [p, m.get(String(gkey(chrom, p))) ?? null]));
  }
  /** db.dbsnp.find_one({id}) per rsID digits → Map(id → { function } | null): the first doc, then its row at its position */
  async function findById(ids) {
    const parts = new Map(); for (const id of new Set(ids)) { const p = Number(BigInt(id) >> 26n); if (!parts.has(p)) parts.set(p, []); parts.get(p).push(BigInt(id)); }
    const docs = new Map();
    for (const [p, keys] of parts) for (const [k, v] of await firstBy(`dbsnp/rsid/part-${String(p).padStart(2, "0")}.parquet`, "rsid", keys, ["rsid", "chrom", "key37", "key38"])) docs.set(k, v);
    const out = new Map(ids.map(id => [String(id), null]));
    for (const [b, col] of [["grch37", "key37"], ["grch38", "key38"]]) { // the doc's function, from whichever build has its position
      const byChrom = new Map();
      for (const [k, d] of docs) if (out.get(k) === null && d[col] != null) { const c = Number(d.chrom); if (!byChrom.has(c)) byChrom.set(c, []); byChrom.get(c).push([k, BigInt(d[col])]); }
      for (const [c, list] of byChrom) {
        const f = await open(`dbsnp/${b}/chr${chromName(c)}.parquet`); if (!f) continue;
        const rows = (await sortedLookup({ file: f.file, metadata: f.metadata, column: "key", keys: list.map(x => x[1]), columns: ["key", "rsid", "function"], compressors })).rows;
        for (const [k, kk] of list) { const r = rows.find(x => x.key === kk && String(x.rsid) === k); if (r) out.set(k, { function: r.function }); }
      }
    }
    return out;
  }
  /** db.dbsnp.find({chromosome, position}) per position → Map(pos → [{ rsid, function }] in dump order) */
  async function findAllAt(build, chrom, positions) {
    const out = new Map(positions.map(p => [p, []])), f = positions.length ? await open(`dbsnp/${BUILDS[build].positions}/chr${chrom}.parquet`) : null;
    if (f) for (const r of (await sortedLookup({ file: f.file, metadata: f.metadata, column: "key", keys: positions.map(p => gkey(chrom, p)), columns: ["key", "rsid", "function"], compressors })).rows)
      out.get(Number(BigInt(r.key) % 1_000_000_000n))?.push(r);
    return out;
  }
  /** rsIDs with more than one dbSNP doc (dbsnp/duplicated_rsids.parquet, read once): Set of digit strings */
  let dups;
  const duplicated = () => (dups ??= (async () => {
    const f = await open("dbsnp/duplicated_rsids.parquet"); if (!f) return new Set();
    const { parquetRead } = await import("hyparquet"); let rows = [];
    await parquetRead({ file: f.file, metadata: f.metadata, columns: ["rsid"], compressors, onComplete: r => { rows = r; } });
    return new Set(rows.map(r => String(r[0])));
  })());
  return { resolveRsids, resolveCoords, resolveInputs, findAt, findById, findAllAt, duplicated };
}
