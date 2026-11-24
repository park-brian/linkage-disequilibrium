// The query planner over a release directory (manifest.json + Parquet files): resolves variant inputs (rsIDs, chr:pos)
// with LDlink parity and, on request, selects each variant's genotype row. Big rsID lists resolve "hybrid": a sample
// finds the regions, the genotype files' rsID labels there answer most of the list in one scan, the rest go direct —
// chosen when the manifest prices the scan below the direct lookups. It fetches no genotype rows itself: the ports read
// the rows they use (select, ctx.row), and the region view reads windows (windowRows). decisions.md: "The planner
// fetches no genotype rows for the ports".
//
//   const pl = await openPlanner(root, openFile)          // openFile(path) → AsyncBuffer
//   const r = await pl.query({ build, inputs, tiers })
//   r.variants    one entry per input: { input, rsid, chrom, pos, found, warning? } (LDlink parity)
//   r.plan        the rsID strategy and its predicted bytes
//   query({ ..., select: "parity" | "allele" }) also sets variants[i].selected = { status, row, parity? }: the genotype
//   row each variant uses — LDlink's last overlapping row, or the allele-aware match
import { parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { sortedLookup } from "../core/parquet.js";
import { createResolver } from "./resolve.js";
import { createSource, required } from "../core/source.js";
import { tierPath } from "../core/genome.js";
import { firstToken } from "./input.js";
import { candidateRows, pickAllele, posRange } from "../core/genotypes.js";
import { pickParity } from "./select.js";
import { groupBy } from "../core/collect.js";

const LABEL_COLS = ["pos", "rsid", "rsid_source", "rsid_pos_delta"];
/** the genotype columns whose per-row-group sizes the manifest carries (what the hybrid strategy prices) */
export const PRICED_COLUMNS = LABEL_COLS;

/**
 * The manifest's row groups as objects { min, max, keyIndex?, payloadOffsetIndex?, colBytes? }. Versions 2 and 3
 * store them column-wise ({ min: [...], max: [...], ..., colBytes: { pos: [...], ... } }), and v3 only for the genotype
 * and rsID index files; version 1 as objects (bucket v4/).
 */
export function expandManifest(m) {
  if ((m.version ?? 1) < 2) return m;
  const files = {};
  for (const [path, f] of Object.entries(m.files)) {
    const R = f.rowGroups;
    if (!R) { files[path] = f; continue; } // v3: no row groups where the planner prices nothing
    const cb = R.colBytes && Object.entries(R.colBytes);
    files[path] = { ...f, rowGroups: R.min.map((min, g) => ({
      // null min / max = no key statistics (v1 leaves them out)
      min: min ?? undefined, max: R.max[g] ?? undefined, keyIndex: R.keyIndex?.[g] ?? null,
      ...(R.payloadOffsetIndex ? { payloadOffsetIndex: R.payloadOffsetIndex[g] } : {}),
      ...(cb ? { colBytes: Object.fromEntries(cb.map(([c, v]) => [c, v[g]])) } : {}),
    })) };
  }
  return { ...m, files };
}
// cost-model constants (bytes), calibrated in the restructure; verify/data/planner.js checks hybrid == direct
const COST = {
  rsidPage: 8_000,         // pages read per rsID key across the 5 rsID-index columns
  clusterGap: 2_000_000,   // positions closer than this form one region
  hybridSample: 32,        // rsIDs resolved directly to discover regions
  directMax: 32,           // lists up to this size always resolve directly
};

/**
 * @param {string} root release root (path or URL prefix)
 * @param {(path: string, byteLength?: number) => Promise<AsyncBuffer>} openFile
 * @param {{ manifest?: object }} [opts] manifest: the parsed manifest.json, when the caller already has it (browser: fetch)
 */
export async function openPlanner(root, openFile, { manifest } = {}) {
  manifest = expandManifest(manifest ?? JSON.parse(new TextDecoder().decode(await (await openFile(`${root}/manifest.json`)).slice(0))));
  // one cached { file, metadata } per path: size and footer from the manifest (no size probe, one exact footer read)
  const size = rel => { const f = manifest.files[rel]; return f && { byteLength: f.byteLength, footer: f.footer.length }; };
  const { file: open } = createSource(root, openFile, { size }), need = required(open);
  const resolver = createResolver({ file: rel => (manifest.files[rel] ? open(rel) : Promise.resolve(null)) });

  const tierFiles = (build, chrom, tiers) => tiers.map(t => tierPath(build, chrom, t)).filter(p => manifest.files[p]);

  // estimated compressed bytes of `cols` for rows with lo <= pos <= hi: the footer, the touched row groups' offset
  // indexes of those columns (approximated by the payload / key offset-index length), and their pages prorated by span
  function scanBytes(path, lo, hi, cols) {
    const f = manifest.files[path];
    let b = f.footer.length;
    for (const rg of f.rowGroups) {
      if (rg.max < lo || rg.min > hi) continue;
      b += (rg.payloadOffsetIndex?.[1] ?? rg.keyIndex?.[3] ?? 0) * cols.length;
      const span = Math.max(1, rg.max - rg.min), frac = Math.min(1, (Math.min(hi, rg.max) - Math.max(lo, rg.min) + 1) / span);
      for (const c of cols) b += (rg.colBytes?.[c] ?? 0) * frac;
    }
    return b;
  }
  // Direct rsID lookups: per partition file touched, its footer; per row group touched, the offset indexes of the
  // 5 read columns (~keyIndex offset-index length x 5) plus the key column index; per key, one page of each column.
  function directRsidBytes(ids) {
    const parts = groupBy(ids.filter(id => id !== null).map(id => [Number(id >> 26n), Number(id)]), ([p]) => p);
    let b = 0;
    for (const [p, xs] of parts) {
      const f = manifest.files[`dbsnp/rsid/part-${String(p).padStart(2, "0")}.parquet`]; if (!f) continue;
      b += f.footer.length;
      const touched = new Set();
      for (const [, k] of xs) { const g = f.rowGroups.findIndex(rg => k >= rg.min && k <= rg.max); if (g >= 0) touched.add(g); }
      for (const g of touched) { const ki = f.rowGroups[g].keyIndex; if (ki) b += ki[1] + ki[3] * 5; }
      b += xs.length * COST.rsidPage;
    }
    return b;
  }
  /** positions → regions [lo, hi] of positions closer than COST.clusterGap */
  function clusters(positions) {
    const out = [];
    for (const p of [...new Set(positions)].sort((a, b) => a - b)) {
      const c = out.at(-1);
      if (c && p - c.hi <= COST.clusterGap) c.hi = p; else out.push({ lo: p, hi: p });
    }
    return out;
  }
  async function labelRows(path, lo, hi) {
    const { file, metadata } = await need(path);
    let rows = [];
    await parquetRead({ file, metadata, compressors, rowFormat: "object", columns: LABEL_COLS, utf8: false,
      filter: { pos: { $gte: lo, $lte: hi } }, onComplete: r => { rows = r; } });
    return rows;
  }

  /** rsIDs → positions with LDlink parity; large clustered lists via region scans of the genotype rsID labels */
  async function resolveRsids(rsids, build, tiers, plan, forceRsid) {
    const ids = rsids.map(s => { const t = firstToken(s).replace(/^[rs]+|[rs]+$/g, ""); return /^\d+$/.test(t) ? BigInt(t) : null; });
    const direct = async idx => {
      const r = await resolver.resolveRsids(idx.map(i => rsids[i]), build);
      return new Map(idx.map((i, k) => [i, r[k]]));
    };
    const all = rsids.map((_, i) => i);
    if (forceRsid === "direct" || (!forceRsid && rsids.length <= COST.directMax)) {
      plan.rsid = { strategy: "direct", n: rsids.length, predicted: Math.round(directRsidBytes(ids)) };
      return direct(all);
    }
    // duplicated rsIDs always go direct (their answer depends on dump order)
    const dupFile = await open("dbsnp/duplicated_rsids.parquet"); // absent: no duplicates
    const dupRows = dupFile ? (await sortedLookup({ ...dupFile, column: "rsid", keys: ids.filter(x => x !== null), columns: ["rsid"],
      compressors })).rows : [];
    const dups = new Set(dupRows.map(r => String(r.rsid)));
    const eligible = all.filter(i => ids[i] !== null && !dups.has(String(ids[i])));
    // discover regions from a sample
    const step = Math.max(1, Math.floor(eligible.length / COST.hybridSample));
    const sample = eligible.filter((_, k) => k % step === 0).slice(0, COST.hybridSample), res = await direct(sample);
    const found = [...res.values()].filter(r => r.found && r.pos != null);
    const regions = [...groupBy(found, r => r.chrom)].flatMap(([chrom, rs]) => clusters(rs.map(r => r.pos))
      .map(c => ({ chrom, lo: c.lo - COST.clusterGap / 2, hi: c.hi + COST.clusterGap / 2 })));
    const scanCost = regions.reduce((a, g) => a + tierFiles(build, g.chrom, tiers)
      .reduce((b, f) => b + scanBytes(f, g.lo, g.hi, LABEL_COLS), 0), 0);
    const sampleCost = directRsidBytes(sample.map(i => ids[i])),
      dupCost = manifest.files["dbsnp/duplicated_rsids.parquet"]?.byteLength ?? 0;
    const directCost = directRsidBytes(eligible.filter(i => !sample.includes(i)).map(i => ids[i]));
    plan.rsid = { n: rsids.length, duplicated: dups.size, sample: sample.length, regions: regions.length,
      scanPredicted: Math.round(scanCost), directPredicted: Math.round(directCost), fixedPredicted: Math.round(sampleCost + dupCost) };
    if (forceRsid === "direct" || (!forceRsid && scanCost >= directCost)) {
      plan.rsid.strategy = "direct"; plan.rsid.directPredicted = Math.round(directRsidBytes(ids) + dupCost);
      return direct(all);
    }
    plan.rsid.strategy = "hybrid";
    // scan labels: accept only labels that guarantee a current dbSNP doc for that rsID at pos + rsid_pos_delta
    // (source 1 vcf, 3 pos+type, 4 allele match; 2 is a merged VCF ID whose current doc may sit elsewhere)
    const label = new Map();
    await Promise.all(regions.flatMap(g => tierFiles(build, g.chrom, tiers).map(async f => {
      for (const r of await labelRows(f, g.lo, g.hi))
        if (r.rsid != null && (r.rsid_source === 1 || r.rsid_source === 3 || r.rsid_source === 4))
          label.set(String(r.rsid), { chrom: g.chrom, pos: r.pos + (r.rsid_pos_delta ?? 0) });
    })));
    const out = new Map(res), left = [];
    for (const i of all) {
      if (out.has(i)) continue;
      const hit = ids[i] !== null && !dups.has(String(ids[i])) ? label.get(String(ids[i])) : undefined;
      // labels are current IDs: the direct path returns ref_rsid null for those (the index stores ref_rsid only for merged docs)
      if (hit) out.set(i, { input: rsids[i], found: true, rsid: `rs${ids[i]}`, ref_rsid: null, chrom: hit.chrom, pos: hit.pos,
        via: "label" });
      else left.push(i);
    }
    plan.rsid.matchedByLabel = all.length - sample.length - left.length; plan.rsid.leftoverDirect = left.length;
    plan.rsid.leftoverPredicted = Math.round(directRsidBytes(left.map(i => ids[i])));
    for (const [i, r] of await direct(left)) out.set(i, r);
    return out;
  }

  /** LDlink's input chain (resolve.js resolveInputs): first token; coordinates → rsID (replace_coord_rsid); then every
   *  rsID → find_one position, through the hybrid/direct rsID resolution → one variant per input */
  async function resolveInputs(build, inputs, tiers, plan, forceRsid) {
    const toks = inputs.map(firstToken), asRs = [...toks], warn = new Map();
    const coIdx = toks.map((t, i) => (t && t.slice(0, 2) !== "rs" ? i : -1)).filter(i => i >= 0);
    if (coIdx.length) {
      plan.coords = { n: coIdx.length };
      const co = await resolver.resolveCoords(coIdx.map(i => toks[i]), build);
      coIdx.forEach((i, k) => { if (co[k].rsid) asRs[i] = co[k].rsid; if (co[k].warning) warn.set(i, co[k].warning); });
    }
    const rsIdx = asRs.map((t, i) => (t && t.slice(0, 2) === "rs" ? i : -1)).filter(i => i >= 0);
    const rs = await resolveRsids(rsIdx.map(i => asRs[i]), build, tiers, plan, forceRsid);
    const variants = inputs.map((input, i) => ({ found: false, input, token: toks[i], rsid: null, ref_rsid: null, chrom: null,
      pos: null }));
    rsIdx.forEach((i, k) => { variants[i] = { ...rs.get(k), input: inputs[i], token: toks[i] }; });
    for (const [i, w] of warn) variants[i].warning = w;
    return variants;
  }
  /** per variant, the row it uses from the overlap window [pos - max_ref_len + 1, pos] over all tiers (LDlink's last
   *  overlapping row, or the allele-aware match) → variants[i].selected */
  async function selectRows(build, variants, tiers, select, plan) {
    plan.select = { mode: select, windows: 0 };
    const SEL = ["pos", "rsid", "ref", "alt", "vcf_order", "multiallelic_source",
      ...(select === "allele" ? ["allele_is_label", "allele_rsids"] : [])];
    for (const [chrom, vs] of groupBy(variants.filter(v => v?.pos != null && v.chrom), v => v.chrom)) {
      const files = await Promise.all(tierFiles(build, chrom, tiers).map(async f => ({ tier: f.split(".").at(-2), ...(await need(f)) })));
      const { rows: cand, reach } = await candidateRows(files, vs.map(v => v.pos), { columns: SEL }), range = posRange(cand);
      plan.select.windows += vs.length;
      for (const v of vs) {
        const w = range(v.pos - reach + 1, v.pos);
        v.selected = select === "allele" ? pickAllele(w, v.pos, v.ref_rsid ?? v.rsid ?? "") : pickParity(w, v.pos);
      }
    }
  }
  async function query({ build, inputs = [], tiers = ["common", "low", "rare"], forceRsid, select }) {
    const plan = { build, tiers };
    const variants = await resolveInputs(build, inputs, tiers, plan, forceRsid);
    if (select) await selectRows(build, variants, tiers, select, plan);
    return { variants, plan };
  }
  // file(rel) → { file, metadata } (cached; footer length from the manifest) and the resolver, shared with callers
  return { query, manifest, file: open, resolver };
}
