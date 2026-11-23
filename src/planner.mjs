// Query planner over a release directory (manifest.json + Parquet files). Given variant inputs (rsIDs, chr:pos) and/or
// regions, it resolves variants with LDlink parity, then fetches genotype rows choosing, per position cluster,
// a region scan or page-level point reads by estimated bytes from the manifest.
//
//   const pl = await openPlanner(root, openFile)          // openFile(path) → AsyncBuffer
//   const r = await pl.query({ build, inputs, regions, tiers })
//   r.rows        genotype rows (pos, rsid, tier, haps, missing, male_haploid, multiallelic_source), sorted by pos
//   r.variants    one entry per input: { input, rsid, chrom, pos, found, warning? } (LDlink parity)
//   r.plan        decisions with predicted bytes (compare with measured bytes in verify_planner.mjs)
//   query({ ..., select: "parity" | "allele" }) also sets variants[i].selected = { status, row, parity? }: the genotype
//   row each variant uses — LDlink's last overlapping row, or the allele-aware match (select_rows.mjs)
import { parquetMetadataAsync, parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { sortedLookup, pool } from "./parquet.mjs";
import { createResolver } from "./dbsnp.mjs";
import { parseCoord, firstToken, tierPath } from "./genome.mjs";
import { candidateRows, pickParity, pickAllele } from "./genotypes.mjs";

const ROW_COLUMNS = ["pos", "rsid", "haps", "missing", "male_haploid", "multiallelic_source"];
// cost-model constants (bytes), calibrated in steps 1-2; verify_planner.mjs checks predictions against measurements
export const COST = {
  rsidPage: 8_000,         // pages read per rsID key across the 5 rsID-index columns
  pointOverhead: 8_000,    // per point read beyond the haps page: pos/other column pages and index slices
  clusterGap: 2_000_000,   // positions closer than this form one cluster
  hybridSample: 32,        // rsIDs resolved directly to discover regions
  directMax: 32,           // lists up to this size always resolve directly
};

export async function openPlanner(root, openFile) {
  const manifest = JSON.parse(new TextDecoder().decode(await (await openFile(`${root}/manifest.json`)).slice(0)));
  const cache = new Map();
  function open(path) { // footer fetched with one exact range request from the manifest
    if (!cache.has(path)) cache.set(path, (async () => {
      const file = await openFile(`${root}/${path}`), m = manifest.files[path];
      return { file, metadata: await parquetMetadataAsync(file, { initialFetchSize: m ? m.footer.length : 64 << 10 }) };
    })());
    return cache.get(path);
  }
  const resolver = createResolver(root, async p => (await open(p.slice(root.length + 1))).file, { exists: p => !!manifest.files[p] });

  const tierFiles = (build, chrom, tiers) => tiers.map(t => tierPath(build, chrom, t)).filter(p => manifest.files[p]);

  // Per-file fixed cost of touching it: footer (exact length from the manifest) + offset indexes of the requested
  // columns for each touched row group (approximated by the payload/key offset-index length x columns).
  function fixedBytes(path, lo, hi, ncols) {
    const f = manifest.files[path];
    let b = f.footer.length;
    for (const rg of f.rowGroups) if (!(rg.max < lo || rg.min > hi)) b += ((rg.payloadOffsetIndex?.[1] ?? rg.keyIndex?.[3] ?? 0) * ncols);
    return b;
  }
  // estimated compressed bytes of `cols` for rows with lo <= pos <= hi (row groups prorated by position span)
  function scanBytes(path, lo, hi, cols) {
    let b = fixedBytes(path, lo, hi, cols.length);
    for (const rg of manifest.files[path].rowGroups) {
      if (rg.max < lo || rg.min > hi) continue;
      const span = Math.max(1, rg.max - rg.min), f = Math.min(1, (Math.min(hi, rg.max) - Math.max(lo, rg.min) + 1) / span);
      for (const c of cols) b += (rg.colBytes?.[c] ?? 0) * f;
    }
    return b;
  }
  // estimated bytes for k points in [lo, hi]: expected distinct haps pages hit, P*(1 - (1 - 1/P)^k), each costing an
  // average haps page plus a share of the other columns' pages; plus the file's fixed cost
  function pointBytes(path, lo, hi, k) {
    const f = manifest.files[path];
    let pages = 0, hapsBytes = 0, otherBytes = 0;
    for (const rg of f.rowGroups) {
      if (rg.max < lo || rg.min > hi) continue;
      const span = Math.max(1, rg.max - rg.min), frac = Math.min(1, (Math.min(hi, rg.max) - Math.max(lo, rg.min) + 1) / span);
      const rgPages = rg.payloadOffsetIndex ? Math.max(1, Math.round(rg.payloadOffsetIndex[1] / 12)) : 1;
      pages += rgPages * frac; hapsBytes += (rg.colBytes?.haps ?? 0) * frac;
      otherBytes += ((rg.colBytes?.pos ?? 0) + (rg.colBytes?.rsid ?? 0) + (rg.colBytes?.multiallelic_source ?? 0) + (rg.colBytes?.male_haploid ?? 0)) * frac;
    }
    const P = Math.max(1, pages), hit = P * (1 - Math.pow(1 - 1 / P, k));
    return fixedBytes(path, lo, hi, ROW_COLUMNS.length) + hit * (hapsBytes / P) + Math.min(otherBytes, hit * COST.pointOverhead);
  }
  // estimated bytes to resolve k chr:pos inputs in [lo, hi] through the dbSNP position index (parity labels)
  function coordBytes(build, chrom, lo, hi, k) {
    const path = `dbsnp/${build === "grch37" ? "grch37" : "grch38"}/chr${chrom}.parquet`, f = manifest.files[path];
    if (!f) return 0;
    const base = BigInt(chrom === "X" ? 23 : chrom === "Y" ? 24 : +chrom) * 1000000000n, klo = Number(base + BigInt(Math.max(0, lo))), khi = Number(base + BigInt(Math.max(0, hi)));
    let pages = 0, bytes = 0, b = f.footer.length;
    for (const rg of f.rowGroups) {
      if (Number(rg.max) < klo || Number(rg.min) > khi) continue;
      const span = Math.max(1, Number(rg.max) - Number(rg.min)), frac = Math.min(1, (Math.min(khi, Number(rg.max)) - Math.max(klo, Number(rg.min)) + 1) / span);
      const rgPages = rg.keyIndex ? Math.max(1, Math.round(rg.keyIndex[3] / 12)) : 1;
      pages += rgPages * frac; bytes += rg.bytes * 3 / 5 * frac; // key, rsid, is_ref of 5 columns
      b += rg.keyIndex ? rg.keyIndex[1] + rg.keyIndex[3] * 3 : 0;  // key column index + offset indexes of 3 columns
    }
    const P = Math.max(1, pages), hit = P * (1 - Math.pow(1 - 1 / P, k));
    return b + hit * (bytes / P);
  }
  // Direct rsID lookups: per partition file touched, its footer; per row group touched, the offset indexes of the
  // 5 read columns (~keyIndex offset-index length x 5) plus the key column index; per key, one page of each column.
  function directRsidBytes(ids) {
    const parts = new Map();
    for (const id of ids) { if (id === null) continue; const p = Number(id >> 26n); if (!parts.has(p)) parts.set(p, []); parts.get(p).push(Number(id)); }
    let b = 0;
    for (const [p, keys] of parts) {
      const f = manifest.files[`dbsnp/rsid/part-${String(p).padStart(2, "0")}.parquet`]; if (!f) continue;
      b += f.footer.length;
      const touched = new Set();
      for (const k of keys) { const g = f.rowGroups.findIndex(rg => k >= rg.min && k <= rg.max); if (g >= 0) touched.add(g); }
      for (const g of touched) { const ki = f.rowGroups[g].keyIndex; if (ki) b += ki[1] + ki[3] * 5; }
      b += keys.length * COST.rsidPage;
    }
    return b;
  }
  function clusters(positions) {
    const s = [...new Set(positions)].sort((a, b) => a - b), out = [];
    for (const p of s) { const c = out.at(-1); if (c && p - c.hi <= COST.clusterGap) { c.hi = p; c.points.push(p); } else out.push({ lo: p, hi: p, points: [p] }); }
    return out;
  }

  async function readRange(path, lo, hi, columns) {
    const { file, metadata } = await open(path);
    let rows = [];
    await parquetRead({ file, metadata, compressors, rowFormat: "object", columns, utf8: false, filter: { pos: { $gte: lo, $lte: hi } }, onComplete: r => { rows = r; } }); // utf8: binary haps/missing
    return rows;
  }
  async function readPoints(path, keys, columns) {
    const { file, metadata } = await open(path);
    return (await sortedLookup({ file, metadata, column: "pos", keys, columns, compressors, payloadColumn: "haps" })).rows;
  }

  /** rsIDs → positions with LDlink parity; large clustered lists via region scans of the genotype rsid labels */
  async function resolveRsids(rsids, build, tiers, plan, forceRsid) {
    const norm = rsids.map(firstToken);
    const ids = norm.map(s => { const t = s.replace(/^[rs]+|[rs]+$/g, ""); return /^\d+$/.test(t) ? BigInt(t) : null; });
    const direct = async idx => { const r = await resolver.resolveRsids(idx.map(i => rsids[i]), build); return new Map(idx.map((i, k) => [i, r[k]])); };
    const all = rsids.map((_, i) => i);
    if (forceRsid === "direct" || (!forceRsid && rsids.length <= COST.directMax)) { plan.rsid = { strategy: "direct", n: rsids.length, predicted: Math.round(directRsidBytes(ids)) }; return direct(all); }
    // duplicated rsIDs always go direct (their answer depends on dump order)
    const dupFile = await open("dbsnp/duplicated_rsids.parquet");
    const dupRows = (await sortedLookup({ ...dupFile, column: "rsid", keys: ids.filter(x => x !== null), columns: ["rsid"], compressors })).rows;
    const dups = new Set(dupRows.map(r => String(r.rsid)));
    const eligible = all.filter(i => ids[i] !== null && !dups.has(String(ids[i])));
    // discover regions from a sample
    const step = Math.max(1, Math.floor(eligible.length / COST.hybridSample)), sample = eligible.filter((_, k) => k % step === 0).slice(0, COST.hybridSample);
    const res = await direct(sample);
    const byChrom = new Map();
    for (const r of res.values()) if (r.found && r.pos != null) { if (!byChrom.has(r.chrom)) byChrom.set(r.chrom, []); byChrom.get(r.chrom).push(r.pos); }
    const regions = [...byChrom].flatMap(([chrom, ps]) => clusters(ps).map(c => ({ chrom, lo: c.lo - COST.clusterGap / 2, hi: c.hi + COST.clusterGap / 2 })));
    const LABEL_COLS = ["pos", "rsid", "rsid_source", "rsid_pos_delta"];
    const scanCost = regions.reduce((a, g) => a + tierFiles(build, g.chrom, tiers).reduce((b, f) => b + scanBytes(f, g.lo, g.hi, LABEL_COLS), 0), 0);
    const sampleCost = directRsidBytes(sample.map(i => ids[i])), dupCost = manifest.files["dbsnp/duplicated_rsids.parquet"]?.byteLength ?? 0;
    const directCost = directRsidBytes(eligible.filter(i => !sample.includes(i)).map(i => ids[i]));
    plan.rsid = { n: rsids.length, duplicated: dups.size, sample: sample.length, regions: regions.length, scanPredicted: Math.round(scanCost), directPredicted: Math.round(directCost),
      fixedPredicted: Math.round(sampleCost + dupCost) };
    if (forceRsid === "direct" || (!forceRsid && scanCost >= directCost)) { plan.rsid.strategy = "direct"; plan.rsid.directPredicted = Math.round(directRsidBytes(ids) + dupCost); return direct(all); }
    plan.rsid.strategy = "hybrid";
    // scan labels: accept only labels that guarantee a current dbSNP doc for that rsID at pos + rsid_pos_delta
    // (source 1 vcf, 3 pos+type, 4 allele match; 2 is a merged VCF ID whose current doc may sit elsewhere)
    const label = new Map();
    await Promise.all(regions.flatMap(g => tierFiles(build, g.chrom, tiers).map(async f => {
      for (const r of await readRange(f, g.lo, g.hi, LABEL_COLS))
        if (r.rsid != null && (r.rsid_source === 1 || r.rsid_source === 3 || r.rsid_source === 4)) label.set(String(r.rsid), { chrom: g.chrom, pos: r.pos + (r.rsid_pos_delta ?? 0) });
    })));
    const out = new Map(res), left = [];
    for (const i of all) {
      if (out.has(i)) continue;
      const hit = ids[i] !== null && !dups.has(String(ids[i])) ? label.get(String(ids[i])) : undefined;
      // labels are current IDs: the direct path returns ref_rsid null for those (the index stores ref_rsid only for merged docs)
      if (hit) out.set(i, { input: rsids[i], found: true, rsid: `rs${ids[i]}`, ref_rsid: null, chrom: hit.chrom, pos: hit.pos, via: "label" });
      else left.push(i);
    }
    plan.rsid.matchedByLabel = all.length - sample.length - left.length; plan.rsid.leftoverDirect = left.length;
    plan.rsid.leftoverPredicted = Math.round(directRsidBytes(left.map(i => ids[i])));
    for (const [i, r] of await direct(left)) out.set(i, r);
    return out;
  }

  async function query({ build, inputs = [], regions = [], tiers = ["common", "low", "rare"], force, forceRsid, select }) {
    const plan = { build, tiers, fetch: [] };
    // LDlink's input chain (resolve.mjs resolveInputs): first token; coordinates → rsID (replace_coord_rsid); then
    // every rsID → find_one position — through the hybrid/direct rsID resolution below
    const toks = inputs.map(firstToken);
    const coIdx = toks.map((t, i) => (t && t.slice(0, 2) !== "rs" ? i : -1)).filter(i => i >= 0);
    const asRs = [...toks], warn = new Map();
    if (coIdx.length) {
      // predicted cost of LDlink-parity resolution through the dbSNP position index
      const byC = new Map();
      for (const i of coIdx) { const p = parseCoord(toks[i]); if (p?.key == null) continue; if (!byC.has(p.chro)) byC.set(p.chro, []); byC.get(p.chro).push(Number(p.pos)); }
      plan.coords = { n: coIdx.length, predicted: Math.round([...byC].reduce((a, [c, ps]) => a + clusters(ps).reduce((b, k) => b + coordBytes(build, c, k.lo, k.hi, k.points.length), 0), 0)) };
      const co = await resolver.resolveCoords(coIdx.map(i => toks[i]), build);
      coIdx.forEach((i, k) => { if (co[k].rsid) asRs[i] = co[k].rsid; if (co[k].warning) warn.set(i, co[k].warning); });
    }
    const rsIdx = asRs.map((t, i) => (t && t.slice(0, 2) === "rs" ? i : -1)).filter(i => i >= 0);
    const rs = await resolveRsids(rsIdx.map(i => asRs[i]), build, tiers, plan, forceRsid);
    const variants = inputs.map((input, i) => ({ found: false, input, token: toks[i], rsid: null, ref_rsid: null, chrom: null, pos: null }));
    rsIdx.forEach((i, k) => { variants[i] = { ...rs.get(k), input: inputs[i], token: toks[i] }; });
    for (const [i, w] of warn) variants[i].warning = w;
    // per chromosome: region scans + point clusters, each choosing scan vs points by estimated bytes
    const want = new Map(); // chrom → { points:Set, regions:[] }
    const get = c => { if (!want.has(c)) want.set(c, { points: new Set(), regions: [] }); return want.get(c); };
    for (const v of variants) if (v?.pos != null && v.chrom) get(v.chrom).points.add(v.pos);
    // regions: chromosome names normalized like LDlink's (chr prefix, case, X/Y, 23/24); overlapping or adjacent
    // regions merged per chromosome so no row is fetched twice; invalid regions reported, not silently dropped
    const CHR = c => { const t = String(c ?? "").trim().replace(/^chr/i, "").toUpperCase(); return t === "23" ? "X" : t === "24" ? "Y" : /^([1-9]|1\d|2[0-2]|X|Y)$/.test(t) ? t : null; };
    const merged = new Map(); plan.invalidRegions = [];
    for (const g of regions) {
      const c = CHR(g.chrom), lo = Number(g.start), hi = Number(g.end);
      if (!c || !Number.isFinite(lo) || !Number.isFinite(hi) || hi < lo) { plan.invalidRegions.push(g); continue; }
      if (!merged.has(c)) merged.set(c, []); merged.get(c).push([lo, hi]);
    }
    for (const [c, iv] of merged) {
      iv.sort((a, b) => a[0] - b[0]);
      const out = []; for (const [lo, hi] of iv) { const last = out.at(-1); if (last && lo <= last[1] + 1) last[1] = Math.max(last[1], hi); else out.push([lo, hi]); }
      get(c).regions.push(...out);
    }
    const rows = [];
    const tasks = [];
    for (const [chrom, w] of want) {
      for (const [lo, hi] of w.regions) for (const f of tierFiles(build, chrom, tiers)) {
        plan.fetch.push({ file: f, kind: "scan", lo, hi, predicted: Math.round(scanBytes(f, lo, hi, ROW_COLUMNS)) });
        tasks.push(async () => { const t = f.split(".").at(-2); for (const r of await readRange(f, lo, hi, ROW_COLUMNS)) rows.push({ ...r, tier: t }); });
      }
      const inRegion = p => w.regions.some(([lo, hi]) => p >= lo && p <= hi);
      for (const c of clusters([...w.points].filter(p => !inRegion(p)))) for (const f of tierFiles(build, chrom, tiers)) {
        const s = scanBytes(f, c.lo, c.hi, ROW_COLUMNS), pt = pointBytes(f, c.lo, c.hi, c.points.length);
        const kind = force ?? (s <= pt ? "scan" : "points");
        plan.fetch.push({ file: f, kind, lo: c.lo, hi: c.hi, points: c.points.length, predicted: Math.round(Math.min(s, pt)), scanPredicted: Math.round(s), pointsPredicted: Math.round(pt) });
        tasks.push(async () => {
          const pts = new Set(c.points);
          const got = kind === "scan" ? (await readRange(f, c.lo, c.hi, ROW_COLUMNS)).filter(r => pts.has(r.pos)) : await readPoints(f, c.points, ROW_COLUMNS);
          const t = f.split(".").at(-2); for (const r of got) rows.push({ ...r, tier: t }); // no spread: results can be large
        });
      }
    }
    await pool(tasks, 8, t => t());
    rows.sort((a, b) => a.pos - b.pos);
    if (select) { // per-variant row choice over the overlap window [pos - max_ref_len + 1, pos], all tiers
      plan.select = { mode: select, windows: 0 };
      const SEL = ["pos", "rsid", "ref", "alt", "vcf_order", "multiallelic_source", ...(select === "allele" ? ["allele_is_label", "allele_rsids"] : [])];
      const byChrom = new Map();
      for (const v of variants) if (v?.pos != null && v.chrom) { if (!byChrom.has(v.chrom)) byChrom.set(v.chrom, []); byChrom.get(v.chrom).push(v); }
      for (const [chrom, vs] of byChrom) {
        const files = await Promise.all(tierFiles(build, chrom, tiers).map(async f => ({ tier: f.split(".").at(-2), ...(await open(f)) })));
        const { rows: cand, reach } = await candidateRows(files, vs.map(v => v.pos), { columns: SEL });
        plan.select.windows += vs.length;
        for (const v of vs) {
          const w = cand.filter(r => r.pos >= v.pos - reach + 1 && r.pos <= v.pos);
          v.selected = select === "allele" ? pickAllele(w, v.pos, v.ref_rsid ?? v.rsid ?? "") : pickParity(w, v.pos);
        }
      }
    }
    plan.predicted = {
      genotypes: plan.fetch.reduce((a, f) => a + f.predicted, 0) + (plan.rsid?.strategy === "hybrid" ? plan.rsid.scanPredicted : 0),
      rsids: plan.rsid?.strategy === "hybrid" ? plan.rsid.fixedPredicted + plan.rsid.leftoverPredicted : plan.rsid?.strategy === "direct" ? (plan.rsid.directPredicted ?? plan.rsid.predicted ?? 0) : 0,
      coords: plan.coords?.predicted ?? 0,
    };
    plan.predictedBytes = plan.predicted.genotypes + plan.predicted.rsids + plan.predicted.coords;
    return { rows, variants, plan };
  }
  // file(rel) → { file, metadata } (cached; footer length from the manifest) and the resolver, shared with callers
  return { query, manifest, scanBytes, pointBytes, file: open, resolver };
}
