// Annotation tracks with LDlink's query semantics (find_one = first row in dump order; find = all, in index order).
// Files: <root>/annotations/{forge_rsid,platforms,gtex_tissues,gtex_genes,gene_names}.parquet and
// annotations/<grch37|grch38>/{regulome,arrays,gwas,gtex_eqtl}.parquet (build_annotations.js), keyed by genome key.
import { parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { rangeScan, sortedLookup } from "./parquet.js";
import { BUILDS, key as gkey, keyPos } from "./genome.js";

/** over a source's file(rel) → { file, metadata } (e.g. planner.file) */
export function createAnnotations(file) {
  /** first row per key in file (= dump) order: Map(String(key) → row) */
  async function firstBy(path, column, keys, columns) {
    const m = new Map(); if (!keys.length) return m;
    const { file: f, metadata } = await file(path);
    for (const r of (await sortedLookup({ file: f, metadata, column, keys, columns, compressors })).rows) { const k = String(r[column]); if (!m.has(k)) m.set(k, r); }
    return m;
  }
  /** FORGEdb find_one({snp_id}) per rsID digits → Map(id → score) */
  const forge = async ids => new Map([...(await firstBy("annotations/forge_rsid.parquet", "rsid", [...new Set(ids)].map(BigInt), ["rsid", "score"]))].map(([k, r]) => [k, r.score]));
  /** FORGEdb by position: the per-build file holds each rsID once, at its first dbSNP doc's position → Map(pos → [{ rsid, score }]) */
  async function forgeAt(build, chrom, positions) {
    const out = new Map(positions.map(p => [p, []])); if (!positions.length) return out;
    const { file: f, metadata } = await file(`annotations/${BUILDS[build].positions}/forge.parquet`);
    for (const r of (await sortedLookup({ file: f, metadata, column: "key", keys: positions.map(p => gkey(chrom, p)), columns: ["key", "rsid", "score"], compressors })).rows) out.get(keyPos(r.key))?.push(r);
    return out;
  }
  /** RegulomeDB find_one({chromosome, position}) → Map(pos → score) */
  async function regulome(build, chrom, positions) {
    const m = await firstBy(`annotations/${BUILDS[build].positions}/regulome.parquet`, "key", positions.map(p => gkey(chrom, p)), ["key", "score"]);
    return new Map(positions.map(p => [p, m.get(String(gkey(chrom, p)))?.score ?? null]));
  }
  /** db.platforms, in collection order: [{ platform, code }] */
  let plat;
  const platforms = () => (plat ??= (async () => { const { file: f, metadata } = await file("annotations/platforms.parquet"); let r = [];
    await parquetRead({ file: f, metadata, compressors, columns: ["platform", "code"], rowFormat: "object", onComplete: x => { r = x; } }); return r; })());
  /** LDlink's population sample lists (annotations/panels.parquet, read once): Map(pop → [sample IDs in file order]) */
  let pan;
  const panels = () => (pan ??= (async () => { const { file: f, metadata } = await file("annotations/panels.parquet"); let r = [];
    await parquetRead({ file: f, metadata, compressors, columns: ["pop", "sample"], rowFormat: "object", onComplete: x => { r = x; } });
    const m = new Map(); for (const { pop, sample } of r) { if (!m.has(pop)) m.set(pop, []); m.get(pop).push(sample); } return m; })());
  /** db.snp_col.find({chromosome, position}) → Map("chrom:pos" → [platforms "a|b" per doc, dump order]) */
  async function arrays(build, points) {
    const { file: f, metadata } = await file(`annotations/${BUILDS[build].positions}/arrays.parquet`), out = new Map();
    const keys = points.map(([c, p]) => gkey(c, p)), names = new Map(points.map(([c, p]) => [String(gkey(c, p)), `${c}:${p}`]));
    if (keys.length) for (const r of (await sortedLookup({ file: f, metadata, column: "key", keys, columns: ["key", "platforms"], compressors })).rows) {
      const k = names.get(String(r.key)); if (!out.has(k)) out.set(k, []); out.get(k).push(r.platforms ?? "");
    }
    return out;
  }
  /** db.gwas_catalog.find({chromosome, position: {$gte: lo, $lte: hi}}) in index order (position, then insertion) */
  async function gwas(build, chrom, lo, hi) {
    const { file: f, metadata } = await file(`annotations/${BUILDS[build].positions}/gwas.parquet`);
    const rows = await rangeScan({ file: f, metadata, column: "key", lo: gkey(chrom, lo), hi: gkey(chrom, hi), compressors,
      columns: ["key", "trait", "pubmedid", "snp_id_current", "risk_allele_frequency", "or_beta", "ci_95_text", "p_value", "chromosome"] });
    return rows.map(r => ({ ...r, pos: keyPos(r.key) }));
  }
  /** the LDexpress $match → $lookup gtex_tissues → $lookup gtex_genes pipeline (both $unwind: inner joins) at one
   *  position, index order; → [{ gene_id, tissue, slope, pval_nominal, tissueSiteDetail, gene_name_v26, ncbi_id }] */
  let gt;
  const gtexRefs = () => (gt ??= (async () => {
    const read = async (path, columns) => { const { file: f, metadata } = await file(path); let r = []; await parquetRead({ file: f, metadata, compressors, columns, rowFormat: "object", onComplete: x => { r = x; } }); return r; };
    const tissues = new Map(), genes = new Map();
    for (const t of await read("annotations/gtex_tissues.parquet", ["tissueSiteDetailId", "tissueSiteDetail"])) { if (!tissues.has(t.tissueSiteDetailId)) tissues.set(t.tissueSiteDetailId, []); tissues.get(t.tissueSiteDetailId).push(t); }
    for (const g of await read("annotations/gtex_genes.parquet", ["gene_id_v26", "gene_name_v26", "ncbi_id"])) { if (!genes.has(g.gene_id_v26)) genes.set(g.gene_id_v26, []); genes.get(g.gene_id_v26).push(g); }
    return { tissues, genes };
  })());
  /** every eQTL in [lo, hi] for the tissues (a track: position, gene, tissue, slope, p), joined like LDexpress */
  async function eqtlRange(build, chrom, lo, hi, tissueIds) {
    const { file: f, metadata } = await file(`annotations/${BUILDS[build].positions}/gtex_eqtl.parquet`), { tissues, genes } = await gtexRefs(), want = new Set(tissueIds);
    const out = [];
    for (const r of await rangeScan({ file: f, metadata, column: "key", lo: gkey(chrom, lo), hi: gkey(chrom, hi), columns: ["key", "gene_id", "tissue", "slope", "pval_nominal"], compressors })) {
      if (!want.has(r.tissue)) continue;
      const g = genes.get(r.gene_id)?.[0], t = tissues.get(r.tissue)?.[0];
      out.push({ pos: keyPos(r.key), gene: g?.gene_name_v26 ?? r.gene_id, tissue: t?.tissueSiteDetail ?? r.tissue, slope: r.slope, p: r.pval_nominal });
    }
    return out;
  }
  async function eqtl(build, chrom, positions, tissueIds) {
    const { file: f, metadata } = await file(`annotations/${BUILDS[build].positions}/gtex_eqtl.parquet`), { tissues, genes } = await gtexRefs(), want = new Set(tissueIds);
    const keys = [...new Set(positions)].map(p => gkey(chrom, p)), out = new Map();
    if (keys.length) for (const r of (await sortedLookup({ file: f, metadata, column: "key", keys, columns: ["key", "gene_id", "tissue", "slope", "pval_nominal"], compressors })).rows) {
      if (!want.has(r.tissue)) continue;
      const pos = keyPos(r.key); if (!out.has(pos)) out.set(pos, []);
      for (const t of tissues.get(r.tissue) ?? []) for (const g of genes.get(r.gene_id) ?? [])
        out.get(pos).push({ gene_id: r.gene_id, slope: r.slope, pval_nominal: r.pval_nominal, tissueSiteDetail: t.tissueSiteDetail, gene_name_v26: g.gene_name_v26, ncbi_id: g.ncbi_id });
    }
    return out;
  }
  /** db.genes_name_coords.find_one({name}) → { chrom, begin, end } for the build (null fields = "NA"), or null */
  let gn;
  async function gene(name, build) {
    gn ??= (async () => { const { file: f, metadata } = await file("annotations/gene_names.parquet"); let r = [];
      await parquetRead({ file: f, metadata, compressors, rowFormat: "object", onComplete: x => { r = x; } }); const m = new Map(); for (const g of r) if (!m.has(g.name)) m.set(g.name, g); return m; })();
    const g = (await gn).get(name), b = BUILDS[build].positions;
    return g ? { chrom: g[`chromosome_${b}`], begin: g[`begin_${b}`], end: g[`end_${b}`] } : null;
  }
  /** refGene transcripts overlapping [lo, hi] (keyed by txStart; the longest human transcripts are < 2.5 Mb) */
  async function genes(build, chrom, lo, hi) {
    const { file: f, metadata } = await file(`annotations/${BUILDS[build].positions}/genes.parquet`);
    const rows = await rangeScan({ file: f, metadata, column: "key", lo: gkey(chrom, Math.max(0, lo - 2_500_000)), hi: gkey(chrom, hi), compressors,
      columns: ["key", "txStart", "txEnd", "name", "name2", "strand", "cdsStart", "cdsEnd", "exonStarts", "exonEnds"] });
    return rows.filter(r => r.txEnd >= lo);
  }
  /** recombination rate points (cM/Mb) at positions in [lo, hi] */
  async function recomb(build, chrom, lo, hi) {
    const { file: f, metadata } = await file(`annotations/${BUILDS[build].positions}/recomb.parquet`);
    return (await rangeScan({ file: f, metadata, column: "key", lo: gkey(chrom, lo), hi: gkey(chrom, hi), columns: ["key", "rate"], compressors })).map(r => ({ pos: keyPos(r.key), rate: r.rate }));
  }
  return { genes, recomb, forge, forgeAt, regulome, platforms, panels, arrays, gwas, eqtl, eqtlRange, gene, gtexTissueIds: async () => [...(await gtexRefs()).tissues.keys()] };
}
