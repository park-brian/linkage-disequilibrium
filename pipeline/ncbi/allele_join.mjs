// Step 8: allele-aware rsID ↔ 1000G join vs. (a) the release's position+type rsID labels and (b) LDlink's
// genotype-row choice, on one chromosome.
// usage: node allele_join.mjs <ncbi_chr.parquet> <genotypes dir/chrN> <out.json>
//   ncbi_chr.parquet: rsid, pos, ref, alt (VCF; comma list), p_dump (LDlink's position from the dump, may be null)
//   genotypes: <dir>/chrN.{common,low,rare}.parquet (columns pos, rsid, rsid_source, rsid_candidates, ref, alt)
// Alleles are compared in minimal form: trim the common suffix, then the common prefix (POS moves right), keeping
// at least one base on each side, so left-aligned VCF records with different anchors/extents compare equal.
import { writeFileSync } from "node:fs";
import { asyncBufferFromFile, parquetMetadataAsync, parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { useDataCwd } from "../../tools/config.mjs";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.mjs)

const [ncbiPath, genoPrefix, outPath] = process.argv.slice(2);
async function read(path, columns) {
  const file = await asyncBufferFromFile(path), metadata = await parquetMetadataAsync(file);
  const out = Object.fromEntries(columns.map(c => [c, []]));
  await parquetRead({ file, metadata, compressors, columns, rowFormat: "array",
    onComplete: rows => { for (const r of rows) columns.forEach((c, i) => out[c].push(r[i])); } });
  return out;
}
export function minimal(pos, ref, alt) {
  let r = ref, a = alt, p = pos;
  while (r.length > 1 && a.length > 1 && r.at(-1) === a.at(-1)) { r = r.slice(0, -1); a = a.slice(0, -1); }
  while (r.length > 1 && a.length > 1 && r[0] === a[0]) { r = r.slice(1); a = a.slice(1); p++; }
  return `${p}:${r}:${a}`;
}
const num = v => (typeof v === "bigint" ? Number(v) : v);

const t0 = performance.now();
const N = await read(ncbiPath, ["rsid", "pos", "ref", "alt", "p_dump"]);
const G = { pos: [], rsid: [], rsid_source: [], rsid_candidates: [], ref: [], alt: [], vcf_order: [], tier: [], row: [] };
for (const tier of ["common", "low", "rare"]) {
  const g = await read(`${genoPrefix}.${tier}.parquet`, ["pos", "rsid", "rsid_source", "rsid_candidates", "ref", "alt", "vcf_order"]);
  for (let i = 0; i < g.pos.length; i++) { for (const k of Object.keys(g)) G[k].push(num(g[k][i])); G.tier.push(tier); G.row.push(i); }
}
// genotype rows in exact source-VCF order (vcf_order, PARQUET_LAYOUT.md §16)
const order = G.pos.map((_, i) => i).sort((a, b) => G.vcf_order[a] - G.vcf_order[b]);
const byAllele = new Map(); // minimal key → genotype row indexes
const posList = order.map(i => G.pos[i]);
for (const i of order) { const k = minimal(G.pos[i], G.ref[i], G.alt[i]); (byAllele.get(k) ?? byAllele.set(k, []).get(k)).push(i); }
const maxRefLen = G.ref.reduce((m, r) => Math.max(m, r.length), 1);
// LDlink: tabix chr:p-p returns rows overlapping p (pos <= p <= pos+len(ref)-1), in VCF order; it keeps the last
function ldlinkRows(p) {
  let lo = 0, hi = posList.length; const from = p - maxRefLen + 1;
  while (lo < hi) { const m = (lo + hi) >> 1; if (posList[m] < from) lo = m + 1; else hi = m; }
  const rows = [];
  for (let j = lo; j < posList.length && posList[j] <= p; j++) { const i = order[j]; if (G.pos[i] + G.ref[i].length - 1 >= p) rows.push(i); }
  return rows;
}

// (a) label check: for each genotype row, which rsIDs carry exactly its alleles?
const ncbiByAllele = new Map(); // only keys that some genotype row has (a full chromosome can exceed Map's 2^24 limit)
for (let i = 0; i < N.rsid.length; i++) for (const a of N.alt[i].split(",")) {
  const k = minimal(num(N.pos[i]), N.ref[i], a); if (!byAllele.has(k)) continue;
  (ncbiByAllele.get(k) ?? ncbiByAllele.set(k, []).get(k)).push(num(N.rsid[i]));
}
const lab = { rows: G.pos.length, by_source: {} };
for (let i = 0; i < G.pos.length; i++) {
  const c = ncbiByAllele.get(minimal(G.pos[i], G.ref[i], G.alt[i])) ?? [];
  const label = G.rsid[i];
  const cls = c.length === 0 ? (label == null ? "no_label_no_allele_match" : "labelled_but_no_allele_match")
    : label == null ? "allele_match_but_unlabelled"
    : c.includes(label) ? (c.length === 1 ? "agree_unique" : "agree_among_several") : (c.length === 1 ? "disagree_unique" : "disagree_several");
  const s = (lab.by_source[G.rsid_source[i]] ??= {}); s[cls] = (s[cls] ?? 0) + 1;
  lab[cls] = (lab[cls] ?? 0) + 1;
}
// (b) LD impact: for each NCBI rsID, LDlink's row (last overlapping row at the dump position) vs the allele-aware row
const imp = { rsids: N.rsid.length, no_dump_position: 0, ldlink_no_row: 0, samples: {} };
const sample = (k, o) => { const s = (imp.samples[k] ??= []); if (s.length < 20) s.push(o); };
for (let i = 0; i < N.rsid.length; i++) {
  const p = num(N.p_dump[i]); if (p == null) { imp.no_dump_position++; continue; }
  const rows = ldlinkRows(p); if (!rows.length) { imp.ldlink_no_row++; continue; }
  const pick = rows.at(-1);
  const match = new Set(); for (const a of N.alt[i].split(",")) for (const j of byAllele.get(minimal(num(N.pos[i]), N.ref[i], a)) ?? []) match.add(j);
  const cls = match.size === 0 ? "ldlink_row_alleles_do_not_match" : match.has(pick) ? (match.size === 1 ? "same_row" : "same_row_among_several")
    : match.size === 1 ? "different_row" : "different_row_several";
  imp[cls] = (imp[cls] ?? 0) + 1;
  if (cls !== "same_row" && cls !== "same_row_among_several") {
    // where LDlink's row comes from, and whether LDlink would warn: its GRCh37 check compares the VCF ID with the query
    // rsID ("Genomic position does not match RS number at 1000G position"); GRCh38 VCFs carry no IDs, so no warning
    const where = G.pos[pick] === p ? "same_pos" : "overlap_only";
    const idMatch = (G.rsid_source[pick] === 1 || G.rsid_source[pick] === 2) && G.rsid[pick] === num(N.rsid[i]);
    const s = ((imp.detail ??= {})[cls] ??= {}); const k = `${where}${idMatch ? ",vcf_id_equals_rsid" : ""}`; s[k] = (s[k] ?? 0) + 1;
  }
  if (cls !== "same_row") sample(cls, { rsid: num(N.rsid[i]), p_dump: p, ncbi: `${num(N.pos[i])}:${N.ref[i]}>${N.alt[i]}`,
    ldlink_row: `${G.pos[pick]}:${G.ref[pick]}>${G.alt[pick]} (${G.tier[pick]}, rs${G.rsid[pick]})`, overlapping_rows: rows.length,
    allele_rows: [...match].map(j => `${G.pos[j]}:${G.ref[j]}>${G.alt[j]}`) });
}
const out = { ncbi: ncbiPath, genotypes: genoPrefix, label_check: lab, ld_impact: imp, ms: Math.round(performance.now() - t0) };
writeFileSync(outPath, JSON.stringify(out, null, 2));
console.log(JSON.stringify({ label_check: lab, ld_impact: { ...imp, samples: undefined } }));
