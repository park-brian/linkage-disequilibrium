// Convert LDlink's ldsc reference panels (/data/ldscore/<group>/<chr>.l2.ldscore.gz + .l2.M[_5_50]) to Parquet,
// reproducing how ldsc reads them (ldscore/parse.py: ldscore(), M()):
//   - chromosomes 1..22 that have any <chr>.* file; always reads <chr>.l2.ldscore.gz (the .gz probe order means stray
//     uncompressed copies are ignored)
//   - concatenated in chromosome order, then a *stable* sort by (CHR, BP) — kept here as the file's row order
//   - MAF and CM columns are dropped by ldsc when both are present; we keep them as optional columns
//   - drop_duplicates(SNP, keep first) is applied by the reader (all rows are kept in the file)
//   - M_5_50 / M per chromosome go into key-value metadata (ldsc sums the first line's floats)
// Also converts w_hm3.snplist (SNP A1 A2) for munge_sumstats --merge-alleles.
// usage: node pipeline/ldscore/convert_ref.js <ldscore dir> <releaseRoot>
import { readdirSync, readFileSync, existsSync, mkdirSync, statSync, writeFileSync } from "node:fs";
import { gunzipSync, zstdCompressSync, constants } from "node:zlib";
import { ParquetWriter } from "hyparquet-writer";
import { fileWriter } from "hyparquet-writer/src/node.js";
import { useDataCwd } from "../../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

const [src, root] = process.argv.slice(2);
const ZSTD = b => zstdCompressSync(b, { params: { [constants.ZSTD_c_compressionLevel]: 3 } });
mkdirSync(`${root}/ldscore`, { recursive: true });
const report = {};

function writeTable(path, schema, columnData, kv) {
  const w = new ParquetWriter({ writer: fileWriter(path), schema, codec: "ZSTD", compressors: { ZSTD }, kvMetadata: kv });
  w.write({ columnData, rowGroupSize: 1_000_000, pageSize: 64 << 10 });
  w.finish();
}
const firstLineFloats = p => readFileSync(p, "utf8").split(/\r?\n/)[0].trim().split(/\s+/).filter(Boolean).map(Number);

const groups = readdirSync(src, { withFileTypes: true }).filter(e => e.isDirectory() && e.name !== "ldscore").map(e => e.name).sort();
for (const g of groups) {
  const dir = `${src}/${g}`, names = readdirSync(dir);
  const chrs = []; for (let c = 1; c <= 22; c++) if (names.some(n => n.startsWith(`${c}.`))) chrs.push(c);
  const rows = []; let header = null, extra = [];
  const M = {}, M550 = {};
  for (const c of chrs) {
    const f = `${dir}/${c}.l2.ldscore.gz`;
    if (!existsSync(f)) throw new Error(`${g}: chromosome ${c} has files but no .l2.ldscore.gz`);
    const lines = gunzipSync(readFileSync(f)).toString("utf8").split(/\r?\n/);
    const h = lines[0].split("\t");
    if (!header) { header = h; extra = h.filter(x => !["CHR", "SNP", "BP", "L2"].includes(x)); }
    else if (h.join("\t") !== header.join("\t")) throw new Error(`${g}: header differs on chr ${c}`);
    const ix = Object.fromEntries(h.map((x, i) => [x, i]));
    for (let i = 1; i < lines.length; i++) {
      if (!lines[i]) continue;
      const t = lines[i].split("\t");
      rows.push({ CHR: Number(t[ix.CHR]), SNP: t[ix.SNP], BP: Number(t[ix.BP]), L2: Number(t[ix.L2]),
        ...Object.fromEntries(extra.map(x => [x, Number(t[ix[x]])])), ord: rows.length });
    }
    if (existsSync(`${dir}/${c}.l2.M_5_50`)) M550[c] = firstLineFloats(`${dir}/${c}.l2.M_5_50`);
    if (existsSync(`${dir}/${c}.l2.M`)) M[c] = firstLineFloats(`${dir}/${c}.l2.M`);
  }
  // stable sort by (CHR, BP); Array.prototype.sort is stable, ord makes it explicit
  rows.sort((a, b) => a.CHR - b.CHR || a.BP - b.BP || a.ord - b.ord);
  const schema = [{ name: "root", num_children: 4 + extra.length },
    { name: "CHR", type: "INT32", repetition_type: "REQUIRED" }, { name: "SNP", type: "BYTE_ARRAY", converted_type: "UTF8",
      repetition_type: "REQUIRED" },
    { name: "BP", type: "INT32", repetition_type: "REQUIRED" }, { name: "L2", type: "DOUBLE", repetition_type: "REQUIRED" },
    ...extra.map(x => ({ name: x, type: "DOUBLE", repetition_type: "OPTIONAL" }))];
  const out = `${root}/ldscore/${g}.parquet`;
  writeTable(out, schema, [
    { name: "CHR", data: Int32Array.from(rows, r => r.CHR) }, { name: "SNP", data: rows.map(r => r.SNP) },
    { name: "BP", data: Int32Array.from(rows, r => r.BP) }, { name: "L2", data: Float64Array.from(rows, r => r.L2) },
    ...extra.map(x => ({ name: x, data: rows.map(r => r[x]) })),
  ], [
    { key: "source", value: `LDlink /data/ldscore/${g} (ldsc reference panel)` },
    { key: "chromosomes", value: chrs.join(",") },
    { key: "M_5_50", value: JSON.stringify(M550) }, { key: "M", value: JSON.stringify(M) },
    { key: "ldsc_order", value: "stable sort by (CHR, BP) of chromosome-ordered files; reader must drop_duplicates(SNP, keep first)" },
  ]);
  // count duplicate SNP ids (sharded sets: a single Set caps at ~16.7M entries)
  const shards = Array.from({ length: 64 }, () => new Set()); let dups = 0;
  for (const r of rows) { let h = 0; for (let i = 0; i < r.SNP.length; i++) h = (h * 31 + r.SNP.charCodeAt(i)) | 0;
    const s = shards[h & 63]; if (s.has(r.SNP)) dups++; else s.add(r.SNP); }
  report[g] = { rows: rows.length, duplicate_snps: dups, chromosomes: chrs.length, extra_columns: extra,
    M_5_50_chrs: Object.keys(M550).length, M_chrs: Object.keys(M).length, bytes: statSync(out).size };
  console.log(JSON.stringify({ group: g, ...report[g] }));
}
// HapMap3 SNP list for --merge-alleles
{
  const lines = readFileSync(`${src}/w_hm3.snplist`, "utf8").split(/\r?\n/).filter(Boolean), h = lines[0].split("\t");
  const t = lines.slice(1).map(l => l.split("\t"));
  const out = `${root}/ldscore/w_hm3.snplist.parquet`;
  writeTable(out, [{ name: "root", num_children: 3 }, ...h.map(x => ({ name: x, type: "BYTE_ARRAY", converted_type: "UTF8",
    repetition_type: "REQUIRED" }))],
    h.map((x, i) => ({ name: x, data: t.map(r => r[i]) })), [{ key: "source", value: "LDlink /data/ldscore/w_hm3.snplist" }]);
  report["w_hm3.snplist"] = { rows: t.length, bytes: statSync(out).size };
  console.log(JSON.stringify({ file: "w_hm3.snplist", ...report["w_hm3.snplist"] }));
}
writeFileSync(`${root}/ldscore/build_report.json`, JSON.stringify(report, null, 2));
