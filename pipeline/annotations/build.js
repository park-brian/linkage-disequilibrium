// Write the annotation Parquet files with hyparquet-writer from the DuckDB intermediates (annotations.sql).
// usage: node pipeline/annotations/build.js <releaseRoot> [name...]
// Every file keeps the intermediate's row order — (key, dump order) — and is checked for it and for its row count.
import { mkdirSync, statSync, writeFileSync, existsSync, readFileSync } from "node:fs";
import { dirname } from "node:path";
import { zstdCompressSync, constants } from "node:zlib";
import { parquetMetadataAsync, parquetRead, asyncBufferFromFile } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { ParquetWriter } from "hyparquet-writer";
import { fileWriter } from "hyparquet-writer/src/node.js";
import { useDataCwd } from "../../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

const [root, ...only] = process.argv.slice(2);
const ZSTD = input => zstdCompressSync(input, { params: { [constants.ZSTD_c_compressionLevel]: 3 } });
const CHUNK = 1_000_000;
const SRC = name => `annot_src/${name}.parquet`;

// column helpers: [name, physical type, options]
const KEY = (name = "key") => ({ name, type: "INT64", required: true, encoding: "DELTA_BINARY_PACKED", columnIndex: true,
  conv: v => BigInt(v) });
const I64 = name => ({ name, type: "INT64", conv: v => (v == null ? null : BigInt(v)) });
const I32 = name => ({ name, type: "INT32", conv: v => (v == null ? null : Number(v)) });
const F64 = name => ({ name, type: "DOUBLE", conv: v => (v == null ? null : Number(v)) });
const STR = name => ({ name, type: "BYTE_ARRAY", converted_type: "UTF8", conv: v => (v == null ? null : String(v)) });
const STRKEY = name => ({ ...STR(name), required: true, columnIndex: true, statistics: true });

const builds = [["grch37", "37"], ["grch38", "38"]];
const SPECS = [
  ...builds.flatMap(([b, n]) => [
    { name: `recomb${n}`, src: SRC(`k_recomb${n}`), out: `annotations/${b}/recomb.parquet`, page: 16 << 10, cols: [KEY(), F64("rate")] },
    { name: `regulome${n}`, src: SRC(`k_regulome${n}`), out: `annotations/${b}/regulome.parquet`, page: 16 << 10, cols: [KEY(),
      STR("score")] },
    { name: `forge${n}`, src: SRC(`k_forge${n}`), out: `annotations/${b}/forge.parquet`, page: 16 << 10, cols: [KEY(), I64("rsid"),
      I32("score")], noOrd: true },
    // genes: (key, txEnd, dump order), as annotations.sql sorts them; LDlink's refGene query has no sort (set semantics)
    { name: `genes${n}`, src: SRC(`k_genes${n}`), out: `annotations/${b}/genes.parquet`, page: 64 << 10, tieBreak: "txEnd",
      cols: [KEY(), STR("chrom"), I32("txStart"), I32("txEnd"), STR("name"), STR("name2"), STR("strand"), I32("cdsStart"), I32("cdsEnd"),
        I32("exonCount"),
             STR("exonStarts"), STR("exonEnds"), STR("exonFrames"), STR("cdsStartStat"), STR("cdsEndStat"), I32("bin"), I32("score")] },
    { name: `gwas${n}`, src: SRC(`k_gwas${n}`), out: `annotations/${b}/gwas.parquet`, page: 64 << 10,
      cols: [KEY(), ...["date_added", "pubmedid", "first_author", "date", "journal", "study", "trait", "mapped_trait", "mapped_gene",
        "strongest_risk_allele", "snps", "merged",
        "snp_id_current", "context", "risk_allele_frequency", "p_value", "pvalue_mlog", "p_value_text", "or_beta", "ci_95_text",
          "study_accession", "chromosome"].map(STR)] },
    { name: `gtex${n}`, src: SRC(`k_gtex${n}`), out: `annotations/${b}/gtex_eqtl.parquet`, page: 64 << 10,
      cols: [KEY(), STR("variant_id"), STR("gene_id"), STR("tissue"), F64("maf"), F64("pval_nominal"), F64("slope"), F64("slope_se"),
        I32("tss_distance")] },
    { name: `arrays${n}`, src: SRC(`k_arrays${n}`), out: `annotations/${b}/arrays.parquet`, page: 16 << 10, cols: [KEY(),
      STR("platforms")] },
  ]),
  { name: "forge_rsid", src: SRC("k_forge_rsid"), out: "annotations/forge_rsid.parquet", page: 8 << 10, cols: [KEY("rsid"), I32("score")] },
  { name: "gene_names", src: SRC("k_gene_names"), out: "annotations/gene_names.parquet", page: 16 << 10,
    cols: [STRKEY("name"), STR("chromosome_grch37"), I32("begin_grch37"), I32("end_grch37"), STR("chromosome_grch38"), I32("begin_grch38"),
      I32("end_grch38")], strKey: true },
  { name: "gtex_tissues", src: SRC("k_gtex_tissues"), out: "annotations/gtex_tissues.parquet", page: 64 << 10, noKey: true,
    cols: ["tissueSiteDetailId", "tissueSiteDetail", "tissueSiteDetailAbbr", "tissueSite", "colorHex", "colorRgb", "datasetId", "hasEGenes",
      "eGeneCount", "rnaSeqAndGenotypeSampleCount", "uberonId"].map(STR) },
  { name: "gtex_genes", src: SRC("k_gtex_genes"), out: "annotations/gtex_genes.parquet", page: 16 << 10, cols: [STRKEY("gene_id_v26"),
    STR("gene_name_v26"), STR("ncbi_id")], strKey: true },
  { name: "platforms", src: SRC("k_platforms"), out: "annotations/platforms.parquet", page: 16 << 10, noKey: true, cols: [STR("platform"),
    STR("code")] },
];

async function readChunk(file, metadata, cols, lo, hi) { // rows [lo, hi), placed by absolute row number
  const out = Object.fromEntries(cols.map(c => [c, new Array(hi - lo)])), seen = Object.fromEntries(cols.map(c => [c,
    new Uint8Array(hi - lo)]));
  await parquetRead({ file, metadata, columns: cols, compressors, rowStart: lo, rowEnd: hi,
    onChunk: ({ columnName, columnData, rowStart }) => {
      for (let i = 0; i < columnData.length; i++) { const r = rowStart + i; if (r < lo || r >= hi) continue;
        out[columnName][r - lo] = columnData[i]; seen[columnName][r - lo] = 1; }
    } });
  for (const c of cols) if (seen[c].indexOf(0) !== -1) throw new Error(`column ${c}: unfilled rows in [${lo}, ${hi})`);
  return out;
}

const report = {};
for (const spec of SPECS.filter(s => !only.length || only.includes(s.name))) {
  if (!existsSync(spec.src)) { console.log(JSON.stringify({ name: spec.name, skipped: "no intermediate" })); continue; }
  const t0 = performance.now();
  const file = await asyncBufferFromFile(spec.src), metadata = await parquetMetadataAsync(file), total = Number(metadata.num_rows);
  const srcCols = metadata.schema.slice(1).map(s => s.name);
  const hasOrd = srcCols.includes("ord") && !spec.noOrd;
  const outPath = `${root}/${spec.out}`; mkdirSync(dirname(outPath), { recursive: true });
  const schema = [{ name: "root", num_children: spec.cols.length }, ...spec.cols.map(c => ({ name: c.name, type: c.type,
    ...(c.converted_type ? { converted_type: c.converted_type } : {}), repetition_type: c.required ? "REQUIRED" : "OPTIONAL" }))];
  const w = new ParquetWriter({ writer: fileWriter(outPath), schema, codec: "ZSTD", compressors: { ZSTD },
    kvMetadata: [{ key: "source",
      value: spec.name.startsWith("gwas") ? "LDlink gwas_catalog backup 2025-10-24" : "LDlink mongodump 2023-04" }, { key: "order",
      value: spec.tieBreak ? `(key, ${spec.tieBreak}, dump order); LDlink's query has no sort` : "(key, dump order) for LDlink parity" }] });
  // the declared order must hold — parity depends on it — so an unsorted intermediate aborts the build (B4)
  let prev = null, prevTb = null, prevOrd = null, firstBad = null;
  for (let s = 0; s < total; s += CHUNK) {
    const e = Math.min(s + CHUNK, total), need = [...new Set([...spec.cols.map(c => c.name), ...(hasOrd ? ["ord"] : [])])];
    const c = await readChunk(file, metadata, need, s, e);
    if (!spec.noKey) {
      const k = spec.cols[0].name;
      for (let i = 0; i < e - s; i++) {
        const v = spec.strKey ? String(c[k][i]) : BigInt(c[k][i]), o = hasOrd ? BigInt(c.ord[i]) : null,
          tb = spec.tieBreak ? Number(c[spec.tieBreak][i]) : 0;
        if (prev !== null && firstBad === null && (v < prev || (v === prev && (tb < prevTb || (tb === prevTb && hasOrd && o <= prevOrd)))))
          firstBad = s + i;
        prev = v; prevTb = tb; prevOrd = o;
      }
    }
    w.write({ rowGroupSize: CHUNK, pageSize: spec.page,
      columnData: spec.cols.map(col => ({ name: col.name, data: c[col.name].map(col.conv),
        ...(col.encoding ? { encoding: col.encoding } : {}),
        columnIndex: !!col.columnIndex, statistics: !!col.columnIndex || !!col.statistics })) });
  }
  w.finish();
  if (firstBad !== null)
    throw new Error(`${outPath}: input not in the declared order at row ${firstBad} (${spec.tieBreak ? `key, ${spec.tieBreak}, ord` : "key, ord"})`);
  const written = Number((await parquetMetadataAsync(await asyncBufferFromFile(outPath))).num_rows);
  if (written !== total) throw new Error(`${outPath}: wrote ${written}, expected ${total}`);
  report[spec.out] = { rows: total, bytes: statSync(outPath).size, bytes_per_row: +(statSync(outPath).size / Math.max(1, total)).toFixed(2),
    sorted: spec.noKey ? null : firstBad === null, seconds: +((performance.now() - t0) / 1000).toFixed(1) };
  console.log(JSON.stringify({ file: spec.out, ...report[spec.out] }));
}
// a partial build (named specs) updates its entries in the existing report instead of replacing the whole report
const reportPath = `${root}/annotations/build_report.json`;
const prior = only.length && existsSync(reportPath) ? JSON.parse(readFileSync(reportPath, "utf8")) : {};
writeFileSync(reportPath, JSON.stringify({ ...prior, ...report }, null, 2));
