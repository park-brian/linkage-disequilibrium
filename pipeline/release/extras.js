// The release's non-genomic inputs, so the app needs one data URL and nothing else:
//  - annotations/panels.parquet: LDlink's population sample lists (<POP>.txt from $LD_PANELS), one row per (pop, sample) in each
//    list's file order; pops sorted. The ports take each population's IDs as a set (LDlink: list(set(ids))).
//  - examples/: the inputs LDlink ships for its example buttons: LDassoc's prostate files, LDscore's 22.bed/bim/fam,
//    and LDherit / LDcorrelation's BBJ sumstats (copied byte for byte).
//  - ldscore/ldsc_counts.parquet: per panel, the SNP count ldsc logs ("Read reference panel LD Scores for N SNPs"):
//    distinct SNP IDs after its sep=r'\s+' read, where a row with an empty SNP field shifts left and its BP becomes the
//    SNP (src/ldsc/panels.js). Counted here once, so a browser never needs a set of up to 17 million IDs.
// Then rebuild the manifest (node pipeline/release/manifest.js <release>).
// usage: node pipeline/release/extras.js <release>
import { readdirSync, readFileSync, copyFileSync, mkdirSync, statSync } from "node:fs";
import { zstdCompressSync } from "node:zlib";
import { parquetWriteFile } from "hyparquet-writer";
import { asyncBufferFromFile, parquetMetadataAsync, parquetReadObjects } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { useDataCwd, PANELS, LDASSOC_EXAMPLES, data } from "../../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

const [root] = process.argv.slice(2);
if (!root) { console.error("usage: extras.js <release>"); process.exit(2); }

const pops = readdirSync(PANELS).filter(f => f.endsWith(".txt")).map(f => f.slice(0, -4)).sort(), pop = [], sample = [];
for (const p of pops) for (const s of readFileSync(`${PANELS}/${p}.txt`, "utf8").split(/\s+/).filter(Boolean)) { pop.push(p); sample.push(s); }
parquetWriteFile({ filename: `${root}/annotations/panels.parquet`, codec: "ZSTD", compressors: { ZSTD: zstdCompressSync }, rowGroupSize: pop.length,
  columnData: [{ name: "pop", data: pop, type: "STRING", nullable: false }, { name: "sample", data: sample, type: "STRING", nullable: false }] });

const LDSCORE = data("ldscore_src/ldscore");
const EXAMPLES = [
  ...["prostate_example_grch37.txt", "prostate_example_grch38.txt"].map(f => [`${LDASSOC_EXAMPLES}/${f}`, f]),
  ...["22.bed", "22.bim", "22.fam", "BBJ_HDLC22.txt", "BBJ_LDLC22.txt"].map(f => [`${LDSCORE}/${f}`, f]),
];
// ldsc's logged panel sizes; V8 caps a Set near 16.7 M entries, so the distinct IDs are counted in 16 shards
const counts = { pop: [], rows: [] };
for (const f of readdirSync(`${root}/ldscore`).filter(f => f.endsWith(".parquet") && !/^(w_hm3|ldsc_counts)/.test(f)).sort()) {
  const file = await asyncBufferFromFile(`${root}/ldscore/${f}`), metadata = await parquetMetadataAsync(file), shards = Array.from({ length: 16 }, () => new Set());
  let rowStart = 0;
  for (const rg of metadata.row_groups) {
    const rowEnd = rowStart + Number(rg.num_rows);
    for (const { SNP, BP } of await parquetReadObjects({ file, metadata, compressors, columns: ["SNP", "BP"], rowStart, rowEnd })) {
      const id = SNP == null || SNP === "" ? String(BP) : SNP; let h = 0;
      for (let i = 0; i < id.length; i++) h = (h * 31 + id.charCodeAt(i)) & 15;
      shards[h].add(id);
    }
    rowStart = rowEnd;
  }
  counts.pop.push(f.replace(/\.parquet$/, "")); counts.rows.push(BigInt(shards.reduce((a, s) => a + s.size, 0)));
}
parquetWriteFile({ filename: `${root}/ldscore/ldsc_counts.parquet`, codec: "ZSTD", compressors: { ZSTD: zstdCompressSync },
  columnData: [{ name: "pop", data: counts.pop, type: "STRING", nullable: false }, { name: "rows", data: counts.rows, type: "INT64", nullable: false }] });

mkdirSync(`${root}/examples`, { recursive: true });
for (const [from, name] of EXAMPLES) copyFileSync(from, `${root}/examples/${name}`);
console.log(JSON.stringify({ ldscRows: Object.fromEntries(counts.pop.map((p, i) => [p, Number(counts.rows[i])])), pops: pops.length, samples: pop.length, panelsBytes: statSync(`${root}/annotations/panels.parquet`).size,
  examples: Object.fromEntries(EXAMPLES.map(([, n]) => [n, statSync(`${root}/examples/${n}`).size])) }));
