// The release's non-genomic inputs, so the app needs one data URL and nothing else:
//  - annotations/panels.parquet: LDlink's population sample lists (<POP>.txt from $LD_PANELS), one row per (pop, sample) in each
//    list's file order; pops sorted. The ports take each population's IDs as a set (LDlink: list(set(ids))).
//  - examples/: the inputs LDlink ships for its example buttons: LDassoc's prostate files, LDscore's 22.bed/bim/fam,
//    and LDherit / LDcorrelation's BBJ sumstats (copied byte for byte).
// Then rebuild the manifest (node pipeline/release/manifest.mjs <release>).
// usage: node pipeline/release/extras.mjs <release>
import { readdirSync, readFileSync, copyFileSync, mkdirSync, statSync } from "node:fs";
import { zstdCompressSync } from "node:zlib";
import { parquetWriteFile } from "hyparquet-writer";
import { useDataCwd, PANELS, LDASSOC_EXAMPLES, data } from "../../tools/config.mjs";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.mjs)

const [root] = process.argv.slice(2);
if (!root) { console.error("usage: extras.mjs <release>"); process.exit(2); }

const pops = readdirSync(PANELS).filter(f => f.endsWith(".txt")).map(f => f.slice(0, -4)).sort(), pop = [], sample = [];
for (const p of pops) for (const s of readFileSync(`${PANELS}/${p}.txt`, "utf8").split(/\s+/).filter(Boolean)) { pop.push(p); sample.push(s); }
parquetWriteFile({ filename: `${root}/annotations/panels.parquet`, codec: "ZSTD", compressors: { ZSTD: zstdCompressSync }, rowGroupSize: pop.length,
  columnData: [{ name: "pop", data: pop, type: "STRING", nullable: false }, { name: "sample", data: sample, type: "STRING", nullable: false }] });

const LDSCORE = data("ldscore_src/ldscore");
const EXAMPLES = [
  ...["prostate_example_grch37.txt", "prostate_example_grch38.txt"].map(f => [`${LDASSOC_EXAMPLES}/${f}`, f]),
  ...["22.bed", "22.bim", "22.fam", "BBJ_HDLC22.txt", "BBJ_LDLC22.txt"].map(f => [`${LDSCORE}/${f}`, f]),
];
mkdirSync(`${root}/examples`, { recursive: true });
for (const [from, name] of EXAMPLES) copyFileSync(from, `${root}/examples/${name}`);
console.log(JSON.stringify({ pops: pops.length, samples: pop.length, panelsBytes: statSync(`${root}/annotations/panels.parquet`).size,
  examples: Object.fromEntries(EXAMPLES.map(([, n]) => [n, statSync(`${root}/examples/${n}`).size])) }));
