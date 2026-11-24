// JS munge vs reference munge_sumstats.py output, then the full JS pipeline (raw sumstats → munge → h2/rg) vs ldsc.
// usage: node verify/ldsc/munge.js ldscore_src/ldscore ldsc_runs/h2rg release_v4/ldscore   (data-relative)
import { readFileSync, writeFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { asyncBufferFromFile, parquetMetadataAsync, parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { mungeSumstats } from "../../src/parity/ldsc/munge.js";
import { estimateH2, estimateRg, dedupPanel } from "../../src/parity/ldsc/regress.js";
import { useDataCwd } from "../../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

const [srcDir, runDir, panelDir] = process.argv.slice(2);
async function readParquet(path, columns) {
  const file = await asyncBufferFromFile(path), metadata = await parquetMetadataAsync(file);
  const out = Object.fromEntries(columns.map(c => [c, []]));
  await parquetRead({ file, metadata, compressors, columns, onComplete: rows => { for (const r of rows) columns.forEach((c,
    i) => out[c].push(r[i])); } });
  return { out, kv: Object.fromEntries((metadata.key_value_metadata ?? []).map(x => [x.key, x.value])) };
}
const hm3 = (await readParquet(`${panelDir}/w_hm3.snplist.parquet`, ["SNP", "A1", "A2"])).out;
const results = {};
const munged = {};
for (const f of ["BBJ_HDLC22", "BBJ_LDLC22"]) {
  const t0 = performance.now();
  const m = mungeSumstats(readFileSync(`${srcDir}/${f}.txt`, "utf8"), hm3);
  const ref = gunzipSync(readFileSync(`${runDir}/${f}.sumstats.gz`)).toString("utf8").replace(/\r\n/g, "\n");
  const a = m.text.split("\n"), b = ref.split("\n");
  let diff = 0, first = null; for (let i = 0; i < Math.max(a.length, b.length); i++) if (a[i] !== b[i]) { diff++; first ??= { line: i,
    js: a[i], ldsc: b[i] }; }
  results[f] = { identical: m.text === ref, differing_lines: diff, first_diff: first, ms: Math.round(performance.now() - t0), log: m.log };
  console.log(JSON.stringify({ file: f, identical: m.text === ref, differing_lines: diff, first_diff: first, ms: results[f].ms,
    nonmissing: m.log.nonmissing, n_min: m.log.n_min, removed_maf: m.log.removed_maf }));
  munged[f] = m.text;
}
// full JS pipeline
const { out, kv } = await readParquet(`${panelDir}/eas.parquet`, ["SNP", "L2"]);
const P = dedupPanel(out.SNP, out.L2), M = Object.values(JSON.parse(kv.M_5_50)).reduce((a, v) => a + v[0], 0);
const h = estimateH2(munged.BBJ_HDLC22, P, M), href = JSON.parse(readFileSync(`${runDir}/ref_h2_eas.json`, "utf8")).hsq;
const r = estimateRg(munged.BBJ_HDLC22, munged.BBJ_LDLC22, P, M), rref = JSON.parse(readFileSync(`${runDir}/ref_rg_eas.json`, "utf8"));
const rel = (x, y) => Math.abs(x - y) / Math.abs(y);
results.pipeline = { h2: [h.tot, href.tot, rel(h.tot, href.tot)], h2_se: [h.tot_se, href.tot_se, rel(h.tot_se, href.tot_se)],
  intercept: [h.intercept, href.intercept],
  rg: [r.rg_ratio, rref.rg_ratio, rel(r.rg_ratio, rref.rg_ratio)], rg_se: [r.rg_se, rref.rg_se, rel(r.rg_se, rref.rg_se)], p: [r.p,
    rref.p] };
console.log(JSON.stringify({ pipeline: results.pipeline }));
writeFileSync(`${runDir}/verify_munge.json`, JSON.stringify(results, null, 2));
const TOL = 1e-10, p = results.pipeline;
if (!results.BBJ_HDLC22.identical || !results.BBJ_LDLC22.identical || [p.h2, p.h2_se, p.rg, p.rg_se].some(([, ,
  r]) => !(r <= TOL))) process.exitCode = 1;
