// Compare the JS LD score regression with the reference ldsc (full-precision dumps from ref_dump.py), using the
// Python-munged .sumstats.gz files and the Parquet reference panels. usage:
//   node verify/ldsc/regress.js ldsc_runs/h2rg release_v4/ldscore   (data-relative)
import { readFileSync, writeFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { asyncBufferFromFile, parquetMetadataAsync, parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { estimateH2, estimateRg, dedupPanel, h2ObsToLiab } from "../../src/ldsc/regress.js";
import { useDataCwd } from "../../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

const [runDir, panelDir] = process.argv.slice(2);
async function panel(group) {
  const file = await asyncBufferFromFile(`${panelDir}/${group}.parquet`), metadata = await parquetMetadataAsync(file);
  const kv = Object.fromEntries(metadata.key_value_metadata.map(x => [x.key, x.value]));
  const snps = [], l2 = [];
  await parquetRead({ file, metadata, compressors, columns: ["SNP", "L2"], onChunk: () => {}, onComplete: rows => { for (const r of rows) { snps.push(r[0]); l2.push(r[1]); } } });
  const M550 = JSON.parse(kv.M_5_50), M = Object.values(M550).reduce((a, v) => a + v[0], 0); // np.sum over chromosomes
  return { ...dedupPanel(snps, l2), M };
}
const read = f => gunzipSync(readFileSync(`${runDir}/${f}`)).toString("utf8");
const ref = f => JSON.parse(readFileSync(`${runDir}/${f}`, "utf8"));
const rel = (a, b) => (a === b ? 0 : Math.abs(a - b) / Math.max(Math.abs(b), 1e-300));
const results = [];
function compare(label, got, exp, keys) {
  const diffs = Object.fromEntries(keys.map(k => [k, { js: got[k], ldsc: exp[k], rel: rel(got[k], exp[k]) }]));
  const worst = Math.max(...Object.values(diffs).map(d => d.rel));
  results.push({ label, worst_rel_diff: worst, diffs }); console.log(JSON.stringify({ label, worst_rel_diff: worst, ...Object.fromEntries(keys.map(k => [k, [got[k], exp[k]]])) }));
}
const hk = ["tot", "tot_se", "intercept", "intercept_se", "mean_chisq", "lambda_gc", "ratio", "ratio_se", "twostep_filtered"];
for (const g of ["eas", "eur"]) {
  const P = await panel(g), t0 = performance.now();
  const h = estimateH2(read("BBJ_HDLC22.sumstats.gz"), P, P.M);
  console.log(`  h2 ${g}: ${h.n_snp} SNPs, ${Math.round(performance.now() - t0)} ms`);
  compare(`h2 HDL chr22 ${g}`, h, ref(`ref_h2_${g}.json`).hsq, hk);
  if (g === "eas") {
    const L = ref("ref_h2_eas_liab.json").hsq, f = h2ObsToLiab(1, 0.3, 0.1);
    compare("h2 liability factor (P=0.3, K=0.1)", { liab_factor: f, liab_tot: f * h.tot }, { liab_factor: L.liab_factor, liab_tot: L.liab_factor * L.tot }, ["liab_factor", "liab_tot"]);
    const t1 = performance.now(), r = estimateRg(read("BBJ_HDLC22.sumstats.gz"), read("BBJ_LDLC22.sumstats.gz"), P, P.M);
    console.log(`  rg ${g}: ${r.n_snp} SNPs, ${Math.round(performance.now() - t1)} ms`);
    const R = ref("ref_rg_eas.json");
    compare("rg HDL/LDL chr22 eas", r, R, ["rg_ratio", "rg_se", "z", "p"]);
    compare("rg: gencov", r.gencov, R.gencov, ["tot", "tot_se", "intercept", "intercept_se", "mean_z1z2", "p", "z"]);
    compare("rg: hsq1", r.hsq1, R.hsq1, hk); compare("rg: hsq2", r.hsq2, R.hsq2, hk);
  }
}
writeFileSync(`${runDir}/verify_regress.json`, JSON.stringify(results, null, 2));
if (results.some(r => !(r.worst_rel_diff <= 1e-10))) process.exitCode = 1;
