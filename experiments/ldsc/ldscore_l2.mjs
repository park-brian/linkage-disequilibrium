// CLI for src/ldsc/l2.mjs (ldsc.py --l2): node ldsc/ldscore_l2.mjs <bfilePrefix> <cm|kb> <window> <outPrefix>
import { readFileSync, writeFileSync } from "node:fs";
import { gzipSync } from "node:zlib";
import { readBim, readFamCount, ldscoreL2, formatOutputs } from "../../src/ldsc/l2.mjs";

const [prefix, unit, win, out] = process.argv.slice(2);
if (!prefix) { console.error("usage: node ldsc/ldscore_l2.mjs <bfilePrefix> <cm|kb> <window> <outPrefix>"); process.exit(1); }
const t0 = performance.now();
const bim = readBim(readFileSync(`${prefix}.bim`, "utf8")), n = readFamCount(readFileSync(`${prefix}.fam`, "utf8"));
const res = ldscoreL2(new Uint8Array(readFileSync(`${prefix}.bed`)), bim, n, unit === "kb" ? { windKb: +win } : { windCm: +win });
const o = formatOutputs(res);
writeFileSync(`${out}.l2.ldscore.gz`, gzipSync(o.ldscore)); writeFileSync(`${out}.l2.M`, o.M); writeFileSync(`${out}.l2.M_5_50`, o.M_5_50);
console.log(JSON.stringify({ snps_read: bim.length, individuals: n, snps_kept: res.M, M_5_50: res.M_5_50, ms: Math.round(performance.now() - t0) }));
