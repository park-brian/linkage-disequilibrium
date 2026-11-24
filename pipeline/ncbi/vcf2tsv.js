// NCBI dbSNP VCF (.gz) → TSV, one row per record on NC_ chromosomes, file order kept.
// Columns: ord, chrom code (1-22, X=23, Y=24, MT=25), pos, rsid (no "rs"), ref, alt (comma list), n_alt,
//   vc (INFO VC), fx (INFO function flags in fixed order, comma list), start (LDlink-style position: POS+1 when REF
//   is longer than 1 base and every ALT shares REF's first base, i.e. the VCF carries an anchor base; else POS).
// usage: node pipeline/ncbi/vcf2tsv.js <vcf.gz> <out.tsv> [maxRecords]      (stats JSON on stdout)
// Reads and decompresses in-process and splits lines on raw buffers (readline over a pipe was ~30 MB/s on Windows).
import { createReadStream, openSync, writeSync, closeSync } from "node:fs";
import { createGunzip } from "node:zlib";
import { useDataCwd } from "../../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

const [src, dst, maxArg] = process.argv.slice(2);
const max = maxArg ? +maxArg : Infinity;
const FLAGS = ["NSF", "NSM", "NSN", "SYN", "U3", "U5", "ASS", "DSS", "INT", "R3", "R5"];
const FLAGSET = new Set(FLAGS);
const chromOf = contig => {
  if (contig.startsWith("NC_0000")) { const n = parseInt(contig.slice(7, 9), 10); return n >= 1 && n <= 24 ? n : 0; }
  return contig.startsWith("NC_012920.") ? 25 : 0;
};
const fd = openSync(dst, "w");
let out = [], outLen = 0, ord = 0, kept = 0, other = 0, done = false;
const flush = () => { if (out.length) { writeSync(fd, out.join("")); out = []; outLen = 0; } };
function handle(line) {
  if (line.charCodeAt(0) === 35) return; // '#'
  if (ord >= max) { done = true; return; }
  const o = ord++;
  const f = line.split("\t", 8), chrom = chromOf(f[0]);
  if (!chrom) { other++; return; }
  const pos = +f[1], ref = f[3], alt = f[4], info = f[7];
  let vc = "", fx = "";
  // INFO scan without splitting the long FREQ/CLN* values: walk ';'-separated keys
  let i = 0;
  while (i < info.length) {
    let j = info.indexOf(";", i); if (j < 0) j = info.length;
    const c = info.charCodeAt(i);
    if (c === 86 && info.startsWith("VC=", i)) vc = info.slice(i + 3, j); // 'V'
    else if (j - i <= 3) { const k = info.slice(i, j); if (FLAGSET.has(k)) fx += (fx ? "," : "") + k; }
    i = j + 1;
  }
  if (fx.includes(",")) fx = fx.split(",").sort((a, b) => FLAGS.indexOf(a) - FLAGS.indexOf(b)).join(",");
  let nAlt = 1, anchored = ref.length > 1 && alt.charCodeAt(0) === ref.charCodeAt(0);
  for (let k = alt.indexOf(","); k >= 0; k = alt.indexOf(",", k + 1)) { nAlt++;
    if (alt.charCodeAt(k + 1) !== ref.charCodeAt(0)) anchored = false; }
  const s = `${o}\t${chrom}\t${pos}\t${f[2].startsWith("rs") ? f[2].slice(2) : ""}\t${ref}\t${alt}\t${nAlt}\t${vc}\t${fx}\t${anchored ? pos + 1 : pos}\n`;
  out.push(s); outLen += s.length; kept++;
  if (outLen > 8 << 20) flush();
}
const t0 = performance.now();
let rest = "";
const stream = createReadStream(src, { highWaterMark: 4 << 20 }).pipe(createGunzip({ chunkSize: 4 << 20 }));
for await (const chunk of stream) {
  const text = rest + chunk.toString("latin1");
  let start = 0;
  for (let nl = text.indexOf("\n"); nl >= 0; nl = text.indexOf("\n", start)) { handle(text.slice(start, nl)); start = nl + 1;
    if (done) break; }
  rest = text.slice(start);
  if (done) { stream.destroy(); break; }
}
if (!done && rest) handle(rest);
flush(); closeSync(fd);
console.log(JSON.stringify({ src, records: ord, kept, other_contigs: other, s: Math.round((performance.now() - t0) / 1000) }));
