// Window cost of each annotation track (cold, local file byte counting). usage: node bench_annotations.mjs <releaseRoot>
import { openSync, readSync, fstatSync, writeFileSync } from "node:fs";
import { parquetMetadataAsync, parquetRead, cachedAsyncBuffer } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { rangeScan } from "../src/parquet.mjs";
import { useDataCwd, repo } from "../tools/config.mjs";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.mjs)

const [root] = process.argv.slice(2);
const MAX_TRANSCRIPT = 2_500_000; // longest refGene transcript (DMD) is ~2.2 Mb: overlap scans start this far left
function counting(path) {
  const fd = openSync(path, "r"), byteLength = fstatSync(fd).size, stats = { requests: 0, bytes: 0 };
  return { stats, file: cachedAsyncBuffer({ byteLength, slice(s, e = byteLength) {
    stats.requests++; stats.bytes += e - s; const b = Buffer.allocUnsafe(e - s); readSync(fd, b, 0, e - s, s);
    return Promise.resolve(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength)); } }, { minSize: 0 }) };
}
const K = (c, p) => BigInt(c) * 1_000_000_000n + BigInt(p);
const results = [];
for (const [span, lo] of [[1e6, 25e6], [1e7, 25e6]]) {
  const hi = lo + span;
  for (const [track, columns, from] of [
    ["recomb", ["key", "rate"], lo], ["regulome", ["key", "score"], lo], ["forge", ["key", "rsid", "score"], lo],
    ["genes", ["key", "name", "name2", "txStart", "txEnd", "exonStarts", "exonEnds", "strand"], lo - MAX_TRANSCRIPT],
    ["gwas", ["key", "trait", "pubmedid", "snp_id_current", "p_value", "or_beta", "risk_allele_frequency", "ci_95_text"], lo],
    ["gtex_eqtl", ["key", "gene_id", "tissue", "slope", "pval_nominal"], lo], ["arrays", ["key", "platforms"], lo],
  ]) {
    const { stats, file } = counting(`${root}/annotations/grch37/${track}.parquet`), t0 = performance.now();
    const metadata = await parquetMetadataAsync(file, { initialFetchSize: 64 << 10 });
    let rows = await rangeScan({ file, metadata, column: "key", lo: K(22, Math.max(0, from)), hi: K(22, hi), columns, compressors });
    if (track === "genes") rows = rows.filter(r => r.txEnd >= lo);
    const o = { track, window: `${span / 1e6} Mb`, rows: rows.length, requests: stats.requests, KB: Math.round(stats.bytes / 1024), ms: Math.round(performance.now() - t0) };
    results.push(o); console.log(JSON.stringify(o));
  }
}
writeFileSync(repo("bench/results/v2_annotations_bench.json"), JSON.stringify(results, null, 2));
