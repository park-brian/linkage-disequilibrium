// Stream LDlink's mongodump dbsnp.bson.gz (gzip on stdin) → TSV on stdout, keeping dump order.
// Columns: ord, id, ref_id, chromosome, position_grch37, position_grch38, type, function  (raw strings as stored)
// `ord` is the 0-based document order in the dump: Mongo returns equal-key matches in this order after a
// mongorestore, which LDlink's find_one / find rely on (parity). Stats go to stderr every 50M docs and at the end.
// usage: aws s3 cp s3://nci-cbiit-dceg-data-nonprod/ldlink/deploy/mongodb/dbsnp.bson.gz - \
//          | node pipeline/dbsnp/bson2tsv.js > dbsnp_dump/dbsnp_dump.tsv
import { createGunzip } from "node:zlib";
import { writeFileSync } from "node:fs";
import { useDataCwd } from "../../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

const FIELDS = ["id", "ref_id", "chromosome", "position_grch37", "position_grch38", "type", "function"];
const stats = { docs: 0, nonString: {}, chromosomes: {}, missingField: {} };
const dec = new TextDecoder();
const clean = s => (s.includes("\t") || s.includes("\n") ? s.replace(/[\t\n\r]/g, " ") : s);

let pending = Buffer.alloc(0), outBuf = [], outLen = 0;
function flushOut(force) {
  if (outLen > (1 << 22) || (force && outLen)) { process.stdout.write(outBuf.join("")); outBuf = []; outLen = 0; }
}
function parseDoc(b, start, end) {
  const rec = {};
  let i = start + 4;
  while (i < end - 1) {
    const t = b[i++]; let j = i; while (b[j] !== 0) j++;
    const k = b.toString("latin1", i, j); i = j + 1;
    switch (t) {
      case 2: { const l = b.readInt32LE(i); rec[k] = dec.decode(b.subarray(i + 4, i + 3 + l)); i += 4 + l; break; }
      case 7: i += 12; break;
      case 1: case 9: case 18: rec[k] = t === 1 ? b.readDoubleLE(i) : Number(b.readBigInt64LE(i)); (stats.nonString[k] ??= 0, stats.nonString[k]++); i += 8; break;
      case 16: rec[k] = b.readInt32LE(i); (stats.nonString[k] ??= 0, stats.nonString[k]++); i += 4; break;
      case 8: rec[k] = !!b[i]; i += 1; break;
      case 10: rec[k] = null; break;
      case 3: case 4: { const l = b.readInt32LE(i); i += l; (stats.nonString[k] ??= 0, stats.nonString[k]++); break; }
      default: throw new Error(`unsupported BSON type ${t} for field ${k} at doc ${stats.docs}`);
    }
  }
  return rec;
}
const gunzip = process.stdin.pipe(createGunzip({ chunkSize: 1 << 20 }));
gunzip.on("data", chunk => {
  const b = pending.length ? Buffer.concat([pending, chunk]) : chunk;
  let off = 0;
  while (off + 4 <= b.length) {
    const len = b.readInt32LE(off);
    if (off + len > b.length) break;
    const r = parseDoc(b, off, off + len);
    const c = r.chromosome; stats.chromosomes[c] = (stats.chromosomes[c] || 0) + 1;
    const row = [stats.docs];
    for (const f of FIELDS) {
      const v = r[f];
      if (v === undefined) { stats.missingField[f] = (stats.missingField[f] || 0) + 1; row.push(""); }
      else row.push(v === null ? "" : clean(String(v)));
    }
    const line = row.join("\t") + "\n"; outBuf.push(line); outLen += line.length;
    stats.docs++; off += len;
    if (stats.docs % 50_000_000 === 0) process.stderr.write(JSON.stringify({ progress: stats.docs, t: new Date().toISOString() }) + "\n");
  }
  pending = b.subarray(off);
  flushOut(false);
});
gunzip.on("end", () => {
  flushOut(true);
  if (pending.length) stats.trailingBytes = pending.length;
  process.stderr.write(JSON.stringify(stats) + "\n");
  writeFileSync("dbsnp_dump/bson2tsv_stats.json", JSON.stringify(stats, null, 2));
  // LDlink queries dbsnp with strings (find_one({"id": str}), position: str(pos)): a non-string value in these fields
  // would be unreachable in Mongo but resolvable here, so fail instead of silently diverging (ROADMAP B5)
  const bad = ["id", "ref_id", "chromosome", "position_grch37", "position_grch38"].filter(f => stats.nonString[f]);
  if (bad.length) { process.stderr.write(`non-string values in queried fields: ${bad.map(f => `${f}=${stats.nonString[f]}`).join(", ")}\n`); process.exitCode = 1; }
});
