// Stream any LDlink mongodump collection (.bson.gz on stdin) → TSV on stdout, keeping dump order.
// usage: aws s3 cp s3://.../<coll>.bson.gz - | node pipeline/dbsnp/bson_stream.js field1 field2 ... > <coll>.tsv
// Output columns: ord, then the requested fields. Values are written as stored (strings stay strings; numbers are
// printed with JS formatting). Arrays of sub-documents become '|'-joined values of their single field
// (e.g. snp_col.data: [{platform: "..."}, ...] → "A|B|C"). Tabs/newlines inside values become spaces.
// A per-field type histogram goes to stderr at the end.
import { createGunzip } from "node:zlib";

const FIELDS = process.argv.slice(2);
// Type fidelity (ROADMAP B5): LDlink's Mongo queries are type-sensitive — numeric range queries ($gte/$lte on
// position fields) never match a string, even "12345". A string that looks numeric is counted as "numeric_string",
// and for the fields named in NUMERIC_QUERY (env, comma list) it is written as "s:<value>", so the build's try_cast
// leaves it un-keyed exactly as LDlink's query leaves it unmatched.
const NUMERIC_QUERY = new Set((process.env.NUMERIC_QUERY ?? "").split(",").filter(Boolean));
const isNumericString = v => typeof v === "string" && /^\s*[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?\s*$/.test(v);
const dec = new TextDecoder();
const types = Object.fromEntries(FIELDS.map(f => [f, {}]));
let docs = 0;
const clean = s => (/[\t\n\r]/.test(s) ? s.replace(/[\t\n\r]/g, " ") : s);

function parse(b, i, end) {
  const rec = {};
  while (i < end) {
    const t = b[i++]; let j = i; while (b[j] !== 0) j++;
    const k = b.toString("latin1", i, j); i = j + 1;
    switch (t) {
      case 1: rec[k] = b.readDoubleLE(i); i += 8; break;
      case 2: { const l = b.readInt32LE(i); rec[k] = dec.decode(b.subarray(i + 4, i + 3 + l)); i += 4 + l; break; }
      case 3: case 4: { const l = b.readInt32LE(i); const sub = parse(b, i + 4, i + l - 1); rec[k] = t === 4 ? Object.values(sub) : sub; i += l; break; }
      case 7: i += 12; break;
      case 8: rec[k] = !!b[i]; i += 1; break;
      case 9: case 18: rec[k] = Number(b.readBigInt64LE(i)); i += 8; break;
      case 10: rec[k] = null; break;
      case 16: rec[k] = b.readInt32LE(i); i += 4; break;
      default: throw new Error(`unsupported BSON type ${t} (${k}) at doc ${docs}`);
    }
  }
  return rec;
}
const fmt = v => {
  if (v === undefined || v === null) return "";
  if (Array.isArray(v)) return v.map(x => (x && typeof x === "object" ? Object.values(x).join("/") : String(x))).join("|");
  if (typeof v === "object") return JSON.stringify(v);
  return String(v);
};
let pending = Buffer.alloc(0), out = [], outLen = 0;
const gz = process.stdin.pipe(createGunzip({ chunkSize: 1 << 20 }));
gz.on("data", chunk => {
  const b = pending.length ? Buffer.concat([pending, chunk]) : chunk;
  let off = 0;
  while (off + 4 <= b.length) {
    const len = b.readInt32LE(off); if (off + len > b.length) break;
    const r = parse(b, off + 4, off + len - 1);
    const row = [docs];
    for (const f of FIELDS) {
      const v = r[f], ty = v === undefined ? "missing" : v === null ? "null" : Array.isArray(v) ? "array" : isNumericString(v) ? "numeric_string" : typeof v;
      types[f][ty] = (types[f][ty] || 0) + 1;
      row.push(ty === "numeric_string" && NUMERIC_QUERY.has(f) ? `s:${clean(v)}` : clean(fmt(v)));
    }
    const line = row.join("\t") + "\n"; out.push(line); outLen += line.length;
    docs++; off += len;
  }
  pending = b.subarray(off);
  if (outLen > 1 << 22) { process.stdout.write(out.join("")); out = []; outLen = 0; }
});
gz.on("end", () => {
  if (outLen) process.stdout.write(out.join(""));
  process.stderr.write(JSON.stringify({ docs, trailingBytes: pending.length, types }) + "\n");
});
// truncated inputs (e.g. a ranged sample of the dump): flush what was parsed, report, and exit non-zero
gz.on("error", e => {
  if (outLen) process.stdout.write(out.join(""));
  process.stderr.write(JSON.stringify({ docs, error: e.message, types }) + "\n");
  process.exitCode = 2;
});
