// Stream a mongoexport --jsonArray file (stdin) → TSV on stdout (ord + requested fields), without loading it all.
// Handles pretty-printed arrays and MongoDB extended JSON scalars ({"$numberInt": "5"}, {"$numberLong": ...},
// {"$numberDouble": ...}, {"$oid": ...}). Tabs/newlines inside values become spaces.
// usage: aws s3 cp s3://.../LDLink.gwas_catalog20251024.json - | node json_array_stream.js F1 F2 ... > gwas.tsv
const FIELDS = process.argv.slice(2);
// Type fidelity (ROADMAP B5): LDlink's Mongo queries are type-sensitive — numeric range queries ($gte/$lte on
// position fields) never match a string, even "12345". A string that looks numeric is counted as "numeric_string",
// and for the fields named in NUMERIC_QUERY (env, comma list) it is written as "s:<value>", so the build's try_cast
// leaves it un-keyed exactly as LDlink's query leaves it unmatched.
const NUMERIC_QUERY = new Set((process.env.NUMERIC_QUERY ?? "").split(",").filter(Boolean));
const isNumericString = v => typeof v === "string" && /^\s*[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?\s*$/.test(v);
const types = Object.fromEntries(FIELDS.map(f => [f, {}]));
let docs = 0, depth = 0, inStr = false, esc = false, buf = [], out = [], outLen = 0;
const WRAPPED_NUMBER = ["$numberInt", "$numberLong", "$numberDouble", "$numberDecimal"];
const scalar = v => {
  if (v && typeof v === "object" && !Array.isArray(v)) {
    for (const k of WRAPPED_NUMBER) if (k in v) return { num: v[k] }; // extended JSON number: a number in Mongo
    for (const k of ["$oid", "$date"]) if (k in v) return v[k];
  }
  return v;
};
const clean = s => (/[\t\n\r]/.test(s) ? s.replace(/[\t\n\r]/g, " ") : s);
function emit(text) {
  const o = JSON.parse(text), row = [docs];
  for (const f of FIELDS) {
    const s = scalar(o[f]), wrapped = s && typeof s === "object" && "num" in s, v = wrapped ? s.num : s;
    const ty = v === undefined ? "missing" : v === null ? "null" : wrapped ? "number" : isNumericString(v) ? "numeric_string" : typeof v;
    types[f][ty] = (types[f][ty] || 0) + 1;
    row.push(v === undefined || v === null ? "" : ty === "numeric_string" && NUMERIC_QUERY.has(f) ? `s:${clean(v)}` : clean(typeof v === "object" ? JSON.stringify(v) : String(v)));
  }
  const line = row.join("\t") + "\n"; out.push(line); outLen += line.length; docs++;
  if (outLen > 1 << 22) { process.stdout.write(out.join("")); out = []; outLen = 0; }
}
process.stdin.setEncoding("utf8");
process.stdin.on("data", chunk => {
  let objStart = depth > 0 ? 0 : -1; // an object still open from the previous chunk continues at index 0
  for (let i = 0; i < chunk.length; i++) {
    const c = chunk[i];
    if (inStr) { if (esc) esc = false; else if (c === "\\") esc = true; else if (c === '"') inStr = false; continue; }
    if (c === '"') { inStr = true; continue; }
    if (c === "{") { if (depth === 0) { buf = []; objStart = i; } depth++; }
    else if (c === "}") { depth--; if (depth === 0) { buf.push(chunk.slice(objStart, i + 1)); emit(buf.join("")); buf = []; objStart = -1; } }
  }
  if (depth > 0) buf.push(chunk.slice(objStart));
});
process.stdin.on("end", () => {
  if (outLen) process.stdout.write(out.join(""));
  process.stderr.write(JSON.stringify({ docs, unterminated: depth !== 0, types }) + "\n");
});
