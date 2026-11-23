// G3s, the parity snapshot: every LDlink comparison case (adapters/*.js cases()) through our ports (ours()), hashed and
// compared with verify/ldlink/snapshot.json. The hashes are taken right after a full G3 (compare.js) passes, so a
// matching hash means "still what LDlink returns" without running the oracle. Outputs are canonicalized (sorted keys)
// and the LDscore family's four time lines masked. Needs no Docker; LD_HTTP=<release URL> (or --http <url>) reads the
// release over HTTP with the app's read path (Range, manifest sizes, coalesce).
// usage: node verify/ldlink/snapshot.js [--update] [--http <url>] [module...]     (npm run verify:snapshot)
import { readFileSync, writeFileSync, existsSync } from "node:fs";
import { createHash } from "node:crypto";
import { nodeContext, canon, CONFIG, MODULES } from "./lib.js";

const args = process.argv.slice(2), update = args.includes("--update");
const http = args.includes("--http") ? args[args.indexOf("--http") + 1] : undefined;
const pick = args.filter((a, i) => !a.startsWith("--") && args[i - 1] !== "--http");
const FILE = CONFIG.repo("verify/ldlink/snapshot.json"), saved = existsSync(FILE) ? JSON.parse(readFileSync(FILE, "utf8")) : {};
const TIME = /^(Conversion finished at|Beginning analysis at|Analysis finished at|Total time elapsed:).*$/gm;
const json = v => JSON.stringify(v, (k, x) => (typeof x === "bigint" ? `${x}n` : typeof x === "string" ? x.replace(TIME, "$1 …") : x));
const hash = v => createHash("sha256").update(json(canon(v))).digest("hex").slice(0, 16);

const ctx = await nodeContext(undefined, { http }), out = { ...saved }, totals = {};
for (const m of pick.length ? pick : MODULES) {
  const A = await import(`./adapters/${m}.js`), t0 = performance.now(), cases = await A.cases(ctx), got = [];
  for (const c of cases) got.push([c.name, hash(await A.ours(ctx, c))]);
  const want = saved[m] ?? [], wantBy = new Map(want), differ = got.filter(([n, h]) => wantBy.get(n) !== h);
  totals[m] = { cases: got.length, same: got.length - differ.length, differ: differ.length, seconds: Math.round((performance.now() - t0) / 1000) };
  for (const [n] of differ.slice(0, 5)) console.log(JSON.stringify({ module: m, case: n, saved: wantBy.get(n) ?? null }));
  console.log(JSON.stringify({ module: m, ...totals[m] }));
  if (update) out[m] = got;
}
// one case per line, so a diff of snapshot.json names the cases that changed
const block = ([m, cs]) => ` ${JSON.stringify(m)}: [\n${cs.map(c => `  ${JSON.stringify(c)}`).join(",\n")}\n ]`;
const write = o => `{\n${Object.entries(o).map(block).join(",\n")}\n}\n`;
if (update) { writeFileSync(FILE, write(out)); console.log(`wrote ${FILE}`); }
else if (Object.values(totals).some(t => t.differ > 0)) process.exitCode = 1;
