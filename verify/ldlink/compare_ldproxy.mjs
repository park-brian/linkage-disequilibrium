// LDproxy: local LDlink (REST GET /LDlinkRest/ldproxy, text; the warning from its proxy<ref>.json in the shared tmp
// directory) vs modules/ldproxy.mjs, line by line. Queries: the first inputs of every LDmatrix case plus edge cases.
// usage: node verify/ldlink/compare_ldproxy.mjs <ldmatrix_cases.json> <out.json> [ldlinkBase] [chunks=NUM_SUBPROCESSES]
import { readFileSync, writeFileSync, existsSync } from "node:fs";
import { randomUUID } from "node:crypto";
import { asyncBufferFromFile, cachedAsyncBuffer } from "hyparquet";
import { loadKernel, KERNEL_WASM } from "../../src/ld/kernel.mjs";
import { createContext } from "../../src/ldlink/context.mjs";
import { ldproxy } from "../../src/ldlink/ldproxy.mjs";
import * as CONFIG from "../../tools/config.mjs";
import { useDataCwd, repo } from "../../tools/config.mjs";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.mjs)

const [casesPath, outPath, base = CONFIG.LDLINK_URL, chunks = "4"] = process.argv.slice(2);
const TMP = CONFIG.LDLINK_TMP, PANELS = CONFIG.PANELS;
const ROOT = process.env.ROOT ?? CONFIG.RELEASE;
const ctx = await createContext({ root: ROOT, open: async p => cachedAsyncBuffer(await asyncBufferFromFile(p)),
  readPanel: async p => readFileSync(`${PANELS}/${p}.txt`, "utf8").split(/\s+/).filter(Boolean), kernel: await loadKernel(readFileSync(KERNEL_WASM)) });

const qs = [];
for (const c of JSON.parse(readFileSync(casesPath, "utf8"))) for (const s of c.inputs.slice(0, 3)) qs.push({ case: c.name, build: c.build, pop: c.pop, var: s.replace(/\r$/, "") });
const edge = (case_, v, pop = "CEU", build = "grch37", extra = {}) => qs.push({ case: case_, build, pop, var: v, ...extra });
edge("r2_d=d", "rs3", "EUR", "grch37", { r2_d: "d" }); edge("window 100 kb", "rs3", "EUR", "grch37", { window: 100000 });
edge("bad r2_d", "rs3", "CEU", "grch37", { r2_d: "x" }); edge("window too large", "rs3", "CEU", "grch37", { window: 2000000 });
edge("not in dbSNP", "rs9999999999"); edge("coordinate", "chr13:32446842"); edge("bad population", "rs3", "FOO");
edge("multi-allelic", "rs62224611"); edge("chrY GRCh37", "rs2032597", "EUR"); edge("chrY GRCh38", "rs2032597", "EUR", "grch38");
edge("near chromosome start", "rs587697622", "EUR"); edge("ALL", "rs3", "ALL");

const report = { base, root: ROOT, queries: [], totals: { queries: 0, identical: 0, differ: 0 } };
for (const q of qs) {
  const ref = randomUUID(), args = { r2_d: "r2", window: 500000, ...q };
  const url = `${base}/LDlinkRest/ldproxy?var=${encodeURIComponent(q.var)}&pop=${encodeURIComponent(q.pop)}&r2_d=${args.r2_d}&window=${args.window}&genome_build=${q.build}&reference=${ref}`;
  const body = await (await fetch(url)).text(), jpath = `${TMP}/proxy${ref}.json`;
  const lj = existsSync(jpath) ? (() => { try { return JSON.parse(readFileSync(jpath, "utf8")); } catch { return {}; } })() : {};
  const l = body.trimStart().startsWith("{") ? { error: JSON.parse(body).error } : { text: body, ...(lj.warning !== undefined ? { warning: lj.warning } : {}) };
  const o = await ldproxy(ctx, { ...args, chunks: +chunks });
  const same = JSON.stringify(l) === JSON.stringify(o);
  let detail = {};
  if (!same && l.text && o.text) {
    const a = l.text.split("\n"), b = o.text.split("\n"), diff = [];
    for (let i = 0; i < Math.max(a.length, b.length); i++) if (a[i] !== b[i]) diff.push({ line: i, ldlink: a[i], ours: b[i] });
    detail = { lines: [a.length, b.length], differing_lines: diff.length, examples: diff.slice(0, 6), warning: [l.warning, o.warning] };
  } else if (!same) detail = { ldlink: l.error ?? l.text?.slice(0, 300), ours: o.error ?? o.text?.slice(0, 300) };
  report.totals.queries++; report.totals[same ? "identical" : "differ"]++;
  report.queries.push({ ...q, verdict: same ? "identical" : "differ", lines: l.text ? l.text.split("\n").length - 2 : 0, ...detail });
  console.log(JSON.stringify(report.queries.at(-1)).slice(0, 1200));
}
writeFileSync(outPath, JSON.stringify(report, null, 2));
console.log(JSON.stringify(report.totals));
