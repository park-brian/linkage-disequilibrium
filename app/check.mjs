// Browser parity check: run module ports in this page (src/ core over HTTP Range reads of the release) and the same
// requests against a local LDlink REST API (CORS: LDlink sends Access-Control-Allow-Origin *), then diff the outputs.
// Open /app/check.html?ldlink=http://127.0.0.1:8090 (LDlink defaults to that address).
import { asyncBufferFromUrl, cachedAsyncBuffer } from "hyparquet";
import { loadKernel, KERNEL_WASM } from "../src/ld/kernel.mjs";
import { createContext } from "../src/ldlink/context.mjs";
import { ldmatrix } from "../src/ldlink/ldmatrix.mjs";
import { ldpair } from "../src/ldlink/ldpair.mjs";
import { ldpop } from "../src/ldlink/ldpop.mjs";
import { ldproxy } from "../src/ldlink/ldproxy.mjs";
import { snpclip } from "../src/ldlink/snpclip.mjs";
import { snpchip } from "../src/ldlink/snpchip.mjs";
import { ldtrait } from "../src/ldlink/ldtrait.mjs";
import { ldexpress } from "../src/ldlink/ldexpress.mjs";

const LDLINK = new URLSearchParams(location.search).get("ldlink") ?? "http://127.0.0.1:8090";
const out = document.querySelector("#results"), status = document.querySelector("#status");
const VECTORS = ["rs13426489", "rs4442975", "rs13387042", "rs13412666", "rs6721996", "rs12465506", "rs12613955", "rs12620095"];
const get = async (path, q) => (await fetch(`${LDLINK}/LDlinkRest/${path}?${new URLSearchParams(q)}`)).text();
const post = async (path, body) => (await fetch(`${LDLINK}/LDlinkRest/${path}`, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) })).text();
const asText = t => (t.trimStart().startsWith("{") ? { error: JSON.parse(t).error } : { text: t });
const canon = o => (Array.isArray(o) ? o.map(canon) : o && typeof o === "object" ? Object.fromEntries(Object.keys(o).sort().map(k => [k, canon(o[k])])) : o);
const tieOnly = (a, b) => Array.isArray(a) && Array.isArray(b) && JSON.stringify({ ...a[0], haplotypes: 0 }) === JSON.stringify({ ...b[0], haplotypes: 0 });
const T = "Adipose_Subcutaneous+Whole_Blood+Testis+Lung";

const CASES = [
  ...["grch37", "grch38", "grch38_high_coverage"].flatMap(build => ["r2", "d"].map(st => ({ module: "LDmatrix", label: `${build} CEU ${st}`,
    ours: async ctx => { const o = await ldmatrix(ctx, { inputs: VECTORS, pop: "CEU", build }); return [["RS_number", ...o.ids].join("\t"), ...o.ids.map((id, i) => [id, ...o[st === "r2" ? "r2" : "d"][i]].join("\t"))].join("\n"); },
    theirs: async () => (await get("ldmatrix", { snps: VECTORS.join("\n"), pop: "CEU", r2_d: st, genome_build: build })).replace(/\r/g, "").trim() }))),
  ...[["rs3", "rs4", "CEU"], ["rs3", "rs148890987", "CEU+YRI"], ["rs12465506", "rs6721996", "ALL"], ["chr13:32446842", "RS4", "EUR"]].map(([v1, v2, pop]) => ({ module: "LDpair", label: `${v1} × ${v2} ${pop}`,
    ours: async ctx => canon(await ldpair(ctx, { var1: v1, var2: v2, pop, build: "grch37" })),
    theirs: async () => canon(JSON.parse(await get("ldpair", { var1: v1, var2: v2, pop, genome_build: "grch37", json_out: "true" })).map?.(({ request, ...r }) => r) ?? {}) })),
  ...[["rs3", "rs4", "ALL"], ["rs12465506", "rs6721996", "EUR+EAS"]].map(([v1, v2, pop]) => ({ module: "LDpop", label: `${v1} × ${v2} ${pop}`,
    ours: async ctx => await ldpop(ctx, { var1: v1, var2: v2, pop, build: "grch37" }), theirs: async () => asText(await get("ldpop", { var1: v1, var2: v2, pop, r2_d: "r2", genome_build: "grch37" })) })),
  ...[["rs3", "CEU", "grch37"], ["rs12465506", "EUR", "grch38"]].map(([v, pop, build]) => ({ module: "LDproxy", label: `${v} ${pop} ${build}`,
    ours: async ctx => { const o = await ldproxy(ctx, { var: v, pop, build, r2_d: "r2", window: 500000 }); return { text: o.text, error: o.error }; },
    theirs: async () => ({ ...asText(await get("ldproxy", { var: v, pop, r2_d: "r2", window: "500000", genome_build: build })), ...{} }) })),
  { module: "SNPclip", label: "test vectors CEU r² ≥ 0.5", ours: ctx => snpclip(ctx, { snps: VECTORS.join("\n"), pop: "CEU", r2_threshold: 0.5, maf_threshold: 0.01, build: "grch37" }),
    theirs: async () => asText(await post("snpclip", { snps: VECTORS.join("\n"), pop: "CEU", r2_threshold: "0.5", maf_threshold: "0.01", genome_build: "grch37" })) },
  { module: "SNPchip", label: "test vectors, all arrays", ours: async ctx => snpchip(ctx, { snps: VECTORS.join("\n"), platforms: (await ctx.platforms()).map(p => p.code).join("+"), build: "grch37" }),
    theirs: async ctx => asText(await post("snpchip", { snps: VECTORS.join("\n"), platforms: (await ctx.platforms()).map(p => p.code).join("+"), genome_build: "grch37" })) },
  { module: "LDtrait", label: "GWAS hits CEU", ours: ctx => ldtrait(ctx, { snps: "rs12465506\nrs3\nrs7837688", pop: "CEU", r2_d: "r2", r2_d_threshold: 0.1, window: 500000, build: "grch37" }),
    theirs: async () => asText(await post("ldtrait", { snps: "rs12465506\nrs3\nrs7837688", pop: "CEU", r2_d: "r2", r2_d_threshold: "0.1", window: "500000", genome_build: "grch37" })) },
  { module: "LDexpress", label: "rs3 + rs12465506, 4 tissues", ours: ctx => ldexpress(ctx, { snps: "rs3\nrs12465506", pop: "CEU", tissues: T, r2_d: "r2", r2_d_threshold: 0.1, p_threshold: 0.1, window: 500000, build: "grch37" }),
    theirs: async () => asText(await post("ldexpress", { snps: "rs3\nrs12465506", pop: "CEU", tissues: T, r2_d: "r2", r2_d_threshold: "0.1", p_threshold: "0.1", window: "500000", genome_build: "grch37" })) },
];

function row(cells) { const tr = document.createElement("tr"); for (const c of cells) { const td = document.createElement("td"); td.textContent = c; tr.append(td); } out.append(tr); return tr; }
const strip = o => (o && typeof o === "object" && !Array.isArray(o) ? Object.fromEntries(Object.entries(o).filter(([k, v]) => v !== undefined && k !== "warning")) : o);

const kernel = await loadKernel(new Uint8Array(await (await fetch(KERNEL_WASM)).arrayBuffer()));
const root = new URL("../release", location.href).href, manifest = await (await fetch(`${root}/manifest.json`)).json();
const ctx = await createContext({ root, open: async (url, byteLength) => cachedAsyncBuffer(await asyncBufferFromUrl({ url, byteLength })),
  kernel, manifest });
let same = 0, ties = 0;
for (const [k, c] of CASES.entries()) {
  status.textContent = `${k + 1} / ${CASES.length}: ${c.module} ${c.label}`;
  const t0 = performance.now(), ours = strip(await c.ours(ctx)), ms = Math.round(performance.now() - t0), theirs = strip(await c.theirs(ctx));
  const a = JSON.stringify(ours), b = JSON.stringify(theirs), verdict = a === b ? "identical" : tieOnly(ours, theirs) ? "tie order only" : "DIFFERENT";
  if (verdict === "identical") same++; else if (verdict === "tie order only") ties++;
  row([c.module, c.label, verdict, `${ms} ms`]).dataset.verdict = verdict;
  if (verdict === "DIFFERENT") console.log(c.module, c.label, { ours, theirs });
}
status.textContent = `done: ${same} identical, ${ties} tie order only, ${CASES.length - same - ties} different, of ${CASES.length} (${kernel.backend} kernel)`;
document.body.dataset.done = "1";
