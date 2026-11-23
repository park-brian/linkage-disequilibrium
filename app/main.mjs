// The browser app: every LDlink module, run in the page on the src/ core against a release served over HTTP
// (Range requests). One form per module, built from MODULES; outputs rendered as LDlink prints them.
import { asyncBufferFromUrl, cachedAsyncBuffer } from "hyparquet";
import { loadKernel, KERNEL_WASM } from "../src/ld/kernel.mjs";
import { createContext } from "../src/ldlink/context.mjs";
import { ldmatrix } from "../src/ldlink/ldmatrix.mjs";
import { ldpair } from "../src/ldlink/ldpair.mjs";
import { ldpop } from "../src/ldlink/ldpop.mjs";
import { ldhap, ldhapText } from "../src/ldlink/ldhap.mjs";
import { ldproxy } from "../src/ldlink/ldproxy.mjs";
import { snpclip } from "../src/ldlink/snpclip.mjs";
import { snpchip } from "../src/ldlink/snpchip.mjs";
import { ldtrait } from "../src/ldlink/ldtrait.mjs";
import { ldexpress } from "../src/ldlink/ldexpress.mjs";
import { ldassoc } from "../src/ldlink/ldassoc.mjs";
import { POPS } from "../src/ldlink/common.mjs";
import { renderMatrix } from "./matrix_view.mjs";
import { renderPyramid } from "./pyramid_view.mjs";
import { pyramid } from "./pyramid_pool.mjs";

const RELEASE = new URL("../release", location.href).href.replace(/\/$/, ""), PANELS = new URL("../panels/", location.href).href;
const EXAMPLE = ["rs13426489", "rs4442975", "rs13387042", "rs13412666", "rs6721996", "rs12465506", "rs12613955", "rs12620095"];
const $ = (sel, el = document) => el.querySelector(sel);
const status = msg => { $("#status").textContent = msg; };

// ---------- the core, in the page ----------
let ctxPromise;
function context() {
  return (ctxPromise ??= (async () => {
    const t0 = performance.now();
    const kernel = await loadKernel(new Uint8Array(await (await fetch(KERNEL_WASM)).arrayBuffer()));
    const open = async url => cachedAsyncBuffer(await asyncBufferFromUrl({ url }));
    const readPanel = async pop => (await (await fetch(`${PANELS}${pop}.txt`)).text()).split(/\s+/).filter(Boolean);
    const ctx = await createContext({ root: RELEASE, open, readPanel, kernel });
    status(`ready (${kernel.backend} kernel, ${Math.round(performance.now() - t0)} ms)`);
    return ctx;
  })());
}

// ---------- region view (pyramid): data on the page, LD in a worker ----------
async function region(ctx, a) {
  let chrom, lo, hi;
  const m = /^\s*(?:chr)?([0-9]{1,2}|X|Y)\s*:\s*([\d,_]+)\s*-\s*([\d,_]+)\s*$/i.exec(a.region);
  if (m) { chrom = m[1].toUpperCase(); lo = +m[2].replace(/[,_]/g, ""); hi = +m[3].replace(/[,_]/g, ""); }
  else {
    const v = (await ctx.planner.query({ build: a.build, inputs: [a.region.trim()] })).variants[0];
    if (!v?.found || v.pos == null) return { error: `${a.region} is not a region (chr:start-end) or a variant in dbSNP 155.` };
    chrom = v.chrom; lo = Math.max(0, v.pos - Math.round(+a.width / 2)); hi = v.pos + Math.round(+a.width / 2);
  }
  if (!(hi > lo)) return { error: "The region's end must be after its start." };
  if (hi - lo > 5_000_000) return { error: "Regions up to 5 Mb, please." };
  if (!a.tiers.length) return { error: "Choose at least one MAF tier." };
  const t0 = performance.now();
  const [rows, genes, recomb, gwas, eqtl, smp, ids] = await Promise.all([
    ctx.windowRows(a.build, chrom, lo, hi, a.tiers), ctx.genes(a.build, chrom, lo, hi), ctx.recomb(a.build, chrom, lo, hi),
    ctx.gwas(a.build, chrom, lo, hi), ctx.eqtlRange(a.build, chrom, lo, hi, a.tissues), ctx.samples(a.build, chrom), ctx.popIds(join(a.pop))]);
  const use = rows.filter(r => r.pos >= lo && r.pos <= hi).sort((x, y) => x.pos - y.pos || x.vcf_order - y.vcf_order);
  if (use.length < 2) return { error: "Fewer than two 1000G variants in this region and tier." };
  if (use.length > 60000) return { error: `${use.length.toLocaleString()} variants: narrow the region or the tiers (up to 60,000).` };
  const tData = Math.round(performance.now() - t0), W = Math.min(1100, Math.max(600, $("#output").clientWidth - 40));
  const res = await pyramid(use.map(r => ({ pos: r.pos, haps: r.haps, missing: r.missing, male_haploid: r.male_haploid })), { samples: smp.ids.length,
    haploidMales: [...smp.haploidMales], hapIndexes: smp.ids.flatMap((s, i) => (ids.has(s) ? [2 * i, 2 * i + 1] : [])),
    maxDist: Math.min(+a.maxDist, hi - lo), threshold: +a.threshold, stat: a.r2_d, W, lo, hi, radius: Math.max(0, Math.min(6, Math.floor(W / use.length / 2))) },
    +a.workers > 0 ? { size: +a.workers } : {});
  return { ...res, chrom, lo, hi, rows: use, labels: use.map(r => (r.rsid != null ? `rs${r.rsid}` : `${chrom}:${r.pos}`)), genes, recomb, gwas: gwas.map(r => ({ ...r })),
    eqtl: { rows: eqtl, tissues: a.tissues.length ? a.tissues.join(", ") : "none chosen" }, tData };
}

// ---------- module forms ----------
const f = {
  snps: (name = "snps", value = EXAMPLE.join("\n"), label = "Variants (rsID or chr:pos, one per line)") => ({ name, type: "textarea", label, value }),
  pop: (multiple = true, value = ["CEU"]) => ({ name: "pop", type: "select", label: "Population", options: POPS, multiple, value }),
  build: () => ({ name: "build", type: "select", label: "Genome build", options: ["grch37", "grch38", "grch38_high_coverage"], value: ["grch37"] }),
  stat: () => ({ name: "r2_d", type: "select", label: "Statistic", options: ["r2", "d"], value: ["r2"] }),
  num: (name, label, value) => ({ name, type: "number", label, value }),
  text: (name, label, value) => ({ name, type: "text", label, value }),
};
const join = v => (Array.isArray(v) ? v.join("+") : v);
const MODULES = {
  Pyramid: { fields: [f.text("region", "Region (chr:start-end) or a variant ± half the width", "chr22:30,000,000-31,000,000"), f.num("width", "Width around a variant (bp)", 500000), f.pop(true, ["EUR"]), f.build(), f.stat(),
      { name: "tiers", type: "select", label: "Variants (MAF tiers)", options: [["common", "common (MAF ≥ 5%)"], ["low", "low (1–5%)"], ["rare", "rare (< 1%)"]], multiple: true, value: ["common"] },
      f.num("maxDist", "Max pair distance (bp)", 500000), f.num("threshold", "Hide pairs below", 0.05), f.num("workers", "Workers (0 = one per core, up to 8)", 0),
      { name: "tissues", type: "select", label: "eQTL tissues", options: [], multiple: true, value: ["Whole_Blood", "Lung"], source: "tissues" }],
    run: (ctx, a) => region(ctx, a), show: (o, a) => showPyramid(o, a) },
  LDmatrix: { fields: [f.snps(), f.pop(), f.build(), f.stat()], run: (ctx, a) => ldmatrix(ctx, { inputs: lines(a.snps), pop: join(a.pop), build: a.build }), show: (o, a) => showMatrix(o, a) },
  LDpair: { fields: [f.text("var1", "Variant 1", "rs3"), f.text("var2", "Variant 2", "rs4"), f.pop(), f.build()], run: (ctx, a) => ldpair(ctx, { ...a, pop: join(a.pop) }), show: o => showPair(o) },
  LDpop: { fields: [f.text("var1", "Variant 1", "rs3"), f.text("var2", "Variant 2", "rs4"), f.pop(true, ["EUR"]), f.build()], run: (ctx, a) => ldpop(ctx, { ...a, pop: join(a.pop) }), show: o => showText(o) },
  LDhap: { fields: [f.snps(), f.pop(), f.build()], run: async (ctx, a) => { const o = await ldhap(ctx, { snps: lines(a.snps), pop: join(a.pop), build: a.build }); return o.error ? o : { text: ldhapText(o, a.build), warning: o.warning }; }, show: o => showText(o) },
  LDproxy: { fields: [f.text("var", "Variant", "rs3"), f.pop(), f.build(), f.stat(), f.num("window", "Window (bp)", 500000)], run: (ctx, a) => ldproxy(ctx, { ...a, pop: join(a.pop), window: +a.window }), show: o => showText(o) },
  SNPclip: { fields: [f.snps(), f.pop(), f.build(), f.num("r2_threshold", "r² threshold", 0.1), f.num("maf_threshold", "MAF threshold", 0.01)], run: (ctx, a) => snpclip(ctx, { ...a, pop: join(a.pop), r2_threshold: +a.r2_threshold, maf_threshold: +a.maf_threshold }), show: o => showText(o) },
  SNPchip: { fields: [f.snps(), f.build(), { name: "platforms", type: "select", label: "Arrays", options: [], multiple: true, value: [], source: "platforms" }], run: (ctx, a) => snpchip(ctx, { ...a, platforms: join(a.platforms) }), show: o => showText(o) },
  LDtrait: { fields: [f.snps(), f.pop(), f.build(), f.stat(), f.num("r2_d_threshold", "Threshold", 0.1), f.num("window", "Window (bp)", 500000)], run: (ctx, a) => ldtrait(ctx, { ...a, pop: join(a.pop), r2_d_threshold: +a.r2_d_threshold, window: +a.window }), show: o => showText(o) },
  LDexpress: { fields: [f.snps("snps", "rs3\nrs12465506"), f.pop(), f.build(), { name: "tissues", type: "select", label: "Tissues", options: [], multiple: true, value: ["Whole_Blood", "Lung"], source: "tissues" }, f.stat(), f.num("r2_d_threshold", "LD threshold", 0.1), f.num("p_threshold", "p threshold", 0.1), f.num("window", "Window (bp)", 500000)],
    run: (ctx, a) => ldexpress(ctx, { ...a, pop: join(a.pop), tissues: join(a.tissues), r2_d_threshold: +a.r2_d_threshold, p_threshold: +a.p_threshold, window: +a.window }), show: o => showText(o) },
  LDassoc: { fields: [{ name: "file", type: "file", label: "Association file (blank: LDlink's prostate example)" }, f.text("columns", "Column names (chr, pos, p)", "chr,pos,p"),
    { name: "region", type: "select", label: "Region", options: ["variant", "gene", "region"], value: ["variant"] }, f.text("origin", "Index variant (rsID; blank = lowest p)", "rs7837688"),
    f.text("name", "Gene (gene mode)", "MYC"), f.text("start", "Start (region mode)", "chr8:128000000"), f.text("end", "End (region mode)", "chr8:128700000"), f.num("window", "Window (bp)", 100000), f.pop(true, ["EUR"]), f.build()],
    run: async (ctx, a) => {
      // no file chosen: LDlink's bundled prostate-cancer example for the build
      const text = a.file ? await a.file.text() : await (await fetch(`../examples/prostate_example_${a.build === "grch37" ? "grch37" : "grch38"}.txt`)).text();
      return ldassoc(ctx, { text, columns: a.columns.split(",").map(s => s.trim()), region: a.region, origin: a.origin || undefined, name: a.name, start: a.start, end: a.end,
        window: a.region === "region" ? undefined : +a.window, pop: join(a.pop), build: a.build });
    }, show: o => showText(o) },
};
const lines = s => s.split(/\r?\n/);

function buildForm(name) {
  const form = $("#form"), m = MODULES[name];
  form.replaceChildren();
  for (const fd of m.fields) {
    const label = document.createElement("label"); label.textContent = fd.label;
    let el;
    if (fd.type === "textarea") { el = document.createElement("textarea"); el.value = fd.value; }
    else if (fd.type === "select") {
      el = document.createElement("select"); el.multiple = !!fd.multiple;
      const fill = opts => { el.replaceChildren(...opts.map(o => { const [v, t] = Array.isArray(o) ? o : [o, o]; const op = new Option(t, v); op.selected = fd.value.includes(v); return op; })); };
      fill(fd.options);
      if (fd.source) context().then(async ctx => fill(fd.source === "platforms" ? (await ctx.platforms()).map(p => [p.code, p.platform]) : (await ctx.gtexTissueIds())));
    } else { el = document.createElement("input"); el.type = fd.type; if (fd.value !== undefined) el.value = fd.value; if (fd.type === "number") el.step = "any"; }
    el.name = fd.name; label.append(el); form.append(label);
  }
  const go = document.createElement("button"); go.type = "submit"; go.textContent = `Run ${name}`; form.append(go);
  form.onsubmit = async e => {
    e.preventDefault(); go.disabled = true;
    const a = {};
    for (const fd of m.fields) { const el = form.elements[fd.name]; a[fd.name] = fd.type === "select" ? (fd.multiple ? [...el.selectedOptions].map(o => o.value) : el.value) : fd.type === "file" ? el.files[0] : el.value; }
    try {
      status("running…"); const ctx = await context(), t0 = performance.now();
      const o = await m.run(ctx, a);
      status(`${name}: ${Math.round(performance.now() - t0)} ms`);
      m.show(o, a);
    } catch (err) { console.error(err); show([notice(String(err?.stack ?? err), true)]); status("error"); }
    finally { go.disabled = false; }
  };
}

// ---------- output ----------
const show = nodes => $("#output").replaceChildren(...nodes);
function notice(msg, isError = false) { const d = document.createElement("div"); d.className = `notice${isError ? " error" : ""}`; d.textContent = msg; return d; }
function card(title, ...children) { const c = document.createElement("div"); c.className = "card"; const h = document.createElement("h2"); h.textContent = title; c.append(h, ...children); return c; }
function download(name, text, label = "Download") { const b = document.createElement("button"); b.type = "button"; b.textContent = label; b.onclick = () => { const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([text], { type: "text/plain" })); a.download = name; a.click(); }; return b; }
/** a tab-separated block → table */
function table(rows) {
  const t = document.createElement("table"), [head, ...body] = rows;
  t.append(Object.assign(document.createElement("thead"), {}), document.createElement("tbody"));
  t.tHead.append(tr(head, "th")); for (const r of body) t.tBodies[0].append(tr(r, "td"));
  return t;
}
function tr(cells, tag) { const r = document.createElement("tr"); for (const c of cells) { const e = document.createElement(tag); e.textContent = c; r.append(e); } return r; }
/** LDlink's text outputs: tab-separated blocks (separated by blank lines / "###" rules) and free-text notes */
function showText(o) {
  if (o.error) return show([notice(o.error, true)]);
  const out = [], blocks = o.text.split(/\n\s*\n|\n#{20,}\n/);
  for (const b of blocks) {
    const ls = b.split("\n").filter(l => l.length);
    if (!ls.length) continue;
    if (ls.length > 1 && ls[0].includes("\t")) { const tb = table(ls.map(l => l.split("\t"))); const bar = document.createElement("div"); bar.className = "toolbar"; bar.append(`${ls.length - 1} rows`, download("output.txt", o.text)); out.push(card("Result", bar, tb)); }
    else out.push(notice(ls.join("\n")));
  }
  if (o.warning) out.push(notice(o.warning));
  show(out);
}
function showPair(o) {
  if (!Array.isArray(o)) return show([notice(o.error, true)]);
  const p = o[0], s1 = p.snp1, s2 = p.snp2, c = p.two_by_two.cells;
  const grid = table([["", `${s2.rsnum} ${s2.allele_1.allele}`, `${s2.rsnum} ${s2.allele_2.allele}`, "total"],
    [`${s1.rsnum} ${s1.allele_1.allele}`, c.c11, c.c12, `${s1.allele_1.count} (${s1.allele_1.frequency})`],
    [`${s1.rsnum} ${s1.allele_2.allele}`, c.c21, c.c22, `${s1.allele_2.count} (${s1.allele_2.frequency})`],
    ["total", `${s2.allele_1.count} (${s2.allele_1.frequency})`, `${s2.allele_2.count} (${s2.allele_2.frequency})`, p.two_by_two.total]]);
  const haps = table([["Haplotype", "Count", "Frequency"], ...Object.values(p.haplotypes).map(h => [h.alleles, h.count, h.frequency])]);
  const st = table([["D′", "R²", "χ²", "p"], [p.statistics.d_prime, p.statistics.r2, p.statistics.chisq, p.statistics.p]]);
  show([card(`${s1.rsnum} (${s1.coord}) × ${s2.rsnum} (${s2.coord})`, grid), card("Haplotypes", haps), card("Statistics", st, ...p.corr_alleles.map(x => Object.assign(document.createElement("p"), { textContent: x, className: "note" }))), ...(p.warning ? [notice(p.warning)] : [])]);
}
function showPyramid(o, a) {
  if (o.error) return show([notice(o.error, true)]);
  window.__lastPyramid = o; // debugging aid: the last region result (grids, rows, timings)
  const bar = document.createElement("div"); bar.className = "toolbar";
  bar.append(`chr${o.chrom}:${o.lo.toLocaleString()}–${o.hi.toLocaleString()} · ${o.rows.length.toLocaleString()} variants · ${o.evaluated.toLocaleString()} pairs evaluated, ${o.pairs.toLocaleString()} ≥ ${a.threshold} · data ${o.tData} ms, LD + pixels ${o.ms.wall} ms on ${o.workers} worker${o.workers > 1 ? "s" : ""} (${o.backend}; slowest share ${o.ms.kernel} + ${o.ms.aggregate} ms), merge ${o.ms.merge} ms`);
  const holder = document.createElement("div");
  show([card(`LD pyramid — ${a.r2_d === "r2" ? "r²" : "D′"}, ${join(a.pop)}, ${a.build} (each pixel: the strongest pair there)`, bar, holder)]);
  renderPyramid(holder, { ...o, stat: a.r2_d });
}
async function showMatrix(o, a) {
  if (!o.ids.length) return show([notice("No variants could be used (not in dbSNP / 1000G, or not biallelic).", true)]);
  const ctx = await context(), v = (await ctx.planner.query({ build: a.build, inputs: o.ids })).variants;
  const chrom = v[0].chrom, pos = v.map(x => x.pos), lo = Math.min(...pos), hi = Math.max(...pos), pad = Math.max(5000, (hi - lo) * 0.025);
  const genes = await ctx.genes(a.build, chrom, Math.max(0, lo - pad), hi + pad);
  const tsv = (m, label) => [["RS_number", ...o.ids].join("\t"), ...o.ids.map((id, i) => [id, ...m[i]].join("\t"))].join("\n") + "\n";
  const bar = document.createElement("div"); bar.className = "toolbar";
  bar.append(`${o.ids.length} variants · chr${chrom}:${lo.toLocaleString()}–${hi.toLocaleString()}`, download("r2.txt", tsv(o.r2), "Download r²"), download("dprime.txt", tsv(o.d), "Download D′"));
  const holder = document.createElement("div");
  show([card(`LDmatrix — D′ below the diagonal, r² above (${join(a.pop)}, ${a.build})`, bar, holder),
    card("r²", table([["RS_number", ...o.ids], ...o.ids.map((id, i) => [id, ...o.r2[i]])])), card("D′", table([["RS_number", ...o.ids], ...o.ids.map((id, i) => [id, ...o.d[i]])]))]);
  renderMatrix(holder, { ids: o.ids, pos, chrom, r2: o.r2, d: o.d, genes, lo: lo - pad, hi: hi + pad });
}

// ---------- navigation ----------
const nav = $("#modules");
for (const name of Object.keys(MODULES)) {
  const b = document.createElement("button"); b.type = "button"; b.textContent = name;
  b.onclick = () => { for (const x of nav.children) x.removeAttribute("aria-current"); b.setAttribute("aria-current", "page"); buildForm(name); history.replaceState(null, "", `#${name}`); show([]); };
  nav.append(b);
}
const open = () => { const name = MODULES[location.hash.slice(1)] ? location.hash.slice(1) : "LDmatrix"; [...nav.children].find(b => b.textContent === name).click(); };
window.addEventListener("hashchange", open);
open();
status("loading core…"); context().catch(err => { console.error(err); status("failed to load the release"); show([notice(String(err), true)]); });
