// The module table: each module's form fields, how a submitted form calls its src/ldlink port, and how the result
// is shown. buildForm(name) draws one module's form into #form.
import { ldmatrix } from "../src/ldlink/ldmatrix.js";
import { ldpair } from "../src/ldlink/ldpair.js";
import { ldpop } from "../src/ldlink/ldpop.js";
import { ldhap, ldhapText } from "../src/ldlink/ldhap.js";
import { ldproxy } from "../src/ldlink/ldproxy.js";
import { snpclip } from "../src/ldlink/snpclip.js";
import { snpchip } from "../src/ldlink/snpchip.js";
import { ldtrait } from "../src/ldlink/ldtrait.js";
import { ldexpress } from "../src/ldlink/ldexpress.js";
import { ldassoc } from "../src/ldlink/ldassoc.js";
import { POPS } from "../src/ldlink/common.js";
import { ldherit } from "../src/ldlink/ldherit.js";
import { ldcorrelation } from "../src/ldlink/ldcorrelation.js";
import { ldscore } from "../src/ldlink/ldscore.js";
import { region } from "./region.js";
import { $, join, show, notice, showText, showPair, showPyramid, showMatrix, showProxy, showAssoc, showLog } from "./render.js";

const EXAMPLE = ["rs13426489", "rs4442975", "rs13387042", "rs13412666", "rs6721996", "rs12465506", "rs12613955", "rs12620095"];
const lines = s => s.split(/\r?\n/);
const f = {
  snps: (name = "snps", value = EXAMPLE.join("\n"), label = "Variants (rsID or chr:pos, one per line)") => ({ name, type: "textarea", label, value }),
  pop: (multiple = true, value = ["CEU"]) => ({ name: "pop", type: "select", label: "Population", options: POPS, multiple, value }),
  build: () => ({ name: "build", type: "select", label: "Genome build", options: ["grch37", "grch38", "grch38_high_coverage"], value: ["grch37"] }),
  stat: () => ({ name: "r2_d", type: "select", label: "Statistic", options: ["r2", "d"], value: ["r2"] }),
  num: (name, label, value) => ({ name, type: "number", label, value }),
  text: (name, label, value) => ({ name, type: "text", label, value }),
  annotate: () => ({ name: "annotate", type: "select", label: "Plot labels", options: [["forge", "FORGEdb score"], ["regulome", "RegulomeDB rank"], ["no", "none"]], value: ["forge"] }),
  transcripts: () => ({ name: "collapse", type: "select", label: "Genes", options: [["yes", "collapse transcripts"], ["no", "all transcripts"]], value: ["yes"] }),
  ldscPop: () => ({ name: "pop", type: "select", label: "LD score population", options: LDSC_POPS, value: ["EUR"] }),
  scale: () => ({ name: "scale", type: "select", label: "Scale", options: ["observed", "liability"], value: ["observed"] }),
  tissues: (label = "Tissues") => ({ name: "tissues", type: "select", label, options: [], multiple: true, value: ["Whole_Blood", "Lung"], source: "tissues" }),
};

// LDlink's LD score populations (client/src/components/select/ldscore-pop-select.tsx)
const LDSC_POPS = ["EUR", "EAS", "AFR", "AMR", "ASJ", "EST", "FIN", "NFE", "NWE", "SEU"];
/** one of LDlink's example inputs from the release's examples/ → a Blob (a failed fetch throws, naming the file) */
async function example(ctx, name) {
  const res = await fetch(`${ctx.root}/examples/${name}`);
  if (!res.ok) throw new Error(`could not read the example ${name} (HTTP ${res.status})`);
  return res.blob();
}
/** a sumstats input: the uploaded file, or LDlink's example from the release → { text, filename, isExample } */
async function sumstats(ctx, file, name) {
  if (file) return { text: await file.text(), filename: file.name, isExample: false };
  return { text: await (await example(ctx, name)).text(), filename: name, isExample: true };
}

export const MODULES = {
  Pyramid: {
    fields: [f.text("region", "Region (chr:start-end) or a variant ± half the width", "chr22:30,000,000-31,000,000"), f.num("width", "Width around a variant (bp)", 500000), f.pop(true, ["EUR"]), f.build(), f.stat(),
      { name: "tiers", type: "select", label: "Variants (MAF tiers)", options: [["common", "common (MAF ≥ 5%)"], ["low", "low (1–5%)"], ["rare", "rare (< 1%)"]], multiple: true, value: ["common"] },
      f.num("maxDist", "Max pair distance (bp)", 500000), f.num("threshold", "Hide pairs below", 0.05), f.num("workers", "Workers (0 = one per core, up to 8)", 0), f.tissues("eQTL tissues"),
      { name: "scores", type: "select", label: "FORGEdb / RegulomeDB track", options: [["yes", "show"], ["no", "hide"]], value: ["yes"] }],
    run: (ctx, a) => region(ctx, a, $("#output").clientWidth), show: showPyramid,
  },
  LDmatrix: { fields: [f.snps(), f.pop(), f.build(), f.stat()], run: (ctx, a) => ldmatrix(ctx, { inputs: lines(a.snps), pop: join(a.pop), build: a.build }), show: showMatrix },
  LDpair: { fields: [f.text("var1", "Variant 1", "rs3"), f.text("var2", "Variant 2", "rs4"), f.pop(), f.build()], run: (ctx, a) => ldpair(ctx, { ...a, pop: join(a.pop) }), show: showPair },
  LDpop: { fields: [f.text("var1", "Variant 1", "rs3"), f.text("var2", "Variant 2", "rs4"), f.pop(true, ["EUR"]), f.build()], run: (ctx, a) => ldpop(ctx, { ...a, pop: join(a.pop) }), show: showText },
  LDhap: {
    fields: [f.snps(), f.pop(), f.build()], show: showText,
    run: async (ctx, a) => { const o = await ldhap(ctx, { snps: lines(a.snps), pop: join(a.pop), build: a.build }); return o.error ? o : { text: ldhapText(o, a.build), warning: o.warning }; },
  },
  LDproxy: {
    fields: [f.text("var", "Variant", "rs3"), f.pop(), f.build(), f.stat(), f.num("window", "Window (bp)", 500000), f.annotate(), f.transcripts()], show: showProxy,
    run: (ctx, a) => ldproxy(ctx, { var: a.var, pop: join(a.pop), build: a.build, r2_d: a.r2_d, window: a.window }, { plot: true }),
  },
  SNPclip: {
    fields: [f.snps(), f.pop(), f.build(), f.num("r2_threshold", "r² threshold", 0.1), f.num("maf_threshold", "MAF threshold", 0.01)], show: showText,
    run: (ctx, a) => snpclip(ctx, { ...a, pop: join(a.pop), r2_threshold: a.r2_threshold, maf_threshold: a.maf_threshold }),
  },
  SNPchip: {
    fields: [f.snps(), f.build(), { name: "platforms", type: "select", label: "Arrays", options: [], multiple: true, value: "all", source: "platforms" }], show: showText,
    run: (ctx, a) => snpchip(ctx, { ...a, platforms: join(a.platforms) }),
  },
  LDtrait: {
    fields: [f.snps(), f.pop(), f.build(), f.stat(), f.num("r2_d_threshold", "Threshold", 0.1), f.num("window", "Window (bp)", 500000)], show: showText,
    run: (ctx, a) => ldtrait(ctx, { ...a, pop: join(a.pop), r2_d_threshold: a.r2_d_threshold, window: a.window }),
  },
  LDexpress: {
    fields: [f.snps("snps", "rs3\nrs12465506"), f.pop(), f.build(), f.tissues(), f.stat(), f.num("r2_d_threshold", "LD threshold", 0.1), f.num("p_threshold", "p threshold", 0.1), f.num("window", "Window (bp)", 500000)], show: showText,
    run: (ctx, a) => ldexpress(ctx, { ...a, pop: join(a.pop), tissues: join(a.tissues), r2_d_threshold: a.r2_d_threshold, p_threshold: a.p_threshold, window: a.window }),
  },
  LDassoc: {
    fields: [{ name: "file", type: "file", label: "Association file" },
      { name: "example", type: "select", label: "No file: use", value: ["no"],
        options: [["no", "nothing (choose a file)"], ["yes", "LDlink's prostate example (17 MB)"]] },
      f.text("columns", "Column names (chr, pos, p)", "chr,pos,p"),
      { name: "region", type: "select", label: "Region", options: ["variant", "gene", "region"], value: ["variant"] }, f.text("origin", "Index variant (rsID; blank = lowest p)", "rs7837688"),
      f.text("name", "Gene (gene mode)", "MYC"), f.text("start", "Start (region mode)", "chr8:128000000"), f.text("end", "End (region mode)", "chr8:128700000"), f.num("window", "Window (bp)", 100000), f.pop(true, ["EUR"]), f.build(), f.stat(), f.annotate(), f.transcripts()],
    show: showAssoc,
    run: async (ctx, a) => {
      // no file chosen: LDlink's bundled prostate-cancer example for the build, only when asked for (a 17 MB download)
      if (!a.file && a.example !== "yes") throw new Error("Choose an association file, or LDlink's prostate example.");
      const ex = `prostate_example_${a.build === "grch37" ? "grch37" : "grch38"}.txt`;
      const text = a.file ? await a.file.text() : await (await example(ctx, ex)).text();
      return ldassoc(ctx, { text, columns: a.columns.split(",").map(s => s.trim()), region: a.region, origin: a.origin || undefined, name: a.name, start: a.start, end: a.end,
        window: a.region === "region" ? undefined : a.window, pop: join(a.pop), build: a.build }, { plot: true, stat: a.r2_d });
    },
  },
  LDherit: {
    fields: [{ name: "file", type: "file", label: "Summary statistics (blank: LDlink's example)" },
      { name: "example", type: "select", label: "Example", options: ["BBJ_HDLC22.txt", "BBJ_LDLC22.txt"], value: ["BBJ_HDLC22.txt"] },
      f.ldscPop(), f.scale(), f.text("samp_prev", "Sample prevalence (liability)", ""), f.text("pop_prev", "Population prevalence (liability)", "")],
    show: showLog,
    run: async (ctx, a) => ldherit(ctx, { ...(await sumstats(ctx, a.file, a.example)), pop: a.pop, scale: a.scale, samp_prev: a.samp_prev, pop_prev: a.pop_prev }),
  },
  LDcorrelation: {
    fields: [{ name: "file", type: "file", label: "Summary statistics 1 (blank: example)" },
      { name: "example", type: "select", label: "Example 1", options: ["BBJ_HDLC22.txt", "BBJ_LDLC22.txt"], value: ["BBJ_HDLC22.txt"] },
      { name: "file2", type: "file", label: "Summary statistics 2 (blank: example)" },
      { name: "example2", type: "select", label: "Example 2", options: ["BBJ_HDLC22.txt", "BBJ_LDLC22.txt"], value: ["BBJ_LDLC22.txt"] },
      f.ldscPop(), f.scale(), f.text("samp_prev", "Sample prevalences (two, comma-separated)", ""), f.text("pop_prev", "Population prevalences (two)", "")],
    show: showLog,
    run: async (ctx, a) => {
      const s1 = await sumstats(ctx, a.file, a.example), s2 = await sumstats(ctx, a.file2, a.example2);
      return ldcorrelation(ctx, { text: s1.text, filename: s1.filename, text2: s2.text, filename2: s2.filename, isExample: s1.isExample && s2.isExample,
        pop: a.pop, scale: a.scale, samp_prev: a.samp_prev, pop_prev: a.pop_prev });
    },
  },
  LDscore: {
    fields: [{ name: "files", type: "file", multiple: true, label: "PLINK .bed, .bim, .fam (blank: LDlink's chromosome 22 example)" },
      { ...f.num("ldwindow", "Window", 1), asText: true }, { name: "windUnit", type: "select", label: "Window unit", options: [["cm", "cM"], ["kb", "kb"]], value: ["cm"] }],
    show: showLog,
    run: async (ctx, a) => {
      const picked = a.files?.length ? a.files : null, names = picked ? picked.map(f => f.name) : ["22.bed", "22.bim", "22.fam"], files = {};
      for (const [k, n] of names.entries()) {
        const blob = picked ? picked[k] : await example(ctx, n);
        files[n] = n.endsWith(".bed") ? new Uint8Array(await blob.arrayBuffer()) : await blob.text();
      }
      return ldscore({ filename: names.join(","), ldwindow: a.ldwindow, windUnit: a.windUnit, isExample: !picked, files });
    },
  },
};

/** a field; `given` (from the hash query) replaces its default: text as is, selections split on "+" */
function field(fd, context, given, loading) {
  if (given != null) fd = { ...fd, value: fd.type === "select" ? given.split("+") : given };
  let el;
  if (fd.type === "textarea") { el = document.createElement("textarea"); el.value = fd.value; }
  else if (fd.type === "select") {
    el = document.createElement("select"); el.multiple = !!fd.multiple;
    // value "all" selects every option (SNPchip: LDlink's form starts with all platforms, client/src/app/snpchip/form.tsx)
    const fill = opts => el.replaceChildren(...opts.map(o => { const [v, t] = Array.isArray(o) ? o : [o, o]; const op = new Option(t, v); op.selected = fd.value === "all" || fd.value.includes(v); return op; }));
    fill(fd.options);
    // a list that fails to load keeps its defaults (a linked run still goes ahead, and reports what it can)
    const options = async ctx => (fd.source === "platforms" ? (await ctx.platforms()).map(p => [p.code, p.platform]) : ctx.gtexTissueIds());
    if (fd.source) loading.push(context().then(async ctx => fill(await options(ctx)))
      .catch(err => console.error(`${fd.label} list:`, err)));
  } else { el = document.createElement("input"); el.type = fd.type; if (fd.value !== undefined) el.value = fd.value; if (fd.type === "number") el.step = "any"; if (fd.multiple) el.multiple = true; }
  el.name = fd.name;
  return el;
}
/** a field's value: selections as arrays (multiple) or strings, files as File(s), number fields as numbers (unless
 *  asText: the port parses the text itself, as LDlink does) */
const value = (fd, el) => (fd.type === "select" ? (fd.multiple ? [...el.selectedOptions].map(o => o.value) : el.value)
  : fd.type === "file" ? (fd.multiple ? [...el.files] : el.files[0]) : fd.type === "number" && !fd.asText ? +el.value : el.value);

/** a module's inputs as the hash query: every field but files, selections joined with "+" */
export const toQuery = (name, a) => new URLSearchParams(Object.entries(a).filter(([k]) => MODULES[name].fields.some(fd => fd.name === k && fd.type !== "file"))
  .map(([k, v]) => [k, Array.isArray(v) ? v.join("+") : String(v ?? "")])).toString();

/**
 * draw module `name`'s form; context() → the core (memoized), status(msg) → the status line; query: initial field values
 * (URLSearchParams from the hash), with run=1 to run at once. Running writes the inputs back to the hash.
 */
// the latest run (or form): a result that arrives after another run started, or after the module changed, is dropped
let current = 0, running = false;
export function buildForm(name, { context, status, query = new URLSearchParams() }) {
  const form = $("#form"), m = MODULES[name];
  current++;
  if (running) { running = false; status(""); } // the run in flight no longer reports
  form.replaceChildren();
  const loading = []; // option lists read from the release (SNPchip arrays, GTEx tissues)
  for (const fd of m.fields) { const label = document.createElement("label"); label.textContent = fd.label; label.append(field(fd, context, query.get(fd.name), loading)); form.append(label); }
  const go = document.createElement("button"); go.type = "submit"; go.textContent = `Run ${name}`; form.append(go);
  form.onsubmit = async e => {
    e.preventDefault(); go.disabled = true;
    const run = ++current, stale = () => run !== current;
    running = true;
    const a = Object.fromEntries(m.fields.map(fd => [fd.name, value(fd, form.elements[fd.name])]));
    history.replaceState(null, "", `${location.search}#${name}?${toQuery(name, a)}`);
    try {
      status("running…"); const ctx = await context(), t0 = performance.now();
      const o = await m.run(ctx, a);
      if (stale()) return;
      running = false;
      status(`${name}: ${Math.round(performance.now() - t0)} ms`);
      await m.show(o, a, ctx);
    } catch (err) {
      console.error(err);
      if (!stale()) { running = false; show([notice(String(err?.message ?? err), true)]); status("error"); }
    }
    finally { go.disabled = false; }
  };
  if (query.get("run") === "1") Promise.all(loading).then(() => form.requestSubmit(go)); // run once the options are there
}
