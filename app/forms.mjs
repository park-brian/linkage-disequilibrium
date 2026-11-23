// The module table: each module's form fields, how a submitted form calls its src/ldlink port, and how the result
// is shown. buildForm(name) draws one module's form into #form.
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
import { region } from "./region.mjs";
import { $, join, show, notice, showText, showPair, showPyramid, showMatrix } from "./render.mjs";

const EXAMPLE = ["rs13426489", "rs4442975", "rs13387042", "rs13412666", "rs6721996", "rs12465506", "rs12613955", "rs12620095"];
const lines = s => s.split(/\r?\n/);
const f = {
  snps: (name = "snps", value = EXAMPLE.join("\n"), label = "Variants (rsID or chr:pos, one per line)") => ({ name, type: "textarea", label, value }),
  pop: (multiple = true, value = ["CEU"]) => ({ name: "pop", type: "select", label: "Population", options: POPS, multiple, value }),
  build: () => ({ name: "build", type: "select", label: "Genome build", options: ["grch37", "grch38", "grch38_high_coverage"], value: ["grch37"] }),
  stat: () => ({ name: "r2_d", type: "select", label: "Statistic", options: ["r2", "d"], value: ["r2"] }),
  num: (name, label, value) => ({ name, type: "number", label, value }),
  text: (name, label, value) => ({ name, type: "text", label, value }),
  tissues: (label = "Tissues") => ({ name: "tissues", type: "select", label, options: [], multiple: true, value: ["Whole_Blood", "Lung"], source: "tissues" }),
};

export const MODULES = {
  Pyramid: {
    fields: [f.text("region", "Region (chr:start-end) or a variant ± half the width", "chr22:30,000,000-31,000,000"), f.num("width", "Width around a variant (bp)", 500000), f.pop(true, ["EUR"]), f.build(), f.stat(),
      { name: "tiers", type: "select", label: "Variants (MAF tiers)", options: [["common", "common (MAF ≥ 5%)"], ["low", "low (1–5%)"], ["rare", "rare (< 1%)"]], multiple: true, value: ["common"] },
      f.num("maxDist", "Max pair distance (bp)", 500000), f.num("threshold", "Hide pairs below", 0.05), f.num("workers", "Workers (0 = one per core, up to 8)", 0), f.tissues("eQTL tissues")],
    run: (ctx, a) => region(ctx, a, $("#output").clientWidth), show: showPyramid,
  },
  LDmatrix: { fields: [f.snps(), f.pop(), f.build(), f.stat()], run: (ctx, a) => ldmatrix(ctx, { inputs: lines(a.snps), pop: join(a.pop), build: a.build }), show: showMatrix },
  LDpair: { fields: [f.text("var1", "Variant 1", "rs3"), f.text("var2", "Variant 2", "rs4"), f.pop(), f.build()], run: (ctx, a) => ldpair(ctx, { ...a, pop: join(a.pop) }), show: showPair },
  LDpop: { fields: [f.text("var1", "Variant 1", "rs3"), f.text("var2", "Variant 2", "rs4"), f.pop(true, ["EUR"]), f.build()], run: (ctx, a) => ldpop(ctx, { ...a, pop: join(a.pop) }), show: showText },
  LDhap: {
    fields: [f.snps(), f.pop(), f.build()], show: showText,
    run: async (ctx, a) => { const o = await ldhap(ctx, { snps: lines(a.snps), pop: join(a.pop), build: a.build }); return o.error ? o : { text: ldhapText(o, a.build), warning: o.warning }; },
  },
  LDproxy: { fields: [f.text("var", "Variant", "rs3"), f.pop(), f.build(), f.stat(), f.num("window", "Window (bp)", 500000)], run: (ctx, a) => ldproxy(ctx, { ...a, pop: join(a.pop), window: +a.window }), show: showText },
  SNPclip: {
    fields: [f.snps(), f.pop(), f.build(), f.num("r2_threshold", "r² threshold", 0.1), f.num("maf_threshold", "MAF threshold", 0.01)], show: showText,
    run: (ctx, a) => snpclip(ctx, { ...a, pop: join(a.pop), r2_threshold: +a.r2_threshold, maf_threshold: +a.maf_threshold }),
  },
  SNPchip: {
    fields: [f.snps(), f.build(), { name: "platforms", type: "select", label: "Arrays", options: [], multiple: true, value: [], source: "platforms" }], show: showText,
    run: (ctx, a) => snpchip(ctx, { ...a, platforms: join(a.platforms) }),
  },
  LDtrait: {
    fields: [f.snps(), f.pop(), f.build(), f.stat(), f.num("r2_d_threshold", "Threshold", 0.1), f.num("window", "Window (bp)", 500000)], show: showText,
    run: (ctx, a) => ldtrait(ctx, { ...a, pop: join(a.pop), r2_d_threshold: +a.r2_d_threshold, window: +a.window }),
  },
  LDexpress: {
    fields: [f.snps("snps", "rs3\nrs12465506"), f.pop(), f.build(), f.tissues(), f.stat(), f.num("r2_d_threshold", "LD threshold", 0.1), f.num("p_threshold", "p threshold", 0.1), f.num("window", "Window (bp)", 500000)], show: showText,
    run: (ctx, a) => ldexpress(ctx, { ...a, pop: join(a.pop), tissues: join(a.tissues), r2_d_threshold: +a.r2_d_threshold, p_threshold: +a.p_threshold, window: +a.window }),
  },
  LDassoc: {
    fields: [{ name: "file", type: "file", label: "Association file (blank: LDlink's prostate example)" }, f.text("columns", "Column names (chr, pos, p)", "chr,pos,p"),
      { name: "region", type: "select", label: "Region", options: ["variant", "gene", "region"], value: ["variant"] }, f.text("origin", "Index variant (rsID; blank = lowest p)", "rs7837688"),
      f.text("name", "Gene (gene mode)", "MYC"), f.text("start", "Start (region mode)", "chr8:128000000"), f.text("end", "End (region mode)", "chr8:128700000"), f.num("window", "Window (bp)", 100000), f.pop(true, ["EUR"]), f.build()],
    show: showText,
    run: async (ctx, a) => {
      // no file chosen: LDlink's bundled prostate-cancer example for the build
      const text = a.file ? await a.file.text() : await (await fetch(`../examples/prostate_example_${a.build === "grch37" ? "grch37" : "grch38"}.txt`)).text();
      return ldassoc(ctx, { text, columns: a.columns.split(",").map(s => s.trim()), region: a.region, origin: a.origin || undefined, name: a.name, start: a.start, end: a.end,
        window: a.region === "region" ? undefined : +a.window, pop: join(a.pop), build: a.build });
    },
  },
};

function field(fd, context) {
  let el;
  if (fd.type === "textarea") { el = document.createElement("textarea"); el.value = fd.value; }
  else if (fd.type === "select") {
    el = document.createElement("select"); el.multiple = !!fd.multiple;
    const fill = opts => el.replaceChildren(...opts.map(o => { const [v, t] = Array.isArray(o) ? o : [o, o]; const op = new Option(t, v); op.selected = fd.value.includes(v); return op; }));
    fill(fd.options);
    if (fd.source) context().then(async ctx => fill(fd.source === "platforms" ? (await ctx.platforms()).map(p => [p.code, p.platform]) : await ctx.gtexTissueIds()));
  } else { el = document.createElement("input"); el.type = fd.type; if (fd.value !== undefined) el.value = fd.value; if (fd.type === "number") el.step = "any"; }
  el.name = fd.name;
  return el;
}
const value = (fd, el) => (fd.type === "select" ? (fd.multiple ? [...el.selectedOptions].map(o => o.value) : el.value) : fd.type === "file" ? el.files[0] : el.value);

/** draw module `name`'s form; context() → the core (memoized), status(msg) → the status line */
export function buildForm(name, { context, status }) {
  const form = $("#form"), m = MODULES[name];
  form.replaceChildren();
  for (const fd of m.fields) { const label = document.createElement("label"); label.textContent = fd.label; label.append(field(fd, context)); form.append(label); }
  const go = document.createElement("button"); go.type = "submit"; go.textContent = `Run ${name}`; form.append(go);
  form.onsubmit = async e => {
    e.preventDefault(); go.disabled = true;
    const a = Object.fromEntries(m.fields.map(fd => [fd.name, value(fd, form.elements[fd.name])]));
    try {
      status("running…"); const ctx = await context(), t0 = performance.now();
      const o = await m.run(ctx, a);
      status(`${name}: ${Math.round(performance.now() - t0)} ms`);
      await m.show(o, a, ctx);
    } catch (err) { console.error(err); show([notice(String(err?.stack ?? err), true)]); status("error"); }
    finally { go.disabled = false; }
  };
}
