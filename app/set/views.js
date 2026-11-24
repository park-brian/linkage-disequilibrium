// The answers. What is held (the ticked variants, or the whole set) decides which one the page opens with: one variant
// → its proxies (LDproxy), two → the pair (LDpair, with LDpop's populations), several → the matrix (LDmatrix), a region
// → the locus; the others that fit are one click away (haplotypes = LDhap, prune = SNPclip). Each view runs LDlink's
// port on the held labels, is drawn for reading, and can be taken away: copied, downloaded as LDlink's exact text, or
// saved as an image. A pair LDlink's own resolution refuses is answered from the two 1000G rows, and says so.
import { ldmatrix } from "../../src/parity/ldlink/ldmatrix.js";
import { ldhap, ldhapText } from "../../src/parity/ldlink/ldhap.js";
import { snpclip } from "../../src/parity/ldlink/snpclip.js";
import { ldproxy } from "../../src/parity/ldlink/ldproxy.js";
import { ldpair } from "../../src/parity/ldlink/ldpair.js";
import { ldpop } from "../../src/parity/ldlink/ldpop.js";
import { ldtrait } from "../../src/parity/ldlink/ldtrait.js";
import { ldexpress } from "../../src/parity/ldlink/ldexpress.js";
import { snpchip } from "../../src/parity/ldlink/snpchip.js";
import { cells, ldStats, vcfAlleles } from "../../src/parity/ldlink/common.js";
import { canvas, tooltip, pointer } from "../views/tracks.js";
import { h, ramp, STAT, fix, download, knob, textTable } from "./dom.js";
import { GROUPS } from "./context.js";
import { locus, known } from "./locus.js";

const pops = S => S.pop.replaceAll("+", " + ");
const rowOf = (env, label) => env.set.ok.find(r => r.label === label)?.g ?? known.get(label);

/** a section that runs when opened, and stays open from one variant to the next (the page state's `open` list) */
function section(env, key, title, load) {
  const body = h("div.body"), d = h("details", { open: env.S.open.split(",").includes(key) }, h("summary", title), body);
  const run = () => { if (d.ran) return; d.ran = true; body.append(h("p.aside", "reading…"));
    load(body).catch(err => { console.error(err); body.replaceChildren(h("p.say.error", String(err?.message ?? err))); }); };
  d.addEventListener("toggle", () => {
    const open = new Set(env.S.open.split(",").filter(Boolean)); d.open ? open.add(key) : open.delete(key);
    env.remember({ open: [...open].join(",") }); if (d.open) run(); });
  if (d.open) run();
  return d;
}

async function proxies(env) {
  const { ctx, S, set, labels: [label] } = env, thr = +(S.thr ?? 0.8), key = S.stat === "r2" ? "r2" : "dprime";
  const o = await ldproxy(ctx, { var: label, pop: S.pop, build: S.build, r2_d: S.stat, window: 500000 }, { plot: true });
  if (env.stale()) return;
  const single = { snps: label, pop: S.pop, build: S.build }, more = [
    section(env, "traits", "Traits in LD with it · LDtrait", async body => {
      const t = await ldtrait(ctx, { ...single, r2_d: S.stat, r2_d_threshold: 0.1, window: 500000 });
      if (t.error) return body.replaceChildren(h("p.note", t.error));
      const v = textTable(t.text, 1);
      body.replaceChildren(h("p.where", `${v.count} GWAS Catalog associations within 500 kb at ${STAT[S.stat]} 0.1 or more`,
        h("span.links", download("LDlink text", "ldtrait.txt", t.text))), v.table, v.notes ? h("p.note", v.notes) : null);
    }),
    section(env, "expression", "Expression in LD with it · LDexpress", async body => {
      const ids = await ctx.gtexTissueIds();
      const run = async tissue => {
        const pick = h("label.knob", "tissue", h("select", { onchange: e => run(e.target.value) },
          ids.map(t => h("option", { value: t, selected: t === tissue }, t.replaceAll("_", " ")))));
        body.classList.add("busy");
        const t = await ldexpress(ctx, { ...single, tissues: tissue, r2_d: S.stat, r2_d_threshold: 0.1, p_threshold: 0.1,
          window: 500000 });
        body.classList.remove("busy");
        if (t.error) return body.replaceChildren(h("p.where", pick), h("p.note", t.error));
        const v = textTable(t.text, 1);
        body.replaceChildren(h("p.where", pick, `${v.count} GTEx eQTLs within 500 kb, ${STAT[S.stat]} and p at 0.1`,
          h("span.links", download("LDlink text", "ldexpress.txt", t.text))), v.table, v.notes ? h("p.note", v.notes) : null);
      };
      await run("Whole_Blood");
    }),
    section(env, "arrays", "Arrays that carry it · SNPchip", async body => {
      const t = await snpchip(ctx, { snps: label, build: S.build, platforms: (await ctx.platforms()).map(p => p.code).join("+") });
      if (t.error) return body.replaceChildren(h("p.note", t.error));
      const names = (t.text.split("\n")[1]?.split("\t")[2] ?? "").split(",").filter(Boolean);
      body.replaceChildren(h("p.where", names.length ? `${names.length} arrays` : "On none of LDlink's arrays",
        h("span.links", download("LDlink text", "snpchip.txt", t.text))), h("p.note", names.join(" · ")));
    })];
  if (o.error) { env.headline(`Proxies of ${label}`); return env.show(h("p.say", o.error), more); }
  const rows = o.plot.rows.filter(r => r.rs !== o.plot.query.rs && r[key] >= thr), have = new Set(set.rows.map(r => r.label));
  const fresh = [...new Set(rows.map(r => r.rs))].filter(l => !have.has(l));
  const base = set.kind === "list" ? set.rows.map(r => r.label) : [label];
  env.headline(`${rows.length.toLocaleString()} proxies of ${label}`, knob(`${STAT[S.stat]} at least`, thr, 0.05,
    v => env.update({ thr: String(v) })),
  fresh.length ? h("button.go", { type: "button", title: "Back undoes it",
    onclick: () => env.update({ q: [...new Set([...base, ...fresh])].join(" "), pick: label },
      { note: `Added ${fresh.length} proxies of ${label}.` }) },
  set.kind === "list" ? `Add ${fresh.length} to the set` : `Make a set of these ${fresh.length + 1}`) : null);
  env.copy("Copy rsIDs", () => rows.map(r => r.rs).join("\n"));
  env.exact("LDlink text", "proxy.txt", o.text);
  env.show(rows.length ? h("div.scroll", h("table", h("thead", h("tr", h("th", "Variant"), h("th.num", "Distance"), h("th.num", "r²"),
    h("th.num", "D′"), h("th.num", "MAF"), h("th", "Correlated alleles"), h("th.num", "Regulome"), h("th", "Function"))),
  h("tbody", rows.slice(0, 500).map(r => h("tr.still", h("td.id", r.rs), h("td.num", r.dist.toLocaleString()),
    h("td.num.heat", { style: `--v:${r.r2}` }, fix(r.r2)), h("td.num", fix(r.dprime)), h("td.num", fix(r.maf)), h("td.pos", r.match),
    h("td.num", r.regulome ?? ""), h("td.traits", r.func === "none" ? "" : r.func))))))
    : h("p.say", `No variant reaches ${STAT[S.stat]} ${thr} with ${label}.`), more);
  env.caption(`Variants within 500 kb in LD with ${label} in ${pops(S)}, strongest first. ${o.warning ?? ""}`,
    rows.length > 500 ? " The first 500 are listed." : "");
}

/** the pair as drawn: statistics, the 2×2 of haplotype counts with allele frequencies in the margins, notes */
function pairView({ ids, alleles, counts: [A, B, C, D], stats, notes }) {
  const N = A + B + C + D, cell = v => h("td.num.heat", { style: `--v:${N ? (v / N).toFixed(3) : 0}` }, String(v)),
    f = v => h("td.num.sum", N ? (v / N).toFixed(3) : "–");
  return [h("div.pairtop", h("dl.stats", stats.flatMap(([k, v]) => [h("dt", k), h("dd", v)])),
    h("table.two", h("thead", h("tr", h("th"), h("th.num", `${ids[1]} ${alleles[1][0]}`), h("th.num", `${ids[1]} ${alleles[1][1]}`),
      h("th.num"))),
    h("tbody", h("tr", h("th", `${ids[0]} ${alleles[0][0]}`), cell(A), cell(B), f(A + B)),
      h("tr", h("th", `${ids[0]} ${alleles[0][1]}`), cell(C), cell(D), f(C + D)),
      h("tr", h("th"), f(A + C), f(B + D), h("td.num.sum", String(N)))))), ...notes.map(n => h("p.note", n))];
}

async function pair(env) {
  const { ctx, S, labels } = env;
  const [o, all] = await Promise.all([ldpair(ctx, { var1: labels[0], var2: labels[1], pop: S.pop, build: S.build }),
    ldpop(ctx, { var1: labels[0], var2: labels[1], pop: "ALL", build: S.build })]);
  if (env.stale()) return;
  env.headline(`${labels[0]} and ${labels[1]}`);
  if (!Array.isArray(o)) { // LDlink refuses: the same numbers from the two rows, if the page holds them
    const rows = labels.map(l => rowOf(env, l));
    if (!rows[0] || !rows[1]) return env.say(o.error);
    const c = await ctx.withPrepared(S.build, rows, S.pop, p => ctx.kernel.counts(p, [0], [1])[0]);
    if (env.stale()) return;
    const t = cells(c), st = ldStats(...t), n4 = v => (st ? v.toFixed(4) : "NA");
    env.show(pairView({ ids: labels, alleles: rows.map(r => vcfAlleles(r.ref, r.alt)), counts: t,
      stats: [["r²", n4(st?.r2)], ["D′", n4(st?.dprime)], ["χ²", n4(st?.chisq)],
        ["p", st ? (st.p >= 0.0001 ? st.p.toFixed(4) : "<0.0001") : "NA"]],
      notes: [`LDlink does not answer for this pair: “${o.error.trim()}” These numbers come from the two 1000 Genomes rows at `
        + "these positions; they are not an LDlink output."] }));
    return env.caption(`Haplotype counts and statistics in ${pops(S)}, computed here.`);
  }
  const p = o[0], s1 = p.snp1, s2 = p.snp2, c = p.two_by_two.cells, st = p.statistics;
  let byPop = null;
  if (all.text) {
    const rows = new Map(all.text.split("\n\n")[0].trimEnd().split("\n").slice(1).map(l => l.split("\t")).map(r => [r[1], r]));
    const val = (stat, v) => h("span.val", { style: `background:${ramp(stat, 0.62 * (+v || 0))}` }, v === "NA" ? "–" : (+v).toFixed(2));
    const bar = (code, cls = "") => { const r = rows.get(code);
      return r && h(`div.poprow${cls}`, { title: `${r[0]} · N ${r[2]}\nr² ${r[5]} · D′ ${r[6]}\n${r[3]}\n${r[4]}` }, h("span.code", code),
        val("r2", r[5]), val("d", r[6])); };
    // the five super-populations side by side, each value on its statistic's ink, so all 26 are read without scrolling
    byPop = h("div.popgrid", [["ALL"], ...Object.entries(GROUPS).map(([g, [, subs]]) => [g, ...Object.keys(subs)])].map(([g, ...subs]) => h(
      "div", h("div.poprow.headrow", h("span.code"), h("span.val", "r²"), h("span.val", "D′")),
        bar(g, ".super"), subs.map(code => bar(code)))));
  }
  env.copy("Copy r² and D′", () => `${s1.rsnum}\t${s2.rsnum}\t${pops(S)}\tr2=${st.r2}\tDprime=${st.d_prime}`);
  env.exact("LDlink JSON", "ldpair.json", JSON.stringify(o, null, 2));
  if (all.text) env.exact("LDlink populations", "ldpop.txt", all.text);
  env.show(h("p.where", `${s1.coord} and ${s2.coord}, in ${pops(S)}`), pairView({ ids: [s1.rsnum, s2.rsnum],
    alleles: [[s1.allele_1.allele, s1.allele_2.allele], [s2.allele_1.allele, s2.allele_2.allele]],
    counts: [c.c11, c.c12, c.c21, c.c22].map(Number), stats: [["r²", st.r2], ["D′", st.d_prime], ["χ²", st.chisq], ["p", st.p]],
    notes: [p.corr_alleles.join(". "), ...(p.warning ? [p.warning] : [])] }),
  byPop ? h("p.sub", "The same pair in every 1000 Genomes population") : null, byPop);
}

async function matrix(env) {
  const { ctx, S, labels } = env;
  if (labels.length < 2) return env.say("A matrix needs at least two variants.");
  if (labels.length > 300) return env.say(`${labels.length.toLocaleString()} variants: tick up to 300 for a matrix, or read the `
    + "region in the locus.");
  const o = await ldmatrix(ctx, { inputs: labels, pop: S.pop, build: S.build }), n = o.ids.length;
  if (env.stale()) return;
  if (n < 2) return env.say("Fewer than two of these variants are usable in 1000 Genomes.");
  const named = n <= 40, L = named ? 92 : 0, W = Math.max(320, env.figure.clientWidth), cell = Math.max(2, Math.min(34,
    Math.floor((Math.min(W, 860) - L) / n))), size = n * cell, { cv, g } = canvas(L + size, L + size);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    if (i === j) continue;
    const v = i < j ? o.r2[i][j] : o.d[i][j], gap = cell >= 8 ? 1 : 0;
    g.fillStyle = v === "NA" ? "#eceae4" : ramp(i < j ? "r2" : "d", +v);
    g.fillRect(L + j * cell, L + i * cell, cell - gap, cell - gap);
    if (cell >= 30 && v !== "NA") { g.fillStyle = +v > 0.55 ? "#fff" : "#6b7078";
      g.font = "10px ui-monospace, monospace"; g.textAlign = "center";
      g.textBaseline = "middle"; g.fillText((+v).toFixed(2), L + j * cell + cell / 2 - 0.5, L + i * cell + cell / 2); }
  }
  if (named) {
    g.font = "11px ui-monospace, monospace"; g.textBaseline = "middle"; g.fillStyle = "#16181d";
    o.ids.forEach((id, i) => {
      g.textAlign = "right"; g.fillText(id, L - 8, L + i * cell + cell / 2);
      g.save(); g.translate(L + i * cell + cell / 2, L - 8); g.rotate(-Math.PI / 2);
      g.textAlign = "left"; g.fillText(id, 0, 0); g.restore();
    });
  }
  const at = (mx, my) => { const i = Math.floor((my - L) / cell), j = Math.floor((mx - L) / cell);
    return i >= 0 && j >= 0 && i < n && j < n && i !== j ? [i, j] : null; };
  const tip = tooltip(cv, (mx, my) => { const c = at(mx, my); cv.style.cursor = c ? "pointer" : "default";
    return c && `${o.ids[c[0]]} × ${o.ids[c[1]]}\nr² ${o.r2[c[0]][c[1]]} · D′ ${o.d[c[0]][c[1]]}\nClick for the pair`; });
  cv.addEventListener("click", e => { const c = at(...pointer(cv, e)); if (c) env.hold([o.ids[c[0]], o.ids[c[1]]]); });
  const tsv = m => [["RS_number", ...o.ids], ...o.ids.map((id, i) => [id, ...m[i]])].map(r => r.join("\t")).join("\n") + "\n";
  env.headline(`r² and D′ among ${n} variants`);
  env.image(cv, "ldmatrix.png");
  env.exact("LDlink r²", "r2.txt", tsv(o.r2)); env.exact("LDlink D′", "dprime.txt", tsv(o.d));
  env.show(cv, tip);
  env.caption(h("span.ramp", { style: `background:linear-gradient(90deg,${ramp("r2", 0)},${ramp("r2", 1)})` }), " r² above the diagonal, ",
    h("span.ramp", { style: `background:linear-gradient(90deg,${ramp("d", 0)},${ramp("d", 1)})` }),
    ` D′ below, in ${pops(S)}, in position order. Click a cell for that pair.`);
}

async function haplotypes(env) {
  const { ctx, S, labels } = env;
  if (labels.length > 30) return env.say(`${labels.length.toLocaleString()} variants: haplotypes are read over up to 30. `
    + "Tick the ones to phase together.");
  const o = await ldhap(ctx, { snps: labels, pop: S.pop, build: S.build });
  if (env.stale()) return;
  if (o.error) return env.say(o.error);
  const snps = Object.keys(o.snps).sort().map(k => o.snps[k]), all = Object.keys(o.haplotypes).sort().map(k => o.haplotypes[k]);
  const haps = all.filter((x, i) => i === 0 || x.Frequency >= 0.01), rest = all.length - haps.length;
  const ref = haps[0].Haplotype.split("_"), max = haps[0].Frequency;
  env.headline(`${all.length} haplotype${all.length === 1 ? "" : "s"} over ${snps.length} variants`);
  env.exact("LDlink text", "haplotypes.txt", ldhapText(o, S.build));
  env.show(h("table.haps",
    h("thead", h("tr", h("th"), haps.map(x => h("th.freq", { title: `${x.Count} chromosomes` },
      h("i", { style: `height:${Math.max(2, Math.round((x.Frequency / max) * 44))}px` }), `${(x.Frequency * 100).toFixed(1)}%`)))),
    h("tbody", snps.map((s, i) => h("tr.still", h("td.id", { title: `${s.Coord} · ${s.Alleles}` }, s.RS),
      haps.map(x => { const a = x.Haplotype.split("_")[i]; return h(`td.allele${a === ref[i] ? "" : ".alt"}`, a); }))))));
  env.caption(`Each column is one haplotype in ${pops(S)}, commonest first; a filled cell differs from the commonest. `
    + (rest ? `${rest} under 1% are not drawn (they are in LDlink's text).` : "") + (o.warning ? ` ${o.warning}` : ""));
}

async function prune(env) {
  const { ctx, S, labels } = env, r2t = +(S.r2t ?? 0.1), maft = +(S.maft ?? 0.01);
  if (labels.length > 5000) return env.say(`${labels.length.toLocaleString()} variants: pruning takes up to 5,000.`);
  const o = await snpclip(ctx, { snps: labels.join("\n"), pop: S.pop, r2_threshold: r2t, maf_threshold: maft, build: S.build });
  if (env.stale()) return;
  if (o.error) return env.say(o.error);
  const rows = o.text.trimEnd().split("\n").slice(1).map(l => l.split("\t")), kept = rows.filter(r => r.at(-1) === "Variant kept.");
  env.headline(`${kept.length} of ${rows.length} independent`, knob("r² below", r2t, 0.05, v => env.update({ r2t: String(v) })),
    knob("MAF at least", maft, 0.01, v => env.update({ maft: String(v) })),
    kept.length < rows.length ? h("button.go", { type: "button", title: "Replace the set with the kept variants (Back undoes it)",
      onclick: () => env.update({ q: kept.map(r => r[0]).join(" "), pick: "" },
        { note: `Pruned ${rows.length} variants to ${kept.length}.` }) }, `Keep these ${kept.length}`) : null);
  env.copy("Copy kept rsIDs", () => kept.map(r => r[0]).join("\n"));
  env.exact("LDlink text", "snpclip.txt", o.text);
  env.show(h("div.scroll", h("table", h("tbody", rows.slice(0, 1000).map(r => { const on = r.at(-1) === "Variant kept.";
    return h(`tr.still${on ? "" : ".out"}`, h("td.mark", on ? "●" : "○"), h("td.id", r[0]), h("td.pos", r[1]), h("td.pos", r[2]),
      h("td.why", on ? "" : r.at(-1).replace(", variant removed.", "").replace("Variant in LD with", "in LD with"))); })))));
  env.caption(`Greedy pruning in the order given, in ${pops(S)}: each kept variant removes every later one with `
    + `r² of ${r2t} or more.`);
}

/** key → [name, LDlink's name for it, draw] */
export const VIEWS = { proxies: ["Proxies", "LDproxy", proxies], pair: ["Pair", "LDpair · LDpop", pair],
  matrix: ["Matrix", "LDmatrix", matrix],
  haplotypes: ["Haplotypes", "LDhap", haplotypes], prune: ["Prune", "SNPclip", prune], locus: ["Locus", "region", locus] };
/** the shape of what is held → the views that answer it, the first being the one the page opens with */
export const FOR = { one: ["proxies", "locus"], two: ["pair", "matrix", "haplotypes", "locus"],
  many: ["matrix", "haplotypes", "prune", "locus"], region: ["locus", "prune"], none: [] };
