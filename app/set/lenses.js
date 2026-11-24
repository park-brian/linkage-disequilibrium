// The lenses other than the locus: each is one LDlink module run on the set (or on the selected rows), drawn for
// reading, with LDlink's exact text one link away. Matrix = LDmatrix, Haplotypes = LDhap, Proxies = LDproxy on the
// anchor (and "add to the set"), Prune = SNPclip (and "keep these"), Pair = LDpair with LDpop's populations, and
// Traits = LDtrait, Expression = LDexpress, Arrays = SNPchip as their tables.
import { ldmatrix } from "../../src/parity/ldlink/ldmatrix.js";
import { ldhap, ldhapText } from "../../src/parity/ldlink/ldhap.js";
import { ldproxy } from "../../src/parity/ldlink/ldproxy.js";
import { snpclip } from "../../src/parity/ldlink/snpclip.js";
import { ldpair } from "../../src/parity/ldlink/ldpair.js";
import { ldpop } from "../../src/parity/ldlink/ldpop.js";
import { ldtrait } from "../../src/parity/ldlink/ldtrait.js";
import { ldexpress } from "../../src/parity/ldlink/ldexpress.js";
import { snpchip } from "../../src/parity/ldlink/snpchip.js";
import { cells, ldStats, vcfAlleles } from "../../src/parity/ldlink/common.js";
import { canvas, tooltip, pointer } from "../views/tracks.js";
import { h, ramp, STAT, fix } from "./dom.js";
import { GROUPS } from "./context.js";
import { locus, known } from "./locus.js";

const DARK = "#16181d";
const using = env => `${env.active.length.toLocaleString()} ${env.picked ? "selected " : ""}variant${env.active.length === 1 ? "" : "s"}`;
const number = (label, value, step, onChange) => h("label.knob", label, h("input", { type: "number", value, step, min: 0, max: 1,
  onchange: e => onChange(+e.target.value) }));

async function matrix(env) {
  const { ctx, S, active } = env;
  if (active.length < 2) return env.say("A matrix needs at least two variants.");
  if (active.length > 300) return env.say(`${using(env)}: select up to 300 for a matrix, or read the region in the locus lens.`);
  const o = await ldmatrix(ctx, { inputs: active.map(r => r.label), pop: S.pop, build: S.build }), n = o.ids.length;
  if (env.stale()) return;
  if (n < 2) return env.say("Fewer than two of these variants are usable in 1000 Genomes.");
  const named = n <= 40, L = named ? 92 : 0, W = Math.max(320, env.figure.clientWidth), cell = Math.max(2, Math.min(34,
    Math.floor((Math.min(W, 860) - L) / n))), size = n * cell, { cv, g } = canvas(L + size, L + size);
  const val = (i, j) => (i < j ? o.r2[i][j] : o.d[i][j]);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    if (i === j) continue;
    const v = val(i, j), stat = i < j ? "r2" : "d";
    g.fillStyle = v === "NA" ? "#eceae4" : ramp(stat, +v);
    g.fillRect(L + j * cell, L + i * cell, cell - (cell >= 8 ? 1 : 0), cell - (cell >= 8 ? 1 : 0));
    if (cell >= 30 && v !== "NA") { g.fillStyle = +v > 0.55 ? "#fff" : "#6b7078";
      g.font = "10px ui-monospace, monospace"; g.textAlign = "center";
      g.textBaseline = "middle"; g.fillText((+v).toFixed(2), L + j * cell + cell / 2 - 0.5, L + i * cell + cell / 2); }
  }
  if (named) {
    g.font = "11px ui-monospace, monospace"; g.fillStyle = DARK; g.textBaseline = "middle";
    o.ids.forEach((id, i) => {
      g.textAlign = "right"; g.fillText(id, L - 8, L + i * cell + cell / 2);
      g.save(); g.translate(L + i * cell + cell / 2, L - 8); g.rotate(-Math.PI / 2);
        g.textAlign = "left"; g.fillText(id, 0, 0); g.restore();
    });
  }
  const at = (mx, my) => { const i = Math.floor((my - L) / cell), j = Math.floor((mx - L) / cell);
    return i >= 0 && j >= 0 && i < n && j < n && i !== j ? [i, j] : null; };
  const tip = tooltip(cv, (mx, my) => { const c = at(mx, my); cv.style.cursor = c ? "pointer" : "default";
    return c && `${o.ids[c[0]]} × ${o.ids[c[1]]}\nr² ${o.r2[c[0]][c[1]]} · D′ ${o.d[c[0]][c[1]]}`; });
  cv.addEventListener("click", e => { const c = at(...pointer(cv, e)); if (c) env.select([o.ids[c[0]], o.ids[c[1]]]); });
  const tsv = m => [["RS_number", ...o.ids], ...o.ids.map((id, i) => [id, ...m[i]])].map(r => r.join("\t")).join("\n") + "\n";
  env.headline(`${using(env)} in position order`);
  env.exact("LDlink r²", "r2.txt", tsv(o.r2)); env.exact("LDlink D′", "dprime.txt", tsv(o.d));
  env.show(cv, tip);
  env.caption(h("span.ramp", { style: `background:linear-gradient(90deg,${ramp("r2", 0)},${ramp("r2", 1)})` }), " r² above the diagonal, ",
    h("span.ramp", { style: `background:linear-gradient(90deg,${ramp("d", 0)},${ramp("d", 1)})` }),
    ` D′ below, in ${S.pop.replaceAll("+", " + ")}. Click a cell to select that pair.`);
}

async function haplotypes(env) {
  const { ctx, S, active } = env;
  if (active.length > 30) return env.say(`${using(env)}: haplotypes are read over up to 30. Select the ones to phase together.`);
  if (!active.length) return env.say("No usable variant in this set.");
  const o = await ldhap(ctx, { snps: active.map(r => r.label), pop: S.pop, build: S.build });
  if (env.stale()) return;
  if (o.error) return env.say(o.error);
  const snps = Object.keys(o.snps).sort().map(k => o.snps[k]), all = Object.keys(o.haplotypes).sort().map(k => o.haplotypes[k]);
  const haps = all.filter((x, i) => i === 0 || x.Frequency >= 0.01), rest = all.length - haps.length;
  const ref = haps[0].Haplotype.split("_"), max = haps[0].Frequency;
  const table = h("table.haps",
    h("thead", h("tr", h("th"), haps.map(x => h("th.freq", { title: `${x.Count} chromosomes` },
      h("i", { style: `height:${Math.max(2, Math.round((x.Frequency / max) * 44))}px` }), `${(x.Frequency * 100).toFixed(1)}%`)))),
    h("tbody", snps.map((s, i) => h("tr", h("td.id", { title: `${s.Coord} · ${s.Alleles}` }, s.RS), haps.map(x => {
      const a = x.Haplotype.split("_")[i];
      return h(`td.allele${a === ref[i] ? "" : ".alt"}`, a);
    })))));
  env.headline(`${using(env)}, ${all.length} haplotype${all.length === 1 ? "" : "s"}`);
  env.exact("LDlink text", "haplotypes.txt", ldhapText(o, S.build));
  env.show(table);
  env.caption(`Each column is one haplotype in ${S.pop.replaceAll("+", " + ")}, commonest first; a filled cell differs from the commonest. `
    + (rest ? `${rest} under 1% are not drawn (they are in LDlink's text).` : "") + (o.warning ? ` ${o.warning}` : ""));
}

async function proxies(env) {
  const { ctx, S, set, sel } = env, anchor = sel[0] ?? set.ok[0]?.label, thr = +(S.thr ?? 0.8), key = S.stat === "r2" ? "r2" : "dprime";
  if (!anchor) return env.say("Select a variant to list its proxies.");
  const o = await ldproxy(ctx, { var: anchor, pop: S.pop, build: S.build, r2_d: S.stat, window: 500000 }, { plot: true });
  if (env.stale()) return;
  if (o.error) return env.say(o.error);
  const rows = o.plot.rows.filter(r => r.rs !== o.plot.query.rs && r[key] >= thr), have = new Set(set.rows.map(r => r.label));
  const fresh = [...new Set(rows.map(r => r.rs))].filter(l => !have.has(l));
  const base = set.kind === "list" ? set.rows.map(r => r.label) : [anchor];
  env.headline(`${rows.length.toLocaleString()} proxies of ${anchor} within 500 kb`, number(`${STAT[S.stat]} at least`, thr, 0.05,
    v => env.update({ thr: String(v) })),
  fresh.length ? h("button", { type: "button", onclick: () => env.update({ q:
    [...new Set([...base, ...fresh])].join("\n"), sel: anchor }) },
    set.kind === "list" ? `Add ${fresh.length} to the set` : `Make a set of these ${fresh.length + 1}`) : null);
  env.exact("LDlink text", "proxy.txt", o.text);
  env.show(rows.length ? h("div.scroll", h("table", h("thead", h("tr", h("th", "Variant"), h("th.num", "Distance"), h("th.num", "r²"),
    h("th.num", "D′"), h("th.num", "MAF"), h("th", "Correlated alleles"), h("th.num", "RDB"), h("th", "Function"))),
  h("tbody", rows.slice(0, 500).map(r => h("tr", { onclick: () => env.select([anchor, r.rs]) }, h("td.id", r.rs),
    h("td.num", r.dist.toLocaleString()), h("td.num.heat", { style: `--v:${r.r2}` }, fix(r.r2)), h("td.num", fix(r.dprime)),
    h("td.num", fix(r.maf)), h("td.pos", r.match), h("td.num", r.regulome ?? ""), h("td.traits", r.func === "none" ? "" : r.func))))))
    : h("p.say", `No variant reaches ${STAT[S.stat]} ${thr} with ${anchor}.`));
  env.caption(`Variants in LD with ${anchor} in ${S.pop.replaceAll("+", " + ")}, strongest first. `
    + (o.warning ?? ""), rows.length > 500 ? " The first 500 are listed." : "");
}

async function prune(env) {
  const { ctx, S, active } = env, r2t = +(S.r2t ?? 0.1), maft = +(S.maft ?? 0.01);
  if (!active.length) return env.say("No usable variant in this set.");
  if (active.length > 5000) return env.say(`${using(env)}: pruning takes up to 5,000.`);
  const o = await snpclip(ctx, { snps: active.map(r => r.label).join("\n"), pop: S.pop, r2_threshold: r2t, maf_threshold: maft,
    build: S.build });
  if (env.stale()) return;
  if (o.error) return env.say(o.error);
  const rows = o.text.trimEnd().split("\n").slice(1).map(l => l.split("\t")), kept = rows.filter(r => r.at(-1) === "Variant kept.");
  env.headline(`${kept.length} of ${rows.length} kept`, number("r² below", r2t, 0.05, v => env.update({ r2t: String(v) })),
    number("MAF at least", maft, 0.01, v => env.update({ maft: String(v) })),
    kept.length < rows.length ? h("button", { type: "button", onclick: () => env.update({ q: kept.map(r => r[0]).join("\n"), sel: "" }) },
      `Keep these ${kept.length}`) : null);
  env.exact("LDlink text", "snpclip.txt", o.text);
  env.show(h("div.scroll", h("table", h("tbody", rows.slice(0, 1000).map(r => { const on = r.at(-1) === "Variant kept.";
    return h(`tr${on ? "" : ".out"}`, h("td.mark", on ? "●" : "○"), h("td.id", r[0]), h("td.pos", r[1]), h("td.pos", r[2]),
      h("td.why", on ? "" : r.at(-1).replace(", variant removed.", "").replace("Variant in LD with", "in LD with"))); })))));
  env.caption(`Greedy pruning in the order given, in ${S.pop.replaceAll("+", " + ")}: each kept variant removes every later one with `
    + `r² of ${r2t} or more.`);
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
      h("tr", h("th"), f(A + C), f(B + D), h("td.num.sum", String(N)))))), ...notes.map(n => h("p.say", n))];
}

/** a pair LDlink's own resolution refuses (a panel label dbSNP places elsewhere, say): the same statistics from the
 *  two 1000G rows the page already holds, with LDlink's answer quoted */
async function pairFromRows(env, rows, refusal) {
  const { ctx, S, sel } = env, c = await ctx.withPrepared(S.build, rows, S.pop, p => ctx.kernel.counts(p, [0], [1])[0]);
  if (env.stale()) return;
  const t = cells(c), st = ldStats(...t), n4 = v => (st ? v.toFixed(4) : "NA");
  env.headline(`${sel[0]} (chr${rows[0].chrom}:${rows[0].pos}) and ${sel[1]} (chr${rows[1].chrom}:${rows[1].pos})`);
  env.show(pairView({ ids: sel, alleles: rows.map(r => vcfAlleles(r.ref, r.alt)), counts: t,
    stats: [["r²", n4(st?.r2)], ["D′", n4(st?.dprime)], ["χ²", n4(st?.chisq)],
      ["p", st ? (st.p >= 0.0001 ? st.p.toFixed(4) : "<0.0001") : "NA"]],
    notes: [`LDlink does not answer for this pair: “${refusal.trim()}” `
      + "These numbers come from the two 1000 Genomes rows at these positions."] }));
  env.caption(`Haplotype counts and statistics in ${S.pop.replaceAll("+", " + ")}, computed here; not an LDlink output.`);
}

async function pair(env) {
  const { ctx, S, sel, set } = env;
  if (sel.length !== 2) return env.say("Select two variants (Ctrl-click the second in the table, or click a pair in the locus or matrix).");
  const [o, pops] = await Promise.all([ldpair(ctx, { var1: sel[0], var2: sel[1], pop: S.pop, build: S.build }),
    ldpop(ctx, { var1: sel[0], var2: sel[1], pop: "ALL", build: S.build })]);
  if (env.stale()) return;
  if (!Array.isArray(o)) {
    const rows = sel.map(l => set.ok.find(r => r.label === l)?.g ?? known.get(l));
    return rows[0] && rows[1] ? pairFromRows(env, rows, o.error) : env.say(o.error);
  }
  const p = o[0], s1 = p.snp1, s2 = p.snp2, c = p.two_by_two.cells, st = p.statistics;
  let byPop = null;
  if (pops.text) {
    const rows = new Map(pops.text.split("\n\n")[0].trimEnd().split("\n").slice(1).map(l => l.split("\t")).map(r => [r[1], r]));
    const fill = (stat, v) => h("span.track", h("i", { style: `width:${(+v || 0) * 100}%;background:${ramp(stat,
      0.25 + 0.75 * (+v || 0))}` }));
    const bar = (code, cls = "") => { const r = rows.get(code);
      return r && h(`div.poprow${cls}`, { title: `${r[0]} · N ${r[2]}\n${r[3]}\n${r[4]}` }, h("span.code", code), fill("r2", r[5]),
        h("span.val", r[5]), fill("d", r[6]), h("span.val", r[6])); };
    byPop = h("div.pops", h("div.poprow.headrow", h("span.code"), h("span.track", "r²"), h("span.val"), h("span.track", "D′"),
      h("span.val")),
    bar("ALL", ".super"), Object.entries(GROUPS).map(([g, [, subs]]) => [bar(g, ".super"), Object.keys(subs).map(code => bar(code))]));
  }
  env.headline(`${s1.rsnum} (${s1.coord}) and ${s2.rsnum} (${s2.coord})`);
  env.exact("LDlink JSON", "ldpair.json", JSON.stringify(o, null, 2));
  if (pops.text) env.exact("LDlink populations", "ldpop.txt", pops.text);
  env.show(pairView({ ids: [s1.rsnum, s2.rsnum],
    alleles: [[s1.allele_1.allele, s1.allele_2.allele], [s2.allele_1.allele, s2.allele_2.allele]],
    counts: [c.c11, c.c12, c.c21, c.c22].map(Number), stats: [["r²", st.r2], ["D′", st.d_prime], ["χ²", st.chisq], ["p", st.p]],
    notes: [p.corr_alleles.join(". "), ...(p.warning ? [p.warning] : [])] }), byPop);
  env.caption(`Haplotype counts and statistics in ${S.pop.replaceAll("+", " + ")}; below, the same pair in every 1000 Genomes population.`);
}

/** a lens whose result is one of LDlink's text tables (LDtrait, LDexpress, SNPchip): the table as it is, the rest as notes */
const textLens = (file, run, describe) => async env => {
  if (!env.active.length) return env.say("No usable variant in this set.");
  const controls = [], o = await run(env, controls);
  if (env.stale()) return;
  env.headline(describe(env), ...controls);
  if (o.error) return env.say(o.error);
  const lines = o.text.split("\n"), cut = lines.findIndex(l => !l.includes("\t")), end = cut < 0 ? lines.length : cut;
  const rows = lines.slice(0, end).map(l => l.split("\t")), notes = lines.slice(end).filter(Boolean).join(" ");
  env.headline(`${(rows.length - 1).toLocaleString()} ${describe(env)}`, ...controls);
  env.exact("LDlink text", file, o.text);
  env.show(h("div.scroll", h("table", h("thead", h("tr", rows[0].map(c => h("th", c)))), h("tbody", rows.slice(1, 501).map(r => h("tr",
    r.map((c, i) => h(i ? "td.cell" : "td.id", { title: c }, c))))))));
  env.caption(rows.length > 501 ? "The first 500 rows are listed. " : "", notes);
};
const labelsOf = env => env.active.map(r => r.label).join("\n");

const traits = textLens("ldtrait.txt", (env, controls) => {
  const { ctx, S } = env, thr = +(S.tthr ?? 0.1);
  controls.push(number(`${STAT[S.stat]} at least`, thr, 0.05, v => env.update({ tthr: String(v) })));
  return ldtrait(ctx, { snps: labelsOf(env), pop: S.pop, build: S.build, r2_d: S.stat, r2_d_threshold: thr, window: 500000 });
}, env => `GWAS Catalog associations in LD with ${using(env)}, within 500 kb`);

const expression = textLens("ldexpress.txt", async (env, controls) => {
  const { ctx, S } = env, tissue = S.tissue ?? "Whole_Blood", ids = await ctx.gtexTissueIds();
  controls.push(h("label.knob", "tissue", h("select", { onchange: e => env.update({ tissue: e.target.value }) },
    ids.map(t => h("option", { value: t, selected: t === tissue }, t.replaceAll("_", " "))))));
  return ldexpress(ctx, { snps: labelsOf(env), pop: S.pop, build: S.build, tissues: tissue, r2_d: S.stat, r2_d_threshold: 0.1,
    p_threshold: 0.1, window: 500000 });
}, env => `GTEx eQTLs in LD with ${using(env)} (${STAT[env.S.stat]} and p at 0.1), within 500 kb`);

const arrays = textLens("snpchip.txt", async env => snpchip(env.ctx, { snps: labelsOf(env), build: env.S.build,
  platforms: (await env.ctx.platforms()).map(p => p.code).join("+") }), () => "variants and the genotyping arrays that carry them");

export const LENSES = { locus: ["Locus", locus], matrix: ["Matrix", matrix], haplotypes: ["Haplotypes", haplotypes],
  proxies: ["Proxies", proxies], prune: ["Prune", prune], pair: ["Pair", pair], traits: ["Traits", traits],
  expression: ["Expression", expression], arrays: ["Arrays", arrays] };
