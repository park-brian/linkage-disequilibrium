// The set as a table: one row per variant, with its position, alleles and frequency in the chosen populations, and
// annotation columns: GWAS Catalog traits at the position, the RegulomeDB rank and the genotyping arrays that carry
// it (SNPchip) are filled in once the table is up; traits in LD (LDtrait) and expression in LD (LDexpress) are
// columns the reader switches on, because they are a port run over the whole set.
// Holding is one state with the usual two gestures: a click on a row holds that variant alone (again: lets go), its
// tick box adds it to or takes it from what is held. What is held decides the answer the page shows (views.js).
import { ldtrait } from "../../src/parity/ldlink/ldtrait.js";
import { ldexpress } from "../../src/parity/ldlink/ldexpress.js";
import { h, mb, fix, download } from "./dom.js";

const SHOWN = 300; // rows drawn (a region can hold thousands; the views still use them all)
export const LIMIT = 50; // LDtrait's and LDexpress's limit on a list

/** what is held after a click on `label`: alone it replaces what was held (or lets go), with `add` it joins or leaves */
export function clickHold(pick, label, add) {
  if (add) return pick.includes(label) ? pick.filter(l => l !== label) : [...pick, label];
  return pick.length === 1 && pick[0] === label ? [] : [label];
}

export function summary(set) {
  if (set.error) return set.error;
  if (!set.rows.length) return "Enter variants or a region above.";
  const n = set.rows.length, chroms = [...new Set(set.ok.map(r => r.chrom))], pos = set.ok.map(r => r.pos);
  const where = chroms.length === 1 ? ` on chromosome ${chroms[0]}`
    + (pos.length > 1 ? `, ${mb(Math.min(...pos))} to ${mb(Math.max(...pos))}`
    : ` at ${pos[0].toLocaleString()}`) : chroms.length > 1 ? ` on ${chroms.length} chromosomes` : "";
  const out = n - set.ok.length;
  return `${n.toLocaleString()} ${set.kind === "region" ? "common variant" : "variant"}${n === 1 ? "" : "s"}${where}`
    + (out ? ` · ${out} not usable` : "");
}

/** the columns a reader can switch on: key → [header, LDlink's name, what a cell counts] */
export const EXTRA = { traits: ["Traits in LD", "LDtrait", "GWAS Catalog traits within 500 kb in LD (r² 0.1 or more)"],
  expression: ["eQTLs in LD", "LDexpress", "GTEx whole-blood eQTLs within 500 kb in LD (r² and p at 0.1)"] };

/** draw the table into `el`; cols: the switched-on EXTRA keys; on = { hold(labels), open(label, key) } → the cells to fill */
export function drawTable(el, set, pick, cols, on) {
  const cells = new Map(), shown = set.rows.slice(0, SHOWN);
  el.pick = pick; // as it is now: markTable moves the marks without a redraw
  const body = shown.map(r => {
    const c = { traits: h("td.traits"), reg: h("td.num"), arrays: h("td.num"), extra: Object.fromEntries(cols.map(k => [k, h("td.num.more",
      { onclick: e => { e.stopPropagation(); on.open(r.label, k); } })])) };
    cells.set(r, c);
    const go = e => on.hold(clickHold(el.pick, r.label, e.ctrlKey || e.metaKey || e.shiftKey));
    const tr = h("tr", { "data-label": r.label, "aria-selected": String(pick.includes(r.label)), tabIndex: r.g ? 0 : -1,
      onclick: r.g ? go : null,
      onkeydown: e => { if (r.g && (e.key === "Enter" || e.key === " ") && e.target === tr) { e.preventDefault(); go(e); } } },
    h("td.tick", r.g ? h("input", { type: "checkbox", checked: pick.includes(r.label), "aria-label": `Hold ${r.label}`,
      onclick: e => e.stopPropagation(), onchange: () => on.hold(clickHold(el.pick, r.label, true)) }) : null), h("td.id", r.label));
    if (!r.g) tr.append(h("td.why", { colSpan: 6 + cols.length }, r.status));
    else tr.append(h("td.pos", `${r.chrom}:${r.pos.toLocaleString()}`), h("td.alleles", { title: r.alleles }, r.alleles),
      h("td.num.maf", { title: `alternate allele frequency ${fix(r.af, 4)}` }, h("i", { style: `width:${Math.round(r.maf * 200)}%` }),
        fix(r.maf, 3)), c.traits, c.reg, c.arrays, ...cols.map(k => c.extra[k]));
    if (!r.g) tr.className = "out";
    return tr;
  });
  el.replaceChildren(...(set.rows.length ? [h("div.scroll", h("table",
    h("thead", h("tr", h("th.tick"), h("th", "Variant"), h("th", "Position"), h("th", "Alleles"),
      h("th.num", { title: "Minor allele frequency in the chosen populations" }, "MAF"),
      h("th", { title: "GWAS Catalog associations at this position" }, "Traits here"),
        h("th.num", { title: "RegulomeDB rank" }, "Regulome"),
      h("th.num", { title: "Genotyping arrays that include this position (SNPchip)" }, "Arrays"),
      cols.map(k => h("th.num", { title: `${EXTRA[k][2]} · ${EXTRA[k][1]}` }, EXTRA[k][0])))),
    h("tbody", body))),
  set.rows.length > SHOWN ? h("p.more", `The first ${SHOWN} of ${set.rows.length.toLocaleString()} are listed; every view uses them all.`)
    : null] : []).filter(Boolean));
  return cells;
}

/** move the marks without redrawing (and bring a newly held row into view) */
export function markTable(el, pick) {
  el.pick = pick;
  for (const tr of el.querySelectorAll("tbody tr")) { tr.setAttribute("aria-selected", String(pick.includes(tr.dataset.label)));
    const box = tr.querySelector("input"); if (box) box.checked = pick.includes(tr.dataset.label); }
  el.querySelector('tr[aria-selected="true"]')?.scrollIntoView({ block: "nearest" });
}

/** fill the annotation cells of the listed rows: keyed lookups for RegulomeDB and arrays, one range read per
 *  chromosome for the GWAS Catalog (skipped when the listed rows span more than 5 Mb of it) */
export async function annotate(ctx, build, cells, stale) {
  const by = new Map();
  for (const r of cells.keys()) if (r.g) by.set(r.chrom, [...(by.get(r.chrom) ?? []), r]);
  await Promise.all([...by].map(async ([chrom, rows]) => {
    const pos = rows.map(r => r.pos), lo = Math.min(...pos), hi = Math.max(...pos);
    const [reg, arr, gwas] = await Promise.all([ctx.regulome(build, chrom, pos), ctx.arrays(build, rows.map(r => [chrom, r.pos])),
      hi - lo <= 5_000_000 ? ctx.gwas(build, chrom, lo, hi) : null]);
    if (stale()) return;
    const traits = new Map();
    for (const g of gwas ?? []) traits.set(g.pos, (traits.get(g.pos) ?? new Set()).add(g.trait));
    for (const r of rows) {
      const c = cells.get(r), t = [...(traits.get(r.pos) ?? [])], a = (arr.get(`${chrom}:${r.pos}`) ?? []).flatMap(d => d.split("|"))
        .filter(Boolean);
      c.traits.textContent = t.length ? t.join("; ") : gwas ? "" : "–"; c.traits.title = t.join("\n");
      c.reg.textContent = reg.get(r.pos) ?? ""; c.arrays.textContent = a.length || ""; c.arrays.title = a.join("\n");
    }
  }));
}

/** fill a switched-on column: LDtrait or LDexpress over the whole set (up to LDlink's 50), counted per query variant;
 *  a cell shows how many distinct traits / genes, lists them on hover, and opens that variant's section on a click.
 *  → a download link for LDlink's text, or a note saying why there is none */
export async function annotateLD(ctx, S, set, cells, key, stale) {
  const rows = set.ok.slice(0, LIMIT), common = { snps: rows.map(r => r.label).join("\n"), pop: S.pop, build: S.build, r2_d: "r2",
    r2_d_threshold: 0.1, window: 500000 };
  const o = key === "traits" ? await ldtrait(ctx, common) : await ldexpress(ctx, { ...common, tissues: "Whole_Blood", p_threshold: 0.1 });
  if (stale()) return null;
  const by = new Map();
  if (o.text) {
    const lines = o.text.split("\n").filter(l => l.includes("\t")).map(l => l.split("\t")), head = lines[0];
    const col = Math.max(1, head.findIndex(c => /Gene Symbol|GWAS Trait/i.test(c)));
    for (const l of lines.slice(1)) by.set(l[0], (by.get(l[0]) ?? new Set()).add(l[col]));
  }
  for (const r of rows) { const c = cells.get(r)?.extra[key], names = [...(by.get(r.label) ?? [])]; if (!c) continue;
    c.textContent = names.length || "0";
      c.title = `${names.slice(0, 40).join("\n")}${names.length > 40 ? "\n…" : ""}\nClick for the table`; }
  return o.text ? download(`${EXTRA[key][1]} text`, `${EXTRA[key][1].toLowerCase()}.txt`, o.text) : h("span.aside", o.error);
}
