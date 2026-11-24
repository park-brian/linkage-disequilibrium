// The set as a table: one row per variant, with its position, alleles and frequency in the chosen populations, and
// annotation columns filled in after the table is up: GWAS Catalog traits at the position (LDtrait's source), the
// RegulomeDB rank (LDproxy's column) and the genotyping arrays that carry it (SNPchip). A click selects a row (one =
// the anchor variant of the locus and proxies lenses, two = a pair); Ctrl / Shift click adds to the selection.
import { h, mb, fix } from "./dom.js";

const SHOWN = 300; // rows drawn (a region can hold thousands; the lenses still use them all)

/** the next selection after a click on `label` */
export function clickSelect(sel, label, add) {
  if (add) return sel.includes(label) ? sel.filter(l => l !== label) : [...sel, label];
  return sel.length === 1 && sel[0] === label ? [] : [label];
}

function summary(set) {
  if (set.error) return set.error;
  if (!set.rows.length) return "Enter variants or a region above.";
  const n = set.rows.length, chroms = [...new Set(set.ok.map(r => r.chrom))], pos = set.ok.map(r => r.pos);
  const where = chroms.length === 1 ? ` on chromosome ${chroms[0]}, ${mb(Math.min(...pos))} to ${mb(Math.max(...pos))}`
    : chroms.length > 1 ? ` on ${chroms.length} chromosomes` : "";
  const out = n - set.ok.length;
  return `${n.toLocaleString()} ${set.kind === "region" ? "common variant" : "variant"}${n === 1 ? "" : "s"}${where}`
    + (out ? ` · ${out} not usable` : "");
}

/** draw the set into `el`; select(labels) changes the selection → { cells: Map(row → { traits, reg, arrays }) } */
export function drawTable(el, set, sel, select) {
  const cells = new Map(), shown = set.rows.slice(0, SHOWN);
  const body = shown.map(r => {
    const c = { traits: h("td.traits"), reg: h("td.num"), arrays: h("td.num") };
    cells.set(r, c);
    const tr = h("tr", { "data-label": r.label, "aria-selected": String(sel.includes(r.label)),
      onclick: e => select(clickSelect(el.sel ?? sel, r.label, e.ctrlKey || e.metaKey || e.shiftKey)) },
    h("td.id", r.label));
    if (!r.g) tr.append(h("td.why", { colSpan: 6 }, r.status));
    else tr.append(h("td.pos", `${r.chrom}:${r.pos.toLocaleString()}`), h("td.alleles", { title: r.alleles }, r.alleles),
      h("td.num.maf", { title: `alternate allele frequency ${fix(r.af, 4)}` }, h("i", { style: `width:${Math.round(r.maf * 200)}%` }),
        fix(r.maf, 3)), c.traits, c.reg, c.arrays);
    if (!r.g) tr.className = "out";
    return tr;
  });
  el.sel = sel;
  el.replaceChildren(...[h("p.summary", summary(set)), set.rows.length ? h("div.scroll", h("table",
    h("thead", h("tr", h("th", "Variant"), h("th", "Position"), h("th", "Alleles"), h("th.num", { title: "Minor allele frequency" }, "MAF"),
      h("th", { title: "GWAS Catalog associations at this position" }, "Traits"), h("th.num", { title: "RegulomeDB rank" }, "RDB"),
      h("th.num", { title: "Genotyping arrays that include this position" }, "Arrays"))),
    h("tbody", body))) : null,
  set.rows.length > SHOWN ? h("p.more", `The first ${SHOWN} of ${set.rows.length.toLocaleString()} are listed; every lens uses them all.`)
    : null].filter(Boolean));
  return { cells };
}

/** move the selection marks without redrawing (and bring a newly selected row into view) */
export function markSelection(el, sel) {
  el.sel = sel;
  for (const tr of el.querySelectorAll("tbody tr")) tr.setAttribute("aria-selected", String(sel.includes(tr.dataset.label)));
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
