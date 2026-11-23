// Output: each module's result rendered as LDlink prints it (text blocks → tables, LDpair's 2×2, the matrix and
// pyramid views).
import { renderMatrix } from "./views/matrix.mjs";
import { renderPyramid } from "./views/pyramid.mjs";

export const $ = (sel, el = document) => el.querySelector(sel);
export const join = v => (Array.isArray(v) ? v.join("+") : v);
export const show = nodes => $("#output").replaceChildren(...nodes);
export function notice(msg, isError = false) { const d = document.createElement("div"); d.className = `notice${isError ? " error" : ""}`; d.textContent = msg; return d; }
function card(title, ...children) { const c = document.createElement("div"); c.className = "card"; const h = document.createElement("h2"); h.textContent = title; c.append(h, ...children); return c; }
function download(name, text, label = "Download") { const b = document.createElement("button"); b.type = "button"; b.textContent = label; b.onclick = () => { const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([text], { type: "text/plain" })); a.download = name; a.click(); }; return b; }
function toolbar(...items) { const bar = document.createElement("div"); bar.className = "toolbar"; bar.append(...items); return bar; }
/** rows of cells → table (the first row is the header) */
function table(rows) {
  const t = document.createElement("table"), [head, ...body] = rows;
  t.append(document.createElement("thead"), document.createElement("tbody"));
  t.tHead.append(tr(head, "th")); for (const r of body) t.tBodies[0].append(tr(r, "td"));
  return t;
}
function tr(cells, tag) { const r = document.createElement("tr"); for (const c of cells) { const e = document.createElement(tag); e.textContent = c; r.append(e); } return r; }

/** LDlink's text outputs: tab-separated blocks (separated by blank lines / "###" rules) and free-text notes */
export function showText(o) {
  if (o.error) return show([notice(o.error, true)]);
  const out = [];
  for (const b of o.text.split(/\n\s*\n|\n#{20,}\n/)) {
    const ls = b.split("\n").filter(l => l.length);
    if (!ls.length) continue;
    if (ls.length > 1 && ls[0].includes("\t")) out.push(card("Result", toolbar(`${ls.length - 1} rows`, download("output.txt", o.text)), table(ls.map(l => l.split("\t")))));
    else out.push(notice(ls.join("\n")));
  }
  if (o.warning) out.push(notice(o.warning));
  show(out);
}

export function showPair(o) {
  if (!Array.isArray(o)) return show([notice(o.error, true)]);
  const p = o[0], s1 = p.snp1, s2 = p.snp2, c = p.two_by_two.cells;
  const grid = table([["", `${s2.rsnum} ${s2.allele_1.allele}`, `${s2.rsnum} ${s2.allele_2.allele}`, "total"],
    [`${s1.rsnum} ${s1.allele_1.allele}`, c.c11, c.c12, `${s1.allele_1.count} (${s1.allele_1.frequency})`],
    [`${s1.rsnum} ${s1.allele_2.allele}`, c.c21, c.c22, `${s1.allele_2.count} (${s1.allele_2.frequency})`],
    ["total", `${s2.allele_1.count} (${s2.allele_1.frequency})`, `${s2.allele_2.count} (${s2.allele_2.frequency})`, p.two_by_two.total]]);
  const haps = table([["Haplotype", "Count", "Frequency"], ...Object.values(p.haplotypes).map(h => [h.alleles, h.count, h.frequency])]);
  const st = table([["D′", "R²", "χ²", "p"], [p.statistics.d_prime, p.statistics.r2, p.statistics.chisq, p.statistics.p]]);
  const notes = p.corr_alleles.map(x => Object.assign(document.createElement("p"), { textContent: x, className: "note" }));
  show([card(`${s1.rsnum} (${s1.coord}) × ${s2.rsnum} (${s2.coord})`, grid), card("Haplotypes", haps), card("Statistics", st, ...notes), ...(p.warning ? [notice(p.warning)] : [])]);
}

export function showPyramid(o, a) {
  if (o.error) return show([notice(o.error, true)]);
  window.__lastPyramid = o; // debugging aid: the last region result (grids, rows, timings)
  const bar = toolbar(`chr${o.chrom}:${o.lo.toLocaleString()}–${o.hi.toLocaleString()} · ${o.rows.length.toLocaleString()} variants · ${o.evaluated.toLocaleString()} pairs evaluated, ${o.pairs.toLocaleString()} ≥ ${a.threshold} · data ${o.tData} ms, LD + pixels ${o.ms.wall} ms on ${o.workers} worker${o.workers > 1 ? "s" : ""} (${o.backend}; slowest share ${o.ms.kernel} + ${o.ms.aggregate} ms), merge ${o.ms.merge} ms`);
  const holder = document.createElement("div");
  show([card(`LD pyramid — ${a.r2_d === "r2" ? "r²" : "D′"}, ${join(a.pop)}, ${a.build} (each pixel: the strongest pair there)`, bar, holder)]);
  renderPyramid(holder, { ...o, stat: a.r2_d });
}

export async function showMatrix(o, a, ctx) {
  if (!o.ids.length) return show([notice("No variants could be used (not in dbSNP / 1000G, or not biallelic).", true)]);
  const v = (await ctx.planner.query({ build: a.build, inputs: o.ids })).variants;
  const chrom = v[0].chrom, pos = v.map(x => x.pos), lo = Math.min(...pos), hi = Math.max(...pos), pad = Math.max(5000, (hi - lo) * 0.025);
  const genes = await ctx.genes(a.build, chrom, Math.max(0, lo - pad), hi + pad);
  const grid = m => [["RS_number", ...o.ids], ...o.ids.map((id, i) => [id, ...m[i]])];
  const tsv = m => grid(m).map(r => r.join("\t")).join("\n") + "\n";
  const holder = document.createElement("div");
  show([card(`LDmatrix — D′ below the diagonal, r² above (${join(a.pop)}, ${a.build})`,
    toolbar(`${o.ids.length} variants · chr${chrom}:${lo.toLocaleString()}–${hi.toLocaleString()}`, download("r2.txt", tsv(o.r2), "Download r²"), download("dprime.txt", tsv(o.d), "Download D′")), holder),
    card("r²", table(grid(o.r2))), card("D′", table(grid(o.d)))]);
  renderMatrix(holder, { ids: o.ids, pos, chrom, r2: o.r2, d: o.d, genes, lo: lo - pad, hi: hi + pad });
}
