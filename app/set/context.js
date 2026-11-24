// The context bar: populations, genome build and statistic, set once for every answer (LDlink asks for them on each
// module's form), and an optional second population every answer is then shown in as well.
// Populations open as a panel grouped by 1000 Genomes super-population; choosing a super-population
// replaces its populations and "All" replaces everything, so the choice always reads as a short union ("CEU + YRI").
import { h } from "./dom.js";

export const GROUPS = {
  AFR: ["African", { YRI: "Yoruba in Ibadan, Nigeria", LWK: "Luhya in Webuye, Kenya", GWD: "Gambian in Western Divisions",
    MSL: "Mende in Sierra Leone", ESN: "Esan in Nigeria", ASW: "African ancestry in SW USA", ACB: "African Caribbean in Barbados" }],
  AMR: ["American", { MXL: "Mexican ancestry in Los Angeles", PUR: "Puerto Rican in Puerto Rico", CLM: "Colombian in Medellín",
    PEL: "Peruvian in Lima" }],
  EAS: ["East Asian", { CHB: "Han Chinese in Beijing", JPT: "Japanese in Tokyo", CHS: "Southern Han Chinese",
    CDX: "Chinese Dai in Xishuangbanna", KHV: "Kinh in Ho Chi Minh City" }],
  EUR: ["European", { CEU: "Utah residents, N and W European ancestry", TSI: "Toscani in Italy", FIN: "Finnish in Finland",
    GBR: "British in England and Scotland", IBS: "Iberian in Spain" }],
  SAS: ["South Asian", { GIH: "Gujarati Indian in Houston", PJL: "Punjabi in Lahore", BEB: "Bengali in Bangladesh",
    STU: "Sri Lankan Tamil in the UK", ITU: "Indian Telugu in the UK" }],
};
const BUILDS = [["grch37", "GRCh37"], ["grch38", "GRCh38"], ["grch38_high_coverage", "GRCh38 HC"]];
const superOf = p => Object.keys(GROUPS).find(g => p in GROUPS[g][1]);

/** the next choice after a click on `code`: toggles it, keeping the union minimal and never empty */
export function togglePop(pops, code) {
  if (code === "ALL") return ["ALL"];
  let next = pops.filter(p => p !== "ALL");
  if (next.includes(code)) next = next.filter(p => p !== code);
  else if (GROUPS[code]) next = [...next.filter(p => superOf(p) !== code), code];
  else next = [...next.filter(p => p !== superOf(code)), code];
  return next.length ? next : pops;
}

function segmented(label, options, value, onPick) {
  return h("div.seg", { role: "group", "aria-label": label }, options.map(([v, text]) => h("button", { type: "button",
    "aria-pressed": String(v === value), onclick: () => onPick(v) }, text)));
}

const keep = {}, closers = {};
/** a population button and its panel. `empty`: the choice may be cleared (the comparison), and is one population or group */
function picker(key, pops, text, title, onPick, empty) {
  const on = p => pops.includes(p) || pops.includes("ALL") || pops.includes(superOf(p));
  const next = code => (empty ? (pops.length === 1 && pops[0] === code ? [] : [code]) : togglePop(pops, code));
  const chip = (code, tip, cls = "") => h(`button.pop${cls}`, { type: "button", title: tip, "aria-pressed": String(on(code)),
    onclick: () => onPick(next(code)) }, code);
  const panel = h("div.panel", { hidden: true },
    h("div.all", chip("ALL", "All 26 populations (2,504 people)"), h("span.hint", empty
      ? "Show every answer a second time, in this population. Click it again to stop comparing."
      : "Choose any mix; the answers use the union of their samples.")),
    h("div.groups", Object.entries(GROUPS).map(([g, [name, subs]]) => h("div.group", chip(g, `${name}: every population below`, ".super"),
      h("span.gname", name), Object.entries(subs).map(([code, tip]) => chip(code, tip))))));
  const button = h(`button.pick${empty && !pops.length ? ".quiet" : ""}`, { type: "button", "aria-expanded": "false", title,
    onclick: e => { e.stopPropagation();
      panel.hidden = !panel.hidden; button.setAttribute("aria-expanded", String(!panel.hidden)); } }, text);
  // a click inside the panel redraws the bar (the choice changed): `keep` reopens the new panel and stops the click that
  // is still bubbling from closing it
  const close = e => { if (!keep[key] && !panel.hidden && !panel.contains(e.target)) { panel.hidden = true;
    button.setAttribute("aria-expanded", "false"); } };
  document.removeEventListener("click", closers[key]); document.addEventListener("click", closers[key] = close);
  if (keep[key]) { panel.hidden = false; button.setAttribute("aria-expanded", "true"); }
  panel.addEventListener("click", () => { keep[key] = true; setTimeout(() => { keep[key] = false; }); }, true);
  return h("div.popwrap", button, panel);
}

/** draw the bar into `el`; S = the page state ({ pop: "CEU+YRI", pop2, build, stat }), update(patch) applies a change */
export function contextBar(el, S, update) {
  const pops = S.pop.split("+"), other = S.pop2 ? [S.pop2] : [];
  el.replaceChildren(picker("pop", pops, pops.join(" + "), "Populations", next => update({ pop: next.join("+") })),
    picker("pop2", other, other.length ? `vs ${other[0]}` : "compare", "Compare with a second population",
      next => update({ pop2: next[0] ?? "" }), true),
    segmented("Genome build", BUILDS, S.build, v => update({ build: v })),
    segmented("Statistic", [["r2", "r²"], ["d", "D′"]], S.stat, v => update({ stat: v })));
}
