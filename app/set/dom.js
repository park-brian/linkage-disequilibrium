// Small DOM helpers for the set page (app/set.html): element building, a text download link, and the shared colours
// the canvases draw with (the page's CSS holds the same values as custom properties).

export const $ = (sel, el = document) => el.querySelector(sel);

/** h("td.num", { title }, "0.31", child) → an element; props are assigned (on* as listeners, data-* / aria-* as attributes) */
export function h(spec, props, ...children) {
  if (props instanceof Node || Array.isArray(props) || typeof props !== "object" || props === null) { children.unshift(props); props = {}; }
  const [tag, ...classes] = spec.split("."), el = document.createElement(tag || "div");
  if (classes.length) el.className = classes.join(" ");
  for (const [k, v] of Object.entries(props)) {
    if (v == null || v === false) continue;
    if (k.startsWith("on")) el.addEventListener(k.slice(2), v);
    else if (k.includes("-")) el.setAttribute(k, v);
    else el[k] = v;
  }
  el.append(...children.flat(Infinity).filter(c => c != null && c !== false));
  return el;
}

/** a link that downloads text as a file (LDlink's exact output for the lens) */
export function download(label, name, text) {
  return h("a.exact.act", { href: "#", title: "The text LDlink's API returns for these inputs, byte for byte",
    onclick: e => { e.preventDefault(); const a = h("a", { href: URL.createObjectURL(new Blob([text], { type: "text/plain" })),
      download: name }); a.click(); URL.revokeObjectURL(a.href); } }, label);
}

/** the LD ramps: paper → the statistic's ink (r² brick, D′ blue) */
export const INK = { r2: [176, 48, 32], d: [38, 70, 140] };
export const ramp = (stat, v, a = 1) => { const [r, g, b] = INK[stat], k = Math.max(0, Math.min(1, v));
  return `rgba(${Math.round(255 - (255 - r) * k)},${Math.round(255 - (255 - g) * k)},${Math.round(255 - (255 - b) * k)},${a})`; };
export const STAT = { r2: "r²", d: "D′" };
export const mb = p => `${(p / 1e6).toFixed(3)} Mb`;
export const fix = (v, n = 3) => (v == null || Number.isNaN(+v) ? "–" : (+v).toFixed(n));

/** a button that copies text (made on the click) to the clipboard, and says so */
export function copyButton(label, make) {
  const b = h("button.act", { type: "button", onclick: async () => { try { await navigator.clipboard.writeText(await make());
    b.textContent = "Copied"; } catch { b.textContent = "Could not copy"; } setTimeout(() => { b.textContent = label; }, 1400); } }, label);
  return b;
}
/** a button that saves a canvas as a PNG */
export function imageButton(cv, name) {
  return h("button.act", { type: "button", onclick: () => cv.toBlob(blob => { const a = h("a", { href: URL.createObjectURL(blob),
    download: name }); a.click(); URL.revokeObjectURL(a.href); }) }, "Save image");
}
/** a small labelled number input (a threshold) */
export const knob = (label, value, step, onChange) => h("label.knob", label, h("input", { type: "number", value, step, min: 0, max: 1,
  onchange: e => onChange(+e.target.value) }));
/** one of LDlink's text tables (tab-separated lines, then notes) → [table, notes]; `skip` leading columns are left out */
export function textTable(text, skip = 0, limit = 500) {
  const lines = text.split("\n"), cut = lines.findIndex(l => !l.includes("\t")), end = cut < 0 ? lines.length : cut;
  const rows = lines.slice(0, end).map(l => l.split("\t").slice(skip)), notes = lines.slice(end).filter(Boolean).join(" ");
  return { count: rows.length - 1, notes, table: h("div.scroll", h("table", h("thead", h("tr", rows[0].map(c => h("th", c)))),
    h("tbody", rows.slice(1, limit + 1).map(r => h("tr", r.map((c, i) => h(i ? "td.cell" : "td.id", { title: c }, c))))))) };
}
