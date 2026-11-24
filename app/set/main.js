// The set page (app/set.html): one input, one context, one variant set, several lenses. A design prototype of the
// "one LD module" in docs/REQUIREMENTS.md beside the module-per-tab app (app/index.html), on the same core and ports.
// There is no Run button: the set is resolved as the input settles, and the lens redraws when the set, the context
// (populations, build, statistic) or the selection changes. The whole state lives in the URL's hash.
import { asyncBufferFromUrl } from "hyparquet";
import { loadKernel, KERNEL_WASM } from "../../src/core/ld/kernel.js";
import { createContext } from "../../src/parity/ldlink/context.js";
import { cached, coalesce } from "../../src/core/parquet.js";
import { releaseUrl, readManifest } from "../site.js";
import { $, h, download } from "./dom.js";
import { resolveSet } from "./set.js";
import { contextBar } from "./context.js";
import { drawTable, markSelection, annotate } from "./table.js";
import { LENSES } from "./lenses.js";

const EXAMPLES = [
  ["a list", "rs13426489 rs4442975 rs13387042 rs13412666 rs6721996 rs12465506 rs12613955 rs12620095"],
  ["a region", "chr22:30,000,000-31,000,000"],
  ["one variant", "rs4442975"],
];
const DEFAULTS = { q: EXAMPLES[0][1], pop: "CEU", build: "grch37", stat: "r2", lens: "locus", sel: "" };
const S = { ...DEFAULTS, ...Object.fromEntries(new URLSearchParams(location.hash.slice(1))) };
if (!LENSES[S.lens]) S.lens = "locus";
const sel = () => S.sel.split(",").filter(Boolean);
const status = msg => { $("#status").textContent = msg; };

const RELEASE = releaseUrl();
let ctxPromise;
function context() {
  return (ctxPromise ??= (async () => {
    const kernel = await loadKernel(new Uint8Array(await (await fetch(KERNEL_WASM)).arrayBuffer()));
    const open = async (url, byteLength) => cached(coalesce(await asyncBufferFromUrl({ url, byteLength })));
    return createContext({ root: RELEASE, open, kernel, manifest: await readManifest(RELEASE) });
  })().catch(err => { ctxPromise = undefined; throw err; }));
}

let set = { kind: "list", rows: [], ok: [] }, setRun = 0, lensRun = 0;

/** apply a change of state; redo only what depends on it */
function update(patch) {
  const before = { ...S };
  Object.assign(S, patch);
  for (const k of Object.keys(S)) if (S[k] == null || S[k] === "") delete S[k];
  S.sel ??= ""; S.q ??= "";
  history.replaceState(null, "", `${location.pathname}${location.search}#${new URLSearchParams(S)}`);
  const changed = k => before[k] !== S[k];
  if (changed("pop") || changed("build") || changed("stat")) contextBar($("#context"), S, update);
  if (changed("q") && $("#q").value !== S.q) { $("#q").value = S.q; fit(); }
  if (changed("lens")) drawNav();
  if (changed("q") || changed("pop") || changed("build")) return refreshSet();
  if (changed("sel")) markSelection($("#set"), sel());
  return drawLens();
}

async function refreshSet() {
  const run = ++setRun, t0 = performance.now();
  status("reading…"); $("#set").classList.add("busy");
  let next;
  try { next = await resolveSet(await context(), S); } catch (err) { console.error(err); next = { kind: "list", rows: [], ok: [],
    error: String(err?.message ?? err) }; }
  if (run !== setRun) return;
  set = next; $("#set").classList.remove("busy");
  const { cells } = drawTable($("#set"), set, sel(), labels => update({ sel: labels.join(",") }));
  status(`${set.rows.length.toLocaleString()} variants in ${Math.round(performance.now() - t0)} ms`);
  context().then(ctx => annotate(ctx, S.build, cells, () => run !== setRun)).catch(err => console.error("annotations:", err));
  drawLens();
}

function drawNav() {
  $("#lenses").replaceChildren(...Object.entries(LENSES).map(([key, [name]]) => h("a", { href: "#", "aria-current": key === S.lens ? "page"
    : null, onclick: e => { e.preventDefault(); update({ lens: key }); } }, name)));
}

async function drawLens() {
  const run = ++lensRun, stale = () => run !== lensRun || setRun !== drawLens.set, picked = sel(), t0 = performance.now();
  drawLens.set = setRun;
  const figure = $("#figure"), head = $("#head"), foot = $("#foot"), links = h("span.links"),
    title = h("span.title"), caption = h("p.caption"),
    strip = h("p.strip");
  const inSet = set.ok.filter(r => picked.includes(r.label)), put = (...nodes) => { if (stale()) return;
    figure.classList.remove("busy"); figure.replaceChildren(...nodes.flat().filter(Boolean)); head.replaceChildren(title, links);
    foot.replaceChildren(strip, caption); };
  const env = { S, set, sel: picked, picked: inSet.length >= 2, active: inSet.length >= 2 ? inSet : set.ok, figure, stale, update,
    select: labels => update({ sel: labels.join(",") }),
    say: msg => put(h("p.say", msg)), show: put,
    headline: (text, ...controls) => title.replaceChildren(text, ...controls.filter(Boolean)),
    exact: (label, name, text) => links.append(download(label, name, text)),
    caption: (...nodes) => caption.replaceChildren(...nodes.filter(Boolean)), strip: (...nodes) => strip.replaceChildren(...nodes) };
  if (set.error || !set.rows.length) return put(h("p.say", set.error ?? "Enter variants or a region to begin."));
  figure.classList.add("busy");
  try {
    env.ctx = await context();
    await LENSES[S.lens][1](env);
    if (!stale()) status(`${LENSES[S.lens][0].toLowerCase()} in ${Math.round(performance.now() - t0)} ms`);
  } catch (err) { console.error(err); put(h("p.say.error", String(err?.message ?? err))); }
}

// the input: grows with its text; the set follows once typing pauses (or at once on blur and on Ctrl+Enter)
const q = $("#q");
const fit = () => { q.style.height = "auto"; q.style.height = `${Math.min(q.scrollHeight, 132)}px`; };
let timer;
const commit = () => { clearTimeout(timer); if (q.value !== S.q) update({ q: q.value, sel: "" }); };
q.value = S.q; fit();
q.addEventListener("input", () => { fit(); clearTimeout(timer); timer = setTimeout(commit, 700); });
q.addEventListener("blur", commit);
q.addEventListener("keydown", e => { if (e.key === "Enter" && (e.ctrlKey || e.metaKey
  || !q.value.includes("\n") && !/\s/.test(q.value.trim()))) {
  e.preventDefault(); commit(); } });
$("#examples").append("Try ", ...EXAMPLES.flatMap(([name, text], i) => [i ? ", " : "", h("a", { href: "#",
  onclick: e => { e.preventDefault(); update({ q: text.replaceAll(" ", "\n"), sel: "", lens: "locus" }); } }, name)]), ".");
let resize;
window.addEventListener("resize", () => { clearTimeout(resize); resize = setTimeout(drawLens, 250); });

window.addEventListener("hashchange", () => location.reload()); // a pasted or edited link: start from its state
contextBar($("#context"), S, update);
drawNav();
refreshSet();
