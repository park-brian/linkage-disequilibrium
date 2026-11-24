// The set page (app/set.html): a design prototype of the "one LD module" in docs/REQUIREMENTS.md beside the
// module-per-tab app (app/index.html), on the same core and ports. One rule runs it: WHAT YOU HOLD DECIDES WHAT YOU
// SEE. What is held is the ticked variants, or the whole set when none is ticked; one variant opens its proxies, two
// the pair, several the matrix, a region the locus (views.js), with the other fitting views one click away and the
// locus beside every answer as context. Holding is one state, changed from the table, the chips or a figure.
// The set is read when the input is committed (Enter, a paste, leaving the box); every change of it is a history
// entry, so Back undoes "add proxies" or "keep these". The whole state lives in the URL's hash.
import { asyncBufferFromUrl } from "hyparquet";
import { loadKernel, KERNEL_WASM } from "../../src/core/ld/kernel.js";
import { createContext } from "../../src/parity/ldlink/context.js";
import { cached, coalesce } from "../../src/core/parquet.js";
import { releaseUrl, readManifest } from "../site.js";
import { $, h, download, copyButton, imageButton } from "./dom.js";
import { resolveSet } from "./set.js";
import { contextBar } from "./context.js";
import { drawTable, markTable, annotate, annotateLD, summary, EXTRA, LIMIT } from "./table.js";
import { VIEWS, FOR } from "./views.js";

const EXAMPLES = [
  ["one variant", "rs4442975"], ["two", "rs4442975 rs13387042"],
  ["a list", "rs13426489 rs4442975 rs13387042 rs13412666 rs6721996 rs12465506 rs12613955 rs12620095"],
  ["a region", "chr22:30,000,000-31,000,000"],
];
// q: the input · pick: what is ticked · view: a view chosen over the default · cols: switched-on table columns ·
// open: opened sections · table: "1" / "0" when the reader opened or closed the set's table · pop2: the comparison
const DEFAULTS = { q: EXAMPLES[2][1], pop: "CEU", pop2: "", build: "grch37", stat: "r2", pick: "", view: "", cols: "", open: "",
  table: "" };
const fromHash = () => ({ ...DEFAULTS, ...Object.fromEntries(new URLSearchParams(location.hash.slice(1))) });
let S = fromHash();
const list = k => S[k].split(",").filter(Boolean);
const status = msg => { $("#status").textContent = msg; };
const url = () => `${location.pathname}${location.search}#${new URLSearchParams(Object.entries(S).filter(([k,
  v]) => v !== DEFAULTS[k]))}`;

const RELEASE = releaseUrl();
let ctxPromise;
function context() {
  return (ctxPromise ??= (async () => {
    const kernel = await loadKernel(new Uint8Array(await (await fetch(KERNEL_WASM)).arrayBuffer()));
    const open = async (path, byteLength) => cached(coalesce(await asyncBufferFromUrl({ url: path, byteLength })));
    return createContext({ root: RELEASE, open, kernel, manifest: await readManifest(RELEASE) });
  })().catch(err => { ctxPromise = undefined; throw err; }));
}

let set = { kind: "list", rows: [], ok: [] }, setRun = 0, answerRun = 0, sideRun = 0;
/** what is held: the ticked labels, or every usable variant of the set */
const held = () => { const p = list("pick"); return p.length ? p : set.ok.map(r => r.label); };
/** its shape, which picks the answer */
function shape() {
  const n = held().length;
  return !n ? "none" : !list("pick").length && set.kind === "region" ? "region" : n === 1 ? "one" : n === 2 ? "two" : "many";
}
const view = () => { const fits = FOR[shape()]; return fits.includes(S.view) ? S.view : fits[0]; };

/** apply a change of state and redo only what depends on it; a change of the set is a history entry (note: what did it).
 *  A new shape of what is held goes back to that shape's own view. */
function update(patch, { note } = {}) {
  const before = { ...S }, was = shape();
  Object.assign(S, patch);
  const changed = (...ks) => ks.some(k => before[k] !== S[k]);
  if (changed("q") && !("pick" in patch)) S.pick = "";
  if (changed("q")) S.table = ""; // a new set folds or unfolds by its own shape
  if ((changed("q") || (changed("pick") && shape() !== was)) && !("view" in patch)) S.view = "";
  if (changed("q")) history.pushState({ note }, "", url()); else history.replaceState(history.state, "", url());
  if (changed("pop", "pop2", "build", "stat")) contextBar($("#context"), S, update);
  if (changed("q") && $("#q").value !== S.q) { $("#q").value = S.q; fit(); }
  hint();
  if (changed("q", "pop", "build")) return refreshSet();
  if (changed("cols")) drawSet();
  if (changed("pick")) { markTable($("#table"), list("pick")); drawHeld(); }
  if (changed("pick", "view")) drawNav();
  if (changed("pick", "view", "stat", "pop2", "thr", "r2t", "maft")) drawAnswer();
  if (changed("pick", "view", "stat")) drawSide();
}
/** keep a detail of the page (an opened section, the table's state) in the URL without redrawing anything */
function remember(patch) { Object.assign(S, patch); history.replaceState(history.state, "", url()); }
/** Back, Forward or an edited link: the state is whatever the URL now says */
function restore() {
  const next = fromHash();
  if (JSON.stringify(next) === JSON.stringify(S)) return;
  S = next; $("#q").value = S.q; fit(); hint();
  contextBar($("#context"), S, update); refreshSet();
}

async function refreshSet() {
  const run = ++setRun, t0 = performance.now();
  status("reading the set…"); $("#setbox").classList.add("busy");
  let next;
  try { next = await resolveSet(await context(), S); } catch (err) { console.error(err); next = { kind: "list", rows: [], ok: [],
    error: String(err?.message ?? err) }; }
  if (run !== setRun) return;
  set = next; $("#setbox").classList.remove("busy");
  status(`${set.rows.length.toLocaleString()} variants read in ${Math.round(performance.now() - t0)} ms`);
  drawSet(); drawHeld(); drawNav(); drawAnswer(); drawSide();
}

/** the set, folded to one line unless it is a list worth reading (or the reader opened it): its table, the columns to
 *  switch on, what the last step did and the way back */
function drawSet() {
  const run = setRun, cols = list("cols").filter(k => EXTRA[k]), table = h("div", { id: "table" }), links = h("span.links");
  const open = S.table ? S.table === "1" : set.kind === "list" && (set.rows.length > 2 || set.rows.length > set.ok.length);
  const did = history.state?.note ? h("span.did", history.state.note, " ", h("a", { href: "#",
    onclick: e => { e.preventDefault(); e.stopPropagation(); history.back(); } }, "Undo")) : null;
  const extras = h("p.cols", "Add a column: ", Object.entries(EXTRA).map(([k, [name, ld, about]]) => { const on = cols.includes(k);
    return h("button.act", { type: "button", "aria-pressed": String(on),
      title: set.ok.length > LIMIT ? `For the first ${LIMIT} variants (${ld}'s limit). ${about}` : about,
      onclick: () => update({ cols: (on ? cols.filter(c => c !== k) : [...cols, k]).join(",") }) }, `${on ? "✓ " : ""}${name} · ${ld}`); }),
  links);
  const details = h("details.setbox", { open }, h("summary", h("span.summary", summary(set)), did), table, set.ok.length ? extras : null);
  details.firstChild.addEventListener("click", () => setTimeout(() => remember({ table: details.open ? "1" : "0" })));
  $("#setbox").replaceChildren(details);
  const cells = drawTable(table, set, list("pick"), cols, { hold: labels => update({ pick: labels.join(",") }),
    // a count in a switched-on column: hold that variant, with the section the count came from open
    open: (label, key) => { remember({ open: [...new Set([...list("open"), key])].join(","), pick: "" }); update({ pick: label }); } });
  context().then(async ctx => {
    annotate(ctx, S.build, cells, () => run !== setRun).catch(err => console.error("annotations:", err));
    for (const k of cols) { const link = await annotateLD(ctx, S, set, cells, k, () => run !== setRun); if (link) links.append(link); }
  }).catch(err => console.error("annotations:", err));
}

/** what is held, when it is less than the set: chips (× lets one go) and the way back to everything */
function drawHeld() {
  const pick = list("pick"), n = set.ok.length;
  const how = set.kind === "region" ? "Click a variant in the figure, or tick rows in the set, to look at one, a pair or a few."
    : "Tick rows, or click a variant, to look at one, a pair or a few of them.";
  $("#held").replaceChildren(...(pick.length ? [h("span.how", "Holding"), ...pick.map(l => h("span.chip", l, h("button", { type: "button",
    "aria-label": `Let go of ${l}`, onclick: () => update({ pick: pick.filter(x => x !== l).join(",") }) }, "×"))),
  h("a", { href: "#", onclick: e => { e.preventDefault(); update({ pick: "" }); } },
    n > 1 ? `back to all ${n.toLocaleString()}` : "let go")]
    : n > 2 ? [h("span.how", how)] : []));
}

function drawNav() {
  const now = view();
  $("main").classList.toggle("wide", now === "locus" || !now);
  $("#views").replaceChildren(...FOR[shape()].map(key => h("a", { href: "#", "aria-current": key === now ? "page" : null,
    onclick: e => { e.preventDefault(); update({ view: key }); } }, VIEWS[key][0], h("small", VIEWS[key][1]))));
}

/** one drawing of a view into `holder`, for population `pop`; what was there stays (dimmed) until the new one is ready,
 *  and stale() says a newer drawing has replaced this one */
async function pane(holder, key, pop, stale, compact = false) {
  const title = h("span.title"), acts = h("span.links"), figure = h("div.figure"), caption = h("p.caption");
  const put = (...nodes) => { if (stale()) return; holder.classList.remove("busy");
    figure.replaceChildren(...nodes.flat(Infinity).filter(Boolean)); holder.replaceChildren(h("div.head", title, acts), figure, caption); };
  const env = { S: { ...S, pop }, set, labels: held(), pick: list("pick"), figure: holder, compact, stale, update, remember,
    // holding from inside the big locus keeps the locus: the figure being read does not give way to another view
    hold: labels => update({ pick: labels.join(","), ...(key === "locus" && !compact ? { view: "locus" } : {}) }),
    say: msg => put(h("p.say", msg)), show: put,
    headline: (...nodes) => title.replaceChildren(...nodes.flat().filter(Boolean)),
    exact: (label, name, text) => acts.append(download(label, name, text)), copy: (label, make) => acts.append(copyButton(label, make)),
    image: (cv, name) => acts.append(imageButton(cv, name)), caption: (...nodes) => caption.replaceChildren(...nodes.filter(Boolean)) };
  holder.classList.add("busy");
  try { env.ctx = await context(); await VIEWS[key][2](env); } catch (err) { console.error(err);
    put(h("p.say.error", String(err?.message ?? err))); }
}

/** the answer: the view for what is held, once per population (two when comparing) */
async function drawAnswer() {
  const run = ++answerRun, at = setRun, stale = () => run !== answerRun || at !== setRun, key = view(), t0 = performance.now();
  const el = $("#answer"), both = S.pop2 && S.pop2 !== S.pop ? [S.pop, S.pop2] : [S.pop];
  if (set.error || !key) { delete el.dataset.sig; return el.replaceChildren(h("p.say", set.error ?? (set.rows.length
    ? "None of these is usable in 1000 Genomes: the set says why." : "Enter variants or a region to begin."))); }
  const sig = `${key}|${both}`; // the same panes are redrawn in place while only their content changes
  if (el.dataset.sig !== sig) { el.dataset.sig = sig; el.replaceChildren(...both.map(p => h("section.pane",
    both.length > 1 ? h("p.inpop", `In ${p.replaceAll("+", " + ")}`) : null, h("div.inner")))); }
  await Promise.all(both.map((p, i) => pane(el.children[i].querySelector(".inner"), key, p, stale)));
  if (!stale()) status(`${VIEWS[key][0].toLowerCase()} in ${Math.round(performance.now() - t0)} ms`);
}

/** the locus beside the answer, as context (not when the locus is the answer, or the set has no one locus) */
async function drawSide() {
  const run = ++sideRun, at = setRun, stale = () => run !== sideRun || at !== setRun, el = $("#side");
  if (view() === "locus" || !view() || set.error) return el.replaceChildren();
  if (!el.firstChild) el.replaceChildren(h("div.inner"));
  await pane(el.firstChild, "locus", S.pop, stale, true);
}

// the input: grows with its text, and is read when it is committed: Enter, a paste, or leaving the box. Until then the
// page says the text is not read yet.
const q = $("#q");
function fit() { q.style.height = "auto"; q.style.height = `${Math.min(q.scrollHeight, 132)}px`; }
const commit = () => { if (q.value !== S.q) update({ q: q.value }); };
function hint() {
  const dirty = q.value !== S.q, many = q.value.includes("\n");
  $("#examples").replaceChildren(...(dirty ? [h("button.read", { type: "button", onclick: commit }, "Read this"),
    ` or press ${many ? "Ctrl+Enter" : "Enter"}`]
    : ["Try ", ...EXAMPLES.flatMap(([name, text], i) => [i ? ", " : "", h("a", { href: "#", onclick: e => { e.preventDefault();
      update({ q: text }); } }, name)]), ". Variants are rsIDs or chr:pos, separated by spaces."]));
}
q.value = S.q; fit();
q.addEventListener("input", () => { fit(); hint(); });
q.addEventListener("paste", () => setTimeout(commit));
q.addEventListener("blur", e => { if (!e.relatedTarget?.classList.contains("read")) commit(); });
q.addEventListener("keydown", e => { if (e.key === "Enter" && (e.ctrlKey || e.metaKey || !q.value.includes("\n"))) { e.preventDefault();
  commit(); } });
$("#link").append(copyButton("Copy link", () => location.href));
let resize;
window.addEventListener("resize", () => { clearTimeout(resize); resize = setTimeout(() => { drawAnswer(); drawSide(); }, 250); });
window.addEventListener("popstate", restore);
window.addEventListener("hashchange", restore);

contextBar($("#context"), S, update);
hint();
refreshSet();
