// The browser app: every LDlink module, run in the page on the src/ core against a release served over HTTP
// (Range requests). forms.js = the module table, render.js = outputs, region.js = the pyramid view's data.
import { asyncBufferFromUrl } from "hyparquet";
import { loadKernel, KERNEL_WASM } from "../src/core/ld/kernel.js";
import { createContext } from "../src/parity/ldlink/context.js";
import { cached, coalesce } from "../src/core/parquet.js";
import { MODULES, buildForm } from "./forms.js";
import { $, show, notice } from "./render.js";
import { releaseUrl, readManifest } from "./site.js";

const RELEASE = releaseUrl();
const status = msg => { $("#status").textContent = msg; };

let ctxPromise;
function context() {
  return (ctxPromise ??= (async () => {
    const t0 = performance.now();
    const kernel = await loadKernel(new Uint8Array(await (await fetch(KERNEL_WASM)).arrayBuffer()));
    // byteLength from the manifest: no HEAD request per file; reads issued together share Range requests
    const open = async (url, byteLength) => cached(coalesce(await asyncBufferFromUrl({ url, byteLength })));
    const manifest = await readManifest(RELEASE);
    const ctx = await createContext({ root: RELEASE, open, kernel, manifest });
    status(`ready (${kernel.backend} kernel, ${Math.round(performance.now() - t0)} ms)`);
    return ctx;
  })().catch(err => { ctxPromise = undefined; throw err; })); // a failed load is retried by the next run
}

// routes: #<Module> or #<Module>?<field>=<value>&…[&run=1] (forms.js toQuery); a module link opens a fresh form, and
// every route is a history entry (back and forward move between modules and runs' inputs)
const nav = $("#modules");
function open(name, query) {
  for (const x of nav.children) if (x.textContent === name) x.setAttribute("aria-current", "page"); else x.removeAttribute("aria-current");
  buildForm(name, { context, status, query });
  show([]);
}
for (const name of Object.keys(MODULES)) {
  nav.append(Object.assign(document.createElement("a"), { href: `#${name}`, textContent: name }));
}
const route = () => {
  const [head, q = ""] = location.hash.slice(1).split("?"), name = MODULES[head] ? head : "LDmatrix";
  if (!MODULES[head]) history.replaceState(null, "", `${location.search}#${name}`);
  open(name, new URLSearchParams(q));
};
window.addEventListener("hashchange", route);
route();
status("loading core…"); context().catch(err => { console.error(err); status("failed to load the release"); show([notice(String(err),
  true)]); });
