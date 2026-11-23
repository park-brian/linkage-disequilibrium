// The browser app: every LDlink module, run in the page on the src/ core against a release served over HTTP
// (Range requests). forms.mjs = the module table, render.mjs = outputs, region.mjs = the pyramid view's data.
import { asyncBufferFromUrl, cachedAsyncBuffer } from "hyparquet";
import { loadKernel, KERNEL_WASM } from "../src/ld/kernel.mjs";
import { createContext } from "../src/ldlink/context.mjs";
import { MODULES, buildForm } from "./forms.mjs";
import { $, show, notice } from "./render.mjs";

const RELEASE = new URL("../release", location.href).href.replace(/\/$/, ""), PANELS = new URL("../panels/", location.href).href;
const status = msg => { $("#status").textContent = msg; };

let ctxPromise;
function context() {
  return (ctxPromise ??= (async () => {
    const t0 = performance.now();
    const kernel = await loadKernel(new Uint8Array(await (await fetch(KERNEL_WASM)).arrayBuffer()));
    const open = async url => cachedAsyncBuffer(await asyncBufferFromUrl({ url }));
    const readPanel = async pop => (await (await fetch(`${PANELS}${pop}.txt`)).text()).split(/\s+/).filter(Boolean);
    const ctx = await createContext({ root: RELEASE, open, readPanel, kernel });
    status(`ready (${kernel.backend} kernel, ${Math.round(performance.now() - t0)} ms)`);
    return ctx;
  })());
}

const nav = $("#modules");
for (const name of Object.keys(MODULES)) {
  const b = document.createElement("button"); b.type = "button"; b.textContent = name;
  b.onclick = () => { for (const x of nav.children) x.removeAttribute("aria-current"); b.setAttribute("aria-current", "page"); buildForm(name, { context, status }); history.replaceState(null, "", `#${name}`); show([]); };
  nav.append(b);
}
const route = () => { const name = MODULES[location.hash.slice(1)] ? location.hash.slice(1) : "LDmatrix"; [...nav.children].find(b => b.textContent === name).click(); };
window.addEventListener("hashchange", route);
route();
status("loading core…"); context().catch(err => { console.error(err); status("failed to load the release"); show([notice(String(err), true)]); });
