// Paste into the app page's console (or run with the browser pane's javascript tool) right after a fresh load:
// runs each module cold (one per page load is the honest number; this runs them in sequence, so later modules reuse
// footers the earlier ones read) and reports requests and transferred bytes per module from Resource Timing.
// usage: await (await import("/bench/app_requests.js")).measure(["LDmatrix"])
export async function measure(names) {
  const out = {}, wait = ms => new Promise(r => setTimeout(r, ms));
  while (!/^ready/.test(document.querySelector("#status").textContent)) await wait(100);
  performance.setResourceTimingBufferSize(100000);
  const entries = () => performance.getEntriesByType("resource").filter(e => /\/release\//.test(e.name));
  out.startup = { requests: entries().length, bytes: entries().reduce((s, e) => s + e.transferSize, 0) };
  for (const n of names) {
    const before = entries().length;
    [...document.querySelectorAll("#modules button")].find(b => b.textContent === n).click(); await wait(200);
    const t0 = performance.now(); document.querySelector("#form button[type=submit]").click();
    while (!document.querySelector("#status").textContent.startsWith(`${n}:`) && document.querySelector("#status").textContent !== "error") await wait(50);
    const es = entries().slice(before);
    out[n] = { ms: Math.round(performance.now() - t0), requests: es.length, bytes: es.reduce((s, e) => s + e.transferSize, 0), status: document.querySelector("#status").textContent };
  }
  return out;
}
