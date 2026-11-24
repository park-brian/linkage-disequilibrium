// The set page in a browser (app/set.html, the design prototype): each lens on the example list, the locus on a region
// and on one variant, a pair LDlink refuses, and a list with unusable inputs, in headless Chromium (Playwright) against
// the dev server. Fails on an uncaught page error, a console error, a lens that reports an error, or a lens whose
// figure or table is missing. It checks that the page works, not its numbers: those are the ports' (verify/ldlink).
// --data <url> reads the release from that URL (the page's ?data=) instead of the dev server's.
// usage: node verify/app/set.js [--port 8098] [--data <url>] [--headed]      (npm run verify:set)
import { spawn } from "node:child_process";
import { join } from "node:path";
import { chromium } from "playwright";
import { REPO } from "../../tools/config.js";

const args = process.argv.slice(2), port = args.includes("--port") ? +args[args.indexOf("--port") + 1] : 8098;
const data = args.includes("--data") ? args[args.indexOf("--data") + 1] : null;
const REGION = "chr22:30,000,000-31,000,000", PAIR = "rs4442975,rs13387042";
// [name, the page state (hash), what the lens must have drawn, text the page must contain]
const CASES = [
  ["list: locus", {}, "canvas"], ["list: locus, anchored", { sel: "rs4442975" }, "canvas"], ["list: matrix", { lens: "matrix" }, "canvas"],
  ["list: haplotypes", { lens: "haplotypes" }, "table.haps"], ["list: proxies", { lens: "proxies", sel: "rs4442975" }, "tbody tr"],
  ["list: prune", { lens: "prune" }, "tbody tr"], ["list: pair", { lens: "pair", sel: PAIR }, ".pops .poprow"],
  ["list: traits", { lens: "traits" }, "tbody tr"], ["list: expression", { lens: "expression" }, "p.say, tbody tr"],
  ["list: arrays", { lens: "arrays" }, "tbody tr"],
  ["region: locus, pair", { q: REGION, pop: "EUR", sel: "rs62227979,rs5753037" }, "canvas"],
  ["region: a pair LDlink refuses", { q: REGION, pop: "EUR+YRI", sel: "rs62227979,rs11412352", lens: "pair" }, "table.two",
    "LDlink does not answer"],
  ["one variant: locus", { q: "rs4442975", stat: "d" }, "canvas"],
  ["mixed inputs", { q: "rs3 rs4 banana rs99999999999 rs4442975", lens: "matrix" }, "canvas", "not an rsID or chr:pos"],
];

const server = spawn(process.execPath, [join(REPO, "tools/serve.js"), String(port), ...(data ? ["--no-release"] : [])],
  { stdio: ["ignore", "pipe", "inherit"] });
await new Promise((ok, bad) => { server.stdout.once("data", ok); server.once("exit", c => bad(new Error(`serve.js exited (${c})`))); });
const base = `http://127.0.0.1:${port}/app/set.html${data ? `?data=${encodeURIComponent(data)}` : ""}`, problems = [];
const browser = await chromium.launch({ headless: !args.includes("--headed") });
try {
  const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
  page.on("pageerror", e => problems.push(`page error: ${e.message}`));
  page.on("console", m => { if (m.type() === "error") problems.push(`console error: ${m.text().slice(0, 200)}`); });
  for (const [name, state, must, text] of CASES) {
    const t0 = Date.now(), lens = state.lens ?? "locus";
    await page.goto("about:blank");
    await page.goto(`${base}#${new URLSearchParams(state)}`);
    const done = l => document.querySelector("#status").textContent.startsWith(`${l} in `) && !document.querySelector("#figure.busy");
    await page.waitForFunction(done, lens, { timeout: 120_000 }).catch(() => problems.push(`${name}: the lens never finished`));
    await page.waitForTimeout(300); // the strip under the locus and the table's annotations follow the figure
    const got = await page.evaluate(sel => ({ drawn: !!document.querySelector(`#figure ${sel.split(", ").join(", #figure ")}`),
      error: document.querySelector("#figure .error")?.textContent ?? null, text: document.body.textContent }), must);
    if (got.error) problems.push(`${name}: ${got.error.slice(0, 200)}`);
    if (!got.drawn) problems.push(`${name}: nothing matching "${must}" was drawn`);
    if (text && !got.text.includes(text)) problems.push(`${name}: the page does not say "${text}"`);
    const cv = await page.$("#figure canvas"), b = cv && await cv.boundingBox(); // hover and click: tooltips, selection
    if (b) { for (const [fx, fy] of [[0.3, 0.2], [0.5, 0.6], [0.5, 0.9]]) await page.mouse.move(b.x + fx * b.width, b.y + fy * b.height);
      await page.mouse.click(b.x + 0.5 * b.width, b.y + 0.9 * b.height); await page.waitForTimeout(600); }
    console.log(JSON.stringify({ case: name, seconds: Math.round((Date.now() - t0) / 100) / 10 }));
  }
} finally {
  await browser.close();
  server.kill();
}
for (const p of problems) console.log(`PROBLEM ${p}`);
console.log(JSON.stringify({ cases: CASES.length, problems: problems.length }));
if (problems.length) process.exitCode = 1;
