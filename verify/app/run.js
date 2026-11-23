// G4, the app in a browser: every view run with its default inputs in headless Chromium (Playwright) against the dev
// server, the pointer moved over each figure, then the switch-mid-run and back-button paths. Fails on an uncaught page
// error, a console error, a view that reports "error", or a figure that is missing. Prints line coverage of app/
// and src/ as the page ran them (tools/coverage.js; module workers such as app/workers/pyramid.js are not measured).
// usage: node verify/app/run.js [--port 8099] [--headed] [view...]      (npm run verify:app)
import { spawn, spawnSync } from "node:child_process";
import { mkdtempSync, writeFileSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { pathToFileURL } from "node:url";
import { chromium } from "playwright";
import { REPO } from "../../tools/config.js";

const args = process.argv.slice(2), port = args.includes("--port") ? +args[args.indexOf("--port") + 1] : 8099;
const pick = args.filter((a, i) => !a.startsWith("--") && args[i - 1] !== "--port");
// each view's default form, run once; extra hash inputs where the default is not a run (LDassoc's example is opt-in)
const VIEWS = { Pyramid: "", LDmatrix: "", LDpair: "", LDpop: "", LDhap: "", LDproxy: "", SNPclip: "", SNPchip: "", LDtrait: "",
  LDexpress: "", LDassoc: "example=yes", LDherit: "", LDcorrelation: "", LDscore: "" };
const FIGURE = new Set(["Pyramid", "LDmatrix", "LDproxy", "LDassoc"]); // views that draw a canvas

const server = spawn(process.execPath, [join(REPO, "tools/serve.js"), String(port)], { stdio: ["ignore", "pipe", "inherit"] });
await new Promise((ok, bad) => { server.stdout.once("data", ok); server.once("exit", c => bad(new Error(`serve.js exited (${c})`))); });
const base = `http://127.0.0.1:${port}`, problems = [], results = [];
const browser = await chromium.launch({ headless: !args.includes("--headed") });
try {
  const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
  page.on("pageerror", e => problems.push(`page error: ${e.message}`));
  page.on("console", m => { if (m.type() === "error") problems.push(`console error: ${m.text().slice(0, 200)}`); });
  await page.coverage.startJSCoverage({ resetOnNavigation: false });
  await page.goto(`${base}/app/`);
  await page.waitForFunction(() => /^ready/.test(document.querySelector("#status").textContent), null, { timeout: 120_000 });

  const run = async (name, extra) => {
    const t0 = Date.now();
    await page.evaluate(([n, q]) => {
      document.querySelector("#output").replaceChildren();
      location.hash = `#${n}?${q}${q ? "&" : ""}run=1`;
    }, [name, extra]);
    const done = n => { const s = document.querySelector("#status").textContent; return s === "error" || s.startsWith(`${n}: `); };
    await page.waitForFunction(done, name, { timeout: 600_000 });
    await page.waitForFunction(() => document.querySelector("#output").children.length > 0, null, { timeout: 60_000 }).catch(() => {});
    await page.waitForTimeout(300); // show() may still await annotation tracks
    const status = await page.textContent("#status"), canvases = await page.$$("#output canvas");
    for (const cv of canvases) { // hover over the figure (tooltips, track lookups)
      const b = await cv.boundingBox();
      const at = [[0.2, 0.3], [0.5, 0.5], [0.8, 0.7], [0.5, 0.95]];
      if (b) for (const [fx, fy] of at) await page.mouse.move(b.x + fx * b.width, b.y + fy * b.height);
    }
    await page.mouse.move(0, 0);
    if (status === "error") problems.push(`${name}: ${(await page.textContent("#output")).slice(0, 200)}`);
    if (FIGURE.has(name) && !canvases.length) problems.push(`${name}: no figure`);
    results.push({ view: name, status, canvases: canvases.length, seconds: Math.round((Date.now() - t0) / 100) / 10 });
    console.log(JSON.stringify(results.at(-1)));
  };
  for (const [name, extra] of Object.entries(VIEWS)) if (!pick.length || pick.includes(name)) await run(name, extra);

  // a run abandoned by a module switch leaves nothing behind; back returns to the previous inputs
  const stale = await page.evaluate(async () => { // switch inside the page, the moment the run starts
    const st = () => document.querySelector("#status").textContent, wait = ms => new Promise(r => setTimeout(r, ms));
    let started = false; // a status observer sees "running…" however short the run (timer polling can miss it)
    const obs = new MutationObserver(() => {
      if (st() !== "running…" || started) return;
      started = true; obs.disconnect();
      document.querySelector('#modules a[href="#LDhap"]').click();
    });
    obs.observe(document.querySelector("#status"), { childList: true, characterData: true, subtree: true });
    location.hash = "#LDpair?var1=rs12465506&var2=rs12613955&run=1"; // uncached inputs: the run reads the release
    for (const t0 = Date.now(); !started; await wait(10)) if (Date.now() - t0 > 30_000) return { status: st(), out: "never ran" };
    await wait(3000);
    return { status: st(), out: document.querySelector("#output").textContent };
  });
  if (stale.out === "never ran" || stale.out.includes("rs12465506") || stale.status.startsWith("LDpair"))
    problems.push(`stale LDpair result shown after switching: ${JSON.stringify(stale)}`);
  await page.goBack(); await page.waitForTimeout(300);
  if (!(await page.evaluate(() => location.hash)).startsWith("#LDpair")) problems.push("back did not return to LDpair");

  // coverage: Chromium's entries, mapped from URLs to repo files, in NODE_V8_COVERAGE's format for tools/coverage.js
  const entries = (await page.coverage.stopJSCoverage()).filter(e => e.url.startsWith(base) && !e.url.includes("/node_modules/"));
  const dir = mkdtempSync(join(tmpdir(), "ldcov-app-"));
  const result = entries.map(e => ({ url: pathToFileURL(join(REPO, new URL(e.url).pathname)).href, functions: e.functions }));
  writeFileSync(join(dir, "app.json"), JSON.stringify({ result }));
  for (const include of ["app/", "src/"])
    spawnSync(process.execPath, [join(REPO, "tools/coverage.js"), dir, "--include", include], { stdio: "inherit" });
  rmSync(dir, { recursive: true, force: true });
} finally {
  await browser.close();
  server.kill();
}
for (const p of problems) console.log(`PROBLEM ${p}`);
console.log(JSON.stringify({ views: results.length, problems: problems.length }));
if (problems.length) process.exitCode = 1;
