// The set page in a browser (app/set.html, the design prototype), in headless Chromium (Playwright) against the dev
// server. Its one rule is checked first: what is held decides the answer (one variant → proxies, two → the pair,
// several → the matrix, a region → the locus), whether it is typed or ticked. Then the other views, a switched-on
// column, a second population, a pair LDlink refuses, a list with unusable inputs, and the history (a set change,
// then Back). Fails on an uncaught page error, a console error, a view that reports an error, or one whose figure or
// table is missing. It checks that the page works, not its numbers: those are the ports' (verify/ldlink).
// --data <url> reads the release from that URL (the page's ?data=) instead of the dev server's.
// usage: node verify/app/set.js [--port 8098] [--data <url>] [--headed]      (npm run verify:set)
import { spawn } from "node:child_process";
import { join } from "node:path";
import { chromium } from "playwright";
import { REPO } from "../../tools/config.js";

const args = process.argv.slice(2), port = args.includes("--port") ? +args[args.indexOf("--port") + 1] : 8098;
const data = args.includes("--data") ? args[args.indexOf("--data") + 1] : null;
const REGION = "chr22:30,000,000-31,000,000", ONE = "rs4442975", TWO = `${ONE} rs13387042`;
// [name, the page state (hash), the view that must be current, what must have been drawn (a selector), text the page must contain]
const CASES = [
  ["one variant: its proxies", { q: ONE }, "Proxies", "#answer tbody tr", "proxies of rs4442975"],
  ["two variants: the pair", { q: TWO }, "Pair", "#answer .popgrid .poprow", "The same pair in every"],
  ["a list: the matrix", {}, "Matrix", "#answer canvas", "among 8 variants"],
  ["a region: the locus", { q: REGION, pop: "EUR" }, "Locus", "#answer canvas", "4,675 common variants"],
  ["a list, one ticked: its proxies", { pick: ONE }, "Proxies", "#answer tbody tr", "back to all 8"],
  ["a list, two ticked: the pair", { pick: TWO.replace(" ", ",") }, "Pair", "#answer table.two"],
  ["a list, three ticked: their matrix", { pick: `${TWO.replace(" ", ",")},rs6721996` }, "Matrix", "#answer canvas", "among 3 variants"],
  ["the locus beside an answer", { q: ONE }, "Proxies", "#side canvas"],
  ["haplotypes", { view: "haplotypes" }, "Haplotypes", "#answer table.haps"],
  ["prune", { view: "prune" }, "Prune", "#answer tbody tr", "Keep these"],
  ["a list read as a locus", { view: "locus" }, "Locus", "#answer canvas"],
  ["a second population", { q: ONE, pop2: "YRI" }, "Proxies", "#answer .pane + .pane .figure", "In YRI"],
  ["a region, a pair LDlink refuses", { q: REGION, pop: "EUR+YRI", pick: "rs62227979,rs11412352" }, "Pair", "#answer table.two",
    "LDlink does not answer"],
  ["mixed inputs", { q: "rs3 rs4 banana rs99999999999 rs4442975" }, "Matrix", "#answer canvas", "not an rsID or chr:pos"],
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
  // the page is settled: the answer has reported, and nothing (the set, an answer, the locus beside it) is still being read
  const settled = async name => { await page.waitForTimeout(150);
    await page.waitForFunction(() => / in \d+ ms$/.test(document.querySelector("#status").textContent) && !document.querySelector(".busy"),
      null, { timeout: 120_000 }).catch(() => problems.push(`${name}: the page never settled`)); };
  const open = async (name, state) => { await page.goto("about:blank"); await page.goto(`${base}#${new URLSearchParams(state)}`);
    await settled(name); };
  const current = () => page.evaluate(() => document.querySelector("#views a[aria-current]")?.firstChild.textContent ?? null);
  const fold = () => page.textContent(".setbox > summary");
  const timed = async (name, run) => { const t0 = Date.now(); await run();
    console.log(JSON.stringify({ case: name, seconds: Math.round((Date.now() - t0) / 100) / 10 })); };

  for (const [name, state, view, must, text] of CASES) await timed(name, async () => {
    await open(name, state);
    const got = await page.evaluate(sel => ({ drawn: !!document.querySelector(sel), error: document.querySelector(".error")?.textContent
      ?? null, text: document.body.textContent }), must);
    if (await current() !== view) problems.push(`${name}: the view is ${await current()}, not ${view}`);
    if (got.error) problems.push(`${name}: ${got.error.slice(0, 200)}`);
    if (!got.drawn) problems.push(`${name}: nothing matching "${must}" was drawn`);
    if (text && !got.text.includes(text)) problems.push(`${name}: the page does not say "${text}"`);
    for (const cv of await page.$$("canvas")) { const b = await cv.boundingBox(); // hover: tooltips
      if (b) for (const [fx, fy] of [[0.3, 0.2], [0.5, 0.6], [0.5, 0.9]]) await page.mouse.move(b.x + fx * b.width, b.y + fy * b.height); }
  });

  // holding is one state: a row click holds that variant alone, a tick adds to it, a chip's × lets go; the figure
  // beside the answer does not move while what is held changes
  await timed("holding", async () => {
    const step = "holding", side = () => page.evaluate(() => { const b = document.querySelector("#side canvas").getBoundingClientRect();
      return `${Math.round(b.top)} ${Math.round(b.height)}`; });
    await open(step, {});
    const at = await side();
    await page.click("#table tbody tr:nth-child(2) td.id"); await settled(step);
    if (await current() !== "Proxies") problems.push(`${step}: a row click did not open its proxies (${await current()})`);
    await page.click("#table tbody tr:nth-child(3) input"); await settled(step);
    if (await current() !== "Pair") problems.push(`${step}: a second tick did not open the pair (${await current()})`);
    if (await side() !== at) problems.push(`${step}: the locus beside the answer moved (${at} -> ${await side()})`);
    await page.click("#held .chip button"); await settled(step);
    if (await current() !== "Proxies") problems.push(`${step}: letting one go did not return to proxies (${await current()})`);
    await page.click("#held a"); await settled(step);
    if (await current() !== "Matrix") problems.push(`${step}: "back to all" did not return to the matrix (${await current()})`);
  });

  // a switched-on column runs LDtrait over the set; a count opens that variant with the section it came from open
  await timed("a column", async () => {
    const step = "a column";
    await open(step, { cols: "traits" });
    await page.waitForFunction(() => /^\d+$/.test(document.querySelector("#table td.more")?.textContent ?? ""), null, { timeout: 60_000 })
      .catch(() => problems.push(`${step}: the Traits in LD column was not filled`));
    await page.click("#table td.more"); await settled(step);
    if (!(await page.$("#answer details[open] table"))) await page.waitForSelector("#answer details[open] table", { timeout: 60_000 })
      .catch(() => problems.push(`${step}: the count did not open the variant's traits`));
  });

  // a set change can be undone: one variant → add its proxies (the set grows, and says so) → Back restores the one
  await timed("history", async () => {
    const step = "history";
    await open(step, { q: ONE });
    const before = await fold();
    await page.click("#answer button.go"); await settled(step);
    if (!(await fold()).includes("Undo") || await fold() === before) problems.push(`${step}: adding proxies did not change the set`);
    await page.goBack(); await settled(step);
    if (await fold() !== before) problems.push(`${step}: Back did not restore the set (${await fold()})`);
  });
} finally {
  await browser.close();
  server.kill();
}
for (const p of problems) console.log(`PROBLEM ${p}`);
console.log(JSON.stringify({ cases: CASES.length + 3, problems: problems.length }));
if (problems.length) process.exitCode = 1;
