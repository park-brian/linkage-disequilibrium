// Build the static site for GitHub Pages: the app and the browser code it imports, nothing else from the repo.
//   dist/index.html      → redirects to app/
//   dist/app/            the app (check.html only with --dev: it calls a local LDlink)
//   dist/src/            the browser core (src/ is browser code by rule; workers resolve src/ld/* relative to app/)
//   dist/vendor/<pkg>/   hyparquet, hyparquet-compressors, fzstd, hysnappy (their browser files + LICENSE); the import
//                        maps in app/*.html are rewritten from ../node_modules/ to ../vendor/
//   dist/.nojekyll       Pages serves files as they are
// The release is not copied (34.8 GB; Pages allows 1 GB): --data names where it is served (Range + CORS), written into
// app/site.js. Without --data the site expects ../release beside app/, like the dev server.
// usage: node tools/build_site.js [--data <url>] [--out dist] [--dev]     (npm run build:site)
import { cpSync, mkdirSync, readFileSync, rmSync, writeFileSync, readdirSync, statSync } from "node:fs";
import { join } from "node:path";
import { REPO, repo } from "./config.js";

const args = process.argv.slice(2), opt = name => { const i = args.indexOf(name); return i < 0 ? undefined : args[i + 1]; };
const data = opt("--data"), out = repo(opt("--out") ?? "dist"), dev = args.includes("--dev");
if (data && !/^https?:\/\//.test(data) && !data.startsWith("/") && !data.startsWith("../")) throw new Error(`--data must be a URL or a path: ${data}`);
if (!out.startsWith(REPO) || out === REPO) throw new Error(`refusing to write outside the repo or over it: ${out}`);
const VENDOR = { hyparquet: "src", "hyparquet-compressors": "src", fzstd: "esm", hysnappy: "js" };

rmSync(out, { recursive: true, force: true });
mkdirSync(out, { recursive: true });
cpSync(repo("src"), join(out, "src"), { recursive: true });
cpSync(repo("app"), join(out, "app"), { recursive: true, filter: f => dev || !/check\.(html|mjs)$/.test(f) });
for (const [pkg, dir] of Object.entries(VENDOR)) {
  cpSync(repo(`node_modules/${pkg}/${dir}`), join(out, "vendor", pkg, dir), { recursive: true });
  cpSync(repo(`node_modules/${pkg}/LICENSE`), join(out, "vendor", pkg, "LICENSE"));
}
for (const f of readdirSync(join(out, "app")).filter(f => f.endsWith(".html"))) {
  const p = join(out, "app", f), html = readFileSync(p, "utf8"), next = html.replaceAll("../node_modules/", "../vendor/");
  if (next === html) throw new Error(`${f}: no import map to rewrite`);
  writeFileSync(p, next);
}
const site = join(out, "app", "site.js"), src = readFileSync(site, "utf8"), set = src.replace("export const DATA = null;", `export const DATA = ${JSON.stringify(data ?? null)};`);
if (set === src && data) throw new Error("app/site.js: DATA line not found");
writeFileSync(site, set);
writeFileSync(join(out, "index.html"), '<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0; url=app/"><a href="app/">LD tools</a>\n');
writeFileSync(join(out, ".nojekyll"), "");

const size = d => readdirSync(d, { withFileTypes: true }).reduce((s, e) => s + (e.isDirectory() ? size(join(d, e.name)) : statSync(join(d, e.name)).size), 0);
console.log(JSON.stringify({ out, data: data ?? "../release", dev, bytes: size(out) }));
