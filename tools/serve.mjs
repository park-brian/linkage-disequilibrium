// Local dev server for the browser app: static files with HTTP Range support (what hyparquet's asyncBufferFromUrl
// needs), mapping two roots: the repo (app/, src/, node_modules/) and the release (panels and example inputs included).
// Local only (binds 127.0.0.1). usage: node tools/serve.mjs [port=8095] [release=$LD_RELEASE] [--root <dir>] [--cors] [--no-release]
//   --root dist --no-release serves a built site (tools/build_site.mjs); --cors lets another origin read the release,
//   as a hosted release must (Access-Control-Allow-Origin *, range headers exposed). npm run serve = the defaults.
import { createServer } from "node:http";
import { createReadStream, statSync } from "node:fs";
import { extname, join, normalize, resolve, sep } from "node:path";
import { REPO, RELEASE } from "./config.mjs";

const argv = process.argv.slice(2), flag = n => argv.includes(n), value = n => (argv.includes(n) ? argv[argv.indexOf(n) + 1] : undefined);
const positional = argv.filter((a, i) => !a.startsWith("--") && !(i && argv[i - 1] === "--root"));
const [port = "8095", release = RELEASE] = positional, here = resolve(REPO, value("--root") ?? "."), cors = flag("--cors");
const MOUNTS = [ // URL prefix → directory (first match wins)
  ...(flag("--no-release") ? [] : [["/release/", resolve(release)]]),
  ["/", here],
];
const CORS = { "Access-Control-Allow-Origin": "*", "Access-Control-Allow-Headers": "Range", "Access-Control-Expose-Headers": "Content-Range, Content-Length, Accept-Ranges",
  "Timing-Allow-Origin": "*" }; // lets the page measure cross-origin bytes (bench/app_requests.js)
const TYPES = { ".html": "text/html; charset=utf-8", ".js": "text/javascript", ".mjs": "text/javascript", ".css": "text/css",
  ".json": "application/json", ".wasm": "application/wasm", ".txt": "text/plain; charset=utf-8", ".parquet": "application/octet-stream", ".svg": "image/svg+xml" };

createServer((req, res) => {
  if (cors && req.method === "OPTIONS") { res.writeHead(204, CORS).end(); return; }
  const url = decodeURIComponent(new URL(req.url, "http://x").pathname);
  const [prefix, dir] = MOUNTS.find(([p]) => url.startsWith(p));
  let path = normalize(join(dir, url.slice(prefix.length)));
  if (!path.startsWith(dir + sep) && path !== dir) { res.writeHead(403, cors ? CORS : {}).end(); return; } // no escaping the mount
  let st;
  try { st = statSync(path); if (st.isDirectory()) { path = join(path, "index.html"); st = statSync(path); } } catch { res.writeHead(404, cors ? CORS : {}).end("not found"); return; }
  const headers = { "Content-Type": TYPES[extname(path)] ?? "application/octet-stream", "Accept-Ranges": "bytes", "Cache-Control": "no-cache", ...(cors ? CORS : {}) };
  const m = /^bytes=(\d*)-(\d*)$/.exec(req.headers.range ?? "");
  if (m) {
    let start = m[1] === "" ? Math.max(0, st.size - +m[2]) : +m[1], end = m[1] === "" || m[2] === "" ? st.size - 1 : Math.min(+m[2], st.size - 1);
    if (start > end || start >= st.size) { res.writeHead(416, { "Content-Range": `bytes */${st.size}`, ...(cors ? CORS : {}) }).end(); return; }
    res.writeHead(206, { ...headers, "Content-Range": `bytes ${start}-${end}/${st.size}`, "Content-Length": end - start + 1 });
    if (req.method === "HEAD") return res.end();
    createReadStream(path, { start, end }).pipe(res);
  } else {
    res.writeHead(200, { ...headers, "Content-Length": st.size });
    if (req.method === "HEAD") return res.end();
    createReadStream(path).pipe(res);
  }
}).listen(+port, "127.0.0.1", () => console.log(`serving ${here} on http://127.0.0.1:${port}/app/${MOUNTS.length > 1 ? ` (release: ${MOUNTS[0][1]})` : ""}${cors ? " with CORS" : ""}`));
