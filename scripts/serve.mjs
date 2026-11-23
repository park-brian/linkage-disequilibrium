// Local dev server for the browser app: static files with HTTP Range support (what hyparquet's asyncBufferFromUrl
// needs), mapping four roots — the repo (app/, src/, node_modules/), the release, the population panels and LDlink's LDassoc
// example files.
// Local only (binds 127.0.0.1). usage: node scripts/serve.mjs [port=8095] [release=experiments/out/release_v4]
import { createServer } from "node:http";
import { createReadStream, statSync } from "node:fs";
import { extname, join, normalize, resolve, sep } from "node:path";

const here = resolve(new URL(".", import.meta.url).pathname.replace(/^\/([A-Za-z]:)/, "$1"), "..");
const [port = "8095", release = join(here, "experiments/out/release_v4")] = process.argv.slice(2);
const MOUNTS = [ // URL prefix → directory (first match wins)
  ["/release/", resolve(release)],
  ["/panels/", resolve(here, "../db/data/1000G_population_samples")],
  ["/examples/", resolve(here, "../ldlink-local/data/ldassoc_example")], // LDlink's bundled LDassoc example files
  ["/", here],
];
const TYPES = { ".html": "text/html; charset=utf-8", ".js": "text/javascript", ".mjs": "text/javascript", ".css": "text/css",
  ".json": "application/json", ".wasm": "application/wasm", ".txt": "text/plain; charset=utf-8", ".parquet": "application/octet-stream", ".svg": "image/svg+xml" };

createServer((req, res) => {
  const url = decodeURIComponent(new URL(req.url, "http://x").pathname);
  const [prefix, dir] = MOUNTS.find(([p]) => url.startsWith(p));
  let path = normalize(join(dir, url.slice(prefix.length)));
  if (!path.startsWith(dir + sep) && path !== dir) { res.writeHead(403).end(); return; } // no escaping the mount
  let st;
  try { st = statSync(path); if (st.isDirectory()) { path = join(path, "index.html"); st = statSync(path); } } catch { res.writeHead(404).end("not found"); return; }
  const headers = { "Content-Type": TYPES[extname(path)] ?? "application/octet-stream", "Accept-Ranges": "bytes", "Cache-Control": "no-cache" };
  const m = /^bytes=(\d*)-(\d*)$/.exec(req.headers.range ?? "");
  if (m) {
    let start = m[1] === "" ? Math.max(0, st.size - +m[2]) : +m[1], end = m[1] === "" || m[2] === "" ? st.size - 1 : Math.min(+m[2], st.size - 1);
    if (start > end || start >= st.size) { res.writeHead(416, { "Content-Range": `bytes */${st.size}` }).end(); return; }
    res.writeHead(206, { ...headers, "Content-Range": `bytes ${start}-${end}/${st.size}`, "Content-Length": end - start + 1 });
    if (req.method === "HEAD") return res.end();
    createReadStream(path, { start, end }).pipe(res);
  } else {
    res.writeHead(200, { ...headers, "Content-Length": st.size });
    if (req.method === "HEAD") return res.end();
    createReadStream(path).pipe(res);
  }
}).listen(+port, "127.0.0.1", () => console.log(`serving on http://127.0.0.1:${port}/app/ (release: ${MOUNTS[0][1]})`));
