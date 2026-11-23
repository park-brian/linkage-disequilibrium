// Where the release lives. In the repo it is the dev server's /release/; tools/build_site.js rewrites DATA for a
// hosted site (the release on its own host, served with Range and CORS). ?data=<url> overrides either.
export const DATA = null; // set by tools/build_site.js --data

/** the release URL (no trailing slash): ?data=, else DATA, else ../release beside app/ */
export function releaseUrl() {
  const q = new URLSearchParams(location.search).get("data") ?? DATA;
  return new URL(q ?? "../release", location.href).href.replace(/\/$/, "");
}

/** the release's manifest, or an error that names the URL (a wrong ?data= is the usual cause) */
export async function readManifest(root) {
  let res;
  try { res = await fetch(`${root}/manifest.json`); } catch (e) { throw new Error(`Could not reach the release at ${root} (${e.message}). Is it served with CORS?`); }
  if (!res.ok) throw new Error(`No release at ${root}: manifest.json returned HTTP ${res.status}.`);
  return res.json();
}
