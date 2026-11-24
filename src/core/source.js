// File access for a release root: one cached { file, metadata } per path, and the only test for "absent".
// With a manifest's sizes the opener skips its size probe (HTTP: a HEAD request) and the footer is one exact read.
import { parquetMetadataAsync } from "hyparquet";

/** a missing file (local ENOENT, HTTP 404) — "not found"; anything else is a real error and must surface */
export const isNotFound = e => e?.code === "ENOENT" || e?.status === 404 || e?.statusCode === 404 || /\b404\b|not found/i.test(String(e?.message ?? ""));

/**
 * memoize f() under key in map; a rejected promise is dropped (when it is still the cached one), so one transient
 * failure (a network error) is retried on the next call instead of poisoning the key for the session
 */
export function memo(map, key, f) {
  if (!map.has(key)) {
    const p = (async () => f())();
    map.set(key, p);
    p.catch(() => { if (map.get(key) === p) map.delete(key); });
  }
  return map.get(key);
}
/** f() once (memo with one key) */
export const once = f => { const m = new Map(); return () => memo(m, 0, f); };
/** a file(rel) that throws, naming the file, where file(rel) gives null (an absent file is reported, not destructured) */
export const required = file => async rel => {
  const f = await file(rel);
  if (!f) throw new Error(`the release has no ${rel}`);
  return f;
};

/**
 * @param {string} root release root (path or URL prefix)
 * @param {(path: string, byteLength?: number) => Promise<AsyncBuffer>} open byteLength: the file's size when known
 * @param {{ exists?: (rel: string) => boolean, size?: (rel: string) => { byteLength: number, footer: number } | undefined }} [opts]
 *   exists: e.g. manifest membership (no request for absent files); size: the file's length and its footer's
 * → { root, file(rel) → { file, metadata } | null (absent) }
 */
export function createSource(root, open, { exists, size } = {}) {
  const cache = new Map();
  function file(rel) {
    if (exists && !exists(rel)) return Promise.resolve(null);
    return memo(cache, rel, async () => {
      const s = size?.(rel);
      let f;
      try { f = await open(`${root}/${rel}`, s?.byteLength); } catch (e) { if (isNotFound(e)) return null; throw e; }
      return { file: f, metadata: await parquetMetadataAsync(f, { initialFetchSize: s?.footer ?? 64 << 10 }) };
    });
  }
  return { root, file };
}

export const kv = metadata => Object.fromEntries((metadata.key_value_metadata ?? []).map(x => [x.key, x.value]));
