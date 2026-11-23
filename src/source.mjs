// File access for a release root: one cached { file, metadata } per path, and the only test for "absent".
import { parquetMetadataAsync } from "hyparquet";

/** a missing file (local ENOENT, HTTP 404) — "not found"; anything else is a real error and must surface */
export const isNotFound = e => e?.code === "ENOENT" || e?.status === 404 || e?.statusCode === 404 || /\b404\b|not found/i.test(String(e?.message ?? ""));

/**
 * @param {string} root release root (path or URL prefix)
 * @param {(path: string) => Promise<AsyncBuffer>} open
 * @param {{ exists?: (relPath: string) => boolean }} [opts] exists: e.g. manifest membership (no request for absent files)
 * → { file(rel) → { file, metadata } | null (absent), kv(metadata) → key-value metadata object }
 */
export function createSource(root, open, { exists } = {}) {
  const cache = new Map();
  function file(rel) {
    if (exists && !exists(rel)) return Promise.resolve(null);
    if (!cache.has(rel)) cache.set(rel, (async () => {
      let f;
      try { f = await open(`${root}/${rel}`); } catch (e) { if (isNotFound(e)) return null; throw e; }
      return { file: f, metadata: await parquetMetadataAsync(f, { initialFetchSize: 64 << 10 }) };
    })());
    return cache.get(rel);
  }
  const kv = metadata => Object.fromEntries((metadata.key_value_metadata ?? []).map(x => [x.key, x.value]));
  return { root, file, kv };
}
