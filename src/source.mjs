// File access for a release root: one cached { file, metadata } per path, and the only test for "absent".
// With a manifest's sizes the opener skips its size probe (HTTP: a HEAD request) and the footer is one exact read.
import { parquetMetadataAsync } from "hyparquet";

/** a missing file (local ENOENT, HTTP 404) — "not found"; anything else is a real error and must surface */
export const isNotFound = e => e?.code === "ENOENT" || e?.status === 404 || e?.statusCode === 404 || /\b404\b|not found/i.test(String(e?.message ?? ""));

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
    if (!cache.has(rel)) cache.set(rel, (async () => {
      const s = size?.(rel);
      let f;
      try { f = await open(`${root}/${rel}`, s?.byteLength); } catch (e) { if (isNotFound(e)) return null; throw e; }
      return { file: f, metadata: await parquetMetadataAsync(f, { initialFetchSize: s?.footer ?? 64 << 10 }) };
    })());
    return cache.get(rel);
  }
  return { root, file };
}

export const kv = metadata => Object.fromEntries((metadata.key_value_metadata ?? []).map(x => [x.key, x.value]));
