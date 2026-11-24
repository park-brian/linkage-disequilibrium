// Reads on Parquet files sorted by a key column: point lookups (sortedLookup) and ranges (rangeScan) through one read
// path (matchRows: row-group statistics, the key column's page index, then each column on its own pages), plus the
// AsyncBuffer wrappers the app reads the release through (cached, coalesce).
import { parquetRead, readColumnIndex, readOffsetIndex } from "hyparquet";
import { memo } from "./source.js";

/** Run async tasks with at most `limit` in flight (bounds memory; browsers cap concurrent requests anyway). */
export async function pool(items, limit, fn) {
  const out = new Array(items.length); let next = 0;
  await Promise.all(Array.from({ length: Math.min(limit, items.length) }, async () => {
    while (next < items.length) { const i = next++; out[i] = await fn(items[i], i); }
  }));
  return out;
}

/** Absolute row number of each returned row (rows are returned in file order). */
export const ROW = Symbol("row");
const cmp = (a, b) => (a < b ? -1 : a > b ? 1 : 0);
// first index i with keys[i] >= lo
function lowerBound(keys, lo) { let a = 0, b = keys.length; while (a < b) { const m = (a + b) >> 1; if (keys[m] < lo) a = m + 1; else b = m;
  } return a; }
const anyInRange = (keys, lo, hi) => { const i = lowerBound(keys, lo); return i < keys.length && keys[i] <= hi; };
const columnAt = (metadata, name) => metadata.row_groups[0].columns.findIndex(c => c.meta_data.path_in_schema.join(".") === name);
const view = async (file, start, length) => new DataView(await file.slice(start, start + length));
/** push row span [s, e) onto spans, merged with the last one when within gap rows */
const addSpan = (spans, s, e, gap) => {
  const l = spans.at(-1);
  if (l && s - l[1] <= gap) l[1] = Math.max(l[1], e); else spans.push([s, e]);
};

/**
 * The read path under sortedLookup and rangeScan, on a file sorted by `column`:
 *   1. row groups whose key-column statistics may hold a match (`may(min, max)`);
 *   2. their key-column page indexes (column + offset index, all groups' fetched together) → the key pages that may;
 *   3. the key column alone over those pages → the exact matching rows (`match(value)`), every one, file order;
 *   4. every other requested column on its OWN pages (its offset index): column pages are sized in bytes, so a key
 *      page can span thousands of rows of a wide column, and only the pages holding matched rows are read; or, with
 *      whole (a range, whose rows are contiguous), all columns at once over the key pages.
 * → { rows: [{ ...columns in the requested order, [ROW] }] in file order, reads: row spans read, rowsRead }
 */
async function matchRows(o) {
  const { file, metadata, column, compressors, may, match, whole = false, maxGapRows = 0, concurrency = 16, utf8 = false } = o;
  let { columns } = o;
  columns = columns.includes(column) ? columns : [column, ...columns];
  const k = columnAt(metadata, column), schemaEl = metadata.schema.find(s => s.name === column), groups = [];
  let start = 0;
  for (const rg of metadata.row_groups) {
    const st = rg.columns[k].meta_data.statistics, n = Number(rg.num_rows);
    if (may(st.min_value ?? st.min, st.max_value ?? st.max)) groups.push({ rg, start, n });
    start += n;
  }
  if (!groups.length) return { rows: [], reads: 0, rowsRead: 0 };
  const indexes = await pool(groups, concurrency, async ({ rg }) => {
    const c = rg.columns[k];
    if (!c.column_index_length || !c.offset_index_length) return null; // read the group whole
    const ciStart = Number(c.column_index_offset), oiStart = Number(c.offset_index_offset);
    const ciEnd = ciStart + c.column_index_length, oiEnd = oiStart + c.offset_index_length;
    // writers store all column indexes, then all offset indexes: one span only when the two are close
    if (Math.max(ciStart, oiStart) - Math.min(ciEnd, oiEnd) <= 64 << 10) {
      const lo = Math.min(ciStart, oiStart), v = await view(file, lo, Math.max(ciEnd, oiEnd) - lo);
      return { ci: readColumnIndex({ view: v, offset: ciStart - lo }, schemaEl), oi: readOffsetIndex({ view: v, offset: oiStart - lo }) };
    }
    const [ci, oi] = await Promise.all([view(file, ciStart, c.column_index_length), view(file, oiStart, c.offset_index_length)]);
    return { ci: readColumnIndex({ view: ci, offset: 0 }, schemaEl), oi: readOffsetIndex({ view: oi, offset: 0 }) };
  });
  const spans = [];
  groups.forEach(({ start, n }, g) => {
    if (!indexes[g]) return addSpan(spans, start, start + n, maxGapRows);
    const { ci, oi } = indexes[g], locs = oi.page_locations;
    for (let p = 0; p < locs.length; p++) {
      if (!may(ci.min_values[p], ci.max_values[p])) continue;
      const end = p + 1 < locs.length ? Number(locs[p + 1].first_row_index) : n;
      addSpan(spans, start + Number(locs[p].first_row_index), start + end, maxGapRows);
    }
  });
  if (whole) { // a range: its rows are contiguous, so every column's pages over the key pages are needed anyway
    const out = [];
    await pool(spans, concurrency, ([rowStart, rowEnd]) => parquetRead({ file, metadata, columns, compressors, rowStart, rowEnd,
      rowFormat: "object", useOffsetIndex: true, utf8,
      onComplete: rows => rows.forEach((r, i) => { if (match(r[column])) out.push({ ...r, [ROW]: rowStart + i }); }) }));
    return { rows: out.sort((a, b) => a[ROW] - b[ROW]), reads: spans.length, rowsRead: spans.reduce((a, [s, e]) => a + e - s, 0) };
  }
  // one column over [rowStart, rowEnd) as decoded arrays (onChunk: no row objects); a chunk may run past the range
  const read = (c, [rowStart, rowEnd], each) => parquetRead({ file, metadata, columns: [c], compressors, rowStart, rowEnd,
    useOffsetIndex: true, utf8, onChunk: ({ columnData, rowStart: s }) => {
      for (let i = Math.max(0, rowStart - s), end = Math.min(columnData.length, rowEnd - s); i < end; i++) each(s + i, columnData[i]);
    } });
  const hits = new Map(); // row → { column: value, … }
  await pool(spans, concurrency, sp => read(column, sp, (row, v) => { if (match(v)) hits.set(row, { [column]: v }); }));
  const rowNums = [...hits.keys()].sort((a, b) => a - b);
  let reads = spans.length, rowsRead = spans.reduce((a, [s, e]) => a + e - s, 0);
  const others = columns.filter(c => c !== column);
  if (rowNums.length && others.length) {
    // the groups holding matches, with their matched rows
    const touched = []; // one pass: groups and matched rows are both in file order
    for (let g = 0, r = 0; g < groups.length && r < rowNums.length; g++) {
      const { start, n } = groups[g], rows = [];
      while (r < rowNums.length && rowNums[r] < start + n) rows.push(rowNums[r++]);
      if (rows.length) touched.push({ ...groups[g], rows });
    }
    await Promise.all(others.map(async c => {
      const ci = columnAt(metadata, c), own = [];
      const bounds = await pool(touched, concurrency, async ({ rg, start, n }) => {
        const m = rg.columns[ci];
        if (!m.offset_index_length) return [start, start + n];
        const oi = readOffsetIndex({ view: await view(file, Number(m.offset_index_offset), m.offset_index_length), offset: 0 });
        return [...oi.page_locations.map(l => start + Number(l.first_row_index)), start + n];
      });
      touched.forEach(({ rows }, t) => {
        const b = bounds[t];
        for (let r = 0, p = 0; r < rows.length; r++) { while (b[p + 1] <= rows[r]) p++; addSpan(own, b[p], b[p + 1], 0); }
      });
      reads += own.length; rowsRead += own.reduce((a, [s, e]) => a + e - s, 0);
      await pool(own, concurrency, sp => read(c, sp, (row, v) => { const h = hits.get(row); if (h) h[c] = v; }));
    }));
  }
  const rows = rowNums.map(r => { const h = hits.get(r), x = {}; for (const c of columns) x[c] = h[c]; x[ROW] = r; return x; });
  return { rows, reads, rowsRead };
}

/**
 * Point lookup on a Parquet file sorted by `column`: every row whose key is one of `keys` (unsorted, duplicates
 * allowed; matched by String value, so BigInt and Number keys compare alike), file order. utf8 = false: hyparquet
 * decodes BYTE_ARRAY columns WITHOUT a UTF8/STRING annotation as text by default, which corrupts binary columns
 * (genotype `haps` in the sparse tier, `missing` bitsets); annotated string columns stay strings.
 * → { rows, reads, rowsRead } (matchRows)
 */
export async function sortedLookup({ file, metadata, column, keys, columns, compressors, maxGapRows = 0, concurrency = 16, utf8 = false }) {
  const sorted = [...new Set(keys)].sort(cmp), want = new Set(sorted.map(String));
  return matchRows({ file, metadata, column, columns, compressors, maxGapRows, concurrency, utf8,
    may: (lo, hi) => anyInRange(sorted, lo, hi), match: v => want.has(String(v)) });
}

/**
 * Range read on a file sorted by `column`: rows with lo <= column <= hi, in file order. (hyparquet's `filter` option
 * prunes row groups but not pages; with large row groups that reads whole column chunks.)
 */
export async function rangeScan({ file, metadata, column, lo, hi, columns, compressors, concurrency = 16, utf8 = false }) {
  return (await matchRows({ file, metadata, column, columns, compressors, concurrency, utf8,
    whole: true,
    may: (mn, mx) => !(mx < lo || mn > hi), match: v => v >= lo && v <= hi })).rows;
}

/**
 * hyparquet's cachedAsyncBuffer (a file under minSize is read whole once; larger files cache each slice by range),
 * except that a failed read is not kept: the next read of that range tries again (one network error would otherwise
 * poison the range for the session).
 */
export function cached(file, { minSize = 1 << 19 } = {}) {
  const m = new Map(), n = file.byteLength;
  if (n < minSize) return { byteLength: n, slice: async (s, e) => (await memo(m, "all", () => file.slice(0, n))).slice(s, e) };
  return { byteLength: n, slice: (s, e = n) => memo(m, `${s}:${e}`, () => file.slice(s, e)) };
}

/**
 * An AsyncBuffer that merges reads issued together: slices requested in the same task are sorted and any within `gap`
 * bytes of each other (up to `maxSpan` per request) become one read of the outer range. Over HTTP this turns the
 * page-by-page reads of a pruned lookup into a few Range requests; the bytes each caller gets are unchanged.
 */
export function coalesce(file, { gap = 32 << 10, maxSpan = 8 << 20 } = {}) {
  let queue = [], timer = null;
  function flush() {
    const q = queue.sort((a, b) => a.start - b.start); queue = []; timer = null;
    for (let i = 0; i < q.length;) {
      const start = q[i].start; let end = q[i].end, j = i + 1;
      while (j < q.length && q[j].start - end <= gap && Math.max(end, q[j].end) - start <= maxSpan) end = Math.max(end, q[j++].end);
      const group = q.slice(i, j);
      Promise.resolve(file.slice(start, end)).then(buf => { for (const r of group) r.resolve(buf.slice(r.start - start, r.end - start)); },
        err => group.forEach(r => r.reject(err)));
      i = j;
    }
  }
  return {
    byteLength: file.byteLength,
    slice(start, end = file.byteLength) {
      return new Promise((resolve, reject) => { queue.push({ start, end, resolve, reject }); timer ??= setTimeout(flush, 0); });
    },
  };
}
