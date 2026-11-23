// Batched point lookup on a Parquet file sorted by `column`, using row-group stats + page index
// (ColumnIndex min/max + OffsetIndex) to read only the pages that can contain the keys.
// This is the access pattern the browser app needs for rsID → position and position-list queries.
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
function lowerBound(keys, lo) { let a = 0, b = keys.length; while (a < b) { const m = (a + b) >> 1; if (keys[m] < lo) a = m + 1; else b = m; } return a; }
const anyInRange = (keys, lo, hi) => { const i = lowerBound(keys, lo); return i < keys.length && keys[i] <= hi; };

/**
 * @param {{file, metadata, column: string, keys: (number|bigint)[], columns: string[], compressors?: object, maxGapRows?: number}} o
 * keys need not be sorted or unique. Returns matching rows (objects) for keys found.
 */
// utf8 = false: hyparquet decodes BYTE_ARRAY columns WITHOUT a UTF8/STRING annotation as text by default, which corrupts
// binary columns (genotype `haps` in the sparse tier, `missing` bitsets). Annotated string columns stay strings.
export async function sortedLookup({ file, metadata, column, keys, columns, compressors, maxGapRows = 0, payloadColumn, concurrency = 16, utf8 = false }) {
  if (!columns.includes(column)) columns = [column, ...columns]; // matching needs the key column (rangeScan does the same)
  keys = [...new Set(keys)].sort(cmp);
  const colIdx = metadata.row_groups[0].columns.findIndex(c => c.meta_data.path_in_schema.join(".") === column);
  const schemaEl = metadata.schema.find(s => s.name === column);
  // 1) row-group pruning with chunk statistics
  let groupStart = 0; const groups = [];
  for (const rg of metadata.row_groups) {
    const st = rg.columns[colIdx].meta_data.statistics;
    const lo = st.min_value ?? st.min, hi = st.max_value ?? st.max;
    if (anyInRange(keys, lo, hi)) groups.push({ rg, groupStart });
    groupStart += Number(rg.num_rows);
  }
  if (!groups.length) return { rows: [], ranges: 0, pagesRows: 0 };
  // 2) fetch column + offset index of the key column for each surviving group (indexes may sit next to each
  //    row group's data, so never coalesce across groups); groups without a column index are read whole
  const indexes = await pool(groups, concurrency, async ({ rg }) => {
    const c = rg.columns[colIdx];
    if (!c.column_index_length || !c.offset_index_length) return null;
    const ciStart = Number(c.column_index_offset), oiStart = Number(c.offset_index_offset);
    const ciEnd = ciStart + c.column_index_length, oiEnd = oiStart + c.offset_index_length;
    // Writers store all column indexes, then all offset indexes, at the end of the file: the two can be far
    // apart, so only fetch them as one span when the gap is small.
    if (Math.max(ciStart, oiStart) - Math.min(ciEnd, oiEnd) <= 64 << 10) {
      const lo = Math.min(ciStart, oiStart), view = new DataView(await file.slice(lo, Math.max(ciEnd, oiEnd)));
      return { ci: readColumnIndex({ view, offset: ciStart - lo }, schemaEl), oi: readOffsetIndex({ view, offset: oiStart - lo }) };
    }
    const [ciBuf, oiBuf] = await Promise.all([file.slice(ciStart, ciEnd), file.slice(oiStart, oiEnd)]);
    return { ci: readColumnIndex({ view: new DataView(ciBuf), offset: 0 }, schemaEl), oi: readOffsetIndex({ view: new DataView(oiBuf), offset: 0 }) };
  });
  // 3) page pruning → absolute row ranges
  const ranges = [];
  for (let g = 0; g < groups.length; g++) {
    const { rg, groupStart } = groups[g], n = Number(rg.num_rows);
    if (!indexes[g]) { ranges.push([groupStart, groupStart + n]); continue; }
    const { ci, oi } = indexes[g];
    const locs = oi.page_locations;
    for (let p = 0; p < locs.length; p++) {
      if (!anyInRange(keys, ci.min_values[p], ci.max_values[p])) continue;
      const s = groupStart + Number(locs[p].first_row_index), e = groupStart + (p + 1 < locs.length ? Number(locs[p + 1].first_row_index) : n);
      const last = ranges.at(-1);
      if (last && s - last[1] <= maxGapRows) last[1] = e; else ranges.push([s, e]);
    }
  }
  const want = new Set(keys.map(String));
  const out = [];
  if (!payloadColumn) {
    // 4) read only those rows (offset index lets hyparquet fetch just the pages), concurrently
    await pool(ranges, concurrency, ([rowStart, rowEnd]) => parquetRead({
      file, metadata, columns, compressors, rowStart, rowEnd, rowFormat: "object", useOffsetIndex: true, utf8,
      onComplete: rows => rows.forEach((r, i) => { if (want.has(String(r[column]))) out.push({ ...r, [ROW]: rowStart + i }); }),
    }));
    out.sort((x, y) => x[ROW] - y[ROW]); // file order (ranges complete concurrently)
    return { rows: out, ranges: ranges.length, pagesRows: ranges.reduce((a, [s, e]) => a + e - s, 0) };
  }
  // 4b) Column pages are sized in bytes, so one key-column page can span thousands of payload rows.
  //     Read just the key column over the candidate ranges to get exact row numbers, then map those rows to
  //     the payload column's own pages (its offset index) and read only those row ranges.
  const hitRows = [];
  await pool(ranges, concurrency, ([rowStart, rowEnd]) => parquetRead({
    file, metadata, columns: [column], compressors, rowStart, rowEnd, rowFormat: "object", useOffsetIndex: true, utf8,
    onComplete: rows => rows.forEach((r, i) => { if (want.has(String(r[column]))) hitRows.push(rowStart + i); }),
  }));
  hitRows.sort((a, b) => a - b);
  const payIdx = metadata.row_groups[0].columns.findIndex(c => c.meta_data.path_in_schema[0] === payloadColumn);
  const payRanges = [];
  let gStart = 0, h = 0;
  for (const rg of metadata.row_groups) {
    const n = Number(rg.num_rows), gEnd = gStart + n;
    if (h < hitRows.length && hitRows[h] < gEnd) {
      const c = rg.columns[payIdx];
      let bounds = [0, n];
      if (c.offset_index_length) {
        const view = new DataView(await file.slice(Number(c.offset_index_offset), Number(c.offset_index_offset) + c.offset_index_length));
        bounds = [...readOffsetIndex({ view, offset: 0 }).page_locations.map(l => Number(l.first_row_index)), n];
      }
      for (; h < hitRows.length && hitRows[h] < gEnd; h++) {
        const r = hitRows[h] - gStart; let p = 0;
        while (bounds[p + 1] <= r) p++;
        const s = gStart + bounds[p], e = gStart + bounds[p + 1], last = payRanges.at(-1);
        if (last && s <= last[1]) last[1] = Math.max(last[1], e); else payRanges.push([s, e]);
      }
    }
    gStart = gEnd;
  }
  await pool(payRanges, concurrency, ([rowStart, rowEnd]) => parquetRead({
    file, metadata, columns, compressors, rowStart, rowEnd, rowFormat: "object", useOffsetIndex: true, utf8,
    onComplete: rows => rows.forEach((r, i) => { if (want.has(String(r[column]))) out.push({ ...r, [ROW]: rowStart + i }); }),
  }));
  out.sort((x, y) => x[ROW] - y[ROW]);
  return { rows: out, ranges: payRanges.length, pagesRows: payRanges.reduce((a, [s, e]) => a + e - s, 0) };
}

/**
 * Range read on a file sorted by `column`: rows with lo <= column <= hi, in file order. Uses row-group statistics,
 * then the key column's page index (ColumnIndex min/max + OffsetIndex) to read only the overlapping pages of the
 * requested columns. (hyparquet's `filter` option prunes row groups but not pages; with large row groups that
 * reads whole column chunks.)
 */
export async function rangeScan({ file, metadata, column, lo, hi, columns, compressors, concurrency = 16, utf8 = false }) {
  const colIdx = metadata.row_groups[0].columns.findIndex(c => c.meta_data.path_in_schema.join(".") === column);
  const schemaEl = metadata.schema.find(s => s.name === column);
  let groupStart = 0; const groups = [];
  for (const rg of metadata.row_groups) {
    const st = rg.columns[colIdx].meta_data.statistics, mn = st.min_value ?? st.min, mx = st.max_value ?? st.max;
    if (!(mx < lo || mn > hi)) groups.push({ rg, groupStart });
    groupStart += Number(rg.num_rows);
  }
  const ranges = [];
  const idx = await pool(groups, concurrency, async ({ rg }) => {
    const c = rg.columns[colIdx];
    if (!c.column_index_length || !c.offset_index_length) return null;
    const [ciBuf, oiBuf] = await Promise.all([file.slice(Number(c.column_index_offset), Number(c.column_index_offset) + c.column_index_length),
      file.slice(Number(c.offset_index_offset), Number(c.offset_index_offset) + c.offset_index_length)]);
    return { ci: readColumnIndex({ view: new DataView(ciBuf), offset: 0 }, schemaEl), oi: readOffsetIndex({ view: new DataView(oiBuf), offset: 0 }) };
  });
  groups.forEach(({ rg, groupStart }, g) => {
    const n = Number(rg.num_rows);
    if (!idx[g]) { ranges.push([groupStart, groupStart + n]); return; }
    const { ci, oi } = idx[g], locs = oi.page_locations;
    for (let p = 0; p < locs.length; p++) {
      if (ci.max_values[p] < lo || ci.min_values[p] > hi) continue;
      const s = groupStart + Number(locs[p].first_row_index), e = groupStart + (p + 1 < locs.length ? Number(locs[p + 1].first_row_index) : n);
      const last = ranges.at(-1);
      if (last && s <= last[1]) last[1] = e; else ranges.push([s, e]);
    }
  });
  const parts = await pool(ranges, concurrency, async ([rowStart, rowEnd]) => {
    let rows = [];
    await parquetRead({ file, metadata, columns: columns.includes(column) ? columns : [column, ...columns], compressors, rowStart, rowEnd, rowFormat: "object", useOffsetIndex: true, utf8, onComplete: r => { rows = r; } });
    return rows.filter(r => r[column] >= lo && r[column] <= hi);
  });
  return parts.flat();
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
      Promise.resolve(file.slice(start, end)).then(buf => { for (const r of group) r.resolve(buf.slice(r.start - start, r.end - start)); }, err => group.forEach(r => r.reject(err)));
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
