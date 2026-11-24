// Small collection helpers the readers share: group a list into a Map, and keep the first row per key (the row a
// Mongo find_one returns, since the release keeps the dump's order; docs/DATA.md "Order").

/** items → Map(key(item) → [items in input order]), keys in first-seen order */
export function groupBy(items, key) {
  const m = new Map();
  for (const x of items) { const k = key(x); if (!m.has(k)) m.set(k, []); m.get(k).push(x); }
  return m;
}

/** rows in file order → Map(String(row[column]) → the first row with that key) */
export function firstRows(rows, column) {
  const m = new Map();
  for (const r of rows) { const k = String(r[column]); if (!m.has(k)) m.set(k, r); }
  return m;
}
