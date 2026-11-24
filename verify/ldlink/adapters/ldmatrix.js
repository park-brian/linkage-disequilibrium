// LDmatrix: the r² and D′ matrices (GET /ldmatrix, text), cell by cell; each case × statistic is one comparison.
import { ldmatrix } from "../../../src/parity/ldlink/ldmatrix.js";
import { matrixCases, get } from "../lib.js";

export const cases = () => matrixCases().flatMap(c => ["r2", "d"].map(stat => ({ ...c, name: `${c.name} [${stat}]`, stat })));
export const ldlink = cs => Promise.all(cs.map(async c => {
  const t = await get("ldmatrix", `snps=${encodeURIComponent(c.inputs.join("\n"))}&pop=${encodeURIComponent(c.pop)}&r2_d=${c.stat}&genome_build=${c.build}`);
  if (t.trim().startsWith("{")) return { error: JSON.parse(t).error ?? t };
  const rows = t.trim().split(/\r?\n/).map(l => l.split("\t"));
  return { ids: rows[0].slice(1), cells: rows.slice(1).map(r => r.slice(1)) };
}));
const memo = new Map();
export async function ours(ctx, c) {
  const key = c.name.replace(/ \[(r2|d)\]$/, "");
  if (!memo.has(key)) memo.set(key, await ldmatrix(ctx, c));
  const o = memo.get(key);
  return { ids: o.ids, cells: c.stat === "r2" ? o.r2 : o.d };
}
/** same SNP list and order, and every cell equal (as text, or as the same number) */
export function verdict(l, o) {
  if (l.error) return { verdict: "ldlink_error", detail: { ldlink_error: l.error } };
  if (JSON.stringify(l.ids) !== JSON.stringify(o.ids)) return { verdict: "differ", detail: { only_ldlink: l.ids.filter(x => !o.ids.includes(x)), only_ours: o.ids.filter(x => !l.ids.includes(x)) } };
  const ex = [];
  l.ids.forEach((_, i) => l.ids.forEach((__, j) => { const a = l.cells[i][j], b = o.cells[i][j]; if (a !== b && !(a !== "NA" && b !== "NA" && Number(a) === Number(b))) ex.push({ i: l.ids[i], j: l.ids[j], ldlink: a, ours: b }); }));
  return ex.length ? { verdict: "differ", detail: { mismatches: ex.length, examples: ex.slice(0, 5) } } : { verdict: "identical", detail: { snps: l.ids.length } };
}
