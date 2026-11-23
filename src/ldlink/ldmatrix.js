// LDmatrix on the core: the SNP labels and the r² and D′ matrices as LDlink prints them (Python str(round(x, 3)),
// "NA" when Ms = 0), with LDmatrix's module rules (server/LDmatrix.py:69-140, 186-372):
//   input gate (len > 2, starts rs/chr, ends in a digit); duplicate inputs dropped (validsnp); a position "NA" for the
//   build or chrY on GRCh38 → dropped; not in 1000G → dropped; source record multi-allelic → dropped ("not a
//   biallelic"); kept variants in position order (stable).
import { pyRoundRepr } from "../pyfmt.js";
import { gate } from "./common.js";

const cell = v => (Number.isNaN(v) ? "NA" : pyRoundRepr(v, 3));

/** → { ids: [label], r2: [[string]], d: [[string]] } */
export async function ldmatrix(ctx, { inputs, pop, build = "grch37" }) {
  const r = await ctx.planner.query({ build, inputs, select: "parity" });
  const seen = new Set(), keep = [];
  r.variants.forEach((v, i) => {
    const tok = v.token; if (!tok || seen.has(tok)) return; seen.add(tok); // validsnp drops duplicate inputs
    const label = v.rsid ?? tok; // replace_coord_rsid: a resolved coordinate becomes its rsID
    if (!gate(label) || !v.found || v.pos == null) return;
    if (v.chrom === "Y" && build !== "grch37") return;
    const s = v.selected; if (!s || !s.row || s.status === "not_biallelic") return;
    keep.push({ label, pos: v.pos, chrom: v.chrom, ref: s.row, order: i });
  });
  keep.sort((a, b) => a.pos - b.pos || a.order - b.order);
  if (!keep.length) return { ids: [], r2: [], d: [] };
  const rows = [];
  for (const k of keep) rows.push({ ...(await ctx.row(build, k.chrom, k.ref)), chrom: k.chrom });
  const p = await ctx.prepare(build, rows, pop), n = rows.length, I = [], J = [];
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { I.push(i); J.push(j); }
  const e = ctx.kernel.exact(p, I, J); ctx.kernel.reset();
  const mat = a => Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => cell(a[i * n + j])));
  return { ids: keep.map(k => k.label), r2: mat(e.r2), d: mat(e.dprime) };
}
