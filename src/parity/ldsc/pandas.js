// The pandas text ldsc prints, under its display options (ldsc.py: display.float_format '{:.4f}', wide display):
// DataFrame.to_string / str() layouts, Series.describe() and DataFrame.corr() numbers.
//   column:  numeric → header with a leading space, '{:.4f}' cells (NaN as "NaN"); string → as is
//            right-justified to the widest entry; columns joined with single spaces
//   index:   left-justified labels, a blank header (str(); to_string(index=False) has none)
import { pyRoundStr } from "../pyfmt.js";
import { npSum } from "./munge.js";

const cell = v => (typeof v !== "number" ? String(v) : Number.isNaN(v) ? "NaN" : Number.isFinite(v) ? pyRoundStr(v,
  4) : v > 0 ? "inf" : "-inf");

/** to_string(index=False) of one row: cols = [[name, value]] */
export function pdRow(cols) {
  const cells = cols.map(([name, v]) => {
    const head = typeof v === "number" ? ` ${name}` : name, text = cell(v), w = Math.max(head.length, text.length);
    return [head.padStart(w), text.padStart(w)];
  });
  return `${cells.map(c => c[0]).join(" ")}\n${cells.map(c => c[1]).join(" ")}`;
}

/** str(DataFrame) of numeric columns with a row index: columns = [name], rows = [[label, ...values]] */
export function pdFrame(columns, rows) {
  const iw = Math.max(...rows.map(r => String(r[0]).length));
  const widths = columns.map((c, j) => Math.max(` ${c}`.length, ...rows.map(r => cell(r[j + 1]).length)));
  const head = [" ".repeat(iw), ...columns.map((c, j) => ` ${c}`.padStart(widths[j]))].join(" ");
  return [head, ...rows.map(r => [String(r[0]).padEnd(iw), ...columns.map((_,
    j) => cell(r[j + 1]).padStart(widths[j]))].join(" "))].join("\n");
}

/** Series.describe() for a float column: count, mean, std (ddof 1), min, 25%, 50%, 75%, max (pandas / numpy formulas) */
export function describe(values) {
  const x = Float64Array.from(values.filter(v => !Number.isNaN(v))), n = x.length, mean = npSum(x) / n;
  const sq = Float64Array.from(x, v => (v - mean) ** 2), std = Math.sqrt(npSum(sq) / (n - 1)); // nanops.nanvar: pairwise sums
  const s = Float64Array.from(x).sort(), q = p => { const h = (n - 1) * p, lo = Math.floor(h), f = h - lo;
    return lo + 1 < n ? s[lo] + (s[lo + 1] - s[lo]) * f : s[lo]; };
  return { count: n, mean, std, min: s[0], "25%": q(0.25), "50%": q(0.5), "75%": q(0.75), max: s[n - 1] };
}

/** Pearson correlation as pandas' nancorr computes it (Welford running means and co-moments) */
export function corr(x, y) {
  let n = 0, mx = 0, my = 0, sxx = 0, syy = 0, sxy = 0;
  for (let i = 0; i < x.length; i++) {
    const a = x[i], b = y[i]; if (Number.isNaN(a) || Number.isNaN(b)) continue;
    n++; const dx = a - mx, dy = b - my; mx += dx / n; my += dy / n;
    sxx += dx * (a - mx); syy += dy * (b - my); sxy += dx * (b - my);
  }
  const d = Math.sqrt(sxx * syy);
  return d === 0 ? NaN : sxy / d;
}
