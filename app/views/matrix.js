// LDmatrix plot (REQUIREMENTS.md §2, LDmatrix.py's Bokeh figure): D′ below the diagonal (blue, opacity = D′), r² above
// (red, opacity = r²), the diagonal purple, NA gray; a connector from the evenly spaced columns to the variants' real
// positions; a position rug in Mb; and the refGene track (transcripts packed into rows, exons as boxes).
import { canvas, tooltip, niceStep } from "./tracks.js";
const COLORS = { r2: [205, 40, 40], d: [36, 87, 197], diag: [120, 60, 170], na: [190, 190, 190] };

export function renderMatrix(holder, { ids, pos, r2, d, genes, lo, hi }) {
  const n = ids.length, W = Math.min(holder.clientWidth || 720, 820), cell = W / n;
  // gene rows: greedy packing with an 80 kb gap (LDlink's packing rule), scaled to the plot's bp → px
  const span = Math.max(1, hi - lo), x = p => ((p - lo) / span) * W, rows = [], placed = [];
  for (const g of [...genes].sort((a, b) => a.txStart - b.txStart)) {
    let r = rows.findIndex(end => g.txStart > end + 80_000);
    if (r < 0) { r = rows.length; rows.push(0); }
    rows[r] = Math.max(rows[r], g.txEnd); placed.push({ g, r });
  }
  const H_MATRIX = W, H_CONN = 46, H_RUG = 26, ROW = 22, H_GENES = Math.max(1, rows.length) * ROW + 10, H = H_MATRIX + H_CONN + H_RUG + H_GENES;
  const { cv, g } = canvas(W, H);
  const fill = ([r, gg, b], a) => `rgba(${r},${gg},${b},${a})`;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    const v = i === j ? null : i > j ? d[i][j] : r2[i][j];
    g.fillStyle = i === j ? fill(COLORS.diag, 1) : v === "NA" ? fill(COLORS.na, 1) : fill(i > j ? COLORS.d : COLORS.r2, Math.max(0.04, +v));
    g.fillRect(j * cell, i * cell, Math.ceil(cell), Math.ceil(cell));
  }
  if (cell >= 6) { g.strokeStyle = "rgba(255,255,255,.6)"; g.lineWidth = 0.5; for (let k = 0; k <= n; k++) { g.beginPath(); g.moveTo(0, k * cell); g.lineTo(W, k * cell); g.moveTo(k * cell, 0); g.lineTo(k * cell, W); g.stroke(); } }
  // connector: column centre → genomic position
  const y0 = H_MATRIX + 2, y1 = H_MATRIX + H_CONN;
  g.strokeStyle = "rgba(60,60,60,.45)"; g.lineWidth = 0.8;
  for (let k = 0; k < n; k++) { g.beginPath(); g.moveTo((k + 0.5) * cell, y0); g.lineTo(x(pos[k]), y1); g.stroke(); }
  // rug + Mb axis
  g.strokeStyle = "#333"; for (let k = 0; k < n; k++) { g.beginPath(); g.moveTo(x(pos[k]), y1); g.lineTo(x(pos[k]), y1 + 8); g.stroke(); }
  g.fillStyle = "#555"; g.font = "11px system-ui"; g.textBaseline = "top";
  const step = niceStep(span / 5);
  let lastEnd = -Infinity;
  for (let t = Math.ceil(lo / step) * step; t <= hi; t += step) {
    const label = `${(t / 1e6).toFixed(step < 1e5 ? 3 : 2)} Mb`, w = g.measureText(label).width, lx = Math.min(W - w, Math.max(0, x(t) - w / 2));
    if (lx < lastEnd + 8) continue; // never overlap labels
    g.fillText(label, lx, y1 + 11); lastEnd = lx + w;
  }
  // genes
  const gy = H_MATRIX + H_CONN + H_RUG;
  for (const { g: t, r } of placed) {
    const yy = gy + r * ROW + 6, xs = Math.max(0, x(t.txStart)), xe = Math.min(W, x(t.txEnd));
    g.strokeStyle = "#2b5d34"; g.fillStyle = "#2b5d34"; g.lineWidth = 1;
    g.beginPath(); g.moveTo(xs, yy + 4); g.lineTo(xe, yy + 4); g.stroke();
    const st = (t.exonStarts ?? "").split(",").filter(Boolean).map(Number), en = (t.exonEnds ?? "").split(",").filter(Boolean).map(Number);
    st.forEach((s, k) => { const a = x(s), b = x(en[k]); if (b >= 0 && a <= W) g.fillRect(Math.max(0, a), yy, Math.max(1, Math.min(W, b) - Math.max(0, a)), 8); });
    g.fillStyle = "#1d2330"; g.fillText(`${t.name2}${t.strand === "-" ? " ←" : " →"}`, Math.max(0, Math.min(W - 80, xs)), yy + 9);
  }
  if (!genes.length) { g.fillStyle = "#888"; g.fillText("no RefSeq genes in this region", 4, gy + 4); }
  // hover: the pair and its values
  const tip = tooltip(cv, (mx, my) => {
    if (my >= W || mx >= W) return null;
    const i = Math.floor(my / cell), j = Math.floor(mx / cell);
    return `${ids[i]} × ${ids[j]}\nr² ${r2[i][j]}   D′ ${d[i][j]}`;
  });
  holder.replaceChildren(cv, tip);
}
