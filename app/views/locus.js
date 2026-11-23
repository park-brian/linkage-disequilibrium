// Locus plots: LDlink's LDproxy and LDassoc figures (server/LDcommon.py ldproxy_figure, server/LDassoc.py), on canvas.
// One scatter over position (Mb), a recombination line on a right-hand "Combined Recombination Rate (cM/Mb)" axis,
// optional point labels (FORGEdb / RegulomeDB), a rug strip of the same colours, then the gene track and the axis
// (tracks.js). Each module maps its rows to points; this file only draws them.
//   spec = { title, xlo, xhi, yLabel, y0, y1, points: [{ pos, y, color, alpha, size, stroke?, label?, hover }],
//            recomb: [{ pos, rate }], recombY: rate → y, hline?: y, genes, collapsed, onPoint?(point) }
import { axisOf, genesTrack, axisTrack, layout, drawTracks, hoverAt, line, pointer } from "./tracks.js";

const PLOT = 300, RUG = 22, PAD = { left: 46, right: 46, top: 22 };

export function renderLocus(holder, spec) {
  const W = Math.min(900, Math.max(600, holder.clientWidth || 900)), inner = W - PAD.left - PAD.right;
  const ax = axisOf(spec.xlo, spec.xhi, inner), X = p => PAD.left + ax.x(p);
  const Y = v => PAD.top + PLOT - ((v - spec.y0) / (spec.y1 - spec.y0)) * PLOT;
  const tracks = [genesTrack(spec.genes, ax, { collapsed: spec.collapsed }), axisTrack(ax)];
  const lay = layout(tracks, PAD.top + PLOT + 8 + RUG + 8), H = lay.end + 4, dpr = window.devicePixelRatio || 1;
  const cv = document.createElement("canvas");
  cv.width = W * dpr; cv.height = H * dpr; cv.style.width = `${W}px`; cv.style.height = `${H}px`;
  const g = cv.getContext("2d"); g.scale(dpr, dpr); g.font = "11px system-ui"; g.textBaseline = "top";
  // title and frame
  g.fillStyle = "#1d2330"; g.font = "600 12px system-ui"; g.textAlign = "center"; g.fillText(spec.title, W / 2, 2); g.textAlign = "left"; g.font = "11px system-ui";
  g.strokeStyle = "#c9ced6"; g.strokeRect(PAD.left, PAD.top, inner, PLOT);
  // y axes: the statistic (left), recombination rate −3..103 (right)
  g.fillStyle = "#555"; g.textBaseline = "middle";
  for (const t of ticks(spec.y0, spec.y1, 5)) { line(g, PAD.left - 4, Y(t), PAD.left, Y(t)); g.textAlign = "right"; g.fillText(fmt(t), PAD.left - 6, Y(t)); }
  for (const r of [0, 20, 40, 60, 80, 100]) { const yy = PAD.top + PLOT - ((r + 3) / 106) * PLOT; line(g, W - PAD.right, yy, W - PAD.right + 4, yy); g.textAlign = "left"; g.fillText(String(r), W - PAD.right + 6, yy); }
  g.save(); g.translate(10, PAD.top + PLOT / 2); g.rotate(-Math.PI / 2); g.textAlign = "center"; g.fillText(spec.yLabel, 0, 0); g.restore();
  g.save(); g.translate(W - 8, PAD.top + PLOT / 2); g.rotate(Math.PI / 2); g.textAlign = "center"; g.fillText("Combined Recombination Rate (cM/Mb)", 0, 0); g.restore();
  g.textAlign = "left"; g.textBaseline = "top";
  // clip the plot area: recombination, the significance line, points, labels
  g.save(); g.beginPath(); g.rect(PAD.left, PAD.top, inner, PLOT); g.clip();
  g.strokeStyle = "rgba(0,0,0,.5)"; g.lineWidth = 1; g.beginPath();
  spec.recomb.forEach((r, k) => { const xx = X(r.pos), yy = Y(spec.recombY(r.rate)); k ? g.lineTo(xx, yy) : g.moveTo(xx, yy); }); g.stroke();
  if (spec.hline != null) { g.strokeStyle = "rgba(0,0,255,.5)"; line(g, PAD.left, Y(spec.hline), W - PAD.right, Y(spec.hline)); } // genome-wide significance
  for (const p of spec.points) {
    g.globalAlpha = p.alpha; g.fillStyle = p.color; g.beginPath(); g.arc(X(p.pos), Y(p.y), p.size / 2, 0, 7); g.fill();
    if (p.stroke) { g.globalAlpha = 1; g.strokeStyle = p.stroke; g.stroke(); }
  }
  g.globalAlpha = 1; g.font = "9px system-ui"; g.textAlign = "center"; g.textBaseline = "middle"; g.fillStyle = "#1d2330";
  for (const p of spec.points) if (p.label != null && p.label !== "") g.fillText(String(p.label), X(p.pos), Y(p.y));
  g.restore(); g.font = "11px system-ui"; g.textAlign = "left"; g.textBaseline = "top";
  // rug strip
  const rugTop = PAD.top + PLOT + 8;
  for (const p of spec.points) { g.strokeStyle = p.color; g.globalAlpha = p.alpha; line(g, X(p.pos), rugTop, X(p.pos), rugTop + RUG); }
  g.globalAlpha = 1;
  // gene track and axis (shifted right by the left padding)
  g.save(); g.translate(PAD.left, 0); drawTracks(g, tracks, lay.tops); g.restore();

  const tip = document.createElement("div"); tip.className = "tip"; tip.hidden = true;
  const hit = (mx, my) => {
    let best = null, bd = Infinity;
    for (const p of spec.points) { const d = Math.hypot(X(p.pos) - mx, Y(p.y) - my); if (d <= Math.max(4, p.size / 2) && d < bd) { bd = d; best = p; } }
    return best;
  };
  cv.addEventListener("mousemove", e => {
    const [mx, my] = pointer(cv, e);
    let text = null;
    if (my < PAD.top + PLOT) text = hit(mx, my)?.hover ?? null;
    else if (my >= lay.tops[0]) text = hoverAt(tracks, lay, ax.lo + ((mx - PAD.left) / inner) * ax.span, my);
    cv.style.cursor = text && my < PAD.top + PLOT && spec.onPoint ? "pointer" : "default";
    if (!text) { tip.hidden = true; return; }
    tip.textContent = text; tip.style.left = `${e.clientX + 12}px`; tip.style.top = `${e.clientY + 12}px`; tip.hidden = false;
  });
  cv.addEventListener("mouseleave", () => { tip.hidden = true; });
  if (spec.onPoint) cv.addEventListener("click", e => { const p = hit(...pointer(cv, e)); if (p) spec.onPoint(p); });
  holder.replaceChildren(cv, tip);
}

function ticks(a, b, n) { const step = nice((b - a) / n), out = []; for (let t = Math.ceil(a / step) * step; t <= b + 1e-9; t += step) out.push(+t.toFixed(10)); return out; }
const nice = s => { const p = 10 ** Math.floor(Math.log10(s)), m = s / p; return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p; };
const fmt = t => (Math.abs(t) < 1e-9 ? "0" : String(+t.toFixed(2)));

/** LDlink's LDproxy figure from ldproxy(…, { plot: true }).plot */
export function proxySpec(plot, { pop, annotate = "forge", genes, recomb, collapsed = true }) {
  const stat = plot.stat === "d" ? "dprime" : "r2", r4 = v => String(Math.round(v * 1e4) / 1e4);
  const points = plot.rows.map((r, i) => ({
    pos: r.pos, y: r[stat], alpha: 0.5, size: 9 + r.maf * 14, row: r,
    color: i === 0 ? "blue" : r.func !== "none" && r.func !== "" ? "red" : "orange",
    label: annotate === "forge" ? r.forge : annotate === "regulome" ? r.regulome : null,
    hover: [`Query Variant: ${plot.query.rs} ${plot.query.alleles}`, `Proxy Variant: ${r.rs} ${r.alleles}`, `Distance (Mb): ${r4(r.dist / 1e6)}`, `MAF (Query,Proxy): ${r4(r.mafq)},${r4(r.maf)}`,
      `R²: ${r4(r.r2)}`, `D': ${r4(r.dprime)}`, `Correlated Alleles: ${r.match}`, `FORGEdb Score: ${r.forge ?? ""}`, `RegulomeDB Rank: ${r.regulome ?? ""}`, `Functional Class: ${r.func}`].join("\n"),
  }));
  const pad = 10000; // LDlink's 0.01 Mb whitespace
  return { title: `Proxies for ${plot.query.rs} in ${pop}`, xlo: plot.lo - pad, xhi: plot.hi + pad, yLabel: stat === "r2" ? "R²" : "D'", y0: -0.03, y1: 1.03,
    points, recomb, recombY: rate => rate / 100, genes, collapsed };
}

// LDassoc's 101-step ramp, #FFCCCC → #FF0000 (server/LDassoc.py), indexed by int(stat · 100)
const REDS = Array.from({ length: 101 }, (_, i) => `#FF${Math.round(0xcc * (1 - i / 100)).toString(16).padStart(2, "0").toUpperCase().repeat(2)}`);

/** LDlink's LDassoc figure from ldassoc(…, { plot: true }).plot */
export function assocSpec(plot, { pop, annotate = "forge", genes, recomb, collapsed = true }) {
  const stat = plot.stat === "d" ? "dprime" : "r2", r4 = v => String(Math.round(v * 1e4) / 1e4), q = plot.index.rs;
  const nlp = p => -Math.log10(p), ymax = Math.max(...plot.rows.map(r => nlp(r.p)), 1);
  const missing = plot.missing.map(m => ({ pos: m.pos, y: nlp(m.p), color: "#FFFFFF", stroke: "gray", alpha: 1, size: 12.5,
    hover: [`Variant: chr${plot.chrom}:${m.pos}`, `P-value: ${m.p}`, `Distance (Mb): ${r4((m.pos - plot.index.pos) / 1e6)}`].join("\n") }));
  const points = plot.rows.map(r => ({
    pos: r.pos, y: nlp(r.p), alpha: 0.7, size: 9 + r.maf * 14, row: r,
    color: r.pos === plot.index.pos ? "#0000FF" : REDS[Math.trunc(r[stat] * 100)],
    label: annotate === "forge" ? r.forge : annotate === "regulome" ? r.regulome : null,
    hover: [`Variant: ${r.rs} ${r.alleles}`, `P-value: ${r.p}`, `Distance (Mb): ${r4(r.dist / 1e6)}`, `MAF: ${r4(r.maf)}`, `R² (${q}): ${r4(r.r2)}`, `D' (${q}): ${r4(r.dprime)}`,
      `Correlated Alleles: ${r.match}`, `FORGEdb Score: ${r.forge ?? ""}`, `RegulomeDB Rank: ${r.regulome ?? ""}`, `Functional Class: ${r.func}`].join("\n"),
  }));
  const pad = 10000;
  return { title: `P-values and Regional LD for ${q} in ${pop}`, xlo: plot.lo - pad, xhi: plot.hi + pad, yLabel: "-log10 P-value", y0: -0.03, y1: ymax * 1.03,
    points: [...missing, ...points], recomb, recombY: rate => (rate / 100) * ymax, hline: nlp(5e-8), genes, collapsed };
}
