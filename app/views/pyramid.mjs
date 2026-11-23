// Region view: annotation tracks over a shared genomic axis, and the LD pyramid below it (pairs at midpoint × half
// distance, pixels = the strongest pair). Tracks: GWAS Catalog (−log10 p), GTEx eQTL (−log10 p, chosen tissues),
// recombination rate (cM/Mb), RefSeq genes. Hover anywhere for the underlying record.
const RAMP = { r2: [205, 40, 40], d: [36, 87, 197] };

/** collapse transcripts into one span per gene symbol (LDlink's "collapse transcripts"): union of extents and exons */
function collapse(genes) {
  const by = new Map();
  for (const t of genes) {
    const s = (t.exonStarts ?? "").split(",").filter(Boolean).map(Number), e = (t.exonEnds ?? "").split(",").filter(Boolean).map(Number);
    const g = by.get(t.name2) ?? { name2: t.name2, name: t.name, strand: t.strand, txStart: t.txStart, txEnd: t.txEnd, exons: [], transcripts: 0 };
    g.txStart = Math.min(g.txStart, t.txStart); g.txEnd = Math.max(g.txEnd, t.txEnd); g.transcripts++; s.forEach((a, k) => g.exons.push([a, e[k]])); by.set(t.name2, g);
  }
  return [...by.values()].sort((a, b) => a.txStart - b.txStart);
}
const g0 = document.createElement("canvas").getContext("2d");

export function renderPyramid(holder, v) {
  const { W, H, grid, bi, bj, rows, labels, lo, hi, stat, genes, recomb, gwas, eqtl } = v;
  const span = Math.max(1, hi - lo), x = p => ((p - lo) / span) * W;
  const rowsG = [], placed = [];
  g0.font = "11px system-ui";
  for (const t of collapse(genes)) { // one row per gene; a label needs its width of free space after the gene
    const end = Math.max(x(t.txEnd), x(t.txStart) + g0.measureText(t.name2).width) + 6;
    let r = rowsG.findIndex(e => x(t.txStart) > e); if (r < 0) { r = rowsG.length; rowsG.push(-Infinity); }
    rowsG[r] = end; placed.push({ g: t, r });
  }
  const T = { gwas: 70, eqtl: 56, recomb: 50, genes: Math.max(1, rowsG.length) * 20 + 6, axis: 26 };
  const top = {}; let y = 0; for (const k of ["gwas", "eqtl", "recomb", "genes", "axis"]) { top[k] = y; y += T[k] + (k === "axis" ? 0 : 8); }
  const Y0 = y, Htot = Y0 + H, dpr = window.devicePixelRatio || 1, cv = document.createElement("canvas");
  cv.width = W * dpr; cv.height = Htot * dpr; cv.style.width = `${W}px`; cv.style.height = `${Htot}px`;
  const g = cv.getContext("2d"); g.scale(dpr, dpr); g.font = "11px system-ui"; g.textBaseline = "top";
  const title = (k, t) => { g.fillStyle = "#8a93a3"; g.fillText(t, 2, top[k]); };
  const nlp = p => (p > 0 ? -Math.log10(p) : 0);
  // GWAS points
  title("gwas", "GWAS Catalog −log₁₀ p");
  const gp = gwas.map(r => ({ r, p: Number(r.p_value) })).filter(o => Number.isFinite(o.p)), gmax = Math.max(8, ...gp.map(o => nlp(o.p)));
  const gy = o => top.gwas + T.gwas - (Math.min(nlp(o.p), gmax) / gmax) * (T.gwas - 14);
  g.strokeStyle = "#e6e8ec"; line(g, 0, top.gwas + T.gwas - (nlp(5e-8) / gmax) * (T.gwas - 14), W, top.gwas + T.gwas - (nlp(5e-8) / gmax) * (T.gwas - 14));
  for (const o of gp) { g.fillStyle = nlp(o.p) >= nlp(5e-8) ? "#7a2e9a" : "#b9a3c6"; g.beginPath(); g.arc(x(o.r.pos), gy(o), 2.2, 0, 7); g.fill(); }
  // eQTL ticks
  title("eqtl", `GTEx eQTL −log₁₀ p (${eqtl.tissues})`);
  const emax = Math.max(8, ...eqtl.rows.map(r => nlp(r.p)));
  for (const r of eqtl.rows) { const h = (Math.min(nlp(r.p), emax) / emax) * (T.eqtl - 14); g.fillStyle = r.slope >= 0 ? "rgba(214,110,20,.55)" : "rgba(30,130,140,.55)"; g.fillRect(x(r.pos), top.eqtl + T.eqtl - h, 1.2, h); }
  // recombination line
  title("recomb", "Recombination rate (cM/Mb)");
  const rmax = Math.max(10, ...recomb.map(r => r.rate));
  g.strokeStyle = "#2f6db5"; g.lineWidth = 1; g.beginPath();
  recomb.forEach((r, k) => { const xx = x(r.pos), yy = top.recomb + T.recomb - (r.rate / rmax) * (T.recomb - 14); k ? g.lineTo(xx, yy) : g.moveTo(xx, yy); }); g.stroke();
  // genes
  title("genes", "RefSeq genes");
  for (const { g: t, r } of placed) {
    const yy = top.genes + 12 + r * 20, xs = Math.max(0, x(t.txStart)), xe = Math.min(W, x(t.txEnd));
    g.fillStyle = g.strokeStyle = "#2b5d34"; line(g, xs, yy + 3, xe, yy + 3);
    for (const [s0, e0] of t.exons) { const a = Math.max(0, x(s0)), b = Math.min(W, x(e0)); if (b >= a) g.fillRect(a, yy, Math.max(1, b - a), 7); }
    g.fillStyle = "#1d2330"; g.fillText(t.name2, Math.max(0, Math.min(W - g.measureText(t.name2).width, xs)), yy + 8);
  }
  // axis
  const step = niceStep(span / 6); g.fillStyle = "#555"; g.strokeStyle = "#999"; line(g, 0, top.axis + 2, W, top.axis + 2);
  let lastEnd = -Infinity;
  for (let t = Math.ceil(lo / step) * step; t <= hi; t += step) { const label = `${(t / 1e6).toFixed(step < 1e5 ? 3 : 2)} Mb`, w = g.measureText(label).width, lx = Math.min(W - w, Math.max(0, x(t) - w / 2)); line(g, x(t), top.axis, x(t), top.axis + 6); if (lx > lastEnd + 8) { g.fillText(label, lx, top.axis + 8); lastEnd = lx + w; } }
  // pyramid
  const img = g.createImageData(W * dpr, H * dpr), [cr, cg, cb] = RAMP[stat];
  for (let py = 0; py < H * dpr; py++) for (let px = 0; px < W * dpr; px++) {
    const vv = grid[Math.floor(py / dpr) * W + Math.floor(px / dpr)], o = (py * W * dpr + px) * 4;
    if (vv < 0) { img.data[o + 3] = 0; continue; }
    img.data[o] = 255 - (255 - cr) * vv; img.data[o + 1] = 255 - (255 - cg) * vv; img.data[o + 2] = 255 - (255 - cb) * vv; img.data[o + 3] = 255;
  }
  g.putImageData(img, 0, Y0 * dpr);
  // hover
  const tip = document.createElement("div"); tip.className = "tip"; tip.hidden = true;
  cv.addEventListener("mousemove", e => {
    const rc = cv.getBoundingClientRect(), mx = e.clientX - rc.left, my = e.clientY - rc.top, gpos = lo + (mx / W) * span;
    let text = null;
    if (my >= Y0) {
      const o = Math.floor(my - Y0) * W + Math.floor(mx);
      if (grid[o] >= 0) { const i = bi[o], j = bj[o]; text = `${labels[i]} × ${labels[j]}\n${stat === "r2" ? "r²" : "D′"} ${grid[o].toFixed(3)}   distance ${(rows[j].pos - rows[i].pos).toLocaleString()} bp`; }
    } else if (my < top.eqtl) {
      const near = nearest(gp, gpos, o => o.r.pos, span / W * 4); if (near) text = `${near.r.trait}\n${near.r.pubmedid ? `PMID ${near.r.pubmedid} · ` : ""}rs${near.r.snp_id_current} · p ${near.r.p_value}`;
    } else if (my < top.recomb) {
      const near = nearest(eqtl.rows, gpos, r => r.pos, span / W * 3); if (near) text = `${near.gene} · ${near.tissue}\nslope ${near.slope} · p ${near.p}`;
    } else if (my < top.genes) {
      const near = nearest(recomb, gpos, r => r.pos, span); if (near) text = `${near.rate.toFixed(2)} cM/Mb at ${near.pos.toLocaleString()}`;
    } else if (my < top.axis) {
      const hit = placed.find(({ g: t, r }) => gpos >= t.txStart && gpos <= t.txEnd && Math.abs(my - (top.genes + 12 + r * 20)) < 12); if (hit) text = `${hit.g.name2} (${hit.g.name}) ${hit.g.strand}\n${hit.g.txStart.toLocaleString()}–${hit.g.txEnd.toLocaleString()}`;
    }
    if (!text) { tip.hidden = true; return; }
    tip.textContent = text; tip.style.left = `${e.clientX + 12}px`; tip.style.top = `${e.clientY + 12}px`; tip.hidden = false;
  });
  cv.addEventListener("mouseleave", () => { tip.hidden = true; });
  holder.replaceChildren(cv, tip);
}
function line(g, a, b, c, d) { g.beginPath(); g.moveTo(a, b); g.lineTo(c, d); g.stroke(); }
function nearest(list, p, at, tol) { let best = null, bd = tol; for (const o of list) { const d = Math.abs(at(o) - p); if (d <= bd) { bd = d; best = o; } } return best; }
const niceStep = s => { const p = 10 ** Math.floor(Math.log10(s)), m = s / p; return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p; };
