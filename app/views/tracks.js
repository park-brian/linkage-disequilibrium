// Annotation tracks over one genomic axis, shared by the region view (pyramid.js) and the locus plots. A track is
// { height, gap, draw(g, top), hover(gpos, my, top) → text | null }; layout() stacks them (each owns the band from its
// top to the next track's top, so hover is decided by band) and hoverAt() asks the band under the pointer.
// Tracks: GWAS Catalog (−log10 p), GTEx eQTL (−log10 p, chosen tissues), recombination rate (cM/Mb), RefSeq genes
// (collapsed per symbol, as LDlink's "collapse transcripts"), and the position axis.

/** the shared axis: positions [lo, hi] over W pixels */
export function axisOf(lo, hi, W) {
  const span = Math.max(1, hi - lo);
  return { lo, hi, W, span, x: p => ((p - lo) / span) * W };
}

export const nlp = p => (p > 0 ? -Math.log10(p) : 0);
const GWS = nlp(5e-8);
function title(g, text, top) { g.fillStyle = "#8a93a3"; g.fillText(text, 2, top); }
/** the pointer in the canvas's drawing coordinates (its style size), also when CSS scales it down (max-width: 100%,
 *  aspect kept: one scale for both axes) */
export function pointer(cv, e) {
  const rc = cv.getBoundingClientRect(), k = parseFloat(cv.style.width) / rc.width;
  return [(e.clientX - rc.left) * k, (e.clientY - rc.top) * k];
}
/** a canvas of W × H CSS pixels drawn at the device pixel ratio → { cv, g (scaled to CSS pixels), dpr } */
export function canvas(W, H) {
  const dpr = window.devicePixelRatio || 1, cv = document.createElement("canvas");
  cv.width = W * dpr; cv.height = H * dpr; cv.style.width = `${W}px`; cv.style.height = `${H}px`;
  const g = cv.getContext("2d"); g.scale(dpr, dpr);
  return { cv, g, dpr };
}
/** a tooltip that follows the pointer over cv: textAt(mx, my) in drawing coordinates → text, or null to hide → the tip */
export function tooltip(cv, textAt) {
  const tip = document.createElement("div"); tip.className = "tip"; tip.hidden = true;
  cv.addEventListener("mousemove", e => {
    const text = textAt(...pointer(cv, e));
    if (!text) { tip.hidden = true; return; }
    tip.textContent = text; tip.style.left = `${e.clientX + 12}px`; tip.style.top = `${e.clientY + 12}px`; tip.hidden = false;
  });
  cv.addEventListener("mouseleave", () => { tip.hidden = true; });
  return tip;
}
export function line(g, a, b, c, d) { g.beginPath(); g.moveTo(a, b); g.lineTo(c, d); g.stroke(); }
export function nearest(list, p, at, tol) { let best = null, bd = tol; for (const o of list) { const d = Math.abs(at(o) - p);
  if (d <= bd) { bd = d; best = o; } } return best; }
export const niceStep = s => { const p = 10 ** Math.floor(Math.log10(s)), m = s / p;
  return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p; };

export function gwasTrack(gwas, ax) {
  const H = 70, gp = gwas.map(r => ({ r, p: Number(r.p_value) })).filter(o => Number.isFinite(o.p)), gmax = Math.max(8,
    ...gp.map(o => nlp(o.p)));
  const yOf = (top, v) => top + H - (Math.min(v, gmax) / gmax) * (H - 14);
  return {
    height: H, gap: 8,
    draw(g, top) {
      title(g, "GWAS Catalog −log₁₀ p", top);
      g.strokeStyle = "#e6e8ec"; line(g, 0, top + H - (GWS / gmax) * (H - 14), ax.W, top + H - (GWS / gmax) * (H - 14));
      for (const o of gp) { g.fillStyle = nlp(o.p) >= GWS ? "#7a2e9a" : "#b9a3c6"; g.beginPath(); g.arc(ax.x(o.r.pos), yOf(top, nlp(o.p)),
        2.2, 0, 7); g.fill(); }
    },
    hover(gpos) {
      const near = nearest(gp, gpos, o => o.r.pos, ax.span / ax.W * 4);
      return near && `${near.r.trait}\n${near.r.pubmedid ? `PMID ${near.r.pubmedid} · ` : ""}rs${near.r.snp_id_current} · p ${near.r.p_value}`;
    },
  };
}

export function eqtlTrack(eqtl, ax) {
  const H = 56, emax = Math.max(8, ...eqtl.rows.map(r => nlp(r.p)));
  return {
    height: H, gap: 8,
    draw(g, top) {
      title(g, `GTEx eQTL −log₁₀ p (${eqtl.tissues})`, top);
      for (const r of eqtl.rows) {
        const h = (Math.min(nlp(r.p), emax) / emax) * (H - 14);
        g.fillStyle = r.slope >= 0 ? "rgba(214,110,20,.55)" : "rgba(30,130,140,.55)"; g.fillRect(ax.x(r.pos), top + H - h, 1.2, h);
      }
    },
    hover(gpos) { const near = nearest(eqtl.rows, gpos, r => r.pos, ax.span / ax.W * 3);
      return near && `${near.gene} · ${near.tissue}\nslope ${near.slope} · p ${near.p}`; },
  };
}

/**
 * FORGEdb score (0–10, bars) and RegulomeDB rank (1a strongest … 7, one dot per position, darker = stronger). Not an
 * LDlink track: LDlink shows these scores only as point labels and hovers on its locus plots.
 */
const REG_CLASS = s => { const n = parseInt(s, 10); return n <= 1 ? 0 : n === 2 ? 1 : n === 3 ? 2 : 3; };
const REG_COLOR = ["#7a1f1f", "#c0504d", "#e59b8e", "#e8d3cf"];
export function scoresTrack(forge, regulome, ax) {
  const H = 52, barH = 30;
  return {
    height: H, gap: 8,
    draw(g, top) {
      title(g, "FORGEdb score · RegulomeDB rank (not in LDlink's region views)", top);
      g.fillStyle = "rgba(42,110,90,.6)";
      for (const r of forge) { const h = (r.score / 10) * (barH - 12); g.fillRect(ax.x(r.pos), top + barH - h, 1, Math.max(1, h)); }
      for (const r of regulome) { g.fillStyle = REG_COLOR[REG_CLASS(r.score)]; g.fillRect(ax.x(r.pos),
        top + barH + 6 + REG_CLASS(r.score) * 3, 1.5, 3); }
    },
    hover(gpos, my, top) {
      const tol = ax.span / ax.W * 3;
      if (my < top + barH) { const n = nearest(forge, gpos, r => r.pos, tol);
        return n && `FORGEdb ${n.score} · rs${n.rsid} at ${n.pos.toLocaleString()}`; }
      const n = nearest(regulome, gpos, r => r.pos, tol); return n && `RegulomeDB ${n.score} at ${n.pos.toLocaleString()}`;
    },
  };
}

export function recombTrack(recomb, ax) {
  const H = 50, rmax = Math.max(10, ...recomb.map(r => r.rate));
  return {
    height: H, gap: 8,
    draw(g, top) {
      title(g, "Recombination rate (cM/Mb)", top);
      g.strokeStyle = "#2f6db5"; g.lineWidth = 1; g.beginPath();
      recomb.forEach((r, k) => { const xx = ax.x(r.pos), yy = top + H - (r.rate / rmax) * (H - 14); k ? g.lineTo(xx, yy) : g.moveTo(xx, yy);
        }); g.stroke();
    },
    hover(gpos) { const near = nearest(recomb, gpos, r => r.pos, ax.span);
      return near && `${near.rate.toFixed(2)} cM/Mb at ${near.pos.toLocaleString()}`; },
  };
}

/** collapse transcripts into one span per gene symbol (LDlink's "collapse transcripts"): union of extents and exons */
export function collapse(genes) {
  const by = new Map();
  for (const t of genes) {
    const s = (t.exonStarts ?? "").split(",").filter(Boolean).map(Number), e = (t.exonEnds ?? "").split(",").filter(Boolean).map(Number);
    const g = by.get(t.name2) ?? { name2: t.name2, name: t.name, strand: t.strand, txStart: t.txStart, txEnd: t.txEnd, exons: [],
      transcripts: 0 };
    g.txStart = Math.min(g.txStart, t.txStart); g.txEnd = Math.max(g.txEnd, t.txEnd); g.transcripts++; s.forEach((a, k) => g.exons.push([a,
      e[k]])); by.set(t.name2, g);
  }
  return [...by.values()].sort((a, b) => a.txStart - b.txStart);
}
const measure = document.createElement("canvas").getContext("2d");

/** transcripts as drawn when not collapsed: one entry per transcript, labelled "symbol (accession)" */
function transcripts(genes) {
  return genes.map(t => {
    const s = (t.exonStarts ?? "").split(",").filter(Boolean).map(Number), e = (t.exonEnds ?? "").split(",").filter(Boolean).map(Number);
    return { name2: t.name2, name: t.name, strand: t.strand, txStart: t.txStart, txEnd: t.txEnd, exons: s.map((a, k) => [a, e[k]]),
      label: `${t.name2} (${t.name})` };
  }).sort((a, b) => a.txStart - b.txStart || a.txEnd - b.txEnd);
}

/** RefSeq genes: collapsed per symbol (default, as the region view), or every transcript (LDlink's collapseTranscript=false) */
export function genesTrack(genes, ax, { collapsed = true } = {}) {
  const rows = [], placed = [];
  measure.font = "11px system-ui";
  for (const t of collapsed ? collapse(genes) : transcripts(genes)) { // one row per gene; a label needs its width of free space after the gene
    const end = Math.max(ax.x(t.txEnd), ax.x(t.txStart) + measure.measureText(t.label ?? t.name2).width) + 6;
    let r = rows.findIndex(e => ax.x(t.txStart) > e); if (r < 0) { r = rows.length; rows.push(-Infinity); }
    rows[r] = end; placed.push({ g: t, r });
  }
  return {
    height: Math.max(1, rows.length) * 20 + 6, gap: 8,
    draw(g, top) {
      title(g, collapsed ? "RefSeq genes" : "RefSeq transcripts", top);
      for (const { g: t, r } of placed) {
        const yy = top + 12 + r * 20, xs = Math.max(0, ax.x(t.txStart)), xe = Math.min(ax.W, ax.x(t.txEnd));
        g.fillStyle = g.strokeStyle = "#2b5d34"; line(g, xs, yy + 3, xe, yy + 3);
        for (const [s0, e0] of t.exons) { const a = Math.max(0, ax.x(s0)), b = Math.min(ax.W, ax.x(e0)); if (b >= a) g.fillRect(a, yy,
          Math.max(1, b - a), 7); }
        const text = t.label ?? t.name2; g.fillStyle = "#1d2330"; g.fillText(text, Math.max(0, Math.min(ax.W - g.measureText(text).width,
          xs)), yy + 8);
      }
    },
    hover(gpos, my, top) {
      const hit = placed.find(({ g: t, r }) => gpos >= t.txStart && gpos <= t.txEnd && Math.abs(my - (top + 12 + r * 20)) < 12);
      return hit && `${hit.g.name2} (${hit.g.name}) ${hit.g.strand}\n${hit.g.txStart.toLocaleString()}–${hit.g.txEnd.toLocaleString()}`;
    },
  };
}

export function axisTrack(ax) {
  return {
    height: 26, gap: 0,
    draw(g, top) {
      const step = niceStep(ax.span / 6); g.fillStyle = "#555"; g.strokeStyle = "#999"; line(g, 0, top + 2, ax.W, top + 2);
      let lastEnd = -Infinity;
      for (let t = Math.ceil(ax.lo / step) * step; t <= ax.hi; t += step) {
        const label = `${(t / 1e6).toFixed(step < 1e5 ? 3 : 2)} Mb`, w = g.measureText(label).width, lx = Math.min(ax.W - w, Math.max(0,
          ax.x(t) - w / 2));
        line(g, ax.x(t), top, ax.x(t), top + 6);
        if (lx > lastEnd + 8) { g.fillText(label, lx, top + 8); lastEnd = lx + w; }
      }
    },
    hover: () => null,
  };
}

/** stack tracks from y0: → { tops, end } (end = where the next element starts) */
export function layout(tracks, y0 = 0) {
  const tops = []; let y = y0;
  for (const t of tracks) { tops.push(y); y += t.height + t.gap; }
  return { tops, end: y };
}
export function drawTracks(g, tracks, tops) { tracks.forEach((t, k) => t.draw(g, tops[k])); }
/** the hover text of the track whose band [tops[k], tops[k+1]) holds my (the last band ends at `end`) */
export function hoverAt(tracks, { tops, end }, gpos, my) {
  for (let k = 0; k < tracks.length;
    k++) if (my >= tops[k] && my < (k + 1 < tracks.length ? tops[k + 1] : end)) return tracks[k].hover(gpos, my, tops[k]);
  return null;
}
