// Region view: annotation tracks (tracks.js) over a shared genomic axis, and the LD pyramid below them (pairs at
// midpoint × half distance, pixels = the strongest pair). Hover anywhere for the underlying record.
import { axisOf, gwasTrack, eqtlTrack, scoresTrack, recombTrack, genesTrack, axisTrack, layout, drawTracks, hoverAt } from "./tracks.js";

const RAMP = { r2: [205, 40, 40], d: [36, 87, 197] };

export function renderPyramid(holder, v) {
  const { W, H, grid, bi, bj, rows, labels, lo, hi, stat, genes, recomb, gwas, eqtl, forge, regulome } = v;
  const ax = axisOf(lo, hi, W);
  const tracks = [gwasTrack(gwas, ax), eqtlTrack(eqtl, ax), ...(forge ? [scoresTrack(forge, regulome, ax)] : []), recombTrack(recomb, ax), genesTrack(genes, ax), axisTrack(ax)];
  const lay = layout(tracks), Y0 = lay.end, Htot = Y0 + H, dpr = window.devicePixelRatio || 1, cv = document.createElement("canvas");
  cv.width = W * dpr; cv.height = Htot * dpr; cv.style.width = `${W}px`; cv.style.height = `${Htot}px`;
  const g = cv.getContext("2d"); g.scale(dpr, dpr); g.font = "11px system-ui"; g.textBaseline = "top";
  drawTracks(g, tracks, lay.tops);
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
    const rc = cv.getBoundingClientRect(), mx = e.clientX - rc.left, my = e.clientY - rc.top, gpos = lo + (mx / W) * ax.span;
    let text = null;
    if (my >= Y0) {
      const o = Math.floor(my - Y0) * W + Math.floor(mx);
      if (grid[o] >= 0) { const i = bi[o], j = bj[o]; text = `${labels[i]} × ${labels[j]}\n${stat === "r2" ? "r²" : "D′"} ${grid[o].toFixed(3)}   distance ${(rows[j].pos - rows[i].pos).toLocaleString()} bp`; }
    } else text = hoverAt(tracks, lay, gpos, my);
    if (!text) { tip.hidden = true; return; }
    tip.textContent = text; tip.style.left = `${e.clientX + 12}px`; tip.style.top = `${e.clientY + 12}px`; tip.hidden = false;
  });
  cv.addEventListener("mouseleave", () => { tip.hidden = true; });
  if (v.onPair) cv.addEventListener("click", e => { // a pixel → its strongest pair
    const rc = cv.getBoundingClientRect(), my = e.clientY - rc.top, o = Math.floor(my - Y0) * W + Math.floor(e.clientX - rc.left);
    if (my >= Y0 && grid[o] >= 0) v.onPair(labels[bi[o]], labels[bj[o]], e);
  });
  holder.replaceChildren(cv, tip);
}
