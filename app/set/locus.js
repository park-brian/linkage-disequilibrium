// The locus lens: the set on its genomic axis. From the top: with one variant selected, every common variant's LD with
// it (LDproxy's question, as a track); GWAS Catalog hits; recombination; genes; the set's variants as ticks on the axis;
// and the LD pyramid of the window hanging from that axis (pairs at midpoint × half distance, each pixel the strongest
// pair, computed in the worker pool as the region view does). A click on a tick or a point selects that variant; a
// click in the pyramid selects that pixel's pair, whose two 45° rails are then drawn down to it.
import { pyramid } from "../workers/pool.js";
import { axisOf, gwasTrack, recombTrack, genesTrack, axisTrack, layout, drawTracks, hoverAt, pointer, canvas, tooltip,
  line } from "../views/tracks.js";
import { h, ramp, INK, STAT, fix, mb } from "./dom.js";
import { windowOf, MAX_REGION_ROWS } from "./set.js";
import { clickSelect } from "./table.js";

const MUTED = "#8a8f98", GRID = "#e9e7e1", DARK = "#16181d", MAX_DIST = 500_000, THRESHOLD = 0.05;

// a window's data (rows, tracks, samples) is kept for the last few windows: a population, statistic or selection
// change recomputes LD only
const windows = new Map();
function windowData(ctx, build, { chrom, lo, hi }) {
  const key = `${build}|${chrom}|${lo}|${hi}`;
  if (!windows.has(key)) {
    if (windows.size > 5) windows.delete(windows.keys().next().value);
    windows.set(key, Promise.all([ctx.windowRows(build, chrom, lo, hi, ["common"]), ctx.genes(build, chrom, lo, hi),
      ctx.recomb(build, chrom, lo, hi), ctx.gwas(build, chrom, lo, hi), ctx.samples(build, chrom)])
      .then(([rows, genes, recomb, gwas, smp]) => ({ genes, recomb, gwas: gwas.map(r => ({ ...r })), smp,
        rows: rows.filter(r => r.pos >= lo && r.pos <= hi).sort((a, b) => a.pos - b.pos || a.vcf_order - b.vcf_order) }))
      .catch(e => { windows.delete(key); throw e; }));
  }
  return windows.get(key);
}
/** the genotype rows of the window last drawn, by label: a pair picked in the pyramid need not be in the set */
export const known = new Map();
// the pool's workers answer one pyramid at a time: runs queue here, and a run that is stale when its turn comes is skipped
let queue = Promise.resolve();
const inTurn = fn => { const next = queue.then(fn, fn); queue = next.catch(() => {}); return next; };

/** LD of every window row with the anchor's genotype row → [{ pos, label, v }] */
async function withAnchor(ctx, S, w, d, labels, g) {
  const rows = [g, ...d.rows.map(r => ({ ...r, chrom: w.chrom }))], J = d.rows.map((_, i) => i + 1), I = J.map(() => 0);
  const e = await ctx.withPrepared(S.build, rows, S.pop, p => ctx.kernel.exact(p, I, J)), val = S.stat === "r2" ? e.r2 : e.dprime;
  return d.rows.map((r, i) => ({ pos: r.pos, label: labels[i], v: val[i] }));
}

function proxyTrack(points, ax, stat, anchor, memberPos) {
  const H = 176, T = 20, y = (top, v) => top + T + (1 - v) * (H - T - 8), px = points.filter(p => !Number.isNaN(p.v));
  const near = (gpos, my, top) => { let best = null, bd = 7;
    for (const p of px) { const d = Math.hypot(((p.pos - gpos) / ax.span) * ax.W, y(top, p.v) - my); if (d < bd) { bd = d; best = p; } }
    return best; };
  return {
    height: H, gap: 12,
    draw(g, top) {
      g.fillStyle = MUTED; g.fillText(`${STAT[stat]} with ${anchor.label}`, 2, top);
      for (const t of [0, 0.5, 1]) { g.strokeStyle = GRID; line(g, 0, y(top, t), ax.W, y(top, t)); g.fillStyle = MUTED;
        g.fillText(String(t), ax.W - 18, y(top, t) - 12); }
      for (const p of px) if (!memberPos.has(p.pos)) { g.fillStyle = ramp(stat, 0.12 + 0.88 * p.v, 0.8); g.beginPath();
        g.arc(ax.x(p.pos), y(top, p.v), 2.2, 0, 7); g.fill(); }
      g.strokeStyle = DARK; g.lineWidth = 1.2;
      for (const p of px) if (memberPos.has(p.pos) && p.pos !== anchor.pos) { g.fillStyle = ramp(stat, 0.12 + 0.88 * p.v); g.beginPath();
        g.arc(ax.x(p.pos), y(top, p.v), 4.2, 0, 7); g.fill(); g.stroke(); }
      const X = ax.x(anchor.pos), Y = y(top, 1); g.fillStyle = DARK; g.beginPath(); g.moveTo(X, Y - 6); g.lineTo(X + 5, Y);
      g.lineTo(X, Y + 6); g.lineTo(X - 5, Y); g.fill(); g.lineWidth = 1;
    },
    hover(gpos, my, top) { const p = near(gpos, my, top);
      return p && `${p.label}\n${STAT[stat]} ${fix(p.v)} with ${anchor.label} · `
        + `${Math.abs(p.pos - anchor.pos).toLocaleString()} bp away`; },
    pick: (gpos, my, top) => near(gpos, my, top)?.label,
  };
}

function setTrack(members, ax, selPos) {
  const dense = members.length > 150, named = members.length <= 16, H = named ? 34 : 14;
  const near = gpos => { let best = null, bd = (ax.span / ax.W) * 5;
    for (const m of members) { const d = Math.abs(m.pos - gpos); if (d <= bd) { bd = d; best = m; } } return best; };
  return {
    height: H, gap: 2,
    draw(g, top) {
      let lastEnd = -Infinity;
      for (const m of members) {
        const X = ax.x(m.pos), on = selPos.has(m.pos);
        g.strokeStyle = on ? DARK : dense ? "rgba(22,24,29,.22)" : "rgba(22,24,29,.75)"; g.lineWidth = on ? 2 : 1;
        line(g, X, top + H - (on ? 12 : 8), X, top + H);
        if (!named) continue;
        const w = g.measureText(m.label).width, lx = Math.max(0, Math.min(ax.W - w, X - w / 2));
        if (lx < lastEnd + 6) continue;
        g.fillStyle = on ? DARK : MUTED; g.fillText(m.label, lx, top + 4); lastEnd = lx + w;
      }
      g.lineWidth = 1;
    },
    hover(gpos) { const m = near(gpos); return m && `${m.label}\n${m.pos.toLocaleString()}`; },
    pick: gpos => near(gpos)?.label,
  };
}

export async function locus(env) {
  const { ctx, S, set, sel } = env, first = set.ok.find(r => sel.includes(r.label));
  // a set on several chromosomes: the locus of the selected variant's chromosome
  let w = windowOf(set);
  if (w.error && first) w = windowOf({ kind: "list", ok: set.ok.filter(r => r.chrom === first.chrom) });
  if (w.error) return env.say(w.error);
  const d = await windowData(ctx, S.build, w), n = d.rows.length;
  if (n < 2) return env.say("Fewer than two common 1000 Genomes variants in this window.");
  if (n > MAX_REGION_ROWS) return env.say(`${n.toLocaleString()} common variants in this window: narrow the set.`);
  const W = Math.max(520, Math.min(1280, env.figure.clientWidth)), ids = await ctx.popIds(S.pop);
  const labels = d.rows.map(r => (r.rsid != null ? `rs${r.rsid}` : `chr${w.chrom}:${r.pos}`));
  // label → position: the window's labels, then the set's own (a list's dbSNP label wins over the panel's)
  const posOf = new Map(labels.map((l, i) => [l, d.rows[i].pos])),
    rowOf = new Map(labels.map((l, i) => [l, { ...d.rows[i], chrom: w.chrom }]));
  for (const r of set.ok) if (r.chrom === w.chrom) { posOf.set(r.label, r.pos); rowOf.set(r.label, r.g); }
  known.clear(); for (const [l, r] of rowOf) known.set(l, r);
  const members = set.ok.filter(r => r.chrom === w.chrom && r.pos >= w.lo && r.pos <= w.hi), chosen = sel.filter(l => posOf.has(l)),
    picked = chosen.length || set.ok.length !== 1 ? chosen : [set.ok[0].label]; // a set of one variant is its own anchor
  const anchor = picked.length === 1 ? { label: picked[0], pos: posOf.get(picked[0]) } : null;
  const [res, points] = await Promise.all([inTurn(() => (env.stale() ? null : pyramid(d.rows.map(r => ({ pos: r.pos, haps: r.haps,
    missing: r.missing, male_haploid: r.male_haploid })), { samples: d.smp.ids.length, haploidMales: [...d.smp.haploidMales],
  hapIndexes: d.smp.ids.flatMap((s, i) => (ids.has(s) ? [2 * i, 2 * i + 1] : [])), maxDist: Math.min(MAX_DIST, w.hi - w.lo),
  threshold: THRESHOLD, stat: S.stat, W, lo: w.lo, hi: w.hi, radius: Math.max(0, Math.min(6, Math.round(W / n))) }))),
  anchor ? withAnchor(ctx, S, w, d, labels, rowOf.get(anchor.label)) : null]);
  if (!res || env.stale()) return;

  const ax = axisOf(w.lo, w.hi, W), selPos = new Set(picked.map(l => posOf.get(l))), memberPos = new Set(members.map(m => m.pos));
  const tracks = [...(anchor ? [proxyTrack(points, ax, S.stat, anchor, memberPos)] : []), ...(d.gwas.length ? [gwasTrack(d.gwas, ax)] : []),
    recombTrack(d.recomb, ax), genesTrack(d.genes, ax), setTrack(members, ax, selPos), axisTrack(ax)];
  const lay = layout(tracks), Y0 = lay.end, { H, grid, bi, bj } = res, { cv, g, dpr } = canvas(W, Y0 + H);
  g.font = "11px system-ui"; g.textBaseline = "top";
  drawTracks(g, tracks, lay.tops);
  const img = g.createImageData(W * dpr, H * dpr), [cr, cg, cb] = INK[S.stat];
  for (let py = 0; py < H * dpr; py++) for (let px = 0; px < W * dpr; px++) {
    const v = grid[Math.floor(py / dpr) * W + Math.floor(px / dpr)], o = (py * W * dpr + px) * 4;
    if (v < 0) continue;
    img.data[o] = 255 - (255 - cr) * v; img.data[o + 1] = 255 - (255 - cg) * v;
      img.data[o + 2] = 255 - (255 - cb) * v; img.data[o + 3] = 255;
  }
  g.putImageData(img, 0, Y0 * dpr);
  // the selection: a hairline up through the tracks, and the 45° rails pairs of that variant lie on; two meet at the pair
  g.save(); g.beginPath(); g.rect(0, 0, W, Y0 + H); g.clip(); g.strokeStyle = "rgba(22,24,29,.45)";
  for (const p of selPos) { const X = ax.x(p); line(g, X, 0, X, Y0); line(g, X, Y0, X + H, Y0 + H); line(g, X, Y0, X - H, Y0 + H); }
  if (selPos.size === 2) { const [a, b] = [...selPos].sort((x, y) => x - y), X = ax.x((a + b) / 2), Y = Y0 + ((b - a) / 2 / ax.span) * W;
    g.strokeStyle = DARK; g.lineWidth = 1.5; g.beginPath(); g.arc(X, Y, 5, 0, 7); g.stroke(); }
  g.restore();

  const at = (mx, my) => { const o = Math.floor(my - Y0) * W + Math.floor(mx); return my >= Y0 && grid[o] >= 0 ? o : -1; };
  const band = my => tracks.findIndex((t, k) => my >= lay.tops[k] && my < (k + 1 < tracks.length ? lay.tops[k + 1] : lay.end));
  const pickAt = (mx, my) => { const k = band(my); return k < 0 ? null : tracks[k].pick?.(w.lo + (mx / W) * ax.span, my, lay.tops[k]); };
  const tip = tooltip(cv, (mx, my) => {
    const o = at(mx, my);
    cv.style.cursor = o >= 0 || (my < Y0 && pickAt(mx, my)) ? "pointer" : "default";
    if (my < Y0) return hoverAt(tracks, lay, w.lo + (mx / W) * ax.span, my);
    return o < 0 ? null : `${labels[bi[o]]} × ${labels[bj[o]]}\n${STAT[S.stat]} ${grid[o].toFixed(3)} · ${(d.rows[bj[o]].pos
      - d.rows[bi[o]].pos).toLocaleString()} bp apart`;
  });
  cv.addEventListener("click", e => {
    const [mx, my] = pointer(cv, e), o = at(mx, my);
    if (o >= 0) return env.select([labels[bi[o]], labels[bj[o]]]);
    const label = my < Y0 ? pickAt(mx, my) : null;
    if (label) env.select(clickSelect(sel, label, e.ctrlKey || e.metaKey || e.shiftKey).slice(-2));
  });
  env.headline(`Chromosome ${w.chrom}, ${mb(w.lo)} to ${mb(w.hi)}`);
  env.show(cv, tip);
  const swatch = h("span.ramp", { style: `background:linear-gradient(90deg,${ramp(S.stat, 0)},${ramp(S.stat, 1)})` });
  env.caption(h("span", { title: `${res.evaluated.toLocaleString()} pairs evaluated, `
    + `${res.pairs.toLocaleString()} at ${THRESHOLD} or above · `
    + `LD and pixels ${res.ms.wall} ms on ${res.workers} worker${res.workers > 1 ? "s" : ""} (${res.backend})` },
  `${STAT[S.stat]} 0 `, swatch, ` 1 among the ${n.toLocaleString()} common variants (MAF ≥ 5%) `
    + `of this window in ${S.pop.replaceAll("+", " + ")}. `
    + `Each pixel is the strongest pair at that midpoint and separation, up to ${Math.min(MAX_DIST, w.hi - w.lo) / 1000} kb apart.`));
  // what the selection says, under the figure
  if (picked.length === 2) { // from the two genotype rows themselves, as the pyramid's pixels are (the pair lens asks LDlink's port)
    const e = await ctx.withPrepared(S.build, picked.map(l => rowOf.get(l)), S.pop, p => ctx.kernel.exact(p, [0], [1]));
    if (env.stale()) return;
    env.strip(h("b", `${picked[0]} × ${picked[1]}`), ` r² ${fix(e.r2[0], 4)} · D′ ${fix(e.dprime[0], 4)} · `
      + `${Math.abs(posOf.get(picked[0]) - posOf.get(picked[1])).toLocaleString()} bp apart `,
    h("a", { href: "#", onclick: ev => { ev.preventDefault(); env.update({ lens: "pair" }); } }, "Open the pair"));
  } else if (anchor) env.strip(h("b", anchor.label), ` is the anchor: the upper track is every variant's ${STAT[S.stat]} with it. `,
    h("a", { href: "#", onclick: e => { e.preventDefault(); env.update({ lens: "proxies" }); } }, "List its proxies"));
  else env.strip("Click a variant to anchor on it, or the pyramid to pick a pair.");
}
