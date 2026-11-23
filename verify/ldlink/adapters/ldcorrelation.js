// LDcorrelation: the REST text (GET /LDlinkRest/ldcorrelation) on LDlink's example sumstats (isExample=true): the pair in
// both orders and the same file twice, for each population LDlink offers; the liability scale; bad inputs. The time
// lines are masked (ldherit.js); errors compare on `error`.
import { readFileSync, existsSync } from "node:fs";
import { ldcorrelation } from "../../../src/ldlink/ldcorrelation.js";
import { get, asText } from "../lib.js";
import { verdict as textVerdict } from "./ldherit.js";

const POPS = ["EUR", "EAS", "AFR", "AMR", "ASJ", "EST", "FIN", "NFE", "NWE", "SEU"];
const H = "BBJ_HDLC22.txt", L = "BBJ_LDLC22.txt";
export function cases() {
  const cs = [];
  for (const pop of POPS) cs.push({ name: `HDL x LDL ${pop}`, pop, filename: H, filename2: L });
  for (const pop of ["EAS", "EUR"]) cs.push({ name: `LDL x HDL ${pop}`, pop, filename: L, filename2: H });
  cs.push({ name: "HDL x HDL EAS (rg = 1)", pop: "EAS", filename: H, filename2: H, degenerate: true });
  cs.push({ name: "liability EAS 0.3,0.4 / 0.1,0.2", pop: "EAS", filename: H, filename2: L, scale: "liability", samp_prev: "0.3,0.4", pop_prev: "0.1,0.2" });
  cs.push({ name: "liability EUR 0.5,0.5 / .01,1e-2", pop: "EUR", filename: H, filename2: L, scale: "liability", samp_prev: "0.5,0.5", pop_prev: ".01,1e-2" });
  cs.push({ name: "liability with one samp_prev", pop: "EAS", filename: H, filename2: L, scale: "liability", samp_prev: "0.3", pop_prev: "0.1,0.2" });
  cs.push({ name: "liability pop_prev out of range", pop: "EAS", filename: H, filename2: L, scale: "liability", samp_prev: "0.3,0.4", pop_prev: "0.1,2" });
  cs.push({ name: "liability without pop_prev", pop: "EAS", filename: H, filename2: L, scale: "liability", samp_prev: "0.3,0.4" });
  cs.push({ name: "bad scale", pop: "EAS", filename: H, filename2: L, scale: "Observed" });
  cs.push({ name: "unknown population", pop: "XYZ", filename: H, filename2: L });
  return cs;
}
const query = c => new URLSearchParams({ pop: c.pop, genome_build: "grch37", filename: c.filename, filename2: c.filename2, isExample: "true",
  ...(c.scale ? { scale: c.scale } : {}), ...(c.samp_prev ? { samp_prev: c.samp_prev } : {}), ...(c.pop_prev ? { pop_prev: c.pop_prev } : {}) });
export const ldlink = cs => Promise.all(cs.map(async c => asText(await get("ldcorrelation", query(c)))));
const texts = {}, read = (ctx, f) => { const p = `${ctx.root}/examples/${f}`; return (texts[f] ??= existsSync(p) ? readFileSync(p, "utf8") : null); };
export async function ours(ctx, c) {
  const o = await ldcorrelation(ctx, { ...c, text: read(ctx, c.filename), text2: read(ctx, c.filename2), isExample: true });
  return typeof o === "string" ? { text: o } : { error: o.error };
}
const FAIL = /^An error occurred while running the (second )?command \(exit \d+\): /;
const zMask = s => s.replace(/^Z-score: .*$/m, "Z-score: …").replace(/^(\S+\.sumstats\.gz +\S+\.sumstats\.gz +\S+ +\S+ +)\S+/m, "$1…");
/**
 * as LDherit's, plus two documented cases:
 *  - a failed ldsc run: LDlink returns the failure text with ldsc.py's Python traceback; compared up to the colon
 *  - rg = 1 (the same file twice): se ≈ 1e-7, so z = rg / se is floating-point noise past ~7 digits; every other line
 *    must match, and the verdict says which case this is
 */
export function verdict(l, o, c) {
  if (l.text && o.text && FAIL.test(l.text) && FAIL.test(o.text)) return l.text.match(FAIL)[0] === o.text.match(FAIL)[0] ? { verdict: "identical" } : { verdict: "differ", detail: { ldlink: l.text.slice(0, 200), ours: o.text.slice(0, 200) } };
  if (c?.degenerate && l.text && o.text) { const r = textVerdict({ text: zMask(l.text) }, { text: zMask(o.text) }); return r.verdict === "identical" ? { verdict: "degenerate_z" } : r; }
  return textVerdict(l, o);
}
