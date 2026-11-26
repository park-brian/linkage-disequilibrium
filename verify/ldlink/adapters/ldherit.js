// LDherit: the REST text (GET /LDlinkRest/ldherit) on LDlink's example sumstats (isExample=true) for each population
// LDlink offers, the liability scale, and bad inputs. Four lines carry the time (munge's "Conversion finished at",
// ldsc's "Beginning analysis at" / "Analysis finished at", "Total time elapsed"): both sides are masked there.
// Errors compare on `error` (the 422 for an ldsc failure also carries a Python traceback in `details`).
// Upload cases (`upload`: a file built from the example in each format the normalizer knows, sumstats.js) go as LDlink's
// web client sends them: POST /upload, validate_sumstats, then the analysis on the normalized file.
import { ldherit } from "../../../src/parity/ldlink/ldherit.js";
import { get, asText, example, uploadSumstats } from "../lib.js";
import { uploadFiles } from "./sumstats.js";

const POPS = ["EUR", "EAS", "AFR", "AMR", "ASJ", "EST", "FIN", "NFE", "NWE", "SEU"];
export async function cases(ctx) {
  const cs = [];
  for (const filename of ["BBJ_HDLC22.txt", "BBJ_LDLC22.txt"]) for (const pop of POPS) cs.push({ name: `${filename} ${pop}`, pop,
    filename });
  cs.push({ name: "liability EAS 0.3 / 0.1", pop: "EAS", filename: "BBJ_HDLC22.txt", scale: "liability", samp_prev: "0.3",
    pop_prev: "0.1" });
  cs.push({ name: "liability EUR 0.5 / 0.01", pop: "EUR", filename: "BBJ_LDLC22.txt", scale: "liability", samp_prev: "0.5",
    pop_prev: "0.01" });
  cs.push({ name: "liability AFR 1e-1 / .05", pop: "AFR", filename: "BBJ_HDLC22.txt", scale: "liability", samp_prev: "1e-1",
    pop_prev: ".05" });
  cs.push({ name: "bad scale", pop: "EUR", filename: "BBJ_HDLC22.txt", scale: "bogus" });
  cs.push({ name: "liability without pop_prev", pop: "EUR", filename: "BBJ_HDLC22.txt", scale: "liability", samp_prev: "0.3" });
  cs.push({ name: "samp_prev out of range", pop: "EUR", filename: "BBJ_HDLC22.txt", scale: "liability", samp_prev: "1.5",
    pop_prev: "0.1" });
  cs.push({ name: "samp_prev not a number", pop: "EUR", filename: "BBJ_HDLC22.txt", scale: "liability", samp_prev: "abc",
    pop_prev: "0.1" });
  cs.push({ name: "unknown population", pop: "XYZ", filename: "BBJ_HDLC22.txt" });
  cs.push({ name: "missing example", pop: "EUR", filename: "nope.txt" });
  const rows = (await example(ctx, "BBJ_HDLC22.txt")).split(/\r?\n/).filter(Boolean).map(l => l.split("\t"));
  const ex = rows.slice(1).map(r => Object.fromEntries(rows[0].map((h, k) => [h, r[k]])));
  const up = uploadFiles(ex);
  for (const k of ["ldsc", "plink", "plink_or", "regenie", "saige"]) cs.push({ name: `upload ${k} EUR`, pop: "EUR", upload: up[k] });
  cs.push({ name: "upload regenie liability EAS 0.3 / 0.1", pop: "EAS", upload: up.regenie, scale: "liability", samp_prev: "0.3",
    pop_prev: "0.1" });
  cs.push({ name: "upload that fails validation", pop: "EUR", upload: { name: "plink_N.txt",
    text: up.plink.text.replace("OBS_CT", "N") } });
  return cs;
}
const query = c => new URLSearchParams({ pop: c.pop, genome_build: "grch37", filename: c.filename, isExample: "true",
  ...(c.scale ? { scale: c.scale } : {}), ...(c.samp_prev ? { samp_prev: c.samp_prev } : {}),
    ...(c.pop_prev ? { pop_prev: c.pop_prev } : {}) });
const opts = c => ({ ...(c.scale ? { scale: c.scale } : {}), ...(c.samp_prev ? { samp_prev: c.samp_prev } : {}),
  ...(c.pop_prev ? { pop_prev: c.pop_prev } : {}) });
export const ldlink = cs => Promise.all(cs.map(async c => {
  if (!c.upload) return asText(await get("ldherit", query(c)));
  const { reference, files: [{ name, fileValid }] } = await uploadSumstats([c.upload]);
  const filename = fileValid.valid ? fileValid.normalized_filename : name; // the web client analyses the normalized file
  return asText(await get("ldherit", new URLSearchParams({ pop: c.pop, genome_build: "grch37", filename, reference, isExample: "false",
    ...opts(c) })));
}));
export async function ours(ctx, c) {
  const o = c.upload ? await ldherit(ctx, { ...c, filename: c.upload.name, text: c.upload.text, isExample: false })
    : await ldherit(ctx, { ...c, text: await example(ctx, c.filename), isExample: true });
  return typeof o === "string" ? { text: o } : { error: o.error };
}
const mask = s => s.replace(/^(Conversion finished at|Beginning analysis at|Analysis finished at|Total time elapsed:).*$/gm, "$1 …");
export function verdict(l, o) {
  if (l.error !== undefined || o.error !== undefined) return l.error === o.error ? { verdict: "identical" } : { verdict: "differ",
    detail: { ldlink: l.error ?? l.text?.slice(0, 300), ours: o.error ?? o.text?.slice(0, 300) } };
  const a = mask(l.text), b = mask(o.text);
  if (a === b) return { verdict: "identical" };
  const x = a.split("\n"), y = b.split("\n"), d = [];
  for (let i = 0; i < Math.max(x.length, y.length); i++) if (x[i] !== y[i]) d.push({ line: i, ldlink: x[i], ours: y[i] });
  return { verdict: "differ", detail: { lines: [x.length, y.length], examples: d.slice(0, 5) } };
}
