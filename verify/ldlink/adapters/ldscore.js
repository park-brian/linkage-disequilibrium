// LDscore: the REST text (GET /LDlinkRest/ldscore) on LDlink's example PLINK fileset (22.bed/bim/fam, isExample=true):
// cM and kb windows, Python int() edge cases of ldwindow, units, and file-list errors. The two time lines are masked
// (ldherit.js); errors compare on `error`.
import { ldscore } from "../../../src/ldlink/ldscore.js";
import { get, asText, example } from "../lib.js";
import { verdict } from "./ldherit.js";

const EX = "22.bed,22.bim,22.fam";
export function cases() {
  const cs = [];
  for (const w of ["1", "2", "5", "10"]) cs.push({ name: `cM ${w}`, filename: EX, ldwindow: w, windUnit: "cm" });
  for (const w of ["50", "100", "500", "1000"]) cs.push({ name: `kb ${w}`, filename: EX, ldwindow: w, windUnit: "kb" });
  for (const w of ["0", "1.5", "+3", " 2 "]) cs.push({ name: `ldwindow ${JSON.stringify(w)}`, filename: EX, ldwindow: w, windUnit: "cm" });
  cs.push({ name: "ldwindow -1", filename: EX, ldwindow: "-1", windUnit: "cm" });
  cs.push({ name: "ldwindow abc", filename: EX, ldwindow: "abc", windUnit: "cm" });
  cs.push({ name: "windUnit cM", filename: EX, ldwindow: "1", windUnit: "cM" });
  cs.push({ name: "windUnit Mb", filename: EX, ldwindow: "1", windUnit: "Mb" });
  cs.push({ name: "no .fam", filename: "22.bed,22.bim", ldwindow: "1", windUnit: "cm" });
  cs.push({ name: "mixed base names", filename: "22.bed,21.bim,22.fam", ldwindow: "1", windUnit: "cm" });
  cs.push({ name: "missing example", filename: "23.bed,23.bim,23.fam", ldwindow: "1", windUnit: "cm" });
  return cs;
}
const query = c => new URLSearchParams({ pop: "EUR", genome_build: "grch37", filename: c.filename, ldwindow: c.ldwindow, windUnit: c.windUnit, isExample: "true" });
export const ldlink = cs => Promise.all(cs.map(async c => asText(await get("ldscore", query(c)))));
export async function ours(ctx, c) {
  const files = {};
  for (const f of c.filename.split(",")) files[f] = await example(ctx, f, { binary: f.endsWith(".bed") });
  const o = ldscore({ ...c, isExample: true, files });
  return typeof o === "string" ? { text: o } : { error: o.error };
}
export { verdict };
