// LDscore: the REST text (GET /LDlinkRest/ldscore) on LDlink's example PLINK fileset (22.bed/bim/fam, isExample=true):
// cM and kb windows, Python int() edge cases of ldwindow, units, and file-list errors. The two time lines are masked
// (ldherit.js); errors compare on `error`. Upload cases (`uploads`: the example renamed or cut to its first SNPs) go through
// POST /upload, then the request names them (isExample=false): names, the chromosome-from-name rule, a subset fileset.
import { ldscore } from "../../../src/parity/ldlink/ldscore.js";
import { secureFilename } from "../../../src/parity/ldsc/normalize.js";
import { get, asText, example, upload } from "../lib.js";
import { verdict } from "./ldherit.js";

const EX = "22.bed,22.bim,22.fam";
export async function cases(ctx) {
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
  // uploads: the example fileset under other names, and its first 4,000 SNPs (the .bed is SNP-major: 3 magic bytes, then
  // ceil(samples / 4) bytes per SNP)
  const bed = await example(ctx, "22.bed", { binary: true }), bim = await example(ctx, "22.bim"), fam = await example(ctx, "22.fam");
  const set = (root, b = bed, m = bim) => [{ name: `${root}.bed`, bytes: b }, { name: `${root}.bim`, text: m }, { name: `${root}.fam`,
    text: fam }];
  const per = Math.ceil(fam.split("\n").filter(Boolean).length / 4), K = 4000;
  const cut = [bed.slice(0, 3 + K * per), bim.split("\n").slice(0, K).join("\n") + "\n"];
  const up = (name, files, extra = {}) => cs.push({ name: `upload ${name}`, uploads: files, filename: files.map(f => f.name).join(","),
    ldwindow: "1", windUnit: "cm", ...extra });
  up("renamed cohort_chr22", set("cohort_chr22"));
  up("first 4000 SNPs, kb 200", set("chr22.first4000", ...cut), { ldwindow: "200", windUnit: "kb" });
  up("name with no chromosome", set("cohort"));
  up("name with two chromosomes", set("cohort.10.2"));
  up("chr22 wins over a bare 1", set("chr22_1"));
  up("a name with spaces", set("my cohort chr22"));
  up("missing .fam", set("cohort_chr22").slice(0, 2));
  return cs;
}
const query = c => new URLSearchParams({ pop: "EUR", genome_build: "grch37", filename: c.filename, ldwindow: c.ldwindow,
  windUnit: c.windUnit, isExample: "true" });
export const ldlink = cs => Promise.all(cs.map(async c => {
  if (!c.uploads) return asText(await get("ldscore", query(c)));
  const { reference } = await upload(c.uploads);
  return asText(await get("ldscore", new URLSearchParams({ pop: "EUR", genome_build: "grch37", filename: c.filename, ldwindow: c.ldwindow,
    windUnit: c.windUnit, reference, isExample: "false" })));
}));
export async function ours(ctx, c) {
  const files = {};
  if (c.uploads) for (const f of c.uploads) files[secureFilename(f.name)] = f.bytes ?? f.text; // as the upload saved them
  else for (const f of c.filename.split(",")) files[f] = await example(ctx, f, { binary: f.endsWith(".bed") });
  const o = ldscore({ ...c, isExample: !c.uploads, files });
  return typeof o === "string" ? { text: o } : { error: o.error };
}
export { verdict };
