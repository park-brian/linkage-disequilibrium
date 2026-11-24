// LDassoc: calculate_assoc in the container (the web route is its only entry point; no API) on LDlink's bundled
// prostate example: the table text, plus the JSON's error and warning.
import { ldassoc } from "../../../src/parity/ldlink/ldassoc.js";
import { drive, example } from "../lib.js";

const FILE = { grch37: "prostate_example_grch37.txt", grch38: "prostate_example_grch38.txt",
  grch38_high_coverage: "prostate_example_grch38.txt" };
export function cases() {
  const cs = [], add = (name, c) => cs.push({ name, columns: ["chr", "pos", "p"], build: "grch37", pop: "EUR", ...c });
  for (const build of ["grch37", "grch38", "grch38_high_coverage"]) {
    add(`variant rs7837688 ±100 kb EUR ${build}`, { build, region: "variant", origin: "rs7837688", window: 100000 });
    add(`gene MYC ±100 kb (lowest p origin) ${build}`, { build, region: "gene", geneName: "MYC", window: 100000 });
    add(`region chr8:128.0–128.7 Mb YRI ${build}`, { build, region: "region", start: "chr8:128000000", end: "chr8:128700000", pop: "YRI" });
  }
  add("variant rs6983267 ±500 kb ALL", { region: "variant", origin: "rs6983267", window: 500000, pop: "ALL" });
  add("gene MYC with origin rs6983267", { region: "gene", geneName: "MYC", window: 100000, origin: "rs6983267" });
  add("region with origin, CEU+YRI", { region: "region", start: "chr8:128400000", end: "chr8:128600000", origin: "rs7837688",
    pop: "CEU+YRI" });
  add("region chrX", { region: "region", start: "chrX:2700000", end: "chrX:3500000" });
  add("region chr22 (no association rows)", { region: "region", start: "chr22:30000000", end: "chr22:30100000" });
  add("region > 3 Mb", { region: "region", start: "chr8:120000000", end: "chr8:124000000" });
  add("region start >= end (string compare)", { region: "region", start: "chr8:9000000", end: "chr8:10000000" });
  add("coordinate origin", { region: "variant", origin: "chr8:128539360", window: 100000 });
  add("origin not in dbSNP", { region: "variant", origin: "rs9999999999", window: 100000 });
  add("gene unknown", { region: "gene", geneName: "NOTAGENE", window: 100000 });
  add("origin missing p-value", { region: "variant", origin: "rs3", window: 100000 });
  add("bad population", { region: "gene", geneName: "MYC", window: 100000, pop: "FOO" });
  return cs;
}
// one shape for both sides: the error wins over text and warning (the JSON has both when LDlink fails late)
const norm = x => ({ text: x.text ?? null, error: x.error ?? null, warning: x.error ? null : x.warning ?? null });
export const ldlink = async cs => drive("ldassoc", cs.map(({ geneName, ...c }) => ({ module: "ldassoc", file: FILE[c.build], ...c,
  name: geneName })))
  .map(r => norm(r.exception ? { error: "Internal server error. Please contact LDlink admin." } : { text: r.text ?? undefined,
    error: r.json.error, warning: r.json.warning }));
export async function ours(ctx, { geneName, ...c }) {
  return norm(await ldassoc(ctx, { ...c, name: geneName, text: await example(ctx, FILE[c.build]) }));
}
