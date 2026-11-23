// LDhap: calculate_hap's full JSON (warnings included) and the API's text, from LDlink's Python in the container.
import { ldhap, ldhapText } from "../../../src/ldlink/ldhap.mjs";
import { matrixCases, clean, drive, canon } from "../lib.mjs";

export function cases() {
  const cs = [];
  for (const c of matrixCases()) {
    const ins = c.inputs.map(clean);
    for (let i = 0; i < ins.length; i += 10) cs.push({ name: `${c.name} [${i}..]`, build: c.build, pop: c.pop, snps: ins.slice(i, i + 10) });
    cs.push({ name: `${c.name} [first 30]`, build: c.build, pop: c.pop, snps: ins.slice(0, 30) });
  }
  const edge = (name, snps, pop = "CEU", build = "grch37") => cs.push({ name, build, pop, snps });
  edge("over 30 inputs", Array.from({ length: 31 }, (_, i) => `rs${i + 3}`));
  edge("different chromosomes", ["rs3", "rs4", "rs12465506"]);
  edge("coordinate and its rsID", ["chr13:32446842", "rs3", "rs4"]);
  edge("duplicates and blanks", ["rs3", "", "rs3", "RS4", " rs4 ", "rs3 extra"]);
  edge("nothing valid", ["foo", "rs9999999999"]);
  edge("single multi-allelic", ["rs62224611"]);
  edge("> 1 Mb range", ["rs3", "rs9999", "rs4"]);
  edge("bad population", ["rs3", "rs4"], "FOO");
  edge("chrY on GRCh38", ["rs2032597", "rs2032599", "rs3"], "EUR", "grch38");
  edge("chrY GRCh37 (haploid)", ["rs2032597", "rs2032599"], "EUR");
  return cs;
}
export const ldlink = async cs => drive("ldhap", cs.map(c => ({ module: "ldhap", ...c })))
  .map(r => (r.exception ? { exception: r.exception } : { json: canon(r.json), text: r.text }));
export async function ours(ctx, c) {
  const o = await ldhap(ctx, c);
  return { json: canon(o), text: o.error !== undefined ? null : ldhapText(o, c.build) };
}
export const verdict = (l, o) => (l.exception ? { verdict: "ldlink_exception", detail: { exception: l.exception } } : null); // null: default comparison
