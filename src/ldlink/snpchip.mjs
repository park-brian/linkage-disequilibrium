// SNPchip on the core: LDlink's API text (details<ref>.txt; server/SNPchip.py calculate_chip + createOutputFile) or
// { error }. Per input, the arrays (snp_col docs at the dbSNP position, dump order) among the requested platform codes.
// LDlink's quirks, reproduced: inputs sorted by (chromosome number, position), stable; Position written "22:123";
// the warnings start ". The following …" (an empty warning is joined with ". "); a platform matches by its name or by
// name.rstrip(' Array') — Python strips the *characters* " Aray" — and is printed with the " Array" suffix removed;
// an empty platform request crashes LDlink (append(...).removesuffix on None).
import { TITLE } from "../genome.mjs";
import { addCoordWarning } from "./common.mjs";

const gate = t => t.length > 2 && (t.startsWith("rs") || t.startsWith("chr")) && /\d$/.test(t);
const CHR = c => (c === "X" ? 23 : c === "Y" ? 24 : +c);

export async function snpchip(ctx, { snps: text, platforms: query, build = "grch37" }) {
  if (!TITLE[build]) return { error: `Invalid genome build. Please specify either ${Object.keys(TITLE).join(", ")}.` };
  const lines = text.toLowerCase().split(/\r\n|\n|\r/);
  if (lines.at(-1) === "") lines.pop(); // readlines()
  if (lines.length > 5000) return { error: `Maximum variant list is 5000  RS numbers or coordinates. Your list contains ${lines.length} entries.` };
  const lists = [];
  for (const l of lines) { const t = l.trim().split(/\s+/).filter(Boolean); if (!lists.some(x => x.join(" ") === t.join(" "))) lists.push(t); }
  const o = {}, toks = lists.filter(t => t.length).map(t => t[0]);
  const vs = (await ctx.planner.query({ build, inputs: toks })).variants;
  vs.forEach(v => addCoordWarning(o, v.warning)); // replace_coords_rsid_list warnings go to a dict that is then reset
  const kept = [], warn = [];
  toks.forEach((t, i) => {
    const v = vs[i], label = t.startsWith("rs") ? t : (v.rsid ?? t);
    if (gate(label) && v.found && v.pos != null) kept.push({ label, chr: CHR(v.chrom), chrom: v.chrom, pos: v.pos }); else warn.push(label);
  });
  let warning = "";
  if (!kept.length) return { error: "Input SNP list does not contain any valid RS numbers or coordinates. " };
  if (warn.length) warning = `. The following RS number(s) or coordinate(s) inputs have warnings: ${warn.join(", ")}`;
  kept.sort((a, b) => a.chr - b.chr || a.pos - b.pos);
  const codes = new Set(String(Array.isArray(query) ? query.join("+") : query).split("+"));
  const wanted = new Set((await ctx.platforms()).filter(p => codes.has(p.code)).map(p => p.platform));
  const at = await ctx.arrays(build, kept.map(k => [k.chrom, k.pos]));
  if (query === "" && kept.some(k => at.has(`${k.chrom}:${k.pos}`))) return { error: "Internal server error. Please contact LDlink admin." };
  const none = [], rows = kept.map(k => {
    const got = [];
    for (const doc of at.get(`${k.chrom}:${k.pos}`) ?? []) for (const p of doc.split("|"))
      if (wanted.has(p) || wanted.has(p.replace(/[ Aray]+$/, ""))) got.push(p.endsWith(" Array") ? p.slice(0, -6) : p);
    if (!got.length) none.push(k.label);
    return [k.label, `${k.chrom}:${k.pos}`, got.join(",")].join("\t");
  });
  if (none.length && none.length !== kept.length) warning += `The following RS number did not have any platforms found: ${none.join(", ")}. `;
  else if (none.length === kept.length) warning += "Entries from the input variant list were not found in any of the selected arrays. ";
  const out = [["RS Number", `Position (${TITLE[build]})`, "Arrays"].join("\t"), ...rows].join("\n") + "\n";
  return { text: out + (warning !== "" ? `\nWARNING: ${warning}\n\n` : "") };
}
