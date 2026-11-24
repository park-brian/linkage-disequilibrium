// SNPchip on the core: LDlink's API text (details<ref>.txt; server/SNPchip.py calculate_chip + createOutputFile) or
// { error }. Per input, the arrays (snp_col docs at the dbSNP position, dump order) among the requested platform codes.
// LDlink's quirks, reproduced: inputs sorted by (chromosome number, position), stable; Position written "22:123";
// the warnings start ". The following …" (an empty warning is joined with ". "); a platform matches by its name or by
// name.rstrip(' Array') — Python strips the *characters* " Aray" — and is printed with the " Array" suffix removed;
// an empty platform request crashes LDlink (append(...).removesuffix on None).
import { TITLE } from "../../core/genome.js";
import { dedupeLists, check, failFirst, screen, WARNED } from "./common.js";

const CHR = c => (c === "X" ? 23 : c === "Y" ? 24 : +c);

export async function snpchip(ctx, { snps: text, platforms: query, build = "grch37" }) {
  const bad = failFirst(check.build(build));
  if (bad) return bad;
  const lines = text.toLowerCase().split(/\r\n|\n|\r/);
  if (lines.at(-1) === "") lines.pop(); // readlines()
  const bad2 = failFirst(check.list(lines.length, 5000));
  if (bad2) return bad2;
  const lists = dedupeLists(lines);
  const o = {}, toks = lists.filter(t => t.length).map(t => t[0]);
  // replace_coords_rsid_list's warnings go to a dict that is then reset (o is dropped); no chromosome Y class
  const { kept: ok, warned: warn } = await screen(ctx, build, toks, o, { select: undefined, chrY: false });
  const kept = ok.map(({ label, v }) => ({ label, chr: CHR(v.chrom), chrom: v.chrom, pos: v.pos }));
  let warning = "";
  if (!kept.length) return { error: "Input SNP list does not contain any valid RS numbers or coordinates. " };
  if (warn.length) warning = `. ${WARNED(warn)}`;
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
  if (none.length && none.length !== kept.length)
    warning += `The following RS number did not have any platforms found: ${none.join(", ")}. `;
  else if (none.length === kept.length) warning += "Entries from the input variant list were not found in any of the selected arrays. ";
  const out = [["RS Number", `Position (${TITLE[build]})`, "Arrays"].join("\t"), ...rows].join("\n") + "\n";
  return { text: out + (warning !== "" ? `\nWARNING: ${warning}\n\n` : "") };
}
