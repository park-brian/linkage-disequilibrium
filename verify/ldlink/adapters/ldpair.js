// LDpair: the API JSON per pair (GET /ldpair, json_out). LDlink orders tied haplotype counts by a hash-randomized
// sample order, so a difference only in that order is "tie_order", not "differ".
import { ldpair } from "../../../src/parity/ldlink/ldpair.js";
import { matrixCases, get, canon } from "../lib.js";

/** consecutive inputs of every LDmatrix case, plus error and edge cases (shared with LDpop) */
export function pairCases() {
  const pairs = [];
  for (const c of matrixCases()) for (let i = 0; i + 1 < c.inputs.length; i++) pairs.push({ name: `${c.name} #${i}`, build: c.build, pop: c.pop, var1: c.inputs[i], var2: c.inputs[i + 1] });
  const edge = (name, var1, var2, pop = "CEU", build = "grch37") => pairs.push({ name, build, pop, var1, var2 });
  edge("different chromosomes", "rs3", "rs12465506");
  edge("coordinate + upper case", "chr13:32446842", "RS4");
  edge("invalid second", "rs3", "xx1"); edge("invalid first", "foo", "rs4"); edge("too short", "rs", "rs4");
  edge("not in dbSNP", "rs9999999999", "rs4"); edge("not in dbSNP (second)", "rs3", "rs9999999999");
  edge("bad population", "rs3", "rs4", "FOO"); edge("two populations", "rs3", "rs4", "CEU+YRI"); edge("ALL", "rs3", "rs4", "ALL");
  edge("coordinate not in dbSNP", "chr13:1", "rs4"); edge("inner whitespace", "rs3 rs5", "rs4");
  edge("chrY on GRCh38", "rs2032597", "rs2032599", "EUR", "grch38"); edge("chrY GRCh37 (haploid)", "rs2032597", "rs2032599", "EUR", "grch37");
  edge("rare pair (monomorphic in CEU)", "rs553970249", "rs536672276");
  return pairs;
}
export const pairQuery = c => `var1=${encodeURIComponent(c.var1)}&var2=${encodeURIComponent(c.var2)}&pop=${encodeURIComponent(c.pop)}&genome_build=${c.build}&json_out=true&r2_d=r2`;

export const cases = pairCases;
export const ldlink = cs => Promise.all(cs.map(async c => {
  const body = await get("ldpair", pairQuery(c));
  const j = JSON.parse(body);
  return canon(Array.isArray(j) ? j.map(({ request, ...r }) => r) : j);
}));
export const ours = async (ctx, c) => canon(await ldpair(ctx, c));
/** same haplotype multiset and count sequence, everything else equal: only tied entries permuted */
export function verdict(l, o) {
  if (JSON.stringify(l) === JSON.stringify(o)) return { verdict: "identical" };
  const h = x => Object.values(x?.[0]?.haplotypes ?? {}), ha = h(l), hb = h(o), k = x => x.map(y => JSON.stringify(y)).sort().join();
  const strip = x => JSON.stringify(Array.isArray(x) ? x.map(({ haplotypes, ...r }) => r) : x);
  if (ha.length && k(ha) === k(hb) && ha.map(y => y.count).join() === hb.map(y => y.count).join() && strip(l) === strip(o)) return { verdict: "tie_order" };
  return { verdict: "differ", detail: { ldlink: l, ours: o } };
}
