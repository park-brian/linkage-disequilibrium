// LDpop: the API's text table per pair (GET /ldpop); LDpair's pairs, plus every fourth pair again with pop=ALL.
import { ldpop } from "../../../src/ldlink/ldpop.js";
import { get, asText } from "../lib.js";
import { pairCases, pairQuery } from "./ldpair.js";

export function cases() {
  const p = pairCases();
  return [...p, ...p.filter((_, i) => i % 4 === 0).map(q => ({ ...q, name: `${q.name} [ALL]`, pop: "ALL" }))];
}
export const ldlink = cs => Promise.all(cs.map(async c => asText(await get("ldpop", pairQuery(c)))));
export const ours = (ctx, c) => ldpop(ctx, c);
