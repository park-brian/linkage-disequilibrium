// LDproxy: the API text plus the warning from LDlink's proxy<ref>.json (GET /ldproxy with our own reference).
import { randomUUID } from "node:crypto";
import { ldproxy } from "../../../src/parity/ldlink/ldproxy.js";
import { matrixCases, clean, get, proxyJson } from "../lib.js";

export function cases() {
  const cs = [];
  for (const c of matrixCases()) for (const s of c.inputs.slice(0, 3)) cs.push({ name: `${c.name} ${clean(s)}`, build: c.build, pop: c.pop, var: clean(s) });
  const edge = (name, v, pop = "CEU", build = "grch37", extra = {}) => cs.push({ name, build, pop, var: v, ...extra });
  edge("r2_d=d", "rs3", "EUR", "grch37", { r2_d: "d" }); edge("window 100 kb", "rs3", "EUR", "grch37", { window: 100000 });
  edge("bad r2_d", "rs3", "CEU", "grch37", { r2_d: "x" }); edge("window too large", "rs3", "CEU", "grch37", { window: 2000000 });
  edge("not in dbSNP", "rs9999999999"); edge("coordinate", "chr13:32446842"); edge("bad population", "rs3", "FOO");
  edge("multi-allelic", "rs62224611"); edge("chrY GRCh37", "rs2032597", "EUR"); edge("chrY GRCh38", "rs2032597", "EUR", "grch38");
  edge("near chromosome start", "rs587697622", "EUR"); edge("ALL", "rs3", "ALL");
  return cs.map(c => ({ r2_d: "r2", window: 500000, ...c }));
}
export const ldlink = cs => Promise.all(cs.map(async c => {
  const ref = randomUUID();
  const body = await get("ldproxy", `var=${encodeURIComponent(c.var)}&pop=${encodeURIComponent(c.pop)}&r2_d=${c.r2_d}&window=${c.window}&genome_build=${c.build}&reference=${ref}`);
  if (body.trimStart().startsWith("{")) return { error: JSON.parse(body).error };
  const w = proxyJson(ref).warning;
  return { text: body, ...(w !== undefined ? { warning: w } : {}) };
}));
export const ours = (ctx, c) => ldproxy(ctx, { ...c, chunks: 4 }); // LDlink's NUM_SUBPROCESSES (tools/ldlink-local/server.env.example)
