// The LDlink comparisons: every module port against a local LDlink ($LDLINK_URL; tools/ldlink-local/README.md).
// Each adapter (adapters/<module>.js) supplies cases(ctx), ldlink(cases) → LDlink's answers, ours(ctx, case) → the
// port's answer in the same shape, and optionally verdict(l, o). The default verdict: equal JSON = "identical", with a
// line diff of the texts otherwise.
// Reports: $LD_DATA/verify/ldlink/<module>.json; prints one line per differing case, then the totals (last line).
// usage: node verify/ldlink/compare.js <module|all>     modules: ldmatrix ldpair ldpop ldhap ldproxy snpclip snpchip ldtrait ldexpress ldassoc
import { writeFileSync, mkdirSync } from "node:fs";
import { nodeContext, canon, textDiff, CONFIG, MODULES } from "./lib.js";

const arg = process.argv[2];
if (!arg || (arg !== "all" && !MODULES.includes(arg))) { console.error(`usage: compare.js <${MODULES.join("|")}|all>`); process.exit(2); }

const ctx = await nodeContext(), out = CONFIG.data("verify/ldlink"), all = {};
mkdirSync(out, { recursive: true });
for (const m of arg === "all" ? MODULES : [arg]) {
  const A = await import(`./adapters/${m}.js`), t0 = performance.now();
  const cases = await A.cases(ctx), theirs = await A.ldlink(cases, ctx), totals = { cases: 0, identical: 0, differ: 0 }, rows = [];
  for (const [i, c] of cases.entries()) {
    const l = theirs[i], o = await A.ours(ctx, c);
    let r = A.verdict?.(l, o, c);
    if (!r) {
      const same = JSON.stringify(canon(l)) === JSON.stringify(canon(o));
      r = same ? { verdict: "identical" } : { verdict: "differ", detail: l?.text && o?.text ? textDiff(l.text, o.text) : { ldlink: l, ours: o } };
    }
    totals.cases++; totals[r.verdict] = (totals[r.verdict] ?? 0) + 1;
    rows.push({ name: c.name, verdict: r.verdict, ...(r.verdict === "identical" ? {} : r.detail ?? {}) });
    if (r.verdict === "differ") console.log(JSON.stringify(rows.at(-1)).slice(0, 1500));
  }
  writeFileSync(`${out}/${m}.json`, JSON.stringify({ module: m, ldlink: CONFIG.LDLINK_URL, release: process.env.ROOT || CONFIG.RELEASE, totals, cases: rows }, null, 2));
  all[m] = totals;
  console.log(JSON.stringify({ module: m, ...totals, seconds: Math.round((performance.now() - t0) / 1000) }));
}
if (arg === "all") console.log(JSON.stringify(all));
if (Object.values(all).some(t => t.differ > 0)) process.exitCode = 1; // tie_order / degenerate_z are documented, not failures
