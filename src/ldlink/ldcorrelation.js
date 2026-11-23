// LDcorrelation on the core: LDlink's REST text (server/LDlink.py ldcorrelation route → ldsc_utils.run_correlation_command):
// munge both sumstats (--merge-alleles w_hm3.snplist), then ldsc.py --rg file1,file2 on the population's panel. The
// route returns the output from the first "Total Observed scale" (or, on the liability scale, "Total Liability scale")
// onward: phenotype 1's h2 summary, phenotype 2's, the genetic covariance and correlation, the pandas summary table, and
// ldsc's two closing time lines (written now; the comparer masks them). Errors come back as { error, details? }.
import { pyRoundStr, npStr4 } from "../pyfmt.js";
import { mungeSumstats } from "../ldsc/munge.js";
import { estimateRg, h2ObsToLiab, parseSumstats } from "../ldsc/regress.js";
import { ctime, secToStr } from "../ldsc/log.js";
import { probability, hsqSummary } from "./ldherit.js";

const LDSC_POP = /^[A-Za-z0-9_+-]{1,64}$/;
export const CORR_FAILED = "The genetic correlation calculation failed. This can happen when the selected LD score source does not overlap with the summary statistics (e.g. different chromosome coverage or genome build).";

/** _parse_probability_pair: two comma-separated probabilities → "a,b" (each re-printed with ".15g") */
function probabilityPair(raw, field) {
  const parts = String(raw).split(",").map(s => s.trim());
  if (parts.length !== 2 || parts.some(p => p === "")) throw new Error(`Invalid ${field}: for liability scale in LDcorrelation, provide two comma-separated values between 0 and 1.`);
  return parts.map((p, i) => probability(p, `${field}[${i + 1}]`));
}

/** regressions.Gencov.summary (one annotation) */
function gencovSummary(g, P, K) {
  const liab = P && K && P.every(x => x != null) && K.every(x => x != null);
  const c = liab ? Math.sqrt(h2ObsToLiab(1, P[0], K[0])) * Math.sqrt(h2ObsToLiab(1, P[1], K[1])) : 1, T = liab ? "Liability" : "Observed";
  return [`Total ${T} scale gencov: ${npStr4(c * g.tot)} (${npStr4(c * g.tot_se)})`, `Mean z1*z2: ${npStr4(g.mean_z1z2)}`,
    `Intercept: ${npStr4(g.intercept)} (${npStr4(g.intercept_se)})`].join("\n");
}

/** regressions.RG.summary(silly=False), no constrained gencov intercept */
function rgSummary(r) {
  if (r.negative_hsq) return ["Genetic Correlation: nan (nan) (h2  out of bounds) ", "Z-score: nan (nan) (h2  out of bounds)", "P: nan (nan) (h2  out of bounds)",
    "WARNING: One of the h2's was out of bounds.", "This usually indicates a data-munging error or that h2 or N is low."].join("\n");
  if (r.rg_ratio > 1.2 || r.rg_ratio < -1.2) return ["Genetic Correlation: nan (nan) (rg out of bounds) ", "Z-score: nan (nan) (rg out of bounds)", "P: nan (nan) (rg out of bounds)",
    "WARNING: rg was out of bounds.", "This often means that h2 is not significantly different from zero."].join("\n");
  return [`Genetic Correlation: ${npStr4(r.rg_ratio)} (${npStr4(r.rg_se)})`, `Z-score: ${npStr4(r.z)}`, `P: ${npStr4(r.p)}`].join("\n");
}

/**
 * pandas DataFrame.to_string(index=False) of one row under ldsc's display.float_format '{:.4f}': numeric columns get a
 * space before the header and '{:.4f}' cells (NaN as "NaN"); string columns (file names, 'NA') as they are; each
 * column right-justified to its widest entry, joined with single spaces.
 */
export function pdRow(cols) {
  const cells = cols.map(([name, v]) => {
    const num = typeof v === "number", text = num ? (Number.isNaN(v) ? "NaN" : Number.isFinite(v) ? pyRoundStr(v, 4) : v > 0 ? "inf" : "-inf") : String(v);
    const head = num ? ` ${name}` : name, w = Math.max(head.length, text.length);
    return [head.padStart(w), text.padStart(w)];
  });
  return `${cells.map(c => c[0]).join(" ")}\n${cells.map(c => c[1]).join(" ")}`;
}

/** _get_rg_table for one pair */
function rgTable(r, files, P, K) {
  const liab = P && K && P.every(x => x != null) && K.every(x => x != null), c = liab ? h2ObsToLiab(1, P[1], K[1]) : 1, h = r.hsq2;
  return pdRow([["p1", files[0]], ["p2", files[1]], ["rg", r.rg_ratio], ["se", r.rg_se], ["z", r.z], ["p", r.p],
    [liab ? "h2_liab" : "h2_obs", c * h.tot], [liab ? "h2_liab_se" : "h2_obs_se", c * h.tot_se], ["h2_int", h.intercept], ["h2_int_se", h.intercept_se],
    ["gcov_int", r.gencov.intercept], ["gcov_int_se", r.gencov.intercept_se]]);
}

/**
 * @param {{ pop, filename, filename2, text, text2, isExample?, scale?, samp_prev?, pop_prev? }} a  text / text2: the two
 *   sumstats files' contents; filename / filename2: their names as LDlink sees them
 * → the REST text, or { error, details? }
 */
export async function ldcorrelation(ctx, { pop, filename, filename2, text, text2, isExample = false, scale = "observed", samp_prev = "", pop_prev = "" }) {
  if (scale != null && !["observed", "liability"].includes(String(scale))) return { error: "Invalid scale parameter." };
  if (pop && !LDSC_POP.test(String(pop).trim())) return { error: "Invalid pop parameter." };
  let P = null, K = null;
  if (String(scale).trim().toLowerCase() === "liability") {
    if (!String(samp_prev).trim() || !String(pop_prev).trim())
      return { error: "When scale=liability for LDcorrelation, both samp_prev and pop_prev are required as two comma-separated values between 0 and 1." };
    try { P = probabilityPair(samp_prev, "samp_prev").map(Number); K = probabilityPair(pop_prev, "pop_prev").map(Number); } catch (e) { return { error: e.message }; }
  }
  for (const [f, t] of [[filename, text], [filename2, text2]]) if (isExample && t == null) return { error: `Example file '${f}' was not found.` };
  const t1 = Date.now(), base = f => String(f).replace(/\.[^.]*$/, ""), lpop = String(pop).trim().toLowerCase();
  const [ma, counts] = await Promise.all([ctx.mergeAlleles(), ctx.ldscCounts()]);
  let m1, m2;
  try { m1 = mungeSumstats(text, ma); m2 = mungeSumstats(text2, ma); } catch (e) { return { error: `An error occurred while running the command (exit 1): ${e.message}` }; }
  // run_correlation_command's failure text has no "Second command failed" marker, so the route returns it as the result
  // (LDlink's details are ldsc.py's Python traceback; ours name the cause)
  if (!counts.has(lpop)) return `An error occurred while running the second command (exit 1): ldsc.py could not read the LD scores in ${String(pop).trim().toLowerCase()}/`;
  const need = new Set(parseSumstats(m1.text, { alleles: true }).map(r => r.SNP)), panel = await ctx.ldscPanel(lpop, need);
  const r = estimateRg(m1.text, m2.text, panel, panel.M);
  const files = [`${base(filename)}.sumstats.gz`, `${base(filename2)}.sumstats.gz`];
  const l = x => x + "-".repeat(x.replace(/\n/g, "").length); // sumstats._print_gencor's underline
  const out = [hsqSummary(r.hsq1, P?.[0], K?.[0]), l("\nHeritability of phenotype 2/2\n"), hsqSummary(r.hsq2, P?.[1], K?.[1]),
    l("\nGenetic Covariance\n"), gencovSummary(r.gencov, P, K), l("\nGenetic Correlation\n"), `${rgSummary(r)}\n`,
    `\nSummary of Genetic Correlation Results\n${rgTable(r, files, P, K)}\n`, `Analysis finished at ${ctime()}`, `Total time elapsed: ${secToStr((Date.now() - t1) / 1000)}`];
  return out.map(x => `${x}\n`).join("");
}
