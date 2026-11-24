// LDherit on the core: LDlink's REST text (server/LDlink.py ldherit route → ldscore.ldsc_utils.run_herit_command),
// i.e. munge_sumstats.py --sumstats <file> --merge-alleles w_hm3.snplist, then ldsc.py --h2 against the population's
// LD score panel (ldsc 3.0.1). The route returns the output from "Total Observed scale" onward; on the liability scale
// that marker is absent, so it returns everything: munge's stdout, "\n--------\n", ldsc's stdout.
// Four lines carry the time ("Conversion finished at", "Beginning analysis at", "Analysis finished at", "Total time
// elapsed"); they are written as ldsc would write them now, and the comparer masks them.
// Errors come back as { error } (LDlink's JSON), with details? for the 422 an ldsc failure produces.
import { cFormatG, pyRound, pyStr, npStr4 } from "../pyfmt.js";
import { mungeSumstats } from "../../ldsc/munge.js";
import { estimateH2, h2ObsToLiab, parseSumstats } from "../../ldsc/regress.js";
import { MASTHEAD, callHeader, ctime, secToStr, pyDictRepr } from "../ldsc/log.js";

const LDSC_POP = /^[A-Za-z0-9_+-]{1,64}$/;
const EXAMPLE_DIR = "/data/ldscore";
const DESCRIBE = {
  SNP: "Variant ID (e.g., rs number)", P: "p-Value", A1: "Allele 1, interpreted as ref allele for signed sumstat.",
  A2: "Allele 2, interpreted as non-ref allele for signed sumstat.", N: "Sample size", N_CAS: "Number of cases", N_CON: "Number of controls",
  Z: "Z-score (0 --> no effect; above 0 --> A1 is trait/risk increasing)", OR: "Odds ratio (1 --> no effect; above 1 --> A1 is risk increasing)",
  BETA: "[linear/logistic] regression coefficient (0 --> no effect; above 0 --> A1 is trait/risk increasing)",
  LOG_ODDS: "Log odds ratio (0 --> no effect; above 0 --> A1 is risk increasing)", INFO: "INFO score (imputation quality; higher --> better imputation)",
  FRQ: "Allele frequency", SIGNED_SUMSTAT: "Directional summary statistic as specified by --signed-sumstats.", NSTUDY: "Number of studies in which the SNP was genotyped.",
};
const HERIT_FAILED = "The heritability calculation failed. This can happen when the selected LD score source does not overlap with the summary statistics (e.g. different chromosome coverage or genome build).";

/** _parse_probability_value: float(...) strictly between 0 and 1, re-printed with format(x, ".15g") */
export function probability(raw, field) {
  const s = String(raw).trim(), x = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/.test(s) || /^[+-]?(inf|infinity|nan)$/i.test(s) ? Number(s.replace(/^([+-]?)inf(inity)?$/i, "$1Infinity")) : NaN;
  if (Number.isNaN(x) && !/^[+-]?nan$/i.test(s)) throw new Error(`Invalid ${field}: value must be a number between 0 and 1 (exclusive).`);
  if (!(x > 0 && x < 1)) throw new Error(`Invalid ${field}: value must be between 0 and 1 (exclusive).`);
  return cFormatG(x, 15);
}

/** munge_sumstats.py's stdout for LDlink's call (sumstatsPath as LDlink passes it; out = the file's base name) */
function mungeLog(m, sumstatsPath, base, t0 = Date.now()) {
  const L = m.log, out = [];
  out.push(callHeader("munge_sumstats.py", [["sumstats", sumstatsPath], ["out", base], ["merge-alleles", `${EXAMPLE_DIR}/w_hm3.snplist`]]));
  out.push(pyDictRepr(L.cnames)); // print(cname_translation)
  out.push("Interpreting column names as follows:");
  out.push(`${Object.entries(L.cnames).map(([k, v]) => `${k}:\t${DESCRIBE[v]}`).join("\n")}\n`);
  out.push(`Reading list of SNPs for allele merge from ${EXAMPLE_DIR}/w_hm3.snplist`);
  out.push(`Read ${L.merge_snps} SNPs for allele merge.`);
  out.push(`Reading sumstats from ${sumstatsPath} into memory 5000000 SNPs at a time.`);
  // filter warnings print inside the chunk loop, before its '.' (one chunk per 5,000,000 rows)
  const warn = [];
  if (L.bad_info > 0) warn.push(`WARNING: ${L.bad_info} SNPs had INFO outside of [0,1.5]. The INFO column may be mislabeled.`);
  if (L.bad_frq > 0) warn.push(`WARNING: ${L.bad_frq} SNPs had FRQ outside of [0,1]. The FRQ column may be mislabeled.`);
  if (L.removed_p > 0) warn.push(`WARNING: ${L.removed_p} SNPs had P outside of (0,1]. The P column may be mislabeled.`);
  const chunks = Math.max(1, Math.ceil(L.read / 5_000_000));
  out.push(`${warn.map(w => `${w}\n`).join("")}${".".repeat(chunks)} done`);
  out.push([`Read ${L.read} SNPs from --sumstats file.`, `Removed ${L.removed_merge} SNPs not in --merge-alleles.`, `Removed ${L.removed_na} SNPs with missing values.`,
    `Removed ${L.removed_info} SNPs with INFO <= 0.9.`, `Removed ${L.removed_maf} SNPs with MAF <= 0.01.`, `Removed ${L.removed_p} SNPs with out-of-bounds p-values.`,
    `Removed ${L.removed_alleles} variants that were not SNPs or were strand-ambiguous.`, `${L.remain} SNPs remain.`].join("\n"));
  const afterDups = L.remain - L.removed_dups;
  out.push(`Removed ${L.removed_dups} SNPs with duplicated rs numbers (${afterDups} SNPs remain).`);
  out.push(`Removed ${L.removed_n} SNPs with N < ${pyStr(L.n_min)} (${afterDups - L.removed_n} SNPs remain).`);
  out.push(`Median value of ${L.sign_name} was ${pyStr(L.median_signed)}, which seems sensible.`);
  out.push(`Removed ${L.removed_mismatch} SNPs whose alleles did not match --merge-alleles (${L.nonmissing} SNPs remain).`);
  out.push(`Writing summary statistics for ${L.written} SNPs (${L.nonmissing} with nonmissing beta) to ${base}.sumstats.gz.`);
  out.push("\nMetadata:");
  out.push(`Mean chi^2 = ${pyStr(pyRound(L.mean_chisq, 3))}`);
  if (L.mean_chisq < 1.02) out.push("WARNING: mean chi^2 may be too small.");
  out.push(`Lambda GC = ${pyStr(pyRound(L.median_chisq / 0.4549, 3))}`);
  out.push(`Max chi^2 = ${pyStr(pyRound(L.max_chisq, 3))}`);
  out.push(`${L.n_gws} Genome-wide significant SNPs (some may have been removed by filtering).`);
  out.push(`\nConversion finished at ${ctime()}`);
  out.push(`Total time elapsed: ${secToStr((Date.now() - t0) / 1000)}`);
  return out.map(l => `${l}\n`).join("");
}

/** regressions.Hsq.summary for one annotation (no overlap), with P / K for the liability scale */
export function hsqSummary(h, P, K) {
  const liab = P != null && K != null, c = liab ? h2ObsToLiab(1, P, K) : 1, T = liab ? "Liability" : "Observed", out = [];
  out.push(`Total ${T} scale h2: ${npStr4(c * h.tot)} (${npStr4(c * h.tot_se)})`);
  out.push(`Lambda GC: ${npStr4(h.lambda_gc)}`);
  out.push(`Mean Chi^2: ${npStr4(h.mean_chisq)}`);
  out.push(`Intercept: ${npStr4(h.intercept)} (${npStr4(h.intercept_se)})`);
  if (h.mean_chisq > 1) out.push(h.ratio < 0 ? "Ratio < 0 (usually indicates GC correction)." : `Ratio: ${npStr4(h.ratio)} (${npStr4(h.ratio_se)})`);
  else out.push("Ratio: NA (mean chi^2 < 1)");
  return out.join("\n");
}

/**
 * @param {{ pop, build?, filename, text, isExample?, scale?, samp_prev?, pop_prev? }} a  text: the sumstats file's
 *   contents (the example's, or the upload); filename: its name as LDlink sees it
 * → the REST text, or { error, details? }
 */
export async function ldherit(ctx, { pop, filename, text, isExample = false, scale = "observed", samp_prev = "", pop_prev = "" }) {
  // structural_input_guard (_validate_choice_value), then the route's own checks
  if (scale != null && !["observed", "liability"].includes(String(scale))) return { error: "Invalid scale parameter." };
  if (pop && !LDSC_POP.test(String(pop).trim())) return { error: "Invalid pop parameter." };
  let P = null, K = null;
  if (String(scale).trim().toLowerCase() === "liability") {
    if (!String(samp_prev).trim() || !String(pop_prev).trim()) return { error: "When scale=liability for LDherit, both samp_prev and pop_prev are required as values between 0 and 1." };
    try { P = probability(samp_prev, "samp_prev"); K = probability(pop_prev, "pop_prev"); } catch (e) { return { error: e.message }; }
  }
  if (isExample && text == null) return { error: `Example file '${filename}' was not found.` };
  const t0 = Date.now(), base = String(filename).replace(/\.[^.]*$/, ""), sumstatsPath = isExample ? `${EXAMPLE_DIR}/${filename}` : filename;
  const lpop = String(pop).trim().toLowerCase(), [ma, counts] = await Promise.all([ctx.mergeAlleles(), ctx.ldscCounts()]);
  let m;
  try { m = mungeSumstats(text, ma); } catch (e) { return { error: `An error occurred while running the command (exit 1): ${e.message}` }; }
  const munged = mungeLog(m, sumstatsPath, base, t0);
  const panel = counts.has(lpop) ? await ctx.ldscPanel(lpop, new Set(parseSumstats(m.text).map(r => r.SNP))) : null;
  if (!panel) return { error: HERIT_FAILED, details: `ldsc.py could not read the LD scores for ${pop}` };
  panel.rows = counts.get(lpop);
  // ldsc.py --h2
  const t1 = Date.now(), dir = `${EXAMPLE_DIR}/${String(pop).trim().toLowerCase()}/`, h = estimateH2(m.text, panel, panel.M);
  const opts = [["out", base], ["h2", `${base}.sumstats.gz`], ["ref-ld-chr", dir], ["w-ld-chr", dir], ...(P != null ? [["samp-prev", P], ["pop-prev", K]] : [])];
  const L = [callHeader("ldsc.py", opts), `Beginning analysis at ${ctime()}`, `Reading summary statistics from ${base}.sumstats.gz ...`,
    `RuntimeWarning: compression has no effect when passing a non-binary object as input for file ${base}.sumstats.gz`,
    `Read summary statistics for ${m.log.nonmissing} SNPs.`, `Reading reference panel LD Score from ${dir}[1-22] ... (ldscore_fromlist)`,
    `Read reference panel LD Scores for ${panel.rows} SNPs.`, "Removing partitioned LD Scores with zero variance.",
    `Reading regression weight LD Score from ${dir}[1-22] ... (ldscore_fromlist)`, `Read regression weight LD Scores for ${panel.rows} SNPs.`,
    `After merging with reference panel LD, ${h.n_snp} SNPs remain.`, `After merging with regression SNP LD, ${h.n_snp} SNPs remain.`];
  if (h.n_snp < 200000) L.push("WARNING: number of SNPs less than 200k; this is almost always bad.");
  L.push("Using two-step estimator with cutoff at 30.", hsqSummary(h, P == null ? null : +P, K == null ? null : +K),
    `Analysis finished at ${ctime()}`, `Total time elapsed: ${secToStr((Date.now() - t1) / 1000)}`);
  const result = `${munged}\n--------\n${L.map(l => `${l}\n`).join("")}`;
  const at = result.indexOf("Total Observed scale");
  return at === -1 ? result : result.slice(at);
}
