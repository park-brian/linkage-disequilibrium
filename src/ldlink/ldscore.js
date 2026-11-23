// LDscore on the core: LDlink's REST text (server/LDlink.py ldscore route → ldsc_utils.run_ldsc_command), i.e.
// ldsc.py --bfile <base> --l2 --ld-wind-cm|--ld-wind-kb <window> --out <base> on an uploaded PLINK fileset, computed
// by src/ldsc/l2.js. The route returns the log from "Summary of LD Scores" onward: pandas' describe() of MAF and L2
// (count dropped), the MAF/L2 correlation matrix, and ldsc's two closing time lines (written now; the comparer masks
// them). Errors come back as { error }.
import { readBim, readFamCount, ldscoreL2 } from "../ldsc/l2.js";
import { pdFrame, describe, corr } from "../ldsc/pandas.js";
import { ctime, secToStr } from "../ldsc/log.js";

const PYINT = /^\s*[+-]?\d+(_\d+)*\s*$/; // what Python's int() accepts

/**
 * @param {{ filename: string, ldwindow?, windUnit?, isExample?, files: { [name]: Uint8Array | string } }} a  filename:
 *   LDlink's comma (or semicolon) separated list, e.g. "22.bed,22.bim,22.fam"; files: their contents by name (null if
 *   absent)
 * → the REST text, or { error }
 */
export function ldscore({ filename = "", ldwindow = "1", windUnit = "cm", isExample = false, files }) {
  // structural_input_guard: _validate_number_value(minimum=0), _validate_choice_value({cm, cM, kb})
  const w = Number(String(ldwindow).trim());
  if (String(ldwindow).trim() === "" || Number.isNaN(w) || w < 0) return { error: "Invalid ldwindow parameter." };
  if (!["cm", "cM", "kb"].includes(String(windUnit))) return { error: "Invalid windUnit parameter." };
  const names = String(filename).replace(/;/g, ",").split(",").map(f => f.trim()).filter(Boolean);
  const base = names.length ? names[0].replace(/\.[^.]*$/, "") : "22";
  // LDlink builds an "All uploaded files must have the same base name" message for a mixed list but only logs it; the
  // request goes on to copy each example (failing on the first absent one) and then checks <base>.bed/.bim/.fam
  for (const f of names) if (isExample && files[f] == null) return { error: `Example file '${f}' was not found.` };
  const missing = [".bed", ".bim", ".fam"].map(e => `${base}${e}`).filter(f => files[f] == null);
  if (missing.length) return { error: `Missing files: ${missing.join(", ")}` };
  // run_ldsc_command: int(ldwindow) or 1, at least 1; 'kb' → --ld-wind-kb, anything else --ld-wind-cm
  const win = PYINT.test(String(ldwindow)) && +String(ldwindow).replace(/_/g, "") > 0 ? +String(ldwindow).replace(/_/g, "") : 1;
  const t0 = Date.now();
  const { rows } = ldscoreL2(files[`${base}.bed`], readBim(files[`${base}.bim`]), readFamCount(files[`${base}.fam`]), windUnit === "kb" ? { windKb: win } : { windCm: win });
  const maf = rows.map(r => r.MAF), l2 = rows.map(r => r.L2), dm = describe(maf), dl = describe(l2);
  const stats = ["mean", "std", "min", "25%", "50%", "75%", "max"];
  const c = corr(maf, l2);
  const out = [`Summary of LD Scores in ${base}.l2.ldscore.gz`, pdFrame(["MAF", "L2"], stats.map(k => [k, dm[k], dl[k]])), "",
    "MAF/LD Score Correlation Matrix", pdFrame(["MAF", "L2"], [["MAF", corr(maf, maf), c], ["L2", c, corr(l2, l2)]]),
    `Analysis finished at ${ctime()}`, `Total time elapsed: ${secToStr((Date.now() - t0) / 1000)}`];
  return out.map(x => `${x}\n`).join("");
}
