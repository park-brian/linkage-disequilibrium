// ldsc's log conventions (ldsc.py 3.0.1, the version LDlink runs): the masthead, the "Call:" header of non-default
// options, time.ctime(), sec_to_str(), and Python's str() of a dict of strings (munge_sumstats prints one).
import { pyStr, pyRound } from "../pyfmt.js";

const VERSION = "3.0.1";
export const MASTHEAD = "*********************************************************************\n" +
  "* LD Score Regression (LDSC)\n" + `* Version ${VERSION}\n` + "* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane\n" +
  "* Broad Institute of MIT and Harvard / MIT Department of Mathematics\n" + "* GNU General Public License v3\n" +
  "*********************************************************************\n";

/**
 * MASTHEAD + "Call: \n./<script> \\\n" + one "--opt value \\" per option, the last backslash dropped (ldsc.py main /
 * munge_sumstats.py): options = [[name, value]] in argparse order.
 */
export function callHeader(script, options) {
  const body = options.map(([k, v]) => `--${k} ${v} \\`).join("\n");
  const h = `${MASTHEAD}Call: \n./${script} \\\n${body}`;
  return `${h.slice(0, -1)}\n`;
}

const DAYS = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"], MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
  "Oct", "Nov", "Dec"];
/** time.ctime(): "Mon Nov  3 07:33:52 2026" (local time; the day is space-padded) */
export function ctime(d = new Date()) {
  const p = n => String(n).padStart(2, "0");
  return `${DAYS[d.getDay()]} ${MONTHS[d.getMonth()]} ${String(d.getDate()).padStart(2,
    " ")} ${p(d.getHours())}:${p(d.getMinutes())}:${p(d.getSeconds())} ${d.getFullYear()}`;
}

/** sec_to_str(round(seconds, 2)): "4.64s", "1m:2.3s", "2h:0m:1.0s" (Python divmod on floats, str() of the remainder) */
export function secToStr(seconds) {
  let t = pyRound(seconds, 2);
  const parts = [];
  for (const b of [60, 60, 24]) { const r = t % b; parts.push(r); t = (t - r) / b; } // Python float divmod: floats throughout
  const [s, m, h] = parts, d = t;
  let f = "";
  if (d > 0) f += `${pyStr(d)}d:`;
  if (h > 0) f += `${pyStr(h)}h:`;
  if (m > 0) f += `${pyStr(m)}m:`;
  return `${f}${pyStr(s)}s`;
}

/** str() of a Python dict whose keys and values are plain strings */
export const pyDictRepr = obj => `{${Object.entries(obj).map(([k, v]) => `'${k}': '${v}'`).join(", ")}}`;
