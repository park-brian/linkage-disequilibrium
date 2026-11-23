// Python number formatting/rounding, bit-for-bit (LDlink and ldsc print and threshold with these).
//   pyRoundStr(x, nd) = '%.{nd}f' % x     pyRound(x, nd) = round(x, nd)     nextUp(x) = math.nextafter(x, inf)
// Both '%.nf' and round(x, n) round the double's EXACT binary value to n decimals, ties to even, and keep the sign of
// values that round to zero. The decision is made on the exact decimal expansion: toFixed(100) is exact for
// |x| < 1e21 whenever the double has ≤ 100 fractional digits (every double ≥ 2^-46 in magnitude), and smaller ones
// round to zero for nd ≤ 30 anyway.

export function pyRoundStr(x, nd) {
  if (!Number.isFinite(x)) return Number.isNaN(x) ? "nan" : x > 0 ? "inf" : "-inf";
  const neg = x < 0 || Object.is(x, -0), a = Math.abs(x);
  if (a >= 1e21) return (neg ? "-" : "") + BigInt(Math.round(a)).toString() + (nd > 0 ? "." + "0".repeat(nd) : ""); // integral here
  const [ip, fp] = a.toFixed(100).split(".");
  const head = BigInt(ip + fp.slice(0, nd)), rest = fp.slice(nd).replace(/0+$/, "");
  // rest is the exact remainder after digit nd: > "5" → up, == "5" → tie → even, < "5" → down
  const up = rest === "" ? false : rest[0] > "5" || (rest[0] === "5" && (rest.length > 1 || head % 2n === 1n));
  const v = (head + (up ? 1n : 0n)).toString().padStart(nd + 1, "0");
  return `${neg ? "-" : ""}${nd > 0 ? `${v.slice(0, -nd)}.${v.slice(-nd)}` : v}`;
}
/** round(x, nd) for nd >= 0: the double nearest the correctly rounded decimal (Python returns exactly that) */
export const pyRound = (x, nd) => Number(pyRoundStr(x, nd));

const f64 = new Float64Array(1), u64 = new BigUint64Array(f64.buffer);
/** math.nextafter(x, inf): the next double above x */
export function nextUp(x) {
  if (Number.isNaN(x) || x === Infinity) return x;
  if (x === 0) return Number.MIN_VALUE;
  f64[0] = x; u64[0] += x > 0 ? 1n : -1n;
  return f64[0];
}

/** str(x) / repr(x) for a Python float: shortest round-trip digits, exponent form below 1e-4 and from 1e16 */
export function pyStr(x) {
  if (!Number.isFinite(x)) return Number.isNaN(x) ? "nan" : x > 0 ? "inf" : "-inf";
  if (Object.is(x, -0)) return "-0.0";
  const a = Math.abs(x);
  if (a !== 0 && (a < 1e-4 || a >= 1e16)) { const [m, e] = x.toExponential().split("e"); return `${m}e${e[0]}${e.slice(1).padStart(2, "0")}`; }
  const s = String(x); return s.includes(".") ? s : `${s}.0`;
}
/** str(round(x, nd)) */
export const pyRoundRepr = (x, nd) => pyStr(pyRound(x, nd));
/** C printf("%.{p}g") (GTEx's text format; every GTEx v8 slope/pval_nominal string round-trips through %.6g) */
export function cFormatG(x, p = 6) {
  if (!Number.isFinite(x)) return Number.isNaN(x) ? "nan" : x > 0 ? "inf" : "-inf";
  if (x === 0) return Object.is(x, -0) ? "-0" : "0";
  const [m, e] = x.toExponential(p - 1).split("e"), ex = +e;
  if (ex < -4 || ex >= p) return `${m.replace(/\.?0+$/, "")}e${ex < 0 ? "-" : "+"}${String(Math.abs(ex)).padStart(2, "0")}`;
  return x.toFixed(Math.max(0, p - 1 - ex)).replace(/(\.\d*?)0+$/, "$1").replace(/\.$/, "");
}

/**
 * numpy's str() of a float under ldsc's np.set_printoptions(precision=4) (floatmode "maxprec", the default), with the
 * brackets removed: ldsc's regressions.s(). At most 4 decimals, trailing zeros dropped but the point kept ("1."),
 * scientific when a nonzero |x| < 1e-4 or ≥ 1e8 ("5.e-05", "1.2346e+08"); nan / inf as numpy writes them.
 * Checked against numpy 2.4.3 (test/fixtures/numpy_s4.json).
 */
export function npStr4(x) {
  if (Number.isNaN(x)) return "nan";
  if (!Number.isFinite(x)) return x > 0 ? "inf" : "-inf";
  const a = Math.abs(x), sign = x < 0 || Object.is(x, -0) ? "-" : "";
  const trim = m => m.replace(/(\.\d*?)0+$/, "$1");
  if (a !== 0 && (a < 1e-4 || a >= 1e8)) { // (an exact tie at 5 significant digits needs a short binary fraction: not reachable here)
    const [m, e] = a.toExponential(4).split("e"), n = +e;
    return `${sign}${trim(m)}e${n < 0 ? "-" : "+"}${String(Math.abs(n)).padStart(2, "0")}`;
  }
  return sign + trim(pyRoundStr(a, 4)); // exact value, ties to even, as numpy's dragon4
}
