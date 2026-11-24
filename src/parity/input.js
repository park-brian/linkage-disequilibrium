// LDlink's input parsing (server/LDcommon.py): the first token of an input, and chr:pos as Python's strip("chr"),
// split(":") and int() read it. The parity layer's front door; the core takes chromosome codes and positions.
import { CHROMS, chromCode } from "../core/genome.js";

// validsnp (LDcommon.py:370): snp_raw.lower().strip().split(), and the modules use snp[0], the first token
export const firstToken = s => String(s ?? "").toLowerCase().trim().split(/\s+/)[0] ?? "";
// LDcommon.py:104-110 — only 1-22, X, Y (and lower-case x, y) are recognized
const LOOKUP = Object.fromEntries([...CHROMS.map(c => [c, c]), ["x", "X"], ["y", "Y"]]);
/**
 * LDlink's chr:pos parsing (validsnp + coord.strip("chr").split(":") + int()): → { chro, pos: BigInt, key } or null.
 * key is null when the position can't be a dbSNP key (≤ 0 or ≥ 1e9).
 */
export function parseCoord(s) {
  const t = firstToken(s).replace(/^[chr]+|[chr]+$/g, "").split(":"); // coord.strip("chr").split(":")
  if (t.length <= 1) return null;
  const chro = LOOKUP[t[0]]; if (!chro) return null;
  if (!/^\s*[+-]?\d+(_\d+)*\s*$/.test(t[1])) return null; // int(temp_coord[1]) (Python allows _ separators)
  const pos = BigInt(t[1].trim().replace(/_/g, ""));
  if (pos <= 0n || pos >= 1_000_000_000n) return { chro, pos, key: null };
  return { chro, pos, key: BigInt(chromCode(chro)) * 1_000_000_000n + pos };
}
