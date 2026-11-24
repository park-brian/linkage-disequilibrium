// LDlink's row choice for a resolved variant (LDcommon.py:590-603): tabix chr:p-p returns every VCF row whose REF
// overlaps p and LDlink keeps the LAST one (the POS check there is dead code); with MAF tiers split into files, "last"
// is the highest vcf_order. A source record that was multi-allelic → LDlink's "Variant is not a biallelic." error.
import { lastOverlapping } from "../core/genotypes.js";

/** LDlink parity: { status: "ok" | "not_biallelic" | "not_in_1000g" ("Variant is not in 1000G reference panel."), row } */
export function pickParity(rows, p) {
  const best = lastOverlapping(rows, p);
  if (!best) return { status: "not_in_1000g", row: null };
  return { status: best.multiallelic_source ? "not_biallelic" : "ok", row: best };
}
