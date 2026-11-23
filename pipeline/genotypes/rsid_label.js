// rsID labels for 1000G genotype rows from LDlink's dbSNP dump (+ optional NCBI allele matches). Used by
// build_genotypes_v2.js; pure functions so they can be tested (test/phaseA.test.js).
//
// Where LDlink's dump stores a variant (PARQUET_LAYOUT.md §15): the first number of the REF allele's HGVS, i.e.
// POS for SNV/MNV and plain insertions, POS+1 for deletions/delins (incl. insertions inside repeats). 1000G rows are
// left-aligned VCF records (anchor base for indels), so a row's rsID can sit at POS or POS+1 in the dump. Rows split
// from multi-allelic records keep the padded REF (REF=AT, ALT=GT is really the SNV A>G at POS), so the variant class
// and the SNV position come from the minimal form.
//
// label() returns { rsid, source, delta, candidates, func, vcfRs, posType, alleles }:
//   source 0 none | 1 VCF rsID (current) | 2 VCF rsID merged → current | 3 position + type | 4 NCBI allele match
//   delta  dump position − row POS for the chosen rsID (0 or 1; null when unlabelled). A reader resolving an rsID via
//          the label must use pos + delta to agree with the dump (LDlink's find_one position).
export const TYPES = ["snv", "ins", "del", "delins", "mnv"]; // dbType codes; 255 = unknown (compatible with nothing)

/** trim the common suffix, then the common prefix (POS moves right); at least one base stays on each side */
export function minimalForm(pos, ref, alt) {
  let r = ref, a = alt, p = pos;
  while (r.length > 1 && a.length > 1 && r.at(-1) === a.at(-1)) { r = r.slice(0, -1); a = a.slice(0, -1); }
  while (r.length > 1 && a.length > 1 && r[0] === a[0]) { r = r.slice(1); a = a.slice(1); p++; }
  return [p, r, a];
}
/** 0 snv | 1 insertion | 2 deletion | 3 delins | 4 mnv, on the minimal form */
export function variantClass(ref, alt) {
  return ref.length === 1 && alt.length === 1 ? 0 : ref.length === 1 ? 1 : alt.length === 1 ? 2 : ref.length === alt.length ? 4 : 3;
}
const compatible = (c, t) => (t === 255 ? false : c === 0 ? t === 0 : t !== 0); // snv ↔ snv; indel/mnv classes ↔ any non-snv type

/**
 * @param db { dbPos: Int32Array (sorted), dbRsid, dbRef: Float64Array, dbType: Uint8Array, dbFunc: (string|null)[] }
 * @param alleleIndex optional { match(pos, ref, alt) → ascending current rsIDs }
 */
export function makeLabeler({ dbPos, dbRsid, dbRef, dbType, dbFunc }, alleleIndex = null) {
  const n = dbPos.length;
  function range(pos) { // [lo, hi) rows at pos
    let a = 0, b = n;
    while (a < b) { const m = (a + b) >> 1; if (dbPos[m] < pos) a = m + 1; else b = m; }
    let e = a; while (e < n && dbPos[e] === pos) e++;
    return [a, e];
  }
  const firstDocOf = (rsid, lo, hi) => { for (let k = lo; k < hi; k++) if (dbRsid[k] === rsid) return k; return -1; };

  return function label(pos, vcfId, ref, alt) {
    const [mp, mr, ma] = minimalForm(pos, ref, alt), c = variantClass(mr, ma);
    // dump positions to search, in LDlink-rule preference order
    //   snv: minimal POS; mnv: minimal POS (NCBI MNVs carry no anchor); insertion: POS (VC=INS) then POS+1 (inside a
    //   repeat NCBI calls it delins → POS+1); deletion/delins: POS+1 then POS
    const places = c === 0 ? [mp] : [...new Set(c === 4 ? [mp, pos, mp + 1] : c === 1 ? [pos, pos + 1, mp, mp + 1] : [pos + 1, pos, mp + 1, mp])];
    const ranges = places.map(p => [p, ...range(p)]);
    const m = /rs(\d+)/.exec(vcfId || ""), vcfRs = m ? +m[1] : null;
    const alleles = alleleIndex ? alleleIndex.match(pos, ref, alt) : null;
    let candidates = 0;
    for (const [, lo, hi] of ranges) for (let i = lo; i < hi; i++) if (dbRsid[i] === dbRef[i] && compatible(c, dbType[i])) candidates++;
    const result = (rsid, source, p, lo, hi, posType) => {
      const fi = rsid === null ? -1 : firstDocOf(rsid, lo, hi);
      return { rsid, source, delta: rsid === null ? null : p - pos, candidates, func: fi >= 0 ? dbFunc[fi] : null, vcfRs, posType, alleles };
    };
    // position + type (lowest current rsID at the first place that has one) — also reported as posType
    let pt = null;
    for (const [p, lo, hi] of ranges) {
      let best = -1;
      for (let i = lo; i < hi; i++) if (dbRsid[i] === dbRef[i] && compatible(c, dbType[i]) && (best < 0 || dbRsid[i] < dbRsid[best])) best = i;
      if (best >= 0) { pt = { rsid: dbRsid[best], p, lo, hi }; break; }
    }
    // 1/2: the VCF's own rsID, wherever LDlink stored it
    if (vcfRs !== null) for (const [p, lo, hi] of ranges) for (let i = lo; i < hi; i++) if (dbRsid[i] === vcfRs) {
      const cur = dbRef[i];
      return result(cur, cur === vcfRs ? 1 : 2, p, lo, hi, pt?.rsid ?? null);
    }
    // 4: an NCBI allele match that the dump has as a current doc at one of the places (lowest rsID, first place)
    if (alleles?.length) for (const [p, lo, hi] of ranges) {
      let best = null;
      for (let i = lo; i < hi; i++) if (dbRsid[i] === dbRef[i] && alleles.includes(dbRsid[i]) && (best === null || dbRsid[i] < best)) best = dbRsid[i];
      if (best !== null) return result(best, 4, p, lo, hi, pt?.rsid ?? null);
    }
    // 3: position + type
    if (pt) return result(pt.rsid, 3, pt.p, pt.lo, pt.hi, pt.rsid);
    return result(null, 0, pos, 0, 0, null);
  };
}
