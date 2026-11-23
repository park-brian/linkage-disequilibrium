// Haplotype storage decoding, with no dependencies (safe to import from a Web Worker).
/** Decode a row's haplotypes (bitset tier, or sparse tier: [listRef, varint deltas…]) into a bitset of hbytes. */
export function decodeHaps(haps, hbytes) {
  if (haps.length === hbytes) return haps;
  const out = new Uint8Array(hbytes), listRef = haps[0] === 1, idx = [];
  for (let k = 0, o = 1; o < haps.length;) { let d = 0, s = 0, b; do { b = haps[o++]; d |= (b & 0x7f) << s; s += 7; } while (b & 0x80); k += d; idx.push(k); }
  if (listRef) { out.fill(0xff); for (const x of idx) out[x >> 3] &= ~(1 << (x & 7)); } else for (const x of idx) out[x >> 3] |= 1 << (x & 7);
  return out;
}
