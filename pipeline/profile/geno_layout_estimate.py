# Estimate storage for candidate per-variant haplotype encodings on a 1000G VCF slice. Read-only.
# usage: python -I pipeline/profile/geno_layout_estimate.py <vcf.gz> [n_records] [block]
import gzip, sys, zlib, struct
try:
    from compression import zstd  # Python 3.14+
except ImportError:
    zstd = None
path = sys.argv[1]; N = int(sys.argv[2]) if len(sys.argv) > 2 else 100000; BLOCK = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
LUT = {b"0|0": (0, 0), b"0|1": (0, 1), b"1|0": (1, 0), b"1|1": (1, 1), b"0": (0, None), b"1": (1, None)}
vcf_bytes = n = skipped = 0; bits_blk = []; sparse_blk = []; tot = dict(bits_raw=0, sparse_raw=0, bits_zlib=0, sparse_zlib=0, bits_zstd=0, sparse_zstd=0)
maf_bins = [0, 0, 0, 0]
def flush():
    for key, blk in (("bits", bits_blk), ("sparse", sparse_blk)):
        raw = b"".join(blk); tot[key + "_raw"] += len(raw); tot[key + "_zlib"] += len(zlib.compress(raw, 6))
        if zstd: tot[key + "_zstd"] += len(zstd.compress(raw, 3))
        blk.clear()
with gzip.open(path, "rb") as f:
    for line in f:
        if line[0] == 35: continue
        vcf_bytes += len(line)
        p = line.rstrip(b"\n").split(b"\t")
        if b"," in p[4] or p[4].startswith(b"<"): skipped += 1; continue
        gts = p[9:]; H = 2 * len(gts); v = bytearray((H + 7) // 8); alt_idx = []
        for s, g in enumerate(gts):
            a = LUT.get(g.split(b":", 1)[0], (None, None))
            for h, x in enumerate(a):
                if x == 1:
                    k = 2 * s + h; v[k >> 3] |= 1 << (k & 7); alt_idx.append(k)
        ac = len(alt_idx); minor = alt_idx if ac <= H - ac else [k for k in range(H) if not (v[k >> 3] >> (k & 7)) & 1]
        maf = min(ac, H - ac) / H; maf_bins[0 if maf < 0.001 else 1 if maf < 0.01 else 2 if maf < 0.05 else 3] += 1
        bits_blk.append(bytes(v))
        # sparse: flag byte (which allele is minor) + uint16 count + delta-encoded uint16 indices
        d = []; prev = 0
        for k in minor: d.append(k - prev); prev = k
        sparse_blk.append(struct.pack("<BH", ac > H - ac, len(minor)) + struct.pack("<%dH" % len(d), *d))
        n += 1
        if n % BLOCK == 0: flush()
        if n >= N: break
flush()
print({"records": n, "skipped_multi_or_sv": skipped, "haplotypes": H, "vcf_uncompressed_bytes": vcf_bytes,
       "maf_bins(<0.1%,<1%,<5%,>=5%)": maf_bins, **{k: v for k, v in tot.items()},
       "bytes_per_variant": {k: round(v / n, 1) for k, v in tot.items()}})
