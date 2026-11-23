# Read-only: for the first N biallelic records of a 1000G VCF, compare max MAF over the 5 super-populations
# with max MAF over the 26 sub-populations, to see what a super-pop-based "common" tier would miss.
# usage: python -I tier_check.py <vcf.gz> <panelDir> [N]
import gzip, sys, os, collections
vcf, pdir = sys.argv[1], sys.argv[2]; N = int(sys.argv[3]) if len(sys.argv) > 3 else 200000
SUPER = ["AFR", "AMR", "EAS", "EUR", "SAS"]
SUB = [f[:-4] for f in os.listdir(pdir) if f.endswith(".txt") and f[:-4] not in SUPER + ["ALL"]]
panel = {p: set(open(os.path.join(pdir, p + ".txt")).read().split()) for p in SUPER + SUB}
c = collections.Counter(); n = 0
with gzip.open(vcf, "rt") as f:
    for line in f:
        if line.startswith("##"): continue
        if line.startswith("#"):
            samples = line.rstrip("\n").split("\t")[9:]
            idx = {p: [i for i, s in enumerate(samples) if s in panel[p]] for p in panel}; continue
        p = line.rstrip("\n").split("\t")
        if "," in p[4] or p[4].startswith("<"): continue
        gts = p[9:]
        def maf(ix):
            a = sum(gts[i].count("1") for i in ix); t = 2 * len(ix)
            return min(a, t - a) / t
        sup = max(maf(idx[q]) for q in SUPER); sub = max(maf(idx[q]) for q in SUB)
        for thr in (0.05, 0.01):
            if sup >= thr: c[f"super>={thr}"] += 1
            if sub >= thr: c[f"sub>={thr}"] += 1
            if sub >= thr and sup < thr: c[f"missed_by_super@{thr}"] += 1
        n += 1
        if n >= N: break
print({"records": n, **c})
