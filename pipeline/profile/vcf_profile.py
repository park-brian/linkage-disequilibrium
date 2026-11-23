# Read-only profile of a 1000G VCF.gz: header, sample order, record classes, GT forms.
# usage: python -I pipeline/profile/vcf_profile.py <vcf.gz> <ALL.txt> [max_records] [gt_sample_records]
import gzip, sys, json, collections, re
path, all_txt = sys.argv[1], sys.argv[2]
max_rec = int(sys.argv[3]) if len(sys.argv) > 3 else 0
gt_n = int(sys.argv[4]) if len(sys.argv) > 4 else 3000
panel = [l.strip() for l in open(all_txt) if l.strip()]
out = {"file": path.replace("\\", "/").split("/")[-1], "meta": collections.Counter(), "info_ids": [], "format_ids": []}
c = collections.Counter(); filt = collections.Counter(); gtf = collections.Counter(); fmt = collections.Counter()
prev = None; dup = 0; first = last = None; idform = collections.Counter(); per_sample_ploidy = None
with gzip.open(path, "rb") as f:
    for line in f:
        if line.startswith(b"##"):
            k = line[2:].split(b"=", 1)[0].decode(); out["meta"][k] += 1
            if k in ("INFO", "FORMAT"):
                m = re.search(rb"ID=([^,>]+)", line); out[k.lower() + "_ids"].append(m.group(1).decode())
            continue
        if line.startswith(b"#CHROM"):
            cols = line.rstrip(b"\n").decode().split("\t"); samples = cols[9:]
            out["n_samples"] = len(samples); out["first_samples"] = samples[:3]
            out["in_panel"] = sum(s in set(panel) for s in samples); out["panel_size"] = len(panel)
            out["panel_missing_from_vcf"] = len(set(panel) - set(samples))
            out["same_order_as_panel"] = samples[:len(panel)] == panel
            continue
        p = line.split(b"\t", 9)
        chrom, pos, vid, ref, alt, _, fl, info, fm = p[:9]
        pos = int(pos); c["records"] += 1
        if first is None: first = (chrom.decode(), pos)
        last = (chrom.decode(), pos)
        if pos == prev: dup += 1
        prev = pos
        filt[fl.decode()] += 1; fmt[fm.decode()] += 1
        if b"," in alt: c["multiallelic"] += 1
        elif alt.startswith(b"<"): c["symbolic_sv"] += 1
        elif len(ref) == 1 and len(alt) == 1: c["snv"] += 1
        else: c["indel_or_mnp"] += 1
        idform["." if vid == b"." else ("rs" if vid.startswith(b"rs") else ("multi;" if b";" in vid else "other"))] += 1
        if b"END=" in info: c["has_END"] += 1
        if c["records"] <= gt_n:
            gts = [g.split(b":", 1)[0] for g in p[9].rstrip(b"\n").split(b"\t")]
            for g in gts: gtf[g.decode()] += 1
            if per_sample_ploidy is None and b"," not in alt:
                per_sample_ploidy = collections.Counter("haploid" if len(g) == 1 else "diploid" for g in gts)
        if max_rec and c["records"] >= max_rec: break
out.update(counts=c, duplicate_pos_adjacent=dup, first=first, last=last, filters=filt.most_common(5),
           formats=fmt.most_common(3), id_forms=idform, gt_forms_first_records=gtf.most_common(12),
           ploidy_first_biallelic=per_sample_ploidy)
print(json.dumps(out, default=dict))
