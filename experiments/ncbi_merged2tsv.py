# refsnp-merged.json.bz2 (NCBI dbSNP) → TSV: one row per (retired rsID, target rsID) in merged_into order.
# Columns: rsid, target, k (index in merged_into), n_targets. Stdlib only; run with python -I.
# usage: python -I ncbi_merged2tsv.py <refsnp-merged.json.bz2> <out.tsv>
import bz2
import json
import sys

src, out = sys.argv[1], sys.argv[2]
n = rows = bad = 0
with bz2.open(src, "rt", encoding="utf-8") as f, open(out, "w", encoding="utf-8", newline="\n") as o:
    for line in f:
        n += 1
        try:
            r = json.loads(line)
            tgt = r["merged_snapshot_data"]["merged_into"]
        except (ValueError, KeyError):
            bad += 1
            continue
        for k, t in enumerate(tgt):
            o.write(f"{r['refsnp_id']}\t{t}\t{k}\t{len(tgt)}\n")
            rows += 1
print(json.dumps({"records": n, "rows": rows, "unparsed": bad}))
