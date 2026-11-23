# Stream a mongodump .bson.gz (stdin) and summarize: doc count, per-field type counts, "NA" counts,
# numeric min/max, distinct values for low-cardinality fields. Read-only; no external deps.
# usage: aws s3 cp s3://.../<coll>.bson.gz - | python -I pipeline/profile/bson_profile.py <coll> [max_docs]
import sys, struct, gzip, json, collections
name = sys.argv[1]; max_docs = int(sys.argv[2]) if len(sys.argv) > 2 else 0

def parse(b, i, end):
    d = {}
    while i < end:
        t = b[i]; i += 1; j = b.index(0, i); k = b[i:j].decode(); i = j + 1
        if t == 1: v = struct.unpack_from('<d', b, i)[0]; i += 8
        elif t == 2: l = struct.unpack_from('<i', b, i)[0]; v = b[i+4:i+3+l].decode('utf8', 'replace'); i += 4 + l
        elif t in (3, 4):
            l = struct.unpack_from('<i', b, i)[0]; v = parse(b, i + 4, i + l - 1); i += l
            if t == 4: v = list(v.values())
        elif t == 7: v = None; i += 12; k = "_id" if k == "_id" else k
        elif t == 8: v = bool(b[i]); i += 1
        elif t == 9: v = struct.unpack_from('<q', b, i)[0]; i += 8
        elif t == 10: v = None
        elif t == 16: v = struct.unpack_from('<i', b, i)[0]; i += 4
        elif t == 18: v = struct.unpack_from('<q', b, i)[0]; i += 8
        else: raise ValueError(t)
        d[k] = v
    return d

types = collections.defaultdict(collections.Counter); na = collections.Counter()
rng = {}; distinct = collections.defaultdict(collections.Counter); n = 0; maxlen = collections.Counter()
f = gzip.open(sys.stdin.buffer)
while True:
    h = f.read(4)
    if len(h) < 4: break
    l = struct.unpack('<i', h)[0]; b = h + f.read(l - 4)
    d = parse(b, 4, l - 1); n += 1
    for k, v in d.items():
        if k == "_id": continue
        types[k][type(v).__name__] += 1
        if v in ("NA", "", None): na[k] += 1
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            lo, hi = rng.get(k, (v, v)); rng[k] = (min(lo, v), max(hi, v))
        elif isinstance(v, str):
            maxlen[k] = max(maxlen[k], len(v))
            try:
                x = float(v); lo, hi = rng.get(k, (x, x)); rng[k] = (min(lo, x), max(hi, x))
            except ValueError: pass
        if isinstance(v, list): distinct[k + "[len]"][min(len(v), 50)] += 1
        elif len(distinct[k]) < 60: distinct[k][str(v)[:40]] += 1
        elif str(v)[:40] in distinct[k]: distinct[k][str(v)[:40]] += 1
        else: distinct[k]["<more>"] += 1
    if max_docs and n >= max_docs: break
low_card = {k: dict(c.most_common(30)) for k, c in distinct.items() if "<more>" not in c}
print(json.dumps({"collection": name, "docs": n, "types": types, "na_or_empty": na, "numeric_range": rng,
                  "max_str_len": maxlen, "low_cardinality": low_card}, default=dict))
