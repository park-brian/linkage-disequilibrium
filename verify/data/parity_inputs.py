# V3 (ROADMAP §2.2): LDlink's input handling, as ground truth for resolve.mjs.
# Re-implements, expression for expression, server/LDcommon.py:
#   validsnp            snp = snp_raw.lower().strip().split()           (one entry per input line; modules use snp[0])
#   replace_coord_rsid  rs* kept; else get_rsnum: coord.strip("chr").split(":"), CHROMOSOME_LOOKUP, int(), find()
#                       → first doc / first id==ref_id doc, "Multiple rsIDs (...) map to genomic coordinates <snp>"
#   get_coords          rsid.strip("rs"); str(int(rsid)); find_one({"id"})  (first doc in dump order)
# and records LDmatrix's extra gate (LDmatrix.py:69-72: len > 2, starts with rs/chr, last char a digit) separately,
# because it is module logic, not resolver logic. Dump lookups go through DuckDB on dbsnp_dump/work.duckdb (relative to LD_DATA: verify/data/parity.mjs runs it there)
# (table d: ord = dump order, Mongo semantics as in parity_expected.sql). Stdlib only; run with python -I.
# usage (cwd = LD_DATA): python -I verify/data/parity_inputs.py <inputs.json> <build grch37|grch38> <out.json>
import json
import subprocess
import sys

inputs_path, build, out_path = sys.argv[1:4]
raw_inputs = json.load(open(inputs_path, encoding="utf8"))
POS = {"grch37": "position_grch37", "grch38": "position_grch38"}[build]
CHROMOSOME_LOOKUP = {name: name for name in [str(i) for i in range(1, 23)] + ["X", "Y"]}
CHROMOSOME_LOOKUP.update({"x": "X", "y": "Y"})


def rsid_key(rsid):  # get_coords, before the query
    rsid = rsid.strip("rs")
    try:
        return str(int(rsid))
    except ValueError:
        return None


def coord_key(coord):  # get_rsnum, before the query
    temp_coord = coord.strip("chr").split(":")
    if len(temp_coord) <= 1:
        return None
    chro = CHROMOSOME_LOOKUP.get(str(temp_coord[0]))
    if chro is None:
        return None
    try:
        pos = int(temp_coord[1])
    except ValueError:
        return None
    return (chro, str(pos))


# pass 1: the lookups LDlink would make
first_tokens = []
for snp_raw in raw_inputs:
    snp = snp_raw.lower().strip().split()
    first_tokens.append(snp[0] if snp else None)
rs_keys, co_keys = set(), set()
for t in first_tokens:
    if t is None:
        continue
    if t[0:2] == "rs":
        k = rsid_key(t)
        if k is not None:
            rs_keys.add(k)
    else:
        k = coord_key(t)
        if k is not None:
            co_keys.add(k)


def duck(sql):
    r = subprocess.run(["duckdb", "-no-agent", "-readonly", "-json", "dbsnp_dump/work.duckdb", "-c", sql], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:])
    return json.loads(r.stdout or "[]")


def sql_list(values):
    return ",".join("'" + v.replace("'", "''") + "'" for v in values) or "NULL"


# pass 2: dump docs (Mongo order = ord); coords may resolve to rsIDs that also need get_coords
co_docs = {}
if co_keys:
    conds = " OR ".join(f"(chromosome = '{c}' AND {POS} = '{p}')" for c, p in co_keys)
    for row in duck(f"SELECT chromosome, {POS} AS pos, id, ref_id, ord FROM d WHERE {conds} ORDER BY ord"):
        co_docs.setdefault((row["chromosome"], row["pos"]), []).append(row)


def replace_coord_rsid(snp):
    warning = None
    if snp[0:2] == "rs":
        return snp, warning
    k = coord_key(snp)
    snp_info_lst = co_docs.get(k, []) if k is not None else None
    if snp_info_lst is None:
        return snp, warning
    if len(snp_info_lst) > 1:
        var_id = "rs" + snp_info_lst[0]["id"]
        ref_variants = [s["id"] for s in snp_info_lst if s["id"] == s["ref_id"]]
        if len(ref_variants) > 1:
            var_id = "rs" + ref_variants[0]
            warning = "Multiple rsIDs (" + ", ".join(["rs" + r for r in ref_variants]) + ") map to genomic coordinates " + snp
        elif len(ref_variants) == 0 and len(snp_info_lst) > 1:
            var_id = "rs" + snp_info_lst[0]["id"]
            warning = "Multiple rsIDs (" + ", ".join(["rs" + r for r in ref_variants]) + ") map to genomic coordinates " + snp
        else:
            var_id = "rs" + ref_variants[0]
        return var_id, warning
    if len(snp_info_lst) == 1:
        return "rs" + snp_info_lst[0]["id"], warning
    return snp, warning


replaced = []
for t in first_tokens:
    replaced.append(None if t is None else replace_coord_rsid(t))
for r in replaced:
    if r is not None and r[0][0:2] == "rs":
        k = rsid_key(r[0])
        if k is not None:
            rs_keys.add(k)
rs_docs = {}
if rs_keys:
    for row in duck(f"SELECT id, arg_min(ref_id, ord) AS ref_id, arg_min(chromosome, ord) AS chromosome, arg_min({POS}, ord) AS pos FROM d WHERE id IN ({sql_list(sorted(rs_keys))}) GROUP BY id"):
        rs_docs[row["id"]] = row

out = []
for raw, t, r in zip(raw_inputs, first_tokens, replaced):
    e = {"input": raw, "token": t, "rsid": None, "ref_rsid": None, "chrom": None, "pos": None, "found": False, "warning": None, "ldmatrix_gate": False}
    if t is not None:
        var, warning = r
        e["warning"] = warning
        e["ldmatrix_gate"] = len(var) > 2 and (var[0:2] == "rs" or var[0:3] == "chr") and var[-1].isdigit()
        k = rsid_key(var) if var[0:2] == "rs" else None
        doc = rs_docs.get(k) if k is not None else None
        if doc is not None:
            e.update(found=True, rsid="rs" + doc["id"], ref_rsid=None if doc["ref_id"] == doc["id"] else "rs" + doc["ref_id"], chrom=doc["chromosome"],
                     pos=None if doc["pos"] == "NA" else int(doc["pos"]))
    out.append(e)
json.dump(out, open(out_path, "w", encoding="utf8"), indent=1)
print(json.dumps({"inputs": len(out), "found": sum(e["found"] for e in out), "warnings": sum(e["warning"] is not None for e in out),
                  "gate_rejects_found": sum(e["found"] and not e["ldmatrix_gate"] for e in out)}))
