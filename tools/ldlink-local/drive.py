# Call LDlink's module functions directly inside the backend container (full JSON, warnings included).
# usage: docker compose exec -T backend python3.13 - /data/tmp/<cases>.json [out.jsonl] < drive.py
#   → one JSON line per case, on stdout or (when the modules print, e.g. LDassoc) into out.jsonl
# case: {"module": "ldhap", "snps": [...], "pop": "CEU", "build": "grch37"}
import json, sys, uuid, os
sys.path.insert(0, "/opt/ldlink")
os.chdir("/opt/ldlink")
os.environ["PATH"] = "/opt/ldlink/apache-bin:" + os.environ["PATH"]  # as the image's CMD does (subprocesses call python3)
from LDhap import calculate_hap
from LDassoc import calculate_assoc
import argparse

TMP = os.environ["TMP_DIR"]
sink = open(sys.argv[2], "w") if len(sys.argv) > 2 else sys.stdout
emit = lambda obj: (sink.write(json.dumps(obj) + "\n"), sink.flush())
for case in json.load(open(sys.argv[1])):
    ref = str(uuid.uuid4())
    try:
        if case["module"] == "ldhap":
            path = os.path.join(TMP, "snps" + ref + ".txt")
            with open(path, "w") as f:
                f.write("\n".join(case["snps"]).lower())  # the API route lower-cases the request
            out = json.loads(calculate_hap(path, case["pop"], ref, False, case["build"]))
            text = None
            if "error" not in out:
                with open(os.path.join(TMP, "snps_" + ref + ".txt")) as a, open(os.path.join(TMP, "haplotypes_" + ref + ".txt")) as b:
                    text = a.read() + "\n" + "#" * 85 + "\n\n" + b.read()
            emit({"json": out, "text": text})
        elif case["module"] == "ldassoc":
            # the /LDlinkRestWeb/ldassoc route's argument mapping (LDlink.py ldassoc), example file
            a = argparse.Namespace(window=None, dprime=False, chr=case["columns"][0], bp=case["columns"][1], pval=case["columns"][2],
                                   transcript=False, annotate="forge", origin=case.get("origin"), name=case.get("name"),
                                   start=case.get("start"), end=case.get("end"))
            if case["region"] in ("variant", "gene"): a.window = int(case["window"])
            path = os.environ["DATA_DIR"] + os.environ["LDASSOC_EXAMPLE_DIR"] + case["file"]
            calculate_assoc(path, case["region"], case["pop"], ref, case["build"], True, a)
            out = json.load(open(os.path.join(TMP, "assoc" + ref + ".json")))
            text = None
            if "error" not in out and os.path.exists(os.path.join(TMP, "assoc" + ref + ".txt")):
                text = open(os.path.join(TMP, "assoc" + ref + ".txt")).read()
            emit({"json": {k: out[k] for k in ("error", "warning") if k in out}, "text": text})
        else:
            emit({"exception": "unknown module"})
    except Exception as e:
        emit({"exception": type(e).__name__ + ": " + str(e)})
