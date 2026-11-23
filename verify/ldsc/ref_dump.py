# Run the reference ldsc (CBIIT/ldsc@ldsc313) h2 / rg exactly as ldsc.py does and dump full-precision results as JSON,
# so the JS port can be compared beyond the 4 printed digits.
# usage: <ldsc venv python> verify/ldsc/ref_dump.py <ldsc-ref dir> <out.json> h2 <sumstats.gz> <ref_ld_dir/> [<samp_prev> <pop_prev>]
#        <ldsc venv python> verify/ldsc/ref_dump.py <ldsc-ref dir> <out.json> rg <a.sumstats.gz,b.sumstats.gz> <ref_ld_dir/>
import json
import sys

ldsc_dir, out, mode, sumstats_arg, ld_dir = sys.argv[1:6]
extra = sys.argv[6:]
sys.path.insert(0, ldsc_dir)
sys.argv = [sys.argv[0]]  # ldsc.py parses argv at import time in some versions
import numpy as np  # noqa: E402
import ldsc  # noqa: E402
from ldscore import sumstats  # noqa: E402

argv = ["--" + mode, sumstats_arg, "--ref-ld-chr", ld_dir, "--w-ld-chr", ld_dir, "--out", out + ".ldsc"]
if extra:
    argv += ["--samp-prev", extra[0], "--pop-prev", extra[1]]
args = ldsc.parser.parse_args(argv)
log = ldsc.Logger(out + ".ldsc.log")
f = lambda x: None if x is None or isinstance(x, str) else float(np.asarray(x).reshape(-1)[0])


def hsq_dict(h, P=None, K=None):
    from ldscore import regressions as reg
    c = reg.h2_obs_to_liab(1, P, K) if P is not None and K is not None else 1
    return {"tot": f(h.tot), "tot_se": f(h.tot_se), "intercept": f(h.intercept), "intercept_se": f(h.intercept_se),
            "mean_chisq": f(getattr(h, "mean_chisq", None)), "lambda_gc": f(getattr(h, "lambda_gc", None)),
            "ratio": f(getattr(h, "ratio", None)), "ratio_se": f(getattr(h, "ratio_se", None)),
            "liab_factor": f(c), "summary": h.summary(None, P=P, K=K) if hasattr(h, "lambda_gc") else None,
            "twostep_filtered": None if getattr(h, "twostep_filtered", None) is None else int(h.twostep_filtered)}


if mode == "h2":
    P = float(extra[0]) if extra else None
    K = float(extra[1]) if extra else None
    h = sumstats.estimate_h2(args, log)
    res = {"mode": "h2", "n_snp": int(h.jknife.delete_values.shape[0]), "hsq": hsq_dict(h, P, K)}
else:
    rg = sumstats.estimate_rg(args, log)[0]
    res = {"mode": "rg", "hsq1": hsq_dict(rg.hsq1), "hsq2": hsq_dict(rg.hsq2),
           "gencov": {"tot": f(rg.gencov.tot), "tot_se": f(rg.gencov.tot_se), "intercept": f(rg.gencov.intercept),
                      "intercept_se": f(rg.gencov.intercept_se), "mean_z1z2": f(rg.gencov.mean_z1z2), "p": f(rg.gencov.p), "z": f(rg.gencov.z)},
           "rg_ratio": f(rg.rg_ratio), "rg_se": f(rg.rg_se), "z": f(rg.z), "p": f(rg.p), "summary": rg.summary()}
json.dump(res, open(out, "w"), indent=2)
print(json.dumps(res)[:600])
