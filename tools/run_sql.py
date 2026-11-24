# Run a DuckDB SQL file one statement per process (state persists in the database file; SET lines are prepended to
# every statement). Stops at the first failing statement. Runs with the working directory at LD_DATA (as every
# pipeline/verify script does), so paths inside the SQL and the <db> argument are data-relative; <sql> is resolved
# from where you call it. A block "@for c in 1..25" … "@end" repeats its statements for each value: {c} is the value and
# {c.chr} its chromosome name (23 → X, 24 → Y, 25 → MT). --print lists the expanded statements without running them.
# usage: python -I tools/run_sql.py <db (data-relative)> <sql file> [--print]
import os, re, subprocess, sys, time
repo = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
data = os.path.abspath(os.environ.get("LD_DATA") or os.path.join(repo, "..", "data"))
db, path, only_print = sys.argv[1], os.path.abspath(sys.argv[2]), "--print" in sys.argv[3:]

def expand(lines):
    out, i = [], 0
    while i < len(lines):
        m = re.fullmatch(r"@for (\w+) in (\d+)\.\.(\d+)", lines[i].strip())
        if not m:
            out.append(lines[i]); i += 1; continue
        var, lo, hi, j = m.group(1), int(m.group(2)), int(m.group(3)), i + 1
        while lines[j].strip() != "@end": j += 1
        for v in range(lo, hi + 1):
            name = {23: "X", 24: "Y", 25: "MT"}.get(v, str(v))
            out += [l.replace("{%s.chr}" % var, name).replace("{%s}" % var, str(v)) for l in lines[i + 1:j]]
        i = j + 1
    return out

lines = [l for l in open(path, encoding="utf8").read().splitlines() if not l.strip().startswith("--")]
text = "\n".join(expand(lines))
stmts = [s.strip() for s in text.split(";\n") if s.strip()]
if only_print:
    for s in stmts: print(s.rstrip(";") + ";")
    sys.exit(0)
os.makedirs(data, exist_ok=True); os.chdir(data)
sets = [s for s in stmts if s.upper().startswith(("SET ", "ATTACH "))]  # session state, repeated in every process
for s in [s for s in stmts if not s.upper().startswith(("SET ", "ATTACH "))]:
    t = time.time()
    r = subprocess.run(["duckdb", "-no-agent", "-bail", db, "-c", ";".join(sets + [s.rstrip(";")]) + ";"], capture_output=True, text=True)
    print(f"[{time.time() - t:6.1f}s exit={r.returncode}] {s[:90].replace(chr(10), ' ')}", flush=True)
    if r.returncode != 0:
        print(r.stdout[-1500:], r.stderr[-1500:]); sys.exit(1)
