# Run a DuckDB SQL file one statement per process (state persists in the database file; SET lines are prepended to
# every statement). Stops at the first failing statement. Runs with the working directory at LD_DATA (as every
# pipeline/verify script does), so paths inside the SQL and the <db> argument are data-relative; <sql> is resolved
# from where you call it. usage: python -I tools/run_sql.py <db (data-relative)> <sql file>
import os, subprocess, sys, time
repo = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
data = os.path.abspath(os.environ.get("LD_DATA") or os.path.join(repo, "..", "data"))
db, path = sys.argv[1], os.path.abspath(sys.argv[2])
os.makedirs(data, exist_ok=True); os.chdir(data)
text = "\n".join(l for l in open(path, encoding="utf8").read().splitlines() if not l.strip().startswith("--"))
stmts = [s.strip() for s in text.split(";\n") if s.strip()]
sets = [s for s in stmts if s.upper().startswith(("SET ", "ATTACH "))]  # session state, repeated in every process
for s in [s for s in stmts if not s.upper().startswith(("SET ", "ATTACH "))]:
    t = time.time()
    r = subprocess.run(["duckdb", "-no-agent", "-bail", db, "-c", ";".join(sets + [s.rstrip(";")]) + ";"], capture_output=True, text=True)
    print(f"[{time.time() - t:6.1f}s exit={r.returncode}] {s[:90].replace(chr(10), ' ')}", flush=True)
    if r.returncode != 0:
        print(r.stdout[-1500:], r.stderr[-1500:]); sys.exit(1)
