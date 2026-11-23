# Run a DuckDB SQL file one statement per process (state persists in the database file; SET lines are
# prepended to every statement). Stops at the first failing statement. usage: python -I run_sql_statements.py <db> <sql>
import subprocess, sys, time
db, path = sys.argv[1], sys.argv[2]
text = "\n".join(l for l in open(path, encoding="utf8").read().splitlines() if not l.strip().startswith("--"))
stmts = [s.strip() for s in text.split(";\n") if s.strip()]
sets = [s for s in stmts if s.upper().startswith(("SET ", "ATTACH "))]  # session state, repeated in every process
for s in [s for s in stmts if not s.upper().startswith(("SET ", "ATTACH "))]:
    t = time.time()
    r = subprocess.run(["duckdb", "-no-agent", "-bail", db, "-c", ";".join(sets + [s.rstrip(";")]) + ";"], capture_output=True, text=True)
    print(f"[{time.time() - t:6.1f}s exit={r.returncode}] {s[:90].replace(chr(10), ' ')}", flush=True)
    if r.returncode != 0:
        print(r.stdout[-1500:], r.stderr[-1500:]); sys.exit(1)
