# One-off: move the runtime core from experiments/ to src/ and rewrite every relative import that pointed at a
# moved file (docs/CORE_REVIEW.md §3). Run from the repo root: python -I scripts/move_to_src.py [--dry-run]
import os, re, sys, shutil

ROOT = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
E, S = os.path.join(ROOT, "experiments"), os.path.join(ROOT, "src")
MOVES = {  # experiments-relative → src-relative
    "lookup.mjs": "parquet.mjs",
    "resolve.mjs": "dbsnp.mjs",
    "planner.mjs": "planner.mjs",
    "select_rows.mjs": "genotypes.mjs",
    "pyfmt.mjs": "pyfmt.mjs",
    "wasm/ld_kernel.mjs": "ld/kernel.mjs",
    "wasm/ld_kernel_js.mjs": "ld/kernel_js.mjs",
    "wasm/ld_kernel.c": "ld/kernel.c",
    "wasm/ld_kernel.wasm": "ld/kernel.wasm",
    "wasm/build.sh": "ld/build.sh",
    "ldsc/ldsc_regress.mjs": "ldsc/regress.mjs",
    "ldsc/munge.mjs": "ldsc/munge.mjs",
    "ldsc/ldscore_l2.mjs": "ldsc/l2.mjs",
    "modules/context.mjs": "ldlink/context.mjs",
    "modules/common.mjs": "ldlink/common.mjs",
    **{f"modules/{m}.mjs": f"ldlink/{m}.mjs" for m in ["ldpair", "ldpop", "ldhap", "ldproxy", "snpclip", "snpchip", "ldtrait", "ldexpress", "ldassoc"]},
}
old2new = {os.path.normpath(os.path.join(E, a)): os.path.normpath(os.path.join(S, b)) for a, b in MOVES.items()}
dry = "--dry-run" in sys.argv

def new_location(path):
    return old2new.get(os.path.normpath(path), os.path.normpath(path))

IMPORT = re.compile(r"""((?:from|import)\s*\(?\s*["'])(\.{1,2}/[^"']+)(["'])""")

def rewrite(src_path_old, text):
    here_new = os.path.dirname(new_location(src_path_old))
    def sub(m):
        target_old = os.path.normpath(os.path.join(os.path.dirname(src_path_old), m.group(2)))
        target_new = new_location(target_old)
        rel = os.path.relpath(target_new, here_new).replace("\\", "/")
        if not rel.startswith("."): rel = "./" + rel
        return m.group(1) + rel + m.group(3)
    return IMPORT.sub(sub, text)

# every JS module under experiments (except out/, node_modules/, prior-work/) gets its imports rewritten
files = []
for d, dirs, fs in os.walk(E):
    dirs[:] = [x for x in dirs if x not in ("out", "node_modules", "prior-work")]
    files += [os.path.join(d, f) for f in fs if f.endswith(".mjs")]
changed = 0
for f in files:
    text = open(f, encoding="utf8").read()
    new = rewrite(f, text)
    dest = new_location(f)
    if new != text or dest != os.path.normpath(f):
        changed += 1
        print(("move " if dest != os.path.normpath(f) else "edit ") + os.path.relpath(f, ROOT) + (" -> " + os.path.relpath(dest, ROOT) if dest != os.path.normpath(f) else ""))
    if dry: continue
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    open(dest, "w", encoding="utf8", newline="").write(new)
    if dest != os.path.normpath(f): os.remove(f)
for a, b in MOVES.items():  # non-JS files
    if a.endswith(".mjs"): continue
    src, dst = os.path.join(E, a), os.path.join(S, b)
    print(f"move {os.path.relpath(src, ROOT)} -> {os.path.relpath(dst, ROOT)}")
    if not dry and os.path.exists(src): os.makedirs(os.path.dirname(dst), exist_ok=True); shutil.move(src, dst)
print(f"{changed} JS files touched{' (dry run)' if dry else ''}")
