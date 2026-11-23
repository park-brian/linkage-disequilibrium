# Move files with `git mv` and rewrite every relative module specifier that pointed at a moved file: `import … from
# "./x"`, `import("./x")`, `new URL("./x", import.meta.url)`, `new Worker(new URL(…))`, and src= / href= in .html.
# Each specifier is resolved against its file's OLD location, mapped through the move table, and re-expressed from the
# file's NEW location. usage (from the repo root): python -I tools/move.py <moves.tsv> [--dry-run]
#   moves.tsv: "<old path>\t<new path>" per line (repo-relative; "#" comments; a trailing "/" moves a directory)
import os, re, subprocess, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SKIP = {"node_modules", ".git", "data", "out"}
dry = "--dry-run" in sys.argv
table = [l.rstrip("\n").split("\t") for l in open(sys.argv[1], encoding="utf8") if l.strip() and not l.startswith("#")]

def files_under(rel):
    full = os.path.join(ROOT, rel)
    if os.path.isfile(full): return [rel]
    out = []
    for d, dirs, fs in os.walk(full):
        dirs[:] = [x for x in dirs if x not in SKIP]
        out += [os.path.relpath(os.path.join(d, f), ROOT).replace("\\", "/") for f in fs]
    return out

moves = {}  # old file → new file (repo-relative, forward slashes)
for old, new in table:
    old, new = old.strip(), new.strip()
    if old.endswith("/"):
        for f in files_under(old.rstrip("/")): moves[f] = new.rstrip("/") + "/" + f[len(old.rstrip("/")) + 1:]
    else: moves[old] = new
loc = lambda rel: moves.get(rel, rel)

SPEC = re.compile(r"""((?:\bfrom|\bimport)\s*\(?\s*|new URL\(\s*|\b(?:src|href)=)(["'])(\.{1,2}/[^"']*)(["'])""")
def rewrite(rel_old, text):
    here_old, here_new = os.path.dirname(rel_old), os.path.dirname(loc(rel_old))
    def sub(m):
        target_old = os.path.normpath(os.path.join(here_old, m.group(3))).replace("\\", "/")
        if target_old.startswith(".."): return m.group(0)  # outside the repo: leave
        trailing = "/" if m.group(3).endswith("/") and not target_old.endswith("/") else ""
        target_new = loc(target_old)
        if target_new == target_old and here_new == here_old: return m.group(0)
        r = os.path.relpath(target_new, here_new or ".").replace("\\", "/") + trailing
        if not r.startswith("."): r = "./" + r
        return m.group(1) + m.group(2) + r + m.group(4)
    return SPEC.sub(sub, text)

code = []
for d, dirs, fs in os.walk(ROOT):
    dirs[:] = [x for x in dirs if x not in SKIP]
    code += [os.path.relpath(os.path.join(d, f), ROOT).replace("\\", "/") for f in fs if f.endswith((".mjs", ".js", ".html"))]
edited = 0
for rel in code:
    path = os.path.join(ROOT, rel); text = open(path, encoding="utf8").read(); new = rewrite(rel, text)
    if new != text:
        edited += 1
        if not dry: open(path, "w", encoding="utf8", newline="").write(new)
for old, new in sorted(moves.items()):
    print(f"{old} -> {new}")
    if dry: continue
    os.makedirs(os.path.dirname(os.path.join(ROOT, new)), exist_ok=True)
    subprocess.run(["git", "mv", "-k", old, new], cwd=ROOT, check=True)
print(f"{len(moves)} moves, {edited} files with rewritten specifiers{' (dry run)' if dry else ''}")
