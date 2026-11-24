// Every repo path the living docs name must exist. Scans `backticked` tokens and [links](targets) in README.md,
// AGENTS.md and docs/*.md (not docs/history/); a brace list (`config.{js,sh}`) must exist in every expansion.
// Placeholders ($VAR, <x>, *, ../), URLs and bare names are skipped. Also every `usage:` line in a code file must run a
// script by its repo path (badUsage). Run by `npm test` (test/docs.test.js).
// usage: node tools/check_docs.js → exits 1 and lists any missing path or bad usage line.
import { readFileSync, readdirSync, existsSync } from "node:fs";
import { resolve, dirname, relative, extname } from "node:path";
import { REPO } from "./config.js";

const TOP = ["src", "app", "pipeline", "verify", "test", "bench", "tools", "docs"];
const ACTIVE = ["PLAN.md", "REVIEW.md"];
  // a plan in progress (and its review) names files it will create or delete; checked once in history
const docs = ["README.md", "AGENTS.md", "src/README.md", ...readdirSync(resolve(REPO,
  "docs")).filter(f => f.endsWith(".md") && !ACTIVE.includes(f)).map(f => `docs/${f}`)];

/** → [{ doc, path }] for every named repo path that does not exist */
export function missingPaths() {
  const missing = [];
  for (const doc of docs) {
    const text = readFileSync(resolve(REPO, doc), "utf8");
    const tokens = [...text.matchAll(/`([^`\s]+)`/g)].map(m => [m[1], REPO])
      .concat([...text.matchAll(/\]\(([^)\s#]+)(?:#[^)]*)?\)/g)].map(m => [m[1], resolve(REPO, dirname(doc))]));
    for (const [raw, base] of tokens) {
      const t = raw.replace(/[.,;:]+$/, "").replace(/:\d+$/, "");
      if (/[$<>*|]|^\.\.|^https?:|^\w+:/.test(t)) continue;
      const first = t.split("/")[0];
      const repoRooted = base === REPO ? TOP.includes(first) || /^[A-Z]+\.md$/.test(t) : true;
      if (!repoRooted || !t.includes("/") && base === REPO && !/\.md$/.test(t)) continue;
      // strip a trailing "(…)" call; a {a,b} list names one path per item
      const brace = /\{([^{}]+)\}/.exec(t), bare = t.replace(/\(.*$/, "");
      const each = brace ? brace[1].split(",").map(x => bare.replace(brace[0], x)) : [bare];
      for (const e of each) if (!existsSync(resolve(base, e))) missing.push({ doc, path: e });
    }
  }
  return missing;
}

const SCRIPT = new Set([".js", ".py", ".sh", ".sql"]);
/** → [{ file, line, script }] for every `usage:` line naming a script that is not a repo path (or not there) */
export function badUsage() {
  const bad = [];
  for (const dir of TOP.filter(d => d !== "docs")) for (const f of readdirSync(resolve(REPO, dir), { recursive: true })) {
    const rel = `${dir}/${String(f).replaceAll("\\", "/")}`;
    if (!SCRIPT.has(extname(rel)) || rel.includes("node_modules/")) continue;
    readFileSync(resolve(REPO, rel), "utf8").split(/\r?\n/).forEach((l, i) => {
      const m = /usage:\s*(.*)$/.exec(l);
      if (!m) return;
      for (const tok of m[1].split(/\s+/)) {
        const t = tok.replace(/^["'(]+|["'),;:]+$/g, "");
        if (!SCRIPT.has(extname(t)) || /[$<>*{}]/.test(t)) continue; // a script name, not a placeholder
        if (!t.includes("/") || !existsSync(resolve(REPO, t))) bad.push({ file: rel, line: i + 1, script: t });
      }
    });
  }
  return bad;
}

if (import.meta.url === `file:///${process.argv[1]?.replace(/\\/g, "/").replace(/^\//, "")}`) {
  const m = missingPaths(), u = badUsage();
  for (const x of m) console.log(`${x.doc}: ${x.path}`);
  for (const x of u) console.log(`${x.file}:${x.line}: usage names ${x.script}`);
  console.log(m.length || u.length ? `${m.length} missing, ${u.length} bad usage` : "all doc paths exist; usage lines name repo scripts");
  process.exit(m.length || u.length ? 1 : 0);
}
