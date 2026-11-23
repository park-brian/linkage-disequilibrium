// Every repo path the living docs name must exist. Scans `backticked` tokens and [links](targets) in README.md,
// AGENTS.md and docs/*.md (not docs/history/). Placeholders ($VAR, <x>, {a,b}, *, ../), URLs and bare names are
// skipped. Run by `npm test` (test/docs.test.js); usage: node tools/check_docs.js → exits 1 and lists any missing.
import { readFileSync, readdirSync, existsSync } from "node:fs";
import { resolve, dirname } from "node:path";
import { REPO } from "./config.js";

const TOP = ["src", "app", "pipeline", "verify", "test", "bench", "tools", "docs"];
const ACTIVE_PLAN = "PLAN.md"; // a plan in progress names files it will create (checked once it moves to history)
const docs = ["README.md", "AGENTS.md", "src/README.md", ...readdirSync(resolve(REPO, "docs")).filter(f => f.endsWith(".md") && f !== ACTIVE_PLAN).map(f => `docs/${f}`)];

/** → [{ doc, path }] for every named repo path that does not exist */
export function missingPaths() {
  const missing = [];
  for (const doc of docs) {
    const text = readFileSync(resolve(REPO, doc), "utf8");
    const tokens = [...text.matchAll(/`([^`\s]+)`/g)].map(m => [m[1], REPO])
      .concat([...text.matchAll(/\]\(([^)\s#]+)(?:#[^)]*)?\)/g)].map(m => [m[1], resolve(REPO, dirname(doc))]));
    for (const [raw, base] of tokens) {
      const t = raw.replace(/[.,;:]+$/, "").replace(/:\d+$/, "");
      if (/[$<>{}*|]|^\.\.|^https?:|^\w+:/.test(t)) continue;
      const first = t.split("/")[0];
      const repoRooted = base === REPO ? TOP.includes(first) || /^[A-Z]+\.md$/.test(t) : true;
      if (!repoRooted || !t.includes("/") && base === REPO && !/\.md$/.test(t)) continue;
      // brace-free shorthand like `config.{mjs,sh}` is skipped above; strip a trailing "(…)" call
      const path = resolve(base, t.replace(/\(.*$/, ""));
      if (!existsSync(path)) missing.push({ doc, path: t });
    }
  }
  return missing;
}

if (import.meta.url === `file:///${process.argv[1]?.replace(/\\/g, "/").replace(/^\//, "")}`) {
  const m = missingPaths();
  for (const x of m) console.log(`${x.doc}: ${x.path}`);
  console.log(m.length ? `${m.length} missing` : "all doc paths exist");
  process.exit(m.length ? 1 : 0);
}
