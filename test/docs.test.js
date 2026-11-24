// The living docs only name repo paths that exist (tools/check_docs.js).
import { test } from "node:test";
import assert from "node:assert/strict";
import { missingPaths, badUsage } from "../tools/check_docs.js";

test("docs: every repo path named in README, AGENTS, src/README and docs/*.md exists", () => {
  const missing = missingPaths();
  assert.deepEqual(missing, [], missing.map(m => `${m.doc}: ${m.path}`).join("\n"));
});

test("usage lines: every script a usage: line runs is named by its repo path", () => {
  const bad = badUsage();
  assert.deepEqual(bad, [], bad.map(b => `${b.file}:${b.line}: ${b.script}`).join("\n"));
});
