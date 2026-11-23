// The living docs only name repo paths that exist (tools/check_docs.js).
import { test } from "node:test";
import assert from "node:assert/strict";
import { missingPaths } from "../tools/check_docs.js";

test("docs: every repo path named in README, AGENTS, src/README and docs/*.md exists", () => {
  const missing = missingPaths();
  assert.deepEqual(missing, [], missing.map(m => `${m.doc}: ${m.path}`).join("\n"));
});
