// The correctness lint (eslint.config.js: undefined names, unused values, unreachable code and the like; no style
// rules) finds nothing in the repo.
import { test } from "node:test";
import assert from "node:assert/strict";
import { relative } from "node:path";
import { ESLint } from "eslint";
import { REPO } from "../tools/config.js";

test("lint: eslint (correctness rules only) finds nothing", async () => {
  const results = await new ESLint({ cwd: REPO }).lintFiles(["."]);
  const at = r => relative(REPO, r.filePath).replaceAll("\\", "/");
  const found = results.flatMap(r => r.messages.map(m => `${at(r)}:${m.line} ${m.ruleId}: ${m.message}`));
  assert.deepEqual(found, []);
});
