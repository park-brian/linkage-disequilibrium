// The layer rule (docs/ARCHITECTURE.md, src/): src/core never imports src/parity, and the kernel chain (core/ld/*,
// core/haps.js), which workers load by path without an import map, imports nothing outside itself.
import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync, readdirSync } from "node:fs";
import { join, posix } from "node:path";
import { REPO } from "../tools/config.js";

const files = dir => readdirSync(join(REPO, dir), { recursive: true })
  .map(f => `${dir}/${String(f).replaceAll("\\", "/")}`).filter(f => f.endsWith(".js"));
const imports = f => [...readFileSync(join(REPO, f), "utf8").matchAll(/(?:from|import)\s*\(?\s*"([^"]+)"/g)]
  .map(m => (m[1].startsWith(".") ? posix.normalize(posix.join(posix.dirname(f), m[1])) : m[1]));

test("layers: src/core imports nothing from src/parity", () => {
  const bad = files("src/core").flatMap(f => imports(f).filter(i => i.startsWith("src/parity/")).map(i => `${f} → ${i}`));
  assert.deepEqual(bad, []);
});

test("layers: the kernel chain imports only itself", () => {
  const chain = [...files("src/core/ld"), "src/core/haps.js"];
  const bad = chain.flatMap(f => imports(f).filter(i => !chain.includes(i)).map(i => `${f} → ${i}`));
  assert.deepEqual(bad, []);
});
