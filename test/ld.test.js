// The core LD API (src/core/ld.js) against the LDlink ports on the release: the same statistics where both pick the
// same rows (the example SNPs: the allele-aware row is LDlink's row), and its typed errors.
import { test } from "node:test";
import assert from "node:assert/strict";
import { existsSync } from "node:fs";
import { RELEASE } from "../tools/config.js";

const needData = { skip: process.env.LD_SKIP_DATA === "1" && "LD_SKIP_DATA=1" };
const EX = ["rs13426489", "rs4442975", "rs13387042", "rs13412666", "rs6721996", "rs12465506", "rs12613955", "rs12620095"];
async function setup() {
  assert.ok(existsSync(`${RELEASE}/manifest.json`), `no release at ${RELEASE}`);
  const { asyncBufferFromFile, cachedAsyncBuffer } = await import("hyparquet");
  const { jsKernel } = await import("../src/core/ld/kernel.js");
  const { createLD } = await import("../src/core/ld.js");
  const { createContext } = await import("../src/parity/ldlink/context.js");
  const open = async p => cachedAsyncBuffer(await asyncBufferFromFile(p));
  const ctx = await createContext({ root: RELEASE, open, kernel: jsKernel() });
  return { ld: createLD({ root: RELEASE, open, kernel: jsKernel() }), ctx };
}

test("core LD API: pair, matrix and proxy agree with LDpair / LDmatrix on the example SNPs", needData, async () => {
  const { ld, ctx } = await setup();
  const { ldpair } = await import("../src/parity/ldlink/ldpair.js");
  const { ldmatrix } = await import("../src/parity/ldlink/ldmatrix.js");
  const p = await ld.pair("rs3", "rs4", { pop: "CEU" });
  const port = (await ldpair(ctx, { var1: "rs3", var2: "rs4", pop: "CEU", build: "grch37" }))[0];
  assert.equal(p.r2.toFixed(4), (+port.statistics.r2).toFixed(4));
  assert.equal(p.dprime.toFixed(4), (+port.statistics.d_prime).toFixed(4));
  assert.equal(Object.values(p.haplotypes).reduce((a, b) => a + b, 0), p.n);
  assert.deepEqual([p.a.pos, p.b.pos], [32446842, 32447222]);
  const m = await ld.matrix(EX, { pop: "CEU" }), pm = await ldmatrix(ctx, { inputs: EX, pop: "CEU", build: "grch37" });
  const ix = pm.ids.map(id => m.variants.findIndex(v => v.input === id)); // the port sorts by position
  for (let i = 0; i < pm.ids.length; i++) for (let j = 0; j < pm.ids.length; j++) {
    if (pm.r2[i][j] !== "NA") assert.ok(Math.abs(+pm.r2[i][j] - m.r2[ix[i]][ix[j]]) <= 5e-4, `${pm.ids[i]} × ${pm.ids[j]}`);
  }
  const { query, proxies } = await ld.proxy("rs3", { pop: "CEU", window: 100_000 });
  assert.equal(query.rsid, "rs3");
  assert.ok(proxies.every((x, k) => x.r2 >= 0.1 && (k === 0 || x.r2 <= proxies[k - 1].r2)), "r² ≥ minR2, descending");
  assert.equal(proxies.find(x => x.variant.rsid === "rs4")?.r2.toFixed(6), p.r2.toFixed(6));
});

test("core LD API: typed errors", needData, async () => {
  const { ld } = await setup();
  const code = async f => { try { await f(); return "none"; } catch (e) { assert.equal(e.name, "LDError"); return e.code; } };
  assert.equal(await code(() => ld.pair("xyz", "rs4", { pop: "CEU" })), "bad_input");
  assert.equal(await code(() => ld.pair("rs3", "rs4", { pop: "NOPE" })), "unknown_population");
  assert.equal(await code(() => ld.pair("rs3", "rs13426489", { pop: "CEU" })), "chromosomes_differ");
  assert.equal(await code(() => ld.pair("rs99999999999", "rs4", { pop: "CEU" })), "not_found");
  assert.equal(await code(() => ld.variant("rs3", "hg19")), "bad_input");
  assert.equal((await ld.variant("chr13:32446842")).rsid, "rs3"); // a coordinate: the row that starts there
});
