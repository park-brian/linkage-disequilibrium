// Unit and regression tests (the Phase A fixes, docs/history/ROADMAP.md §2.1, and later kernel work). npm test
// Release-backed tests read tools/config.mjs RELEASE and FAIL when it is missing; LD_SKIP_DATA=1 skips them.
// Each test reproduces the failing input from the code review; it failed before its fix.
import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync, existsSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { loadKernel, jsKernel, KERNEL_WASM } from "../src/ld/kernel.mjs";
import { makeLabeler, minimalForm, variantClass, TYPES } from "../pipeline/genotypes/rsid_label.mjs";
import { mungeSumstats } from "../src/ldsc/munge.mjs";
import { getSeparators } from "../src/ldsc/regress.mjs";
import { pyFixed3 } from "../src/ldsc/l2.mjs";
import { parseCoord } from "../src/genome.mjs";
import { isNotFound } from "../src/source.mjs";
import { RELEASE } from "../tools/config.mjs";

const here = new URL(".", import.meta.url);
// release-backed tests: skipped only on request (LD_SKIP_DATA=1); a configured but missing release is a failure
const needData = { skip: process.env.LD_SKIP_DATA === "1" && "LD_SKIP_DATA=1" };
const haveRelease = () => assert.ok(existsSync(`${RELEASE}/manifest.json`), `no release at ${RELEASE} (set LD_RELEASE / LD_DATA, or LD_SKIP_DATA=1)`);
const fixture = f => new URL(`fixtures/${f}`, here);

// ---------- R1: WASM kernel, identity-copy path keeps bits beyond nh ----------
function rng(seed) { let s = seed >>> 0; return () => ((s = (s * 1664525 + 1013904223) >>> 0) / 2 ** 32); }
const popPrefix = (bits, nh) => { let c = 0; for (let k = 0; k < nh; k++) c += (bits[k >> 3] >> (k & 7)) & 1; return c; };
test("R1: prepare() counts only the first nh haplotypes on the identity path (WASM and JS backends)", async () => {
  for (const k of [await loadKernel(readFileSync(KERNEL_WASM)), jsKernel()]) { // both backends: same results
  const samples = 100, hbytes = Math.ceil((2 * samples) / 8), r = rng(42);
  const rows = Array.from({ length: 30 }, (_, i) => ({ pos: 1000 + i, haps: Uint8Array.from({ length: hbytes }, () => (r() * 256) | 0) }));
  const allAlt = [{ pos: 1, haps: new Uint8Array(hbytes).fill(255) }];
  for (const nh of [40, 63, 64, 96, 128, 150, 200]) {
    const idx = Array.from({ length: nh }, (_, j) => j);
    assert.deepEqual([...k.prepare(allAlt, { samples, hapIndexes: idx }).views().pa], [nh], `all-alt row, nh=${nh}`);
    const p = k.prepare(rows, { samples, hapIndexes: idx });
    assert.deepEqual([...p.views().pa], rows.map(x => popPrefix(x.haps, nh)), `allele counts, nh=${nh}`);
    // LD is invariant to haplotype order: the identity path must equal the general (permuted) path
    const a = k.band(p, { maxDist: 1e9, maxOut: 1 << 12 });
    const q = k.prepare(rows, { samples, hapIndexes: [...idx].reverse() });
    const b = k.band(q, { maxDist: 1e9, maxOut: 1 << 12 });
    assert.equal(a.count, b.count, `pair count, nh=${nh}`);
    assert.deepEqual([...a.r2], [...b.r2], `r2, nh=${nh}`);
    assert.deepEqual([...a.dprime], [...b.dprime], `dprime, nh=${nh}`);
  }
  }
});

// ---------- B1/B2/B8: rsID labels for indels and split multi-allelic rows ----------
function db(docs) { // docs: [pos, rsid, ref_rsid, type, function], any order → sorted typed arrays
  const d = [...docs].sort((x, y) => x[0] - y[0]);
  return { dbPos: Int32Array.from(d, x => x[0]), dbRsid: Float64Array.from(d, x => x[1]), dbRef: Float64Array.from(d, x => x[2]),
    dbType: Uint8Array.from(d, x => { const t = TYPES.indexOf(x[3]); return t < 0 ? 255 : t; }), dbFunc: d.map(x => x[4] ?? null) };
}
test("B1: a deletion's rsID is found at POS+1 (LDlink stores deletions there) with delta 1", () => {
  const label = makeLabeler(db([[100, 11, 11, "snv"], [101, 10, 10, "del", "frameshift"]]));
  const l = label(100, null, "AT", "A");
  assert.equal(l.rsid, 10); assert.equal(l.source, 3); assert.equal(l.delta, 1); assert.equal(l.func, "frameshift");
  const v = label(100, "rs10", "AT", "A"); // the VCF's own rsID is found at POS+1 too
  assert.equal(v.rsid, 10); assert.equal(v.source, 1); assert.equal(v.delta, 1);
});
test("B1: insertions prefer POS (VC=INS) and fall back to POS+1 (delins inside a repeat)", () => {
  assert.deepEqual([makeLabeler(db([[100, 20, 20, "ins"]]))(100, null, "A", "AT")].map(l => [l.rsid, l.delta]), [[20, 0]]);
  assert.deepEqual([makeLabeler(db([[101, 21, 21, "delins"]]))(100, null, "A", "AT")].map(l => [l.rsid, l.delta]), [[21, 1]]);
});
test("B1: an NCBI allele match that the dump has wins over position + type (source 4)", () => {
  const label = makeLabeler(db([[101, 30, 30, "del"], [101, 31, 31, "del"]]), { match: () => [31] });
  const l = label(100, null, "ATT", "A");
  assert.deepEqual([l.rsid, l.source, l.delta, l.alleles], [31, 4, 1, [31]]);
});
test("B2: a split multi-allelic row REF=AT ALT=GT is labelled as the SNV it is", () => {
  assert.deepEqual(minimalForm(100, "AT", "GT"), [100, "A", "G"]);
  assert.equal(variantClass("A", "G"), 0);
  const l = makeLabeler(db([[100, 40, 40, "snv"], [101, 41, 41, "del"]]))(100, null, "AT", "GT");
  assert.deepEqual([l.rsid, l.source, l.delta], [40, 3, 0]);
  const shifted = makeLabeler(db([[101, 42, 42, "snv"]]))(100, null, "CAT", "CGT"); // SNV after prefix trim is at 101
  assert.deepEqual([shifted.rsid, shifted.delta], [42, 1]);
});
test("B8: a dbSNP doc with an unknown type is compatible with nothing", () => {
  const l = makeLabeler(db([[101, 50, 50, "weird"]]))(100, null, "AT", "A");
  assert.equal(l.rsid, null); assert.equal(l.source, 0);
});

// ---------- L1/L2/L4: munge vs reference munge_sumstats.py (fixtures/ref_*.sumstats.gz) ----------
const snplist = (() => {
  const [, ...rows] = readFileSync(fixture("mini.snplist"), "utf8").trim().split("\n").map(l => l.split("\t"));
  return { SNP: rows.map(r => r[0]), A1: rows.map(r => r[1]), A2: rows.map(r => r[2]) };
})();
for (const [name, what] of [["cc", "L1: N derived from N_CAS/N_CON is written"], ["cc_n", "L2: N_CAS/N_CON override an N column"], ["neff", "L4: %.3f of N rounds the exact value"]]) {
  test(`${what} (fixtures/${name}.txt == ref_${name}.sumstats.gz)`, () => {
    const got = mungeSumstats(readFileSync(fixture(`${name}.txt`), "utf8"), snplist).text;
    const want = gunzipSync(readFileSync(fixture(`ref_${name}.sumstats.gz`))).toString("utf8").replace(/\r\n/g, "\n");
    assert.equal(got, want);
  });
}

// ---------- L3: jackknife separators = np.floor(np.linspace(0, N, nb + 1)) ----------
test("L3: block separators match numpy (fixtures from numpy 2.4.3)", () => {
  const f = JSON.parse(readFileSync(new URL("numpy_separators.json", here), "utf8"));
  for (const [k, v] of Object.entries(f)) { const [N, nb] = k.split(",").map(Number); assert.deepEqual(getSeparators(N, nb), v, k); }
});

// ---------- L4: pyFixed3 = Python '%.3f' ----------
test("L4: pyFixed3 matches Python '%.3f' on ties and near-ties", () => {
  const cases = { 0.0005: "0.001", 0.0025: "0.003", 31.6765: "31.677", 110.8665: "110.867", 0.125: "0.125", 0.0625: "0.062",
    2.5: "2.500", "-0.0005": "-0.001", 1e-4: "0.000", 47104.6665: "47104.666", 0: "0.000" };
  for (const [x, want] of Object.entries(cases)) assert.equal(pyFixed3(Number(x)), want, x);
  assert.equal(pyFixed3(-0), "-0.000"); assert.equal(pyFixed3(-0.0004), "-0.000");
});

// ---------- P1: malformed coordinates are "not found", never an exception ----------
test("P1: parseCoord follows LDlink's int() parsing", () => {
  assert.equal(parseCoord("chr22:abc"), null);
  assert.equal(parseCoord("chr22:1.5"), null);
  assert.equal(parseCoord("chr22:12_345").pos, 12345n);
  assert.equal(parseCoord("X:5").chro, "X");
  assert.equal(parseCoord("chr22:0").key, null);
});
test("P1: planner.query survives malformed coordinate inputs", needData, async () => {
  haveRelease();
  const { openPlanner } = await import("../src/planner.mjs");
  const { asyncBufferFromFile, cachedAsyncBuffer } = await import("hyparquet");
  const root = RELEASE;
  const pl = await openPlanner(root, async p => cachedAsyncBuffer(await asyncBufferFromFile(p)));
  const r = await pl.query({ build: "grch37", inputs: ["chr22:abc", "chr22:1.5", "chr22:16051453", "chr22:0"] });
  assert.deepEqual(r.variants.map(v => [v.found, v.pos]), [[false, null], [false, null], [true, 16051453], [false, null]]);
});

// ---------- V1 finding: binary BYTE_ARRAY columns must not be decoded as UTF-8 ----------
const rare = `${RELEASE}/genotypes/grch37/chr22.rare.parquet`;
test("readers return binary haps/missing as bytes, not UTF-8 strings (sparse tier)", needData, async () => {
  haveRelease();
  const { asyncBufferFromFile, parquetMetadataAsync } = await import("hyparquet");
  const { compressors } = await import("hyparquet-compressors");
  const { sortedLookup, rangeScan } = await import("../src/parquet.mjs");
  const p = rare, file = await asyncBufferFromFile(p), metadata = await parquetMetadataAsync(file);
  const scan = await rangeScan({ file, metadata, column: "pos", lo: 16050000, hi: 16060000, columns: ["pos", "haps", "missing"], compressors });
  assert.ok(scan.length > 0);
  for (const r of scan) assert.ok(r.haps instanceof Uint8Array, "rangeScan haps");
  const look = await sortedLookup({ file, metadata, column: "pos", keys: [scan[0].pos], columns: ["pos", "haps"], compressors, payloadColumn: "haps" });
  assert.ok(look.rows[0].haps instanceof Uint8Array, "sortedLookup haps");
  assert.deepEqual([...look.rows[0].haps], [...scan.find(r => r.pos === scan[0].pos).haps]);
});

// ---------- K1: kernel doubles == LDlink's LD_calcs, thresholds on doubles, Python round(x, 3) ----------
test("K1: exact() is bit-identical to Python; band thresholds compare doubles (>= / strict >); pyRound = round()", async () => {
  const { pyRound } = await import("../src/pyfmt.mjs");
  for (const k of [await loadKernel(readFileSync(KERNEL_WASM)), jsKernel()]) { // both backends: same results
  const cases = JSON.parse(readFileSync(fixture("ld_counts.json"), "utf8"));
  const hex = h => { const m = /^(-?)0x([01])\.([0-9a-f]+)p([+-]\d+)$/.exec(h); if (!m) return parseFloat(h); return (m[1] ? -1 : 1) * (parseInt(m[2], 16) + parseInt(m[3], 16) / 16 ** m[3].length) * 2 ** +m[4]; };
  let n = 0;
  for (const c of cases) {
    const samples = Math.ceil(c.n / 2), hb = Math.ceil((2 * samples) / 8), a = new Uint8Array(hb), b = new Uint8Array(hb);
    const set = (arr, h) => { arr[h >> 3] |= 1 << (h & 7); };
    for (let h = 0; h < c.x; h++) { set(a, h); set(b, h); }
    for (let h = c.x; h < c.pa; h++) set(a, h);
    for (let h = c.pa; h < c.pa + c.pb - c.x; h++) set(b, h);
    const p = k.prepare([{ pos: 1, haps: a }, { pos: 2, haps: b }], { samples, hapIndexes: Array.from({ length: c.n }, (_, i) => i) });
    const e = k.exact(p, [0], [1]), r2 = hex(c.r2), dp = hex(c.dp);
    assert.equal(e.r2[0], r2, `r2 ${JSON.stringify(c)}`); assert.equal(e.dprime[0], dp, `dprime ${JSON.stringify(c)}`);
    assert.equal(k.band(p, { maxDist: 10, threshold: r2, maxOut: 4 }).count, 1, "r2 >= r2 keeps the pair");
    assert.equal(k.band(p, { maxDist: 10, threshold: r2, strict: true, maxOut: 4 }).count, 0, "r2 > r2 drops it");
    assert.equal(k.band(p, { maxDist: 10, threshold: dp, stat: "dprime", maxOut: 4 }).count, 1, "D' >= D' keeps the pair");
    assert.equal(String(pyRound(r2, 3)), c.r2_3.replace(/\.0$/, "").replace(/^-0$/, "0"), `round(r2,3) ${c.r2}`);
    assert.equal(String(pyRound(dp, 3)), c.dp_3.replace(/\.0$/, "").replace(/^-0$/, "0"), `round(dp,3) ${c.dp}`);
    n++;
  }
  assert.equal(n, cases.length);
  }
});

// ---------- P5: only a missing file is "not found"; other I/O errors surface ----------
test("P5: resolver reports missing files as not found but rethrows other I/O errors", needData, async () => {
  haveRelease();
  const { createResolver } = await import("../src/dbsnp.mjs");
  const { asyncBufferFromFile, cachedAsyncBuffer } = await import("hyparquet");
  const root = RELEASE;
  const ok = async p => cachedAsyncBuffer(await asyncBufferFromFile(p));
  // missing partition (rsID beyond the last partition): ENOENT → not found
  assert.deepEqual((await createResolver(root, ok).resolveRsids(["rs99999999999999"], "grch37"))[0].found, false);
  // a network-style failure on an existing partition must not look like "not in dbSNP"
  const flaky = async p => { if (/part-00/.test(p)) throw Object.assign(new Error("socket hang up"), { code: "ECONNRESET" }); return ok(p); };
  await assert.rejects(createResolver(root, flaky).resolveRsids(["rs123"], "grch37"), /socket hang up/);
  const http403 = async p => { if (/chr22/.test(p)) throw Object.assign(new Error("HTTP 403 Forbidden"), { status: 403 }); return ok(p); };
  await assert.rejects(createResolver(root, http403).resolveCoords(["chr22:16051453"], "grch37"), /403/);
  // with exists(): an absent file is not found without any request
  let opened = 0; const counting = async p => { opened++; return ok(p); };
  assert.equal((await createResolver(root, counting, { exists: () => false }).resolveRsids(["rs123"], "grch37"))[0].found, false);
  assert.equal(opened, 0);
  assert.ok(isNotFound({ code: "ENOENT" }) && isNotFound({ status: 404 }) && isNotFound(new Error("HTTP 404 Not Found")));
  assert.ok(!isNotFound(new Error("HTTP 4040 bytes")) && !isNotFound({ status: 403 }));
});

// ---------- K2/K3: band() output sizing, truncation, row ranges ----------
test("K2/K3: band sizes output exactly, flags truncation, and handles row ranges and block edge cases", async () => {
  for (const k of [await loadKernel(readFileSync(KERNEL_WASM)), jsKernel()]) { // both backends: same results
  const samples = 50, hb = Math.ceil((2 * samples) / 8), r = rng(9);
  const mk = n => Array.from({ length: n }, (_, i) => ({ pos: 1 + i, haps: Uint8Array.from({ length: hb }, () => (r() * 256) | 0) }));
  const idx = Array.from({ length: 2 * samples }, (_, j) => j);
  // > 4M passing pairs: the default buffer (4M) overflows and band() re-runs with exactly `written`
  const big = k.prepare(mk(3000), { samples, hapIndexes: idx });
  const all = k.band(big, { maxDist: 1e9 });
  assert.equal(all.truncated, false); assert.equal(all.count, all.written); assert.ok(all.count > 1 << 22, `count ${all.count}`);
  const capped = k.band(big, { maxDist: 1e9, maxOut: 1000 });
  assert.equal(capped.count, 1000); assert.equal(capped.truncated, true); assert.equal(capped.written, all.count);
  // row ranges: an empty range is empty (it used to mean "all rows"); a split covers every pair exactly once
  const p = k.prepare(mk(200), { samples, hapIndexes: idx });
  const full = k.band(p, { maxDist: 50, block: 4 });
  assert.equal(k.band(p, { maxDist: 50, block: 4, rows: [0, 0] }).count, 0);
  assert.equal(k.band(p, { maxDist: 50, block: 4, rows: [120, 80] }).count, 0);
  const seen = new Map();
  for (const rows of [[0, 77], [77, 150], [150, 200]]) { const b = k.band(p, { maxDist: 50, block: 4, rows }); for (let t = 0; t < b.count; t++) { const key = `${b.i[t]},${b.j[t]}`; assert.ok(!seen.has(key), `duplicate ${key}`); seen.set(key, b.r2[t]); } }
  assert.equal(seen.size, full.count);
  for (let t = 0; t < full.count; t++) assert.equal(seen.get(`${full.i[t]},${full.j[t]}`), full.r2[t]);
  // block <= 1 with a row range no longer loops forever (B is clamped to 1)
  assert.equal(k.band(p, { maxDist: 50, block: -3, rows: [0, 200] }).count, full.count);
  assert.equal(k.band(p, { maxDist: 50, block: 0, rows: [0, 200] }).count, full.count);
  }
});

// ---------- P4/P6: planner regions and lookup columns ----------
test("P6: sortedLookup matches keys even when columns omit the key column", needData, async () => {
  haveRelease();
  const { asyncBufferFromFile, parquetMetadataAsync } = await import("hyparquet");
  const { compressors } = await import("hyparquet-compressors");
  const { sortedLookup, rangeScan } = await import("../src/parquet.mjs");
  const p = rare, file = await asyncBufferFromFile(p), metadata = await parquetMetadataAsync(file);
  const some = (await rangeScan({ file, metadata, column: "pos", lo: 16050000, hi: 16060000, columns: ["pos"], compressors })).slice(0, 5).map(r => r.pos);
  const got = await sortedLookup({ file, metadata, column: "pos", keys: some, columns: ["ref"], compressors });
  assert.deepEqual(got.rows.map(r => r.pos).sort((a, b) => a - b), [...new Set(some)].sort((a, b) => a - b));
});
test("P4: overlapping regions return each row once; chromosome names are normalized", needData, async () => {
  haveRelease();
  const { openPlanner } = await import("../src/planner.mjs");
  const { asyncBufferFromFile, cachedAsyncBuffer } = await import("hyparquet");
  const root = RELEASE;
  const pl = await openPlanner(root, async p => cachedAsyncBuffer(await asyncBufferFromFile(p)));
  const { createHash } = await import("node:crypto"); // planner rows carry no ref/alt: split multi-allelic rows share pos
  const key = r => `${r.tier}|${r.pos}|${createHash("sha1").update(r.haps).digest("hex")}`;
  const one = await pl.query({ build: "grch37", regions: [{ chrom: "22", start: 30_000_000, end: 30_300_000 }], tiers: ["common"] });
  const two = await pl.query({ build: "grch37", regions: [{ chrom: "chr22", start: 30_000_000, end: 30_200_000 }, { chrom: 22, start: 30_100_000, end: 30_300_000 }, { chrom: "chr22", start: 30_150_000, end: 30_160_000 }], tiers: ["common"] });
  assert.ok(one.rows.length > 0);
  assert.equal(new Set(two.rows.map(key)).size, two.rows.length, "no duplicate rows");
  assert.deepEqual(two.rows.map(key), one.rows.map(key));
  const x = await pl.query({ build: "grch37", regions: [{ chrom: "x", start: 100_000_000, end: 100_050_000 }, { chrom: 23, start: 100_000_000, end: 100_050_000 }], tiers: ["common"] });
  assert.ok(x.rows.length > 0); assert.equal(new Set(x.rows.map(key)).size, x.rows.length);
  const bad = await pl.query({ build: "grch37", regions: [{ chrom: "chr99", start: 1, end: 2 }, { chrom: "22", start: 5, end: 1 }] });
  assert.equal(bad.plan.invalidRegions.length, 2);
});

// ---------- L5–L8: LDscore edge cases vs the reference ldsc ----------
test("L5: --ld-wind of 0 is rejected like ldsc (np.array(..., dtype=bool)); the set window wins", async () => {
  const { ldscoreL2 } = await import("../src/ldsc/l2.mjs");
  // 3 SNPs 1 Mb apart, 4 individuals (1 byte per SNP): a 500 kb / 0.5 cM window is not whole-chromosome
  const bim = [0, 1, 2].map(i => ({ CHR: 22, SNP: `rs${i + 1}`, CM: i, BP: 1 + i * 1_000_000, A1: "A", A2: "G" })), bed = new Uint8Array([0x6c, 0x1b, 0x01, 0b11100100, 0b11011000, 0b01101100]);
  for (const o of [{ windCm: 0 }, { windKb: 0 }, {}, { windCm: 1, windKb: 500 }]) assert.throws(() => ldscoreL2(bed, bim, 4, o), /exactly one --ld-wind/, JSON.stringify(o));
  assert.doesNotThrow(() => ldscoreL2(bed, bim, 4, { windCm: 0, windKb: 500 })); // cm 0 is "unset", kb is used
  assert.doesNotThrow(() => ldscoreL2(bed, bim, 4, { windCm: 0.5 }));
});
test("L6: more jackknife blocks than data points is an error (jackknife.py:104)", async () => {
  const { hsq } = await import("../src/ldsc/regress.mjs");
  const n = 151, r = rng(5), ld = Float64Array.from({ length: n }, () => 1 + 50 * r()), N = new Float64Array(n).fill(1000);
  const chisq = Float64Array.from({ length: n }, (_, i) => (i < 5 ? 40 : 0.5 + r()));
  assert.throws(() => hsq(chisq, ld, ld, N, 1e6, { nBlocks: n, twostep: 30 }), /More blocks than data points/); // step 1 has 146 SNPs
  assert.doesNotThrow(() => hsq(chisq, ld, ld, N, 1e6, { nBlocks: 100, twostep: 30 }));
});
test("L7: erfc(±inf), and chi2.isf(P, 1) matches scipy down to the smallest normal P", async () => {
  const { erfc } = await import("../src/ldsc/stats.mjs");
  const { chi2Isf1 } = await import("../src/ldsc/munge.mjs");
  assert.equal(erfc(Infinity), 0); assert.equal(erfc(-Infinity), 2); assert.ok(Number.isNaN(erfc(NaN)));
  const cases = JSON.parse(readFileSync(fixture("chi2isf.json"), "utf8"));
  let n = 0;
  for (const [ps, want, z3] of cases) {
    const p = parseFloat(ps); if (p < 2.2250738585072014e-308) continue; // subnormal P: scipy/cephes itself loses precision
    const got = chi2Isf1(p), w = parseFloat(want);
    assert.equal(pyFixed3(Math.sqrt(got)), z3, `Z for P=${ps}`);
    if (w > 1e-10) assert.ok(Math.abs(got - w) / w < 1e-9, `chi2 for P=${ps}: ${got} vs ${w}`);
    n++;
  }
  assert.ok(n > 4500, `only ${n} normal-range cases`); // ~4.6% of the random P are subnormal and skipped
});
test("L8: munge edge cases (BOM, hex, inf, extra fields, duplicate columns, info, whitespace) == reference munge_sumstats.py", () => {
  const d = "l8/", sl = { SNP: ["rs1", "rs2"], A1: ["A", "C"], A2: ["G", "T"] };
  for (const f of ["bom", "hex", "inf_n", "extra", "dup_first", "dup_later", "two_p_last", "info_lower", "inf_beta", "tabs_spaces"]) {
    let got; try { got = mungeSumstats(readFileSync(fixture(`${d}${f}.txt`), "utf8"), sl).text; } catch (e) { got = `ERR ${e.message}`; }
    const gz = fixture(`${d}ref_${f}.sumstats.gz`);
    const want = existsSync(gz) ? gunzipSync(readFileSync(gz)).toString().replace(/\r\n/g, "\n")
      : `ERR ${readFileSync(fixture(`${d}ref_${f}.log`), "utf8").split(/\r?\n/).filter(l => /^ValueError: /.test(l)).at(-1).slice(12)}`;
    assert.equal(got, want, f);
  }
});

// ---------- B7 follow-up: GTs with multi-digit allele indexes (1000G chr1:207237233 has 12 ALTs) ----------
test("gtCodes parses single- and multi-digit allele indexes, missing calls and haploid GTs identically", async () => {
  const { gtCodes } = await import("../pipeline/genotypes/vcf_gt.mjs");
  const rec = gts => ["1", "100", ".", "A", "C", ".", ".", ".", "GT", ...gts];
  // fast path (<= 9 ALTs) and the general path agree on single-digit records
  const one = ["0|1", "1|0", ".|0", "2", ".", "1|1:35", "0/2"];
  assert.deepEqual([...gtCodes(rec(one), one.length, 2)], [...gtCodes(rec(one), one.length, 12)]);
  assert.deepEqual([...gtCodes(rec(one), one.length, 2)], [0, 1, 1, 0, -1, 0, 2, -2, -1, -2, 1, 1, 0, 2]);
  // multi-digit indexes: 10, 12 are alleles 10 and 12, not '1' then garbage
  const many = ["10|0", "0|12", "11", "12|.:7", "9|10"];
  assert.deepEqual([...gtCodes(rec(many), many.length, 12)], [10, 0, 0, 12, 11, -2, 12, -1, 9, 10]);
});

test("K4: masked band row ranges partition the full masked band (both backends, any split)", async () => {
  // synthetic chrX-like rows: 300 variants, 40 samples (8 haploid males), some uncalled haplotypes, close positions
  let seed = 7; const rnd = () => ((seed = (seed * 1103515245 + 12345) >>> 0) / 2 ** 32);
  const samples = 40, hb = Math.ceil((2 * samples) / 8), males = new Set([1, 4, 9, 15, 22, 27, 33, 38]);
  const rows = Array.from({ length: 300 }, (_, k) => {
    const haps = new Uint8Array(hb), missing = k % 7 === 0 ? new Uint8Array(hb) : null, f = 0.1 + 0.8 * rnd();
    for (let h = 0; h < 2 * samples; h++) { if (rnd() < f) haps[h >> 3] |= 1 << (h & 7); if (missing && rnd() < 0.05) missing[h >> 3] |= 1 << (h & 7); }
    return { pos: 1000 + k * 37 + Math.floor(rnd() * 30), haps, missing, male_haploid: k % 3 !== 0 };
  }).sort((a, b) => a.pos - b.pos);
  const opts = { samples, haploidMales: males, hapIndexes: Array.from({ length: 2 * samples }, (_, i) => i) };
  const key = b => Array.from({ length: b.count }, (_, k) => `${b.i[k]},${b.j[k]},${b.r2[k]},${b.dprime[k]}`);
  for (const k of [await loadKernel(readFileSync(KERNEL_WASM)), jsKernel()]) {
    const p = k.prepare(rows, opts); assert.ok(p.masked, "rows with male_haploid / missing must use the masked path");
    const full = k.band(p, { maxDist: 2000, threshold: 0 });
    for (const cuts of [[0, 300], [0, 1, 300], [0, 150, 300], [0, 37, 120, 121, 299, 300]]) {
      const parts = [];
      for (let s = 0; s + 1 < cuts.length; s++) parts.push(...key(k.band(p, { maxDist: 2000, threshold: 0, rows: [cuts[s], cuts[s + 1]] })));
      assert.deepEqual(parts, key(full), `${k.backend}: split ${cuts.join("/")}`);
    }
    k.reset();
  }
});

test("K5: masked rows with identical masks (one popcount) match the general path, JS = WASM, counts() too", async () => {
  // chrX-like rows: most share the haploid-male mask (same class); every 5th row has its own uncalled haplotypes
  let seed = 11; const rnd = () => ((seed = (seed * 1103515245 + 12345) >>> 0) / 2 ** 32);
  const samples = 40, hb = Math.ceil((2 * samples) / 8), males = new Set([0, 3, 8, 14, 21, 26, 30, 37]);
  const rows = Array.from({ length: 200 }, (_, k) => {
    const haps = new Uint8Array(hb), missing = k % 5 === 0 ? new Uint8Array(hb) : null, f = 0.05 + 0.9 * rnd();
    for (let h = 0; h < 2 * samples; h++) { if (rnd() < f) haps[h >> 3] |= 1 << (h & 7); if (missing && rnd() < 0.08) missing[h >> 3] |= 1 << (h & 7); }
    return { pos: 5000 + k * 41, haps, missing, male_haploid: k % 17 !== 0 };
  });
  const opts = { samples, haploidMales: males, hapIndexes: Array.from({ length: 2 * samples }, (_, i) => i).filter(i => i % 9 !== 4) };
  const key = b => Array.from({ length: b.count }, (_, k) => `${b.i[k]},${b.j[k]},${b.r2[k]},${b.dprime[k]}`);
  const I = Array.from({ length: 199 }, (_, i) => i), J = I.map(i => (i * 7 + 3) % 200);
  const results = [];
  for (const k of [await loadKernel(readFileSync(KERNEL_WASM)), jsKernel()]) {
    const p = k.prepare(rows, opts), v = p.views();
    const classes = new Set(v.cls), shared = [...v.cls].filter((c, i) => v.cls.indexOf(c) !== i).length;
    assert.ok(classes.size > 1 && shared > 100, `${k.backend}: the fixture must mix shared and distinct masks (${classes.size} classes)`);
    const fast = k.band(p, { maxDist: 3000, threshold: 0 }), fastCounts = k.counts(p, I, J);
    v.cls.set(Int32Array.from({ length: p.n }, (_, i) => i)); // every row its own class: the general path everywhere
    const general = k.band(p, { maxDist: 3000, threshold: 0 }), generalCounts = k.counts(p, I, J);
    assert.ok(fast.count > 1000);
    assert.deepEqual(key(fast), key(general), `${k.backend}: band`);
    assert.deepEqual(fastCounts, generalCounts, `${k.backend}: counts`);
    results.push(key(fast));
    k.reset();
  }
  assert.deepEqual(results[0], results[1], "wasm = js");
});
