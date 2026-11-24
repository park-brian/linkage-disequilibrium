// Verify the Parquet parity resolver (resolve.js) against answers computed from the dump with Mongo semantics
// (parity_expected.sql), plus the reference test-vector SNPs. usage: node verify/data/parity.js <releaseRoot>
import { readFileSync, writeFileSync } from "node:fs";
import { asyncBufferFromFile, cachedAsyncBuffer } from "hyparquet";
import { createResolver } from "../../src/parity/resolve.js";
import { useDataCwd, repo } from "../../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)

const [root] = process.argv.slice(2);
const resolver = createResolver(root, async p => cachedAsyncBuffer(await asyncBufferFromFile(p)));
const load = f => JSON.parse(readFileSync(`dbsnp_dump/${f}`, "utf8"));
const report = {};
const posOf = k => (k == null ? null : Number(BigInt(k) % 1_000_000_000n));
const chromName = c => (c == null ? null : c === 23 ? "X" : c === 24 ? "Y" : c === 25 ? "MT" : String(c));

// 1) rsID → position (find_one semantics), both builds
const exp = load("parity_rsid.json");
for (const build of ["grch37", "grch38"]) {
  const got = await resolver.resolveRsids(exp.map(e => `rs${e.rsid}`), build);
  const bad = [];
  exp.forEach((e, i) => {
    const g = got[i], key = build === "grch37" ? e.key37 : e.key38;
    const ok = e.found ? g.found && g.pos === posOf(key) && g.chrom === chromName(e.chrom) && (g.ref_rsid ?? null) === (e.ref_rsid == null ? null : `rs${e.ref_rsid}`) : !g.found;
    if (!ok) bad.push({ e, g });
  });
  report[`rsid_${build}`] = { cases: exp.length, found: exp.filter(e => e.found).length, mismatches: bad.length, examples: bad.slice(0,
    5) };
}
// 2) chr:pos → rsID (replace_coord_rsid semantics)
for (const [build, f] of [["grch37", "parity_pos37.json"], ["grch38", "parity_pos38.json"]]) {
  const cases = load(f);
  const inputs = cases.map(c => { const k = BigInt(c.key), ch = chromName(Number(k / 1_000_000_000n));
    return `chr${ch}:${k % 1_000_000_000n}`; });
  const got = await resolver.resolveCoords(inputs, build);
  const bad = [];
  cases.forEach((c, i) => {
    const want = c.expected_rsid == null ? null : `rs${c.expected_rsid}`;
    if (got[i].rsid !== want || !!got[i].warning !== c.expect_warning) bad.push({ input: inputs[i], want, warn: c.expect_warning,
      got: got[i] });
  });
  report[`coord_${build}`] = { cases: cases.length, multi_doc: cases.filter(c => c.ids?.length > 1).length,
    warnings_expected: cases.filter(c => c.expect_warning).length, mismatches: bad.length, examples: bad.slice(0, 5) };
}
// 3) reference test vectors (REQUIREMENTS.md §5): same positions the reference LDmatrix run used
const TEST_VECTORS = { grch37: repo("verify/reference/results/tv_grch37_ceu.json"),
  grch38: repo("verify/reference/results/tv_grch38_ceu.json") };
for (const [build, f] of Object.entries(TEST_VECTORS)) {
  const snps = JSON.parse(readFileSync(f, "utf8")).snps;
  const got = await resolver.resolveRsids(snps.map(s => s.rs), build);
  report[`test_vectors_${build}`] = { snps: snps.length, mismatches: snps.filter((s,
    i) => !(got[i].found && got[i].pos === s.pos && got[i].chrom === "2")).map(s => s.rs) };
}
// 4) V3: LDlink's whole input chain on non-canonical inputs (validsnp first token, replace_coord_rsid incl. warning
//    TEXT, get_coords strip("rs") + int()), ground truth from parity_inputs.py (LDlink's Python expressions + the
//    dump), compared with resolver.resolveInputs and the planner (auto and forced-direct rsID resolution)
{
  const { execFileSync } = await import("node:child_process");
  const { openPlanner } = await import("../../src/parity/planner.js");
  const rs = exp.slice(0, 120).map(e => String(e.rsid)), merged = exp.filter(e => e.ref_rsid != null && e.ref_rsid !== e.rsid).slice(0,
    30).map(e => String(e.rsid));
  const multi = load("parity_pos37.json").filter(c => c.ids?.length > 1).slice(0, 60),
    single = load("parity_pos37.json").filter(c => c.ids?.length === 1).slice(0, 60);
  const coord = c => { const k = BigInt(c.key), ch = chromName(Number(k / 1_000_000_000n)); return [ch, k % 1_000_000_000n]; };
  const inputs = [
    ...rs.flatMap(id => [`rs${id}`, `RS${id}`, `  rs${id}\t`, `rs${id}rs`, `rsrs${id}`, `rs0${id}`, `rs+${id}`, `rs${id} rs1`, `${id}`,
      `r${id}`, `rs${id.slice(0, 2)}_${id.slice(2)}`]),
    ...merged.map(id => `rs${id}`),
    ...[...multi, ...single].flatMap(c => { const [ch, p] = coord(c); return [`chr${ch}:${p}`, `CHR${ch}:${p}`, `${ch}:${p}`,
      `chr${ch}:${p} extra`, `chr${ch}:0${p}`, `chr${ch}:+${p}`, `chrchr${ch}:${p}`, `chr${ch.toLowerCase()}:${p}:x`]; }),
    "", "   ", "rs", "rs abc", "chr22:", "chr22", ":123", "chr22:abc", "chr22:1.5", "chr22:-5", "chr23:100", "chrm:100", "chrx:2700000",
      "x:2700000", "Y:2700000", "rs-5", "rs99999999999999",
  ];
  writeFileSync("dbsnp_dump/parity_inputs_in.json", JSON.stringify(inputs));
  for (const build of ["grch37", "grch38"]) {
    execFileSync("python", ["-I", repo("verify/data/parity_inputs.py"), "dbsnp_dump/parity_inputs_in.json", build,
      `dbsnp_dump/parity_inputs_${build}.json`], { stdio: ["ignore", "inherit", "inherit"] });
    const want = JSON.parse(readFileSync(`dbsnp_dump/parity_inputs_${build}.json`, "utf8"));
    const k = v => JSON.stringify({ found: !!v.found, rsid: v.found ? v.rsid : null, ref_rsid: v.found ? v.ref_rsid ?? null : null,
      chrom: v.found ? v.chrom : null, pos: v.found ? v.pos ?? null : null, warning: v.warning ?? null });
    const got = { resolver: await resolver.resolveInputs(inputs, build) };
    for (const forceRsid of [undefined, "direct"]) {
      const pl = await openPlanner(root, async p => cachedAsyncBuffer(await asyncBufferFromFile(p)));
      got[`planner_${forceRsid ?? "auto"}`] = (await pl.query({ build, inputs, forceRsid })).variants;
    }
    for (const [who, g] of Object.entries(got)) {
      const bad = inputs.map((x, i) => [x, k(g[i]), k(want[i])]).filter(([, a, b]) => a !== b).map(([input, got, want]) => ({ input, got,
        want }));
      report[`inputs_${build}_${who}`] = { cases: inputs.length, found: want.filter(w => w.found).length,
        warnings: want.filter(w => w.warning).length, mismatches: bad.length, examples: bad.slice(0, 5) };
    }
    // LDmatrix's own gate (module logic, ROADMAP Phase C): found inputs it would still reject
    report[`inputs_${build}_ldmatrix_gate_rejects`] = want.filter(w => w.found && !w.ldmatrix_gate).map(w => w.input).slice(0, 20);
  }
}
writeFileSync("verify/parity_report.json", JSON.stringify(report, null, 2)); // $LD_DATA/verify, not the release
if (Object.values(report).some(r => r && (Array.isArray(r.mismatches) ? r.mismatches.length : r.mismatches) > 0)) process.exitCode = 1;
console.log(JSON.stringify(report, (k, v) => (k === "examples" && !v.length ? undefined : v), 2));
