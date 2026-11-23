// Verify select_rows.mjs on a rebuilt genotype set against independent ground truth:
//   parity — a direct pass over the source 1000G VCF: for each query position, the LAST VCF record whose REF overlaps
//            it (tabix region semantics), compared by (pos, ref, alt list) with the picked row's source record
//   allele — rows whose alleles match the rsID's NCBI alleles, recomputed here from the NCBI slice (not from the
//            build's allele_rsids column)
// usage: node verify_select.mjs <vcf.gz> <genoPrefix e.g. out/v3_grch37_chr22> <ncbi_slice.parquet> <out.json> [nSample]
//   The NCBI slice must carry p_dump (LDlink's position from the dump, see compare_ncbi.sql / step 8).
import { createReadStream, writeFileSync } from "node:fs";
import { createGunzip } from "node:zlib";
import { createInterface } from "node:readline";
import { asyncBufferFromFile, parquetMetadataAsync, parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";
import { candidateRows, pickParity, pickAllele } from "../src/genotypes.mjs";

const [vcfPath, prefix, ncbiPath, outPath, nArg] = process.argv.slice(2);
const nSample = +(nArg ?? 20000);
const num = v => (typeof v === "bigint" ? Number(v) : v);
function minimal(pos, ref, alt) {
  let r = ref, a = alt, p = pos;
  while (r.length > 1 && a.length > 1 && r.at(-1) === a.at(-1)) { r = r.slice(0, -1); a = a.slice(0, -1); }
  while (r.length > 1 && a.length > 1 && r[0] === a[0]) { r = r.slice(1); a = a.slice(1); p++; }
  return `${p}:${r}:${a}`;
}
// NCBI slice → deterministic sample of rsIDs with a dump position (every k-th; plus all that fall near 1000G rows)
const nf = await asyncBufferFromFile(ncbiPath), nm = await parquetMetadataAsync(nf);
const N = { rsid: [], pos: [], ref: [], alt: [], p_dump: [] };
await parquetRead({ file: nf, metadata: nm, compressors, columns: Object.keys(N), rowFormat: "array",
  onComplete: rows => { for (const r of rows) Object.keys(N).forEach((c, i) => N[c].push(num(r[i]))); } });
// genotype files
const files = await Promise.all(["common", "low", "rare"].map(async tier => {
  const file = await asyncBufferFromFile(`${prefix}/${tier === "rare" ? "rare_sparse_p16k" : `${tier}_bits_p64k`}.parquet`);
  return { tier, file, metadata: await parquetMetadataAsync(file) };
}));
const cols = ["pos", "ref", "alt", "vcf_order", "allele_is_label", "allele_rsids", "multiallelic_source", "rsid"];
// sample: rsIDs whose dump position lies within 1 kb of a 1000G row are the interesting ones; take every k-th of those
// only rsIDs whose dump position overlaps some 1000G row (the cases where parity and allele mode choose a row)
const spans = []; for (const f of files) await parquetRead({ file: f.file, metadata: f.metadata, columns: ["pos", "ref"], compressors, onComplete: r => { for (const [p, ref] of r) spans.push([p, p + ref.length - 1]); } });
spans.sort((a, b) => a[0] - b[0]);
const maxLen = spans.reduce((m, [a, b]) => Math.max(m, b - a + 1), 1);
const near = p => { let a = 0, b = spans.length; while (a < b) { const h = (a + b) >> 1; if (spans[h][0] < p - maxLen + 1) a = h + 1; else b = h; } for (let j = a; j < spans.length && spans[j][0] <= p; j++) if (spans[j][1] >= p) return true; return false; };
const idx = []; for (let i = 0; i < N.rsid.length; i++) if (N.p_dump[i] != null && near(N.p_dump[i])) idx.push(i);
const step = Math.max(1, Math.floor(idx.length / nSample)), sample = idx.filter((_, k) => k % step === 0).slice(0, nSample);
// plus rsIDs whose dump position overlaps NO 1000G row (the not_in_1000g path), 10% of the sample size
const far = []; for (let i = 0; i < N.rsid.length; i++) if (N.p_dump[i] != null && !near(N.p_dump[i])) far.push(i);
const fstep = Math.max(1, Math.floor(far.length / Math.ceil(nSample / 10)));
const farSample = far.filter((_, k) => k % fstep === 0).slice(0, Math.ceil(nSample / 10));
const farSet = new Set(farSample);
sample.push(...farSample);
const t0 = performance.now();
const { rows, reach } = await candidateRows(files, sample.map(i => N.p_dump[i]), { columns: cols });
rows.sort((a, b) => a.pos - b.pos);
const fetchMs = performance.now() - t0;
const windowRows = p => { let a = 0, b = rows.length; while (a < b) { const h = (a + b) >> 1; if (rows[h].pos < p - reach + 1) a = h + 1; else b = h; } const o = []; for (let j = a; j < rows.length && rows[j].pos <= p; j++) o.push(rows[j]); return o; };

// ground truth (parity): stream the VCF once; for each sampled position, keep the last record overlapping it
const want = new Map(sample.map(i => [N.p_dump[i], null])), qs = [...want.keys()].sort((a, b) => a - b);
const rl = createInterface({ input: createReadStream(vcfPath).pipe(createGunzip()), crlfDelay: Infinity });
let qi = 0;
const lastIsSv = new Set(); // positions where LDlink's last overlapping record is a symbolic SV
for await (const line of rl) {
  if (line.charCodeAt(0) === 35) continue;
  const t = line.split("\t", 6), pos = +t[1], ref = t[3], alt = t[4];
  const sv = alt.includes("<"), end = pos + ref.length - 1; // symbolic SVs are not in the release (build skips them)
  while (qi < qs.length && qs[qi] < pos - 1000) qi++;
  for (let k = qi; k < qs.length && qs[k] <= end; k++) if (qs[k] >= pos) { if (sv) lastIsSv.add(qs[k]); else { lastIsSv.delete(qs[k]); want.set(qs[k], { pos, ref, alt }); } }
}
const res = { vcf: vcfPath.split(/[\\/]/).pop(), genotypes: prefix, ncbi: ncbiPath, sampled_rsids: sample.length, reach, candidate_rows: rows.length, fetch_ms: Math.round(fetchMs), ldlink_last_record_is_sv: lastIsSv.size, not_near_sampled: farSample.length,
  parity: { agree: 0, disagree: 0, examples: [] }, allele: { agree: 0, disagree: 0, examples: [], status: {} } };
for (const i of sample) {
  const p = N.p_dump[i], w = windowRows(p);
  // parity
  const got = pickParity(w, p), exp = want.get(p);
  const ok = exp == null ? got.row == null : got.row != null && got.row.pos === exp.pos && got.row.ref === exp.ref && exp.alt.split(",").includes(got.row.alt) && exp.alt.split(",").at(-1) === got.row.alt;
  // LDlink's last overlapping record is a symbolic SV, which the build skips: a known, separate class (not agreement)
  if (lastIsSv.has(p)) { res.parity.ldlink_picks_sv = (res.parity.ldlink_picks_sv ?? 0) + 1; if ((res.parity.sv_examples ??= []).length < 5) res.parity.sv_examples.push({ p, got: got.row && `${got.row.pos}:${got.row.ref}>${got.row.alt}` }); }
  else if (farSet.has(i) && !(exp == null && got.row == null)) { res.parity.disagree++; res.parity.examples.push({ p, why: "expected not_in_1000g", exp, got: got.row }); }
  else if (ok) res.parity.agree++; else { res.parity.disagree++; if (res.parity.examples.length < 10) res.parity.examples.push({ p, exp, got: got.row && { pos: got.row.pos, ref: got.row.ref, alt: got.row.alt, vcf_order: got.row.vcf_order } }); }
  // allele
  const keys = new Set(N.alt[i].split(",").map(a => minimal(N.pos[i], N.ref[i], a)));
  const expRows = w.filter(r => keys.has(minimal(r.pos, r.ref, r.alt))).map(r => r.vcf_order).sort((a, b) => a - b);
  const a = pickAllele(w, p, N.rsid[i]);
  res.allele.status[a.status] = (res.allele.status[a.status] ?? 0) + 1;
  const gotRows = a.candidates ? a.candidates.map(r => r.vcf_order) : a.row ? [a.row.vcf_order] : [];
  if (JSON.stringify(gotRows) === JSON.stringify(expRows)) res.allele.agree++;
  else { res.allele.disagree++; if (res.allele.examples.length < 10) res.allele.examples.push({ rsid: N.rsid[i], ncbi: `${N.pos[i]}:${N.ref[i]}>${N.alt[i]}`, p, expRows, gotRows, status: a.status }); }
}
writeFileSync(outPath, JSON.stringify(res, null, 2));
console.log(JSON.stringify(res, (k, v) => (k === "examples" && !v.length ? undefined : v)));
