// LDlink's upload normalizer (server/sumstats_normalizer.py via POST /upload + GET /validate_sumstats) against ours
// (src/parity/ldsc/normalize.js): per uploaded file, the validation answer (valid, errors, warnings, detected format,
// normalized name) and the normalized file byte for byte. Inputs are built from LDlink's example sumstats: each format
// the normalizer knows, its transforms, and the ways a file fails (formats, rows, quoting, names, selections).
import { normalizeSumstats, secureFilename } from "../../../src/parity/ldsc/normalize.js";
import { uploadSumstats, uploadedFile, example } from "../lib.js";

/** the example as rows of named fields (SNP CHR POS A2 A1 Frq Rsq BETA SE P LOG10P N) */
const parse = text => {
  const [head, ...rows] = text.split(/\r?\n/).filter(Boolean).map(l => l.split("\t"));
  return rows.map(r => Object.fromEntries(head.map((h, k) => [h, r[k]])));
};
const table = (cols, rows, { sep = "\t", eol = "\n" } = {}) => [cols.join(sep), ...rows.map(r => r.join(sep))].join(eol) + eol;
const or = beta => String(Math.exp(+beta));

/** named upload files from the example's rows (n rows of it), one per format and failure the normalizer has */
export function uploadFiles(rows) {
  const swap = (r, k) => (k % 3 === 0 ? [r.A1, r.A2] : [r.A2, r.A1]); // [REF, ALT]: every third row has A1 as REF
  const plink = rows.map((r, k) => [r.SNP, r.A1, ...swap(r, k), r.N, r.BETA, r.P]);
  return {
    ldsc: { name: "ldsc_ready.txt", text: table(["SNP", "A1", "A2", "N", "P", "BETA"], rows.map(r => [r.SNP, r.A1, r.A2, r.N, r.P,
      r.BETA])) },
    plink: { name: "plink_ref_alt.txt", text: table(["ID", "A1", "REF", "ALT", "OBS_CT", "BETA", "P"], plink) },
    plink_or: { name: "plink_or.glm", text: table(["#CHROM", "ID", "A1", "REF", "ALT", "OBS_CT", "OR", "P"],
      rows.map((r, k) => [r.CHR, r.SNP, r.A1, ...swap(r, k), r.N, or(r.BETA), r.P])) },
    plink_a2: { name: "plink_a2.txt", text: table(["ID", "A1", "A2", "OBS_CT", "BETA", "P"], rows.map(r => [r.SNP, r.A1, r.A2, r.N, r.BETA,
      r.P])) },
    regenie: { name: "regenie.regenie", text: table(["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "N", "BETA", "LOG10P"],
      rows.map(r => [r.CHR, r.POS, r.SNP, r.A2, r.A1, r.N, r.BETA, r.LOG10P]), { sep: " " }) },
    saige: { name: "saige.saige", text: table(["CHR", "POS", "MarkerID", "Allele1", "Allele2", "N", "OR", "p.value"],
      rows.map(r => [r.CHR, r.POS, r.SNP, r.A2, r.A1, r.N, or(r.BETA), r.P])) },
  };
}

export async function cases(ctx) {
  const ex = parse(await example(ctx, "BBJ_HDLC22.txt")), few = ex.slice(0, 300), f = uploadFiles(few), cs = [];
  for (const [k, v] of Object.entries(f)) cs.push({ name: `format ${k}`, file: v });
  const bad = (name, text, extra = {}) => cs.push({ name, file: { name: `${name.replace(/\W+/g, "_")}.txt`, text }, ...extra });
  const rowsOf = n => few.slice(0, n);
  // PLINK quirks
  bad("plink N not OBS_CT", table(["ID", "A1", "REF", "ALT", "N", "BETA", "P"], rowsOf(20).map(r => [r.SNP, r.A1, r.A2, r.A1, r.N, r.BETA,
    r.P])));
  bad("plink A1 matches neither", table(["ID", "A1", "REF", "ALT", "OBS_CT", "BETA", "P"],
    rowsOf(15).map(r => [r.SNP, "Z", r.A2, r.A1, r.N, r.BETA, r.P])));
  bad("plink REF only, A1 is REF", table(["ID", "A1", "REF", "OBS_CT", "BETA", "P"], rowsOf(5).map(r => [r.SNP, r.A1, r.A1, r.N, r.BETA,
    r.P])));
  bad("plink ALT only", table(["ID", "A1", "ALT", "OBS_CT", "BETA", "P"], rowsOf(5).map(r => [r.SNP, r.A1, r.A2, r.N, r.BETA, r.P])));
  bad("plink lower-case alleles", table(["ID", "A1", "REF", "ALT", "OBS_CT", "BETA", "P"],
    rowsOf(5).map(r => [r.SNP, r.A1.toLowerCase(), r.A2.toLowerCase(), r.A1, r.N, r.BETA, r.P])));
  // values
  const pv = ["0.5", "1", "1e-300", "0.000_5", " 0.25 ", "nan", "inf", "0", "1.5", "-0.1", "abc", "", "5e-324"];
  bad("p values", table(["ID", "A1", "A2", "OBS_CT", "BETA", "P"], pv.map((p, k) => [`rs${k + 1}`, "A", "G", "100", "0.1", p])));
  bad("regenie LOG10P edges", table(["ID", "ALLELE0", "ALLELE1", "N", "BETA", "LOG10P"],
    ["0", "1", "300", "-1", "nan", "inf", "2.5"].map((p, k) => [`rs${k + 1}`, "A", "G", "100", "0.1", p])));
  bad("odds ratios", table(["MarkerID", "Allele1", "Allele2", "N", "OR", "p.value"],
    ["1", "2.5", "0", "-1", "0.001", "1e300"].map((o, k) => [`rs${k + 1}`, "A", "G", "100", o, "0.5"])));
  bad("too many bad rows (10 kept)", table(["ID", "A1", "A2", "OBS_CT", "BETA", "P"],
    Array.from({ length: 14 }, (_, k) => [`rs${k}`, "A", "G", "100", "0.1", "2"])));
  // reading
  bad("csv with quotes and CRLF", table(["ID", "A1", "A2", "OBS_CT", "BETA", "P"],
    [['"rs1"', "A", "G", "100", '"0.1"', "0.5"], ['"rs2,x"', "C", "T", "100", "0.2", '"0.25"'], ['"rs""3"', "A", "C", "100", "0.3", "0.1"]],
    { sep: ",", eol: "\r\n" }));
  bad("csv header with spaces", "ID, A1, A2, OBS_CT, BETA, P\nrs1,A,G,100,0.1,0.5\n");
  bad("short and long rows", "ID\tA1\tA2\tOBS_CT\tBETA\tP\nrs1\tA\tG\t100\t0.1\nrs2\tA\tG\t100\t0.1\t0.5\textra\n");
  bad("whitespace, blank lines", "ID   A1 A2\tOBS_CT BETA P\n\nrs1  A  G 100 0.1 0.5\n   \nrs2 C T 100 0.2 0.25\n");
  bad("duplicate columns", "ID\tA1\tA2\tOBS_CT\tBETA\tP\tP\nrs1\tA\tG\t100\t0.1\t0.9\t0.5\n");
  bad("mixed-case LDSC headers", "rsid\tPval\tn\tBeta\nrs1\t0.5\t100\t0.1\n");
  bad("BOM in header", "﻿ID\tA1\tA2\tOBS_CT\tBETA\tP\nrs1\tA\tG\t100\t0.1\t0.5\n");
  bad("header only", "ID\tA1\tA2\tOBS_CT\tBETA\tP\n");
  bad("empty file", "");
  bad("no format matches", "rsid\tchr\tpos\nrs1\t1\t100\n");
  // selections and names
  cs.push({ name: "selected plink_raw on a REGENIE file", file: f.regenie, format: "plink_raw" });
  cs.push({ name: "selected regenie_raw", file: f.regenie, format: "regenie_raw" });
  cs.push({ name: "selected pre_munged on a PLINK file", file: f.plink, format: "pre_munged" });
  cs.push({ name: "selected unknown format", file: f.plink, format: "bogus" });
  cs.push({ name: "a name with spaces and accents", file: { name: "my sümstats (v2).txt", text: f.plink.text } });
  return cs;
}

const answer = fv => ({ valid: fv.valid, errors: fv.errors, warnings: fv.warnings, detected_format: fv.detected_format ?? null,
  normalized_filename: fv.normalized_filename });
export async function ldlink(cs) {
  return Promise.all(cs.map(async c => {
    const { reference, files: [{ name, fileValid }] } = await uploadSumstats([c.file], c.format);
    const renamed = fileValid.valid && fileValid.normalized_filename !== name;
    return { saved: name, ...answer(fileValid), normalized: renamed ? await uploadedFile(reference, fileValid.normalized_filename) : null };
  }));
}
export async function ours(ctx, c) {
  const name = secureFilename(c.file.name), n = normalizeSumstats(c.file.text, name, c.format);
  return { saved: name, ...answer(n), normalized: n.valid && n.normalized_filename !== name ? n.text : null };
}
