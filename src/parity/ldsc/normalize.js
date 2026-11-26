// LDlink's upload normalizer (server/sumstats_normalizer.py, "sumstats-normalizer-v1"), run by its validate_sumstats
// route on every uploaded summary statistics file before LDherit / LDcorrelation: detect the format (LDSC-ready,
// PLINK, REGENIE, SAIGE, in that order, or the one the user selected), and unless it is LDSC-ready, rewrite it as
// <name>.ldsc.tsv with the columns SNP A1 A2 N P BETA (p from LOG10P, BETA = log(OR), every P as format(p, ".15g"));
// the analysis then runs on that file. An LDSC-ready file is used as uploaded.
// Reproduced as Python reads and writes it: csv.reader / DictReader (excel dialect) when the header has a comma or a
// tab, else str.split(); the DictReader's row keys are the raw header fields while columns are looked up by their
// stripped names; csv.DictWriter's minimal quoting and "\r\n" line ends; str.strip() / split() whitespace; float().
// Bugs kept: a PLINK file's N comes only from OBS_CT (an N column satisfies detection, then every row fails).
import { cFormatG, pyFloat } from "../pyfmt.js";

const OUT = ["SNP", "A1", "A2", "N", "P", "BETA"];
const FORMATS = [
  { name: "LDSC-ready", required: { SNP: ["SNP", "MARKERNAME", "SNPID", "RS", "RSID", "RS_NUMBER", "RS_NUMBERS"],
    P: ["P", "PVALUE", "P_VALUE", "PVAL", "P_VAL", "GC_PVALUE"], N: ["N", "N_CAS", "N_CON"],
    SIGNED: ["Z", "ZSCORE", "Z_SCORE", "OR", "B", "BETA", "LOG_ODDS"] }, canonical: [], transforms: {} },
  { name: "PLINK", required: { SNP: ["ID"], N: ["OBS_CT", "N"], A1: ["A1"], A2: ["REF", "A2", "ALT"], P: ["P"], BETA: ["BETA", "OR"] },
    canonical: [["ID", "SNP"], ["OBS_CT", "N"], ["A1", "A1"], ["REF", "A2"], ["A2", "A2"], ["ALT", "A2"], ["P", "P"], ["BETA", "BETA"],
      ["OR", "BETA"]], transforms: { OR: "log_or" } },
  { name: "REGENIE", required: { SNP: ["ID"], N: ["N"], A1: ["ALLELE1", "A1"], A2: ["ALLELE0", "A2"], P: ["P", "LOG10P"],
    BETA: ["BETA", "OR"] }, canonical: [["ID", "SNP"], ["N", "N"], ["ALLELE1", "A1"], ["A1", "A1"], ["ALLELE0", "A2"], ["A2", "A2"],
    ["P", "P"], ["LOG10P", "P"], ["BETA", "BETA"], ["OR", "BETA"]], transforms: { LOG10P: "log10p", OR: "log_or" } },
  { name: "SAIGE", required: { SNP: ["MARKERID"], N: ["N"], A1: ["ALLELE2"], A2: ["ALLELE1"], P: ["P_VALUE", "P"], BETA: ["BETA", "OR"] },
    canonical: [["MARKERID", "SNP"], ["N", "N"], ["ALLELE2", "A1"], ["ALLELE1", "A2"], ["P_VALUE", "P"], ["P", "P"], ["BETA", "BETA"],
      ["OR", "BETA"]], transforms: { OR: "log_or" } },
];
const SELECTED = { pre_munged: "LDSC-ready", plink_raw: "PLINK", regenie_raw: "REGENIE", saige_raw: "SAIGE" };

// str.strip() / str.split() whitespace (str.isspace() code points; not the BOM)
const WS = "\\t-\\r\\x1c-\\x20\\x85\\xa0\\u1680\\u2000-\\u200a\\u2028\\u2029\\u202f\\u205f\\u3000";
const STRIP = new RegExp(`^[${WS}]+|[${WS}]+$`, "gu"), SPLIT = new RegExp(`[${WS}]+`, "u");
const strip = s => String(s).replace(STRIP, "");
const split = s => strip(s).split(SPLIT).filter(Boolean);
const clean = h => strip(h).toUpperCase().replaceAll("-", "_").replaceAll(".", "_");

/** the text's lines as Python's text mode with newline="" yields them: split at \r\n, \r or \n, endings kept */
function lines(text) {
  const out = [];
  for (let i = 0, s = 0; i <= text.length; i++) {
    const c = text[i];
    if (i === text.length) { if (s < i) out.push(text.slice(s)); break; }
    if (c === "\n" || c === "\r") { const e = c === "\r" && text[i + 1] === "\n" ? i + 2 : i + 1; out.push(text.slice(s, e)); s = e;
      i = e - 1; }
  }
  return out;
}
/** csv.reader (excel dialect: quotechar ", doublequote, not strict) over the text → records (arrays of fields) */
function csvRecords(text, delim) {
  const recs = [];
  let rec = [], f = "", state = "start", i = 0;
  const endField = () => { rec.push(f); f = ""; }, endRec = () => { endField(); recs.push(rec); rec = []; };
  for (; i < text.length; i++) {
    const c = text[i];
    if (state === "quoted") {
      if (c === '"') state = "quote";
      else f += c;
    } else if (state === "quote") { // after a quote inside a quoted field
      if (c === '"') { f += '"'; state = "quoted"; }
      else if (c === delim) { endField(); state = "start"; }
      else if (c === "\n" || c === "\r") { endRec(); state = "eol"; if (c === "\r" && text[i + 1] === "\n") i++; }
      else { f += c; state = "field"; }
    } else if (c === "\n" || c === "\r") {
      if (state === "start" && rec.length === 0) { if (c === "\r" && text[i + 1] === "\n") i++; recs.push([]); continue; } // a blank line
      endRec(); state = "start"; if (c === "\r" && text[i + 1] === "\n") i++;
    } else if (c === delim) { endField(); state = "start"; }
    else if (c === '"' && state === "start") state = "quoted";
    else { f += c; state = "field"; }
    if (state === "eol") state = "start";
  }
  if (state !== "start" || rec.length || f) endRec();
  return recs;
}
/** csv.DictWriter row (tab, minimal quoting) */
const csvRow = vals => vals.map(v => (/[\t"\r\n]/.test(v) ? `"${v.replaceAll('"', '""')}"` : v)).join("\t") + "\r\n";

/** the file as the normalizer reads it → { columns (stripped), rows: [{ key → value }] } */
function read(text) {
  const ls = lines(text);
  if (!ls.length || ls[0] === "") return { columns: [], rows: [] };
  const head = ls[0], delim = head.includes(",") && !head.includes("\t") ? "," : head.includes("\t") ? "\t" : null;
  if (!delim) {
    const header = split(head), rows = [];
    for (const l of ls.slice(1)) {
      const v = split(l);
      if (v.length) rows.push(Object.fromEntries(header.slice(0, v.length).map((h, k) => [h, v[k]]))); // dict(zip(header, values))
    }
    return { columns: header, rows };
  }
  const [fields = [], ...recs] = csvRecords(text, delim);
  const columns = csvRecords(head, delim)[0]?.map(strip).filter(Boolean) ?? [];
  const rows = recs.filter(r => r.length).map(r => Object.fromEntries(fields.map((h, k) => [h, r[k] ?? ""]))); // DictReader
  return { columns, rows };
}
const byClean = columns => Object.fromEntries(columns.map(c => [clean(c), c]));
function matches(fmt, cols) {
  const mapped = {}, missing = [];
  for (const [canon, aliases] of Object.entries(fmt.required)) {
    const col = aliases.map(a => cols[a]).find(Boolean);
    if (col) mapped[col] = canon; else missing.push(`Missing required ${canon} column. Expected one of: ${aliases.join(", ")}`);
  }
  return { ok: !missing.length, mapped, missing };
}

const num = (value, column) => {
  const x = pyFloat(value);
  if (x === null) throw new Error(`Column ${column} contains a non-numeric value: ${value}`);
  return x;
};
const probability = (value, column) => {
  const x = num(value, column);
  if (!(0 < x && x <= 1)) throw new Error(`Column ${column} contains a p-value outside (0, 1]: ${value}`);
  return cFormatG(x, 15);
};
function transform(value, how, column) {
  if (how === "log10p") {
    const x = num(value, column);
    if (x < 0) throw new Error(`Column ${column} contains a negative -log10(p) value: ${value}`);
    return cFormatG(Math.pow(10, -x), 15);
  }
  if (how === "log_or") {
    const x = num(value, column);
    if (x <= 0) throw new Error(`Column ${column} contains a non-positive odds ratio: ${value}`);
    return cFormatG(Math.log(x), 15);
  }
  return strip(value);
}
/** PLINK's A2 per row (REF / ALT are unordered with respect to A1) */
function plinkA2(row, a1, refCol, altCol, a2Col) {
  if (a2Col) return strip(row[a2Col] ?? "");
  a1 = strip(a1).toUpperCase();
  const ref = refCol ? strip(row[refCol] ?? "").toUpperCase() || null : null,
    alt = altCol ? strip(row[altCol] ?? "").toUpperCase() || null : null;
  if (ref !== null && alt !== null) {
    if (a1 === ref) return alt;
    if (a1 === alt) return ref;
    throw new Error(`A1 allele ${a1} matches neither REF (${ref}) nor ALT (${alt})`);
  }
  if (ref !== null) { if (a1 === ref) throw new Error(`A1 allele ${a1} matches REF (${ref}); cannot determine A2 without an ALT column`);
    return ref; }
  if (alt !== null) { if (a1 === alt) throw new Error(`A1 allele ${a1} matches ALT (${alt}); cannot determine A2 without a REF column`);
    return alt; }
  throw new Error("No REF, ALT, or A2 column found to determine A2");
}
/** os.path.splitext(basename without .gz)[0] + ".ldsc.tsv" */
function normalizedName(filename) {
  let f = String(filename).split("/").at(-1);
  if (f.endsWith(".gz")) f = f.slice(0, -3);
  const dot = f.lastIndexOf("."), lead = f.match(/^\.*/)[0].length; // splitext ignores leading dots
  return `${dot > lead - 1 && dot >= lead ? f.slice(0, dot) : f}.ldsc.tsv`;
}

/**
 * normalize_sumstats_for_ldsc on an uploaded file's text → { valid, errors, warnings, detected_format,
 * normalized_filename, text (what the analysis reads: the rewritten file, or the upload) }
 */
export function normalizeSumstats(text, filename, selectedFormat) {
  const out = { valid: true, errors: [], warnings: [], detected_format: null, normalized_filename: filename, text };
  const fail = (...errors) => ({ ...out, valid: false, errors: [...out.errors, ...errors] });
  const { columns, rows } = read(text);
  if (!columns.length) return fail("File is empty or missing a header row.");
  const cols = byClean(columns);
  let fmt, mapped;
  if (selectedFormat) {
    fmt = FORMATS.find(f => f.name === SELECTED[selectedFormat]);
    if (!fmt) return fail("Unsupported summary statistics format selection. Select PLINK raw, REGENIE raw, SAIGE raw, or Pre-munged.");
    const m = matches(fmt, cols);
    if (!m.ok) {
      out.detected_format = fmt.name;
      const msg = `Uploaded file does not match the selected ${fmt.name} format. Select the correct format or upload a pre-munged`;
      return fail(`${msg} LDSC-ready file.`, ...m.missing);
    }
    mapped = m.mapped;
  } else {
    const tried = [];
    for (const f of FORMATS) { const m = matches(f, cols); if (m.ok) { fmt = f; mapped = m.mapped; break;
      } tried.push(`${f.name}: ${m.missing.join("; ")}`); }
    if (!fmt) return fail("Could not detect a supported summary statistics format.", ...tried);
  }
  out.detected_format = fmt.name;
  if (fmt.name === "LDSC-ready") return out;
  // the output mapping: per canonical column, the first source column present (canonical_columns order)
  const source = {};
  for (const [src, canon] of fmt.canonical) if (cols[src] && !(canon in source)) source[canon] = cols[src];
  const plink = fmt.name === "PLINK", refCol = plink ? cols.REF : undefined, altCol = plink ? cols.ALT : undefined,
    a2Col = plink ? cols.A2 : undefined;
  let body = csvRow(OUT), n = 0;
  const errors = [];
  rows.forEach((row, k) => {
    try {
      const o = {};
      for (const c of OUT) {
        if (plink && c === "A2") { o.A2 = plinkA2(row, o.A1, refCol, altCol, a2Col); continue; }
        if (!source[c]) throw new Error(`No source column mapped for ${c}`);
        o[c] = transform(row[source[c]] ?? "", fmt.transforms[clean(source[c])], source[c]);
      }
      o.P = probability(o.P, "P");
      body += csvRow(OUT.map(c => o[c])); n++;
    } catch (e) { if (errors.length < 10) errors.push(`Row ${k + 2}: ${e.message}`); }
  });
  if (errors.length) return fail(...errors);
  if (!n) return fail("No data rows were found in the uploaded summary statistics file.");
  const warnings = [];
  if (fmt.name === "REGENIE" && "LOG10P" in cols) warnings.push("LOG10P was converted to P for LDSC munging.");
  if (Object.keys(mapped).some(c => clean(c) === "OR")) warnings.push("OR was converted to BETA using log(OR) for LDSC munging.");
  return { ...out, warnings, normalized_filename: normalizedName(filename), text: body };
}

/**
 * An upload as LDlink's web client analyses it: saved under secure_filename, normalized (validate_sumstats), then the
 * analysis runs on the normalized file → { filename, text } or { failed: "Summary statistics validation failed for <name>.",
 * errors } (validate_sumstats_preanalysis; the validation's own errors are what the upload form shows)
 */
export function prepareUpload(text, filename, selectedFormat) {
  const name = secureFilename(filename), n = normalizeSumstats(text, name, selectedFormat);
  if (!n.valid) return { failed: `Summary statistics validation failed for ${name}.`, errors: n.errors };
  return { filename: n.normalized_filename, text: n.text, warnings: n.warnings };
}

/** werkzeug's secure_filename (posix): the name LDlink saves an upload under */
export function secureFilename(name) {
  let s = String(name).normalize("NFKD").replace(/[^\x00-\x7f]/g, "").replaceAll("/", " ");
  s = s.split(/[\t-\r\x1c-\x20]+/).filter(Boolean).join("_").replace(/[^A-Za-z0-9_.-]/g, "");
  return s.replace(/^[._]+|[._]+$/g, "");
}
