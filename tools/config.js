// The one place that knows where things are on this machine. Everything else asks here; environment variables
// override every default. Paths are absolute.
//
//   LD_DATA        generated data root (release/, build intermediates, verify reports)   default <repo>/../data
//   LD_RELEASE     the release the app, tests and verifiers read                         default $LD_DATA/release_v4
//   LD_PANELS      1000G population sample lists (<POP>.txt)                              default <repo>/../db/data/1000G_population_samples
//   LD_VCF         1000G VCFs (GRCh37/, GRCh38/, GRCh38_High_Coverage/)                  default <repo>/../db/1000G
//   LDLINK_URL     a local LDlink (the parity oracle)                                    default http://127.0.0.1:8090
//   LDLINK_HOME    that LDlink's working directory (data/tmp exchange, LDassoc examples)  default <repo>/../ldlink-local
//   LD_REF_NODE_MODULES  node_modules with @gmod/tabix for verify/reference/ldmatrix_reference.js
//                                                                                        default <repo>/../convert-vcf-json/node_modules
//
// Pipeline and verify scripts run with the working directory at LD_DATA (`useDataCwd()`), so data paths in JS, SQL
// and Python are written relative to it; repo files are found from import.meta.url.
import { fileURLToPath } from "node:url";
import { resolve } from "node:path";
import { mkdirSync } from "node:fs";

export const REPO = resolve(fileURLToPath(new URL("..", import.meta.url)));
const env = (name, dflt) => resolve(process.env[name] || dflt);
export const DATA = env("LD_DATA", resolve(REPO, "../data"));
export const RELEASE = env("LD_RELEASE", resolve(DATA, "release_v4"));
export const PANELS = env("LD_PANELS", resolve(REPO, "../db/data/1000G_population_samples"));
export const VCF = env("LD_VCF", resolve(REPO, "../db/1000G"));
export const LDLINK_URL = process.env.LDLINK_URL || "http://127.0.0.1:8090";
export const LDLINK_HOME = env("LDLINK_HOME", resolve(REPO, "../ldlink-local"));
export const LDLINK_TMP = resolve(LDLINK_HOME, "data/tmp");
export const LDASSOC_EXAMPLES = resolve(LDLINK_HOME, "data/ldassoc_example");
export const REF_NODE_MODULES = env("LD_REF_NODE_MODULES", resolve(REPO, "../convert-vcf-json/node_modules"));

/** a repo file as an absolute path: repo("verify/data/parity_inputs.py") */
export const repo = rel => resolve(REPO, rel);
/** a data file as an absolute path: data("dbsnp_dump/work.duckdb") */
export const data = rel => resolve(DATA, rel);
/** run with the working directory at LD_DATA (data-relative paths in scripts, SQL and arguments) */
export function useDataCwd() { mkdirSync(DATA, { recursive: true }); process.chdir(DATA); return DATA; }
