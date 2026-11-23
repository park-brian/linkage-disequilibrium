// The data an LDlink module port needs, composed from the core: one planner (variant resolution + its file cache),
// genotype rows, dbSNP find_one lookups, annotation tracks, population panels and kernel preparation.
// Browser-ready apart from `open`, which the caller supplies (node: local files; browser: Range requests). Population
// panels come from the release (annotations/panels.parquet) unless the caller passes `readPanel`.
import { openPlanner } from "../planner.js";
import { createGenotypes } from "../genotypes.js";
import { createAnnotations } from "../annotations.js";
import { ldscPanel, ldscCounts, mergeAlleles } from "../ldsc/panels.js";
import { memo } from "../source.js";

/**
 * @param {{ root: string, open: (path: string, byteLength?: number) => Promise<AsyncBuffer>, readPanel?: (pop: string) => Promise<string[]>,
 *   kernel, manifest?: object }} o manifest: the parsed manifest.json when the caller has it (else read through open)
 */
export async function createContext({ root, open, readPanel, kernel, manifest }) {
  const planner = await openPlanner(root, open, { manifest }), panels = new Map();
  const genotypes = createGenotypes(planner.file), annotations = createAnnotations(planner.file);
  readPanel ??= async pop => (await annotations.panels()).get(pop) ?? []; // an unknown population has no samples
  /** population sample IDs: the union of the panels in "CEU+YRI" (LDlink: list(set(ids))) */
  async function popIds(pop) {
    return memo(panels, pop, async () => new Set((await Promise.all(pop.split("+").map(readPanel))).flat()));
  }
  /** the kernel options for rows (possibly from several chromosomes, same sample order) over a population */
  async function kernelOptions(build, rows, pop) {
    const chroms = [...new Set(rows.map(r => r.chrom))], ss = await Promise.all(chroms.map(c => genotypes.samples(build, c)));
    if (ss.some(s => s.ids.join() !== ss[0].ids.join())) throw new Error("sample order differs between chromosomes");
    const ids = await popIds(pop), haploidMales = new Set(ss.flatMap(s => [...s.haploidMales]));
    const hapIndexes = ss[0].ids.flatMap((s, i) => (ids.has(s) ? [2 * i, 2 * i + 1] : []));
    return { samples: ss[0].ids.length, haploidMales, hapIndexes };
  }
  /**
   * rows prepared for the kernel over a population → use(p, opts), whose result is returned; the kernel is reset after.
   * use must be synchronous: the WASM kernel's memory is one arena that reset() frees, so nothing else may prepare
   * between this prepare and its reset (two runs at once would reuse each other's buffers). Reading the samples and
   * the panel (the async part) comes first; prepare, use and reset then run without yielding.
   */
  async function withPrepared(build, rows, pop, use) {
    const opts = await kernelOptions(build, rows, pop), p = kernel.prepare(rows, opts);
    try { return use(p, opts); } finally { kernel.reset(); }
  }
  // LD score panels and the HapMap3 list for the LDscore family, read once each (src/ldsc/panels.js)
  const ldsc = new Map();
  return { root, planner, kernel, ldscPanel: (pop, need) => ldscPanel(planner.file, pop, need), ldscCounts: () => memo(ldsc, "counts", () => ldscCounts(planner.file)),
    mergeAlleles: () => memo(ldsc, "w_hm3", () => mergeAlleles(planner.file)), ...genotypes, ...annotations, popIds, withPrepared,
    dbsnpAt: planner.resolver.findAt, dbsnpById: planner.resolver.findById, dbsnpAllAt: planner.resolver.findAllAt, duplicatedRsids: planner.resolver.duplicated };
}
