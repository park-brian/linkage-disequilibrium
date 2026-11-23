// The data an LDlink module port needs, composed from the core: one planner (variant resolution + its file cache),
// genotype rows, dbSNP find_one lookups, annotation tracks, population panels and kernel preparation.
// Browser-ready apart from `open`, which the caller supplies (node: local files; browser: Range requests). Population
// panels come from the release (annotations/panels.parquet) unless the caller passes `readPanel`.
import { openPlanner } from "../planner.mjs";
import { createGenotypes } from "../genotypes.mjs";
import { createAnnotations } from "../annotations.mjs";

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
    if (!panels.has(pop)) panels.set(pop, (async () => new Set((await Promise.all(pop.split("+").map(readPanel))).flat()))());
    return panels.get(pop);
  }
  /** prepare rows (possibly from several chromosomes, same sample order) for the kernel over a population */
  async function prepare(build, rows, pop) {
    const chroms = [...new Set(rows.map(r => r.chrom))], ss = await Promise.all(chroms.map(c => genotypes.samples(build, c)));
    if (ss.some(s => s.ids.join() !== ss[0].ids.join())) throw new Error("sample order differs between chromosomes");
    const ids = await popIds(pop), haploidMales = new Set(ss.flatMap(s => [...s.haploidMales]));
    const hapIndexes = ss[0].ids.flatMap((s, i) => (ids.has(s) ? [2 * i, 2 * i + 1] : []));
    return kernel.prepare(rows, { samples: ss[0].ids.length, haploidMales, hapIndexes });
  }
  return { root, planner, kernel, ...genotypes, ...annotations, popIds, prepare,
    dbsnpAt: planner.resolver.findAt, dbsnpById: planner.resolver.findById, dbsnpAllAt: planner.resolver.findAllAt, duplicatedRsids: planner.resolver.duplicated };
}
