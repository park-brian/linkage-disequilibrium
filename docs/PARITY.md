# Parity with LDlink: what the ports reproduce

Every module in `src/ldlink/` returns LDlink's own output, byte for byte, including the behaviours that are bugs.

**How this was established:** each module was run against a local LDlink (`tools/ldlink-local/`), built from NCI's
source (`../nci-webtools-dceg-linkage/server`), on all three builds (`npm run verify:ldlink`; current totals in
[STATUS.md](STATUS.md)). The investigation and the per-case reports are in `docs/history/CORE_REVIEW.md` §6–10.

**Changing a rule below:** find it in the LDlink source (`file:line` habit), change the port, and re-run that
module's comparison.

## Shared rules (all modules)

**Input chain.** `validsnp`, `replace_coord_rsid` and `get_coords`, ported as `dbsnp.resolveInputs`:
- **Tokens:** the first whitespace token is lower-cased.
- **Coordinates → rsID:** the first doc at the position with `id == ref_id` (a "Multiple rsIDs" warning when
  several qualify, including LDlink's empty-list text), else the first doc.
- **rsID → position:** `find_one({id})`, the first doc in dump order. A position of "NA" means not in dbSNP.

**Genotype row.** tabix `chr:p-p` returns every row overlapping p; LDlink keeps the **last**, which is the highest
`vcf_order` (`genotypes.pickParity`):
- **Multi-allelic:** a multi-allelic source record → "not a biallelic".
- **SVs:** lines containing "END" are dropped.
- **Indels** whose dbSNP position is POS + 1 still resolve to the overlapping row.

**Statistics:**
- **Pair LD:** the 2×2 haplotype table, with D′ by the sign-dependent Dmax and r² = Δ²/Ms; NA when Ms = 0.
- **Haploid males** contribute one haplotype; a haplotype uncalled in either variant is dropped from the pair.
- **Thresholds** compare unrounded doubles.
- **Formatting:** Python's own: `round`, `str(float)`, `'%.nf'` (`pyfmt.mjs`); χ² uses an exact integer numerator.

**Raw VCF IDs.** LDlink prints the VCF's ID in "does not match RS number" warnings and in LDproxy's RS_Number.
Stored as `vcf_id` (GRCh37).

## Per module

| module | output compared | rules worth knowing |
|---|---|---|
| **LDmatrix** | r² and D′ text matrices | input gate (len > 2, starts rs / chr, ends in a digit); duplicates dropped; NA position, chrY on GRCh38, not in 1000G, multi-allelic → dropped; position order; 3 decimals. Only the GT strings `0\|0 0\|1 1\|0 1\|1 0 1` count |
| **LDpair** | API JSON | A..D follow the *sorted allele strings*; p = `2(1 − ½(1 + erf(√χ²/√2)))`, `"<0.0001"`; only the error is returned for a failed pair; the biallelic error appends the *warning* text. **Nondeterministic in LDlink:** tied haplotype counts follow a hash-randomized sample order (per mod_wsgi process); the comparison reports these as `tie_order` |
| **LDpop** | API text | a haploid call counts as one REF and one ALT allele over two; SNP2's frequencies divide by SNP1's allele total; rows in `pop_order`, super-populations expanded |
| **LDhap** | `calculate_hap` JSON + API text | tokens deduped *before* coordinate replacement; same-chromosome error raised after `parse_vcf`; a haplotype with any haploid / uncalled allele dropped whole; count descending, ties by string |
| **LDproxy** | API text + JSON warning | window cut into `NUM_SUBPROCESSES` (4) chunks (`chunkWindow`), each a tabix overlap read: boundary-spanning variants listed twice, the last chunk ends 3 bp past the window; stable sort \|distance\| then statistic; the first line takes RS / statistics from the top row but Coord / Alleles from the query; proxies need R² > 0.01; RS_Number = VCF rsID else the first dbSNP doc at the VCF POS; MAF over haplotypes called in both; the REST layer validates `r2_d` and `window` first |
| **SNPclip** | API text | route keeps only rs / chr lines; MAF from genotype strings; greedy pruning in input order (r² ≥ threshold); an individual diploid at one variant and haploid at the other is skipped for that pair; `maf ≥ 0` with a monomorphic survivor → LDlink's internal error (compares `None ≥ t`) |
| **SNPchip** | API text | sorted by (chromosome number, position), `22:123`; warnings start with ". "; platform match by name or `name.rstrip(' Array')` (strips *characters*); an empty platform request crashes LDlink |
| **LDtrait** | API text | GWAS records in index order (position, then insertion); pair LD from one tabix row per position: NA when multi-allelic or POS ≠ the queried position; R2 / D′ unrounded; columns shifted against their headers ("Risk Allele" = RAF, "Effect Size" = float(OR or BETA), "Beta or OR" = the CI with brackets, letters and spaces removed) |
| **LDexpress** | API text | LDproxy's chunked window; eQTLs inner-joined to tissues and genes (`$lookup` / `$unwind`), p < threshold, `%.6g`; a query not in 1000G deletes the request's population file, so every later query returns nothing; a multi-allelic query is computed on its GT digits; blank entries count toward the 10-entry limit and are then dropped |
| **LDassoc** | `calculate_assoc` text + JSON error / warning (no API) | the origin genotype is the *first* overlapping row (the checks use the *last*); one tabix region per association position (a row is listed under every position it overlaps, and kept if its own POS is in the file); sort \|distance\| then p; correlated alleles proxy-first, indexes crossed; region mode compares coordinate *strings*; a coordinate origin raises KeyError; High Coverage lowest-p search asks tabix for `chrchr8:…` and always fails |

## Deliberate differences

- **The REST shells:** the ports implement module behaviour, plus the REST layer's parameter checks (`Invalid r2_d
  parameter.`, `Invalid window parameter.`). They do not implement tokens, rate limits or the web-only routes, such
  as LDtrait's "too many pairs" prompt.
- **Allele-aware selection:** `select: "allele"` is opt-in. Instead of LDlink's overlapping row, it picks the row
  whose NCBI alleles match the rsID; the default is parity.
