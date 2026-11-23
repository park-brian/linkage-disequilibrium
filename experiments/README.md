# Parquet experiments

Size/speed experiments behind [docs/PARQUET_LAYOUT.md](../docs/PARQUET_LAYOUT.md). Run `npm install` here first. Large outputs go to `out/` (git-ignored); summaries are in `results/`.

| script | purpose |
|---|---|
| `build_genotypes.mjs <vcf.gz> <outDir> <panelDir>` | one pass over a 1000G VCF → candidate genotype layouts (codecs, page sizes, bitset vs sparse, MAF tiers) |
| `bench_genotypes.mjs <outDir> <panelDir>` | region scans, SNP lists, banded LD throughput per population |
| `build_dbsnp.mjs <outDir>` | rsid-sorted and composite-key-sorted dbSNP index candidates (needs `out/dbsnp_src/*` exported by DuckDB from `../../sources/dbsnp/dbsnp.parquet`) |
| `bench_dbsnp.mjs <outDir>`, `bench_dbsnp_cached.mjs` | batched rsID lookups and position-window scans |
| `lookup.mjs` | page-index batched lookup on a sorted column (row-group stats → ColumnIndex/OffsetIndex → only matching pages) |
| `build_genotypes_v2.mjs <vcf.gz> <outDir> <panelDir> <dbsnpPos.parquet> <chromCode>` | v2 build: sub-population MAF tiers, rsID labels from the 2023 dump (VCF rsID, else position + type), dbSNP function, parity flags; one file per tier × page size |
| `bench_v2.mjs <outDir>` | v2 window scans, scattered point reads (payload-page lookup), per-page sizes |
| `dbsnp_bson2tsv.mjs` | stream LDlink's `dbsnp.bson.gz` mongodump (stdin) → TSV keeping dump order (`ord`) |
| `dbsnp_sort.sql` | DuckDB external sort into intermediates (rsID partitions; per-build per-chromosome position files), anomalies, summary |
| `build_dbsnp_v2.mjs <rsid\|pos37\|pos38> <releaseRoot>` | final dbSNP indexes with hyparquet-writer; checks (key, dump order) and row counts |
| `build_manifest.mjs <releaseRoot>` | `manifest.json`: footer ranges, row-group key ranges, page-index locations per file |
| `resolve.mjs` | LDlink-parity resolver: `get_coords` (first doc in dump order) and `replace_coord_rsid` |
| `parity_expected.sql` + `run_sql_statements.py` | expected answers from the dump with Mongo semantics (one DuckDB process per statement) |
| `verify_parity.mjs <releaseRoot>` | resolver vs expected answers and reference test vectors → `parity_report.json` |
| `ld_check.mjs <v2 dir> <vcf.gz> <tbi> <vcfChrom> <panel> <start>` | end-to-end check: r²/D′ from Parquet bitsets (incl. haploid-male mask) vs reference LDlink math on the VCF |
| `wasm/ld_kernel.c`, `wasm/build.sh` | LD kernel in C → standalone WASM (no libc): `ld_kernel_scalar.wasm` (i64.popcnt), `ld_kernel_simd.wasm` (i8x16.popcnt); band, blocked band (4-row register tiling, row ranges for workers), masked band, tiles, pair counts; r² bound pruning |
| `wasm/ld_kernel.mjs` | JS wrapper: packs population haplotypes from v2 Parquet rows (bitset or sparse, missing, haploid males) into WASM memory; `band`, `tile`, `pairCounts` |
| `wasm/bench_kernel.mjs`, `wasm/bench_blocked.mjs`, `wasm/bench_threads.mjs` | exactness vs an independent JS reference, pruning losslessness, blocked = unblocked, speed vs JS, worker-thread scaling |
| `build_dup_rsids.mjs <releaseRoot>` | `dbsnp/duplicated_rsids.parquet`: the 296,308 rsIDs with more than one doc (always resolved directly for parity) |
| `planner.mjs` | query planner: LDlink-parity resolution (direct or hybrid label scan for rsIDs, dbSNP position index for coords), then per-cluster scan vs page reads by manifest-estimated bytes |
| `verify_planner.mjs <releaseRoot>` | hybrid == direct resolution, predicted vs measured bytes per file group, and the alternative strategy measured |
| `presign.sh <releaseRoot> <s3 prefix> <profile>` | presigned GET URLs (12 h) for every release file → `out/presigned.json` (git-ignored; contains signatures) |
| `net_bench.mjs <presigned.json> <releaseRoot> <panelDir>` | planner scenarios over HTTPS Range requests (cold), optional `MAX_INFLIGHT` cap to model browser connection limits |
| `bson_stream.mjs F...` / `json_array_stream.mjs F...` | stream a mongodump `.bson.gz` / a mongoexport JSON array (stdin) → TSV with dump order + chosen fields |
| `fetch_annotations.sh` | fetch LDlink's annotation collections + the 2025 GWAS backup from S3 to `out/annot_src/*.tsv` |
| `annotations.sql` | key and sort annotations per build (DuckDB intermediates), counts of keyed vs. not keyable docs |
| `build_annotations.mjs <releaseRoot>` | write `annotations/` Parquet files with hyparquet-writer (order and row-count checks) |
| `parity_annotations.sql` + `verify_annotations.mjs <releaseRoot>` | expected LDlink answers from the raw dump tables vs. the Parquet files |
| `bench_annotations.mjs <releaseRoot>` | 1/10 Mb window cost per annotation track (uses `rangeScan` page-level range reads) |
| `ldsc/convert_ldscore_ref.mjs` | LDlink's `/data/ldscore/<group>/*.l2.ldscore.gz` + M files → `ldscore/<group>.parquet` (+ `w_hm3.snplist.parquet`) |
| `ldsc/ldscore_l2.mjs <prefix> <cm\|kb> <window> <out>` | ldsc `--l2` port for PLINK bed/bim/fam (byte-identical to LDlink's example) |
| `ldsc/munge.mjs` | `munge_sumstats.py` port (LDlink flags: `--merge-alleles w_hm3.snplist`) |
| `ldsc/ldsc_regress.mjs` | h² / rg LD score regression port (IRWLS, two-step, block jackknife) |
| `ldsc/ref_dump.py` | run the reference ldsc (venv) and dump full-precision h²/rg JSON |
| `ldsc/verify_regress.mjs`, `ldsc/verify_munge.mjs` | JS vs. reference: regressions on Python-munged files; JS munge output + full JS pipeline |
| `fetch_ncbi_dbsnp.sh <outDir>` | download NCBI dbSNP b155 + b157 VCFs (+tbi) and `refsnp-merged.json.bz2` (resumable; md5 where NCBI provides it) |
| `ncbi_vcf2tsv.mjs <vcf.gz> <out.tsv>` | NCBI VCF → TSV (chrom code, pos, rsid, ref, alt, VC, function flags), ~35 min per file; then DuckDB `COPY` to Parquet |
| `ncbi_merged2tsv.py <merged.json.bz2> <out.tsv>` | merge history → (retired rsID, target) pairs |
| `compare_ldlink_parse.sql` | LDlink's own `rsjson_parse_data_mongo.py` output (b155 chr22 JSON) vs. the dump's chr22 docs |
| `compare_ncbi.sql` (`sed s/{B}/b155/`) | genome-wide NCBI VCF vs. dump: rsIDs, chromosome, LDlink-rule positions, type, function (as sets), merges |
| `allele_join.mjs <ncbi_chr.parquet> <genotypes/chrN> <out.json>` | allele-aware rsID ↔ 1000G join: release labels and LDlink's genotype-row choice vs. allele match |
| `build_genotypes_v2.mjs` with `ALLELES=<ncbi slice.parquet>` | adds `vcf_order` (always) and `allele_is_label` / `allele_rsids` (NCBI allele matches) to the genotype files |
| `select_rows.mjs` | genotype-row choice per variant: `pickParity` (LDlink: last overlapping row by `vcf_order`), `pickAllele` (opt-in), `candidateRows` (overlap windows via `rangeScan`) |
| `verify_select.mjs <vcf.gz> <genoDir> <ncbi slice> <out.json> [n]` | both modes vs. independent ground truth (direct VCF pass; recomputed allele matches) |
| `rsid_label.mjs` | rsID labels for genotype rows: VCF rsID, else NCBI allele match, else position + type, searched at the dump position LDlink uses (POS+1 for deletions); `rsid_pos_delta` |
| `test/phaseA.test.mjs` | regression tests for the ROADMAP Phase A fixes (`node --test test/phaseA.test.mjs`); fixtures from reference munge_sumstats.py and numpy in `test/` |
| `kernel_check.mjs <genoPrefix> <vcf.gz> <tbi> <vcfChrom> <start> <out.json>` | the WASM kernel itself (all paths, scalar + SIMD) vs LDlink's LD math on the VCF, incl. prefix-shaped populations, sparse tier, chrX masks |
| `parity_inputs.py <inputs.json> <build> <out.json>` | LDlink's input chain (validsnp → replace_coord_rsid → get_coords) in LDlink's own Python expressions, against the dump; ground truth for `verify_parity.mjs` part 4 |
| `pyfmt.mjs` | Python-exact `'%.nf'` / `round(x, n)` (exact binary value, ties to even, signed zero) and `nextUp`; used by `pyFixed3` and the kernel API |
| `wasm/ld_kernel.mjs` | the LD kernel API: `loadKernel(bytes)` (WASM, JS fallback) / `jsKernel()`; prepare, band, tile, pairCounts, exact, reset |
| `wasm/ld_kernel.c` → `ld_kernel.wasm` (`wasm/build.sh`), `wasm/ld_kernel_js.mjs` | the two backends: scalar WASM and plain JS, same algorithm and identical results |
| `wasm/bench_js_vs_wasm.mjs <common.parquet> <panelDir> <lo> <hi> [maxDist] [R] [out.json]`, `wasm/bench_ldmatrix_size.mjs` | WASM vs JS backend through the same API: outputs checked identical, Mpairs/s; LDmatrix-sized requests (docs/CORE_REVIEW.md §2, §5) |
| `build_all.sh [P] [build ...]`, `phaseB_slices.sql`, `check_release.mjs <root>` | Phase B: per-chromosome inputs, all-chromosome genotype build (resumable), release integrity check |
