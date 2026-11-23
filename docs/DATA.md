# LDlink source data inventory

What the reference LDlink (`../nci-webtools-dceg-linkage`) reads, where it lives, what it looks like, and how we
intend to repackage it as static Parquet for a browser-only, hyparquet-based tool.

Profiling outputs referenced below are in [`scripts/profile/results/`](../scripts/profile/results/); the scripts
that produced them are in [`scripts/profile/`](../scripts/profile/). All profiling was read-only.

> Snapshot date: 2026-11-19. AWS account 201624419975 (us-east-1), the `nci-cbiit-dceg-data-nonprod` bucket.

## 1. Where the data lives

| What | Location | Notes |
|---|---|---|
| App config | SSM `/analysistools/{dev,qa}/ldlink/*` | `s3_bucket=nci-cbiit-dceg-data-nonprod`, `s3_bucket_data_dir=ldlink`, `genotypes_dir=1000G/`, `population_samples_dir=1000G_population_samples/`, `dbsnp_version=155`, `mongodb_database=LDLink`, `ldtrait_src=https://www.ebi.ac.uk/gwas/api/search/downloads/alternative` |
| 1000G genotypes | `s3://nci-cbiit-dceg-data-nonprod/ldlink/1000G/{GRCh37,GRCh38,GRCh38_High_Coverage}/*.vcf.gz` + `.tbi` (71 indexes) | 17.0 / 14.3 / 37.4 GB. Read by `tabix` straight from S3 (`server/LDcommon.py:147-200`). Identical byte sizes to the local copy in `../db/1000G/` |
| Population panels | `s3://…/ldlink/deploy/data/1000G_population_samples/<POP>.txt` (also `../db/data/…`) | one sample ID per line; copied to EFS `/data/` at deploy |
| LDassoc example | `s3://…/ldlink/deploy/data/ldassoc_example/prostate_example_grch3{7,8}.txt` | 16 MB each |
| MongoDB `LDLink` | **live hosts are gone**; mongodump in `s3://…/ldlink/deploy/mongodb/*.bson.gz` (2023-04-07) | see §3 |
| Newer GWAS Catalog | `s3://…/ldlink/backups/LDLink.gwas_catalog20251024.json` (1.3 GB, mongoexport JSON array) | 835,632 docs, latest `DATE ADDED TO CATALOG` 2025-08-31 |

### Live MongoDB status (probe, 2026-11-19)
A temporary Fargate task (`dev-ldlink-api` task definition, ldlink subnets/SG, ECS Exec) was started, used, and
stopped (`STOPPED: ldlink mongo probe done`; no tasks left running with `--started-by ldlink-mongo-probe`).
From inside the VPC, DNS works (`ldlink-dev.nih.gov → 10.208.2.161`) but **both configured Mongo hosts fail to
resolve**: `ncias-d2869-c.nci.nih.gov` (dev, SSM last modified 2023-04-20) and `ncias-q2882-c.nci.nih.gov` (qa).
There is no DocumentDB cluster in the account, and no ldlink ECS service is currently running (task definitions
only). **The S3 mongodump is therefore the authoritative copy of the annotation data**, plus the 2025 GWAS backup.

## 2. 1000 Genomes genotypes

File templates (from `server/LDcommon.py:24-65`; `%s` = chromosome):

| build key | dir | file | contig names |
|---|---|---|---|
| `grch37` | `GRCh37` | `ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz` | `1`…`22`,`X`,`Y` |
| `grch38` | `GRCh38` | `ALL.chr%s.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz` | `1`…`22`,`X` |
| `grch38_high_coverage` | `GRCh38_High_Coverage` | `CCDG_14151_B01_GRM_WGS_2020-08-05_chr%s.filtered.shapeit2-duohmm-phased.vcf.gz` | `chr1`…`chr22`,`chrX` |

(The old browser prototype `sources.js` points at `…v5b…` on the EBI FTP; LDlink uses `…v5…`.)

### Profile (chr22 full; chrX/Y partial) — `results/g37_22.json`, `g38_22.json`, `hc_22.json`, `*_X.json`, `g37_Y.json`

| | GRCh37 | GRCh38 | GRCh38 High Coverage |
|---|---|---|---|
| samples in VCF | 2504 | 2548 | 3202 |
| vs `ALL.txt` (2504) | identical, same order | 2503 present (**NA18498 missing**), +45 extra, different order | all 2504 present (+698 related), different order |
| chr22 records | 1,103,547 | 1,059,079 | 1,070,401 |
| SNV / indel+MNP | 1,055,454 / 41,031 | 993,880 / 65,199 | 929,834 / 140,567 |
| multi-allelic records (`,` in ALT) | 6,348 (GT values up to `4`) | 0 (split into biallelic rows) | 0 (split) |
| symbolic SV (`<…>`) | 714 (655 with `END=`) | 0 | 0 |
| adjacent duplicate POS | 1,166 | 3,260 | **61,774** |
| ID column | 47% `rs…`, 53% `.` | **all `.`** | `chr:pos:ref:alt` (e.g. `22:10519265:CA:C`) on chr22; `.` on chrX |
| FORMAT | `GT` (chrX adds `GL:GP`; chrY some CNV fields) | `GT` | `GT` |
| FILTER | `PASS` (chrX PAR region has `.`) | `PASS` | `PASS` |
| chrX non-PAR | males haploid (`0`/`1`) | **only 106,963 records on all of chrX, all diploid** — looks like a partial/odd release; treat GRCh38 chrX as suspect | males haploid |
| chrY | 1,233 males, haploid, `.` missing calls | n/a (LDlink rejects chrY for GRCh38 builds) | n/a |

Implications:
- **rsIDs never come from the GRCh38 VCFs**; LDlink always resolves rsID → (chr, pos) via dbSNP and then matches
  by position only (`parse_vcf`, `LDcommon.py:577-634`). With duplicate positions it keeps the record whose POS equals
  the query, else the *last* one seen — at a split multi-allelic site that is an arbitrary allele. The new tool should
  match on position **and** alleles where it can, and surface ambiguity.
- High Coverage genotypes agree with phase 3 for ~99.9% of 2504 shared samples at the test SNPs, with ~10–20% of
  samples phase-swapped as whole haplotypes (`scripts/profile/gt_compare_builds.mjs`); 2-locus counts were identical for
  the CEU test vectors. LDlink uses only the 2504 panel samples for High Coverage, too.
- Sample order differs per build → population masks must be built per build from the VCF header.

### Population panels (`1000G_population_samples/*.txt`)
ALL=2504; super-populations AFR=661, AMR=347, EAS=504, EUR=503, SAS=489; 26 populations:
ACB 96, ASW 61, BEB 86, CDX 93, CEU 99, CHB 103, CHS 105, CLM 94, ESN 99, FIN 99, GBR 91, GIH 103, GWD 113,
IBS 107, ITU 102, JPT 104, KHV 99, LWK 99, MSL 85, MXL 64, PEL 85, PJL 96, PUR 104, STU 102, TSI 107, YRI 108.
LDpop additionally has a hard-coded lat/long + display-name table (`server/LDpop.py:343-518`).

## 3. MongoDB `LDLink` collections (from the 2023-04 dump)

Counts are full scans of the dump (`results/m_*.json`) unless marked.

| collection | docs | key fields (types as stored) | indexes | used by |
|---|---|---|---|---|
| `dbsnp` (b155) | ~1.07 B (from prior full Parquet export, 1,073,486,098) | `id` str (digits, no `rs`), `ref_id` str, `chromosome` str (`1`…`22`,`X`,`Y`), `position_grch37`/`position_grch38` **str** or `"NA"`, `function` str (comma list / `"NA"`), `type` str | `id`; `(chromosome, position_grch37)`; `(chromosome, position_grch38)` | every module (rsID↔coord, LDproxy `function`) |
| `refGene_grch37` / `_grch38` | 78,106 / 78,363 | UCSC refGene: `name`, `chrom` (`chr1`…, incl. `chrM`), `strand`, `txStart/txEnd/cdsStart/cdsEnd` int, `exonCount`, `exonStarts/exonEnds` comma-terminated str, `name2` (symbol) | `(chrom, txStart, txEnd)` | gene tracks |
| `recomb` | 3,283,362 | `rate` **str** (cM/Mb, 0–107.8), `chromosome_grch37/38` str (no `chr`, 1–22 + X), `position_grch37` int, `position_grch38` int or `"NA"` (766) | chr+pos per build | LDproxy / LDassoc recombination line |
| `regulome` | 16,357,032 | `score` str categorical (`1a`…`1f`,`2a`…`2c`,`3a`,`3b`,`4`…`7`), `chromosome_grch37/38` **with `chr` prefix**, `position_grch37` int, `position_grch38` int or `"NA"` (22,050) | chr+pos per build | LDproxy / LDassoc |
| `forge_score` | 37,962,826 | `snp_id` str (`rs…`), `score` str `0`–`10` | `snp_id` | LDproxy / LDassoc / LDmatrix diagonal |
| `genes_name_coords` | 26,195 | `name` (upper-case symbol), `chromosome_grch37/38`, `begin_/end_grch37/38` int (`"NA"` for 279 in 38) | `name` | LDassoc gene mode |
| `gwas_catalog` | 462,053 (2025 backup: 835,632) | 38 GWAS Catalog "alternative" TSV columns verbatim (all str, e.g. `DISEASE/TRAIT`, `PUBMEDID`, `SNPS`, `RISK ALLELE FREQUENCY`, `OR or BETA`, `95% CI (TEXT)`, `P-VALUE`) + `chromosome` str, `position_grch37/38` int or `"NA"` | `(chromosome, position_grch37)`, `(chromosome, position_grch38)` | LDtrait |
| `gtex_tissue_eqtl` | 71,478,479 (523,233 lack a GRCh37 position) | GTEx v8 `signif_variant_gene_pairs` columns as **str** (`variant_id`, `gene_id`, `tss_distance`, `maf`, `pval_nominal`, `slope`, `slope_se`, …) + `chromosome_grch37/38`, `position_grch37/38` int, `tissueSiteDetailId` | chr+pos per build | LDexpress |
| `gtex_tissues` | 54 | GTEx v8 tissue info (`tissueSiteDetailId`, `tissueSiteDetail`, `colorHex`, sample counts…) | — | LDexpress tissue picker |
| `gtex_genes` | 56,200 | `gene_id_v26`, `gene_name_v26`, `ncbi_id` (`"NA"` for 17,693) | `gene_id_v26` | LDexpress `$lookup` |
| `snp_col` | 10,411,324 | `chromosome_grch37/38` str (incl. `0`, `XY`, `MT`), `position_grch37/38` int (0 / `"NA"` present), `data: [{platform}]` | chr+pos per build | SNPchip |
| `platforms` | 83 | `platform`, `code` | — | SNPchip |
| `api_users`, `api_log`, `ldscore_runs` | — | token users / usage / LDSC run state | — | **not needed** (server concerns) |

Data-quality quirks to normalize at Parquet build time:
- Numbers stored as strings (`dbsnp` positions, `recomb.rate`, every GTEx numeric, all GWAS fields) and `"NA"` as a
  null sentinel in otherwise-int columns.
- Chromosome spelling differs: `"7"` (dbsnp, recomb, gwas, snp_col, genes) vs `"chr7"` (regulome, refGene) vs
  `chr7` contigs in High Coverage VCFs.
- `gwas_catalog.position_grch37` max is 3,411,805,015 and a `dbsnp.position_grch37` of 69,716,399,231 exists (both impossible) → validate against chromosome lengths.
- The dbSNP dump is grouped by chromosome (first 3M docs are all chr11), not by id; ~97% of `function` values are `"NA"`.
- `snp_col` contains `0`/`XY`/`MT` chromosomes and position 0.
- `dbsnp.ref_id != id` marks merged rsIDs (e.g. rs52815562 → rs2179); LDlink prefers `id == ref_id` records.

## 4. Prior local experiments (context)

`../sources/` holds earlier Parquet work: dbSNP sorted by `id` with row-group
sweeps (best candidates: ~10–20k rows/row group with page index), a position-sorted copy
(`sources/dbsnp/layout/dbsnp.chrpos.parquet`), GRCh37-only haplotype Parquet storing each variant's GTs as one
base64(deflate) string, and `sources/dbsnp/dbsnp2.html` (hyparquet + page index + multipart Range batching).

## 5. Proposed Parquet layout (to validate in the next phase)

> Superseded by the measured design in [PARQUET_LAYOUT.md](PARQUET_LAYOUT.md) (MAF tiers, hybrid dbSNP lookup).

Built with hyparquet-writer from the sources above; read with hyparquet over HTTP Range. One file per
build × chromosome where the query is positional.

| dataset | key / sort | columns | notes |
|---|---|---|---|
| `dbsnp/rsid.parquet` | `rsid` (int) | `rsid`, `ref_rsid`, `chrom`, `pos37`, `pos38`, `type`, `function` | small row groups + page index; replaces Mongo `dbsnp` id lookup |
| `dbsnp/<build>/chr<N>.parquet` | `pos` | `pos`, `rsid`, `function` | coordinate → rsID (LDlink `get_rsnum`), and rsID labels for VCF rows lacking IDs |
| `genotypes/<build>/chr<N>.parquet` | `pos` | `pos`, `rsid` (joined from dbSNP), `ref`, `alt`, `ac`/`an`, `haps` (packed bitset, 1 bit per haplotype, sample order in file metadata), `missing` (bitset, only when non-empty: haploid/`.`) | one row per biallelic allele (multi-allelic split like GRCh38 does); population = precomputed bit mask; LD counts = AND + popcount |
| `annotations/<build>/genes.parquet` | `chrom`, `txStart` | refGene columns, exons as int lists | whole file is small (~6 MB as BSON) |
| `annotations/<build>/recomb/chr<N>.parquet` | `pos` | `pos`, `rate` float | |
| `annotations/<build>/regulome/chr<N>.parquet` | `pos` | `pos`, `score` (dictionary string) | |
| `annotations/forge.parquet` | `rsid` | `rsid`, `score` int8 | |
| `annotations/<build>/gwas/chr<N>.parquet` | `pos` | trimmed GWAS columns (trait, PMID, risk allele, effect, CI, p) | from the 2025 backup / EBI download |
| `annotations/<build>/gtex/chr<N>.parquet` | `pos` | `variant_id`, `gene_id`, `tissue`, `maf`, `pval`, `slope`, `slope_se` (typed) | + small `gtex_tissues`/`gtex_genes` lookup files |
| `annotations/<build>/arrays/chr<N>.parquet` | `pos` | `pos`, `platform_codes` | + `platforms.parquet` |
| `populations.json` | — | per build: sample order + population → sample indexes | tiny; drives the bit masks |

### Genotype encoding estimate — `scripts/profile/geno_layout_estimate.py`, first 100k chr22 records

| | GRCh37 (5008 haps) | High Coverage (6404 haps) |
|---|---|---|
| VCF text (uncompressed) | 10,224 B/variant | 14,009 B/variant |
| packed bitset, raw | 626 B | 801 B |
| bitset, deflate (1000-variant blocks) | 48.6 B | 52.6 B |
| bitset, zstd-3 | **37.5 B** | **42.2 B** |
| sparse minor-allele indexes, deflate | 50.2 B | 48.5 B |
| MAF < 0.1% / < 1% / < 5% / ≥ 5% | 63% / 19% / 7% / 11% | 53% / 31% / 10% / 6% |

→ roughly 40–45 MB per chr22 and ~3–4 GB per build genome-wide (vs 14–37 GB of VCF.gz). A 300-SNP LDmatrix
touches at most 300 pages of ~40 KB (≤ ~12 MB, usually far less because query SNPs cluster), plus dbSNP lookups.
Open questions for the layout phase: which codecs hyparquet-writer can emit vs. what hyparquet can read without
extra compressors (snappy is built in; zstd needs `hyparquet-compressors`), page size vs. request count, and whether
to drop very rare variants (MAF < 0.1%) into a separate file since LD for them is rarely requested.
