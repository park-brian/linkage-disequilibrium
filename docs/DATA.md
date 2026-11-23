# The release: what the app reads

A release is a directory of Parquet files plus `manifest.json`, served as static files with HTTP Range support (local:
`npm run serve`; hosted: a bucket behind HTTP/2). Current: `release_v4` (`$LD_RELEASE`), 324 files + manifest,
34.8 GB. A private copy is at `s3://eaap-ac-test-bucket/ldlink-parquet/v4/`.

**Conventions:**
- **Genome key** = chromosome code · 10⁹ + position (INT64). Codes: 1–22, X = 23, Y = 24, MT = 25. Positions are in
  the file's build.
- **Order** = sorted by key, then by the source's dump order (`ord`). Wherever LDlink reads Mongo with `find_one`,
  the first row for a key is the row LDlink returns.
- **Compression and indexes:** ZSTD; page indexes on the sort column.
- **Metadata:** `source`, `order` and other provenance keys in the key-value metadata.
- **Builds:** `grch37`, `grch38` and `grch38_high_coverage`. High Coverage uses GRCh38 positions, so its dbSNP and
  annotation files are the `grch38` ones (`src/genome.js` `BUILDS`).

## Hosted copy

`https://d2x2yo644ey89o.cloudfront.net` serves `release_v4` as uploaded to `s3://eaap-ac-test-bucket/ldlink-parquet/v5/`
(AWS account 285397596138, us-east-1). The bucket stays private (Block Public Access on); CloudFront distribution
`E1S1X441YI89BW` reads it through an origin access control, and the bucket policy lets only that distribution read
`v5/`. Responses allow any origin (CORS with `Range`, `Timing-Allow-Origin`), are compressed (the manifest travels as
0.6 MB), cached a day at the edge (`Cache-Control: public, max-age=86400`), over HTTP/2 and 3 (price class 100).
A new release goes to a new prefix (`v6/`) and a new origin path, not over `v5/`.

## manifest.json

Version 2: `{ version: 2, created, files: { <path>: { byteLength, numRows, footer: { offset, length }, key, columns,
kv, rowGroups } } }`. `rowGroups` is stored column-wise, one array entry per row group: `{ min: [], max: [], bytes: [],
keyIndex?: [], payloadOffsetIndex?: [], colBytes?: { <column>: [] } }`.
- `min` / `max`: the key column's statistics (`null`: none). `keyIndex`: `[columnIndexOffset, length, offsetIndexOffset,
  length]` of the key column. `payloadOffsetIndex`: `[offset, length]` of the `haps` offset index.
- `colBytes`: genotype files only, for the columns the planner prices (`PRICED_COLUMNS` in `src/planner.js`).
- `src/planner.js` `expandManifest` turns it back into one object per row group, and also reads version 1 (one object
  per row group, `colBytes` for every column), which is what the bucket's `v4/` still holds.
- The planner prices reads from it. Each file opens with its `byteLength` (no size probe) and its footer in one exact
  Range request.
- release_v4: 1.92 MB (0.61 MB gzipped), down from 7.15 MB in version 1 (1.49 MB gzipped). It is the app's whole startup
  download. The version 1 file is kept at `$LD_DATA/release_v4_manifest.v1.json`.
- A path that is not in the manifest is absent (`exists`).

Built by `pipeline/release/manifest.js`.

## genotypes/<build>/chr<N>.<tier>.parquet: 1000 Genomes haplotypes

| | |
|---|---|
| files | 70 chromosomes × 3 MAF tiers = 210 files, 10.6 GB. GRCh37: 1–22, X, Y; GRCh38 and High Coverage: 1–22, X |
| tiers | by `maf_max` over the 26 sub-populations: `common` ≥ 0.05 > `low` ≥ 0.01 > `rare` |
| sort | `pos`, then `vcf_order` |
| columns | `pos` INT32 · `rsid` INT64? · `rsid_pos_delta` INT32? · `rsid_candidates` · `rsid_source` · `ref`, `alt` str · `multiallelic_source` bool · `male_haploid` bool · `ac` · `haps` · `missing` bytes? · `maf_max` float · `function` str? · `vcf_order` INT32 · `vcf_id` str? (GRCh37) · `allele_is_label` bool · `allele_rsids` str? |
| key-value | `samples` (VCF sample order) · `male_haploid_samples` (sample indexes) · `max_ref_len` · `tier_rule` · `rsid_rule` · `vcf_order_rule` · `vcf_id_rule` · `allele_rule` |

**Column semantics:**
- **`haps`:** one bit per haplotype. In `common`/`low` it is a fixed-length bitset (2 × samples bits, 626 bytes for
  2,504 samples). In `rare` it is a sparse varint list: byte 0 is `listRef`, then index deltas
  (`src/haps.js` `decodeHaps`).
- **`missing`:** only true `.` calls.
- **`male_haploid` + `male_haploid_samples`:** they mark which samples' second haplotype is absent (chrX
  non-PAR, chrY).
- **Multi-allelic records** are split into one row per ALT (`multiallelic_source`).
- **`vcf_order`:** the source VCF row index, so LDlink's "last overlapping row" = the highest `vcf_order` across tiers.
- **`rsid`:** the label LDlink would show (`pipeline/genotypes/rsid_label.js`):
  - **`rsid_source`:** 1 = VCF rsID, 2 = merged → current, 3 = position + type, 4 = NCBI allele match.
  - **`rsid_pos_delta`:** the dump position − POS, where deletions sit at POS + 1.
- **`vcf_id`:** the raw VCF ID when it contains "rs" and is not `rs<rsid>`. LDlink prints it in warnings and in
  LDproxy's RS_Number.
- **`allele_rsids` / `allele_is_label`:** the NCBI b155 rsIDs whose alleles match the row (allele-aware mode).

## dbsnp/: LDlink's 2023 dbSNP dump (b155)

| file | rows | columns | for |
|---|---|---|---|
| `rsid/part-NN.parquet` (32, by rsid ≫ 26) | 0.64 B | `rsid`, `ref_rsid`?, `chrom`?, `key37`?, `key38`? | `find_one({id})`: first row per rsID |
| `<grch37\|grch38>/chr<N>.parquet` (48) | | `key`, `rsid`, `is_ref`, `type`?, `function`? | `find({chromosome, position})`, all rows per key in dump order |
| `duplicated_rsids.parquet` | 296,308 | `rsid` | rsIDs with more than one dump doc |

Built by `pipeline/dbsnp/` from the mongodump (`bson2tsv` → `sort.sql` → `build.js`).

## annotations/: LDlink's other collections

| file | key | columns | LDlink query reproduced |
|---|---|---|---|
| `<b>/genes.parquet` | txStart | refGene columns | `getRefGene` (overlap via `annotations.genes`) |
| `<b>/recomb.parquet` | position | `rate` | `getRecomb` |
| `<b>/regulome.parquet` | position | `score` | `get_regDB`: `find_one({chromosome, position})` |
| `<b>/forge.parquet` | the rsID's first dbSNP doc position | `rsid`, `score` | the FORGEdb fast path: each rsID once |
| `forge_rsid.parquet` | rsid | `score` | `get_forgeDB`: `find_one({snp_id})` |
| `<b>/arrays.parquet` | position | `platforms` (`a\|b` per doc) | SNPchip `snp_col.find` |
| `platforms.parquet` | — | `platform`, `code` | SNPchip platform names |
| `<b>/gwas.parquet` | position | GWAS Catalog columns (strings, as stored) | LDtrait `gwas_catalog.find` |
| `<b>/gtex_eqtl.parquet` | position | `gene_id`, `tissue`, `slope`, `pval_nominal`, … (doubles; `%.6g` reproduces GTEx text) | LDexpress `gtex_tissue_eqtl` |
| `gtex_tissues.parquet`, `gtex_genes.parquet` | id | joined like LDexpress's `$lookup` / `$unwind` | |
| `gene_names.parquet` | name | per-build chromosome / begin / end | `genes_name_coords.find_one` (LDassoc gene mode) |
| `panels.parquet` | — | `pop`, `sample` (each list in file order; 32 pops, 7,512 rows) | LDlink's `<POP>.txt` sample lists |

Built by `pipeline/annotations/` (`fetch.sh` → `annotations.sql` → `build.js`); `panels.parquet` by
`pipeline/release/extras.js`.

## ldscore/: LD score reference panels

One file per population (`afr`, `eur`, …): `CHR`, `SNP`, `BP`, `L2`, plus `M` / `M_5_50` in the key-value
metadata. `w_hm3.snplist.parquet` holds HapMap3 SNPs and alleles. Used by `src/ldsc/` (heritability, genetic
correlation).

`ldsc_counts.parquet` (`pop`, `rows`) holds, per panel, the SNP count ldsc logs ("Read reference panel LD Scores for
N SNPs"): distinct SNP IDs after ldsc's `sep=r'\s+'` read, in which a row with an empty SNP field shifts left and its
BP becomes the SNP (eas: 2,079 such rows; 8,854,920 counted against 8,852,842 distinct IDs). Written by
`pipeline/release/extras.js`, so a browser never holds a set of up to 17 million IDs (`src/ldsc/panels.js` keeps
only the SNPs an analysis needs).

## examples/: LDlink's example inputs

Plain files, copied byte for byte by `pipeline/release/extras.js`, for the example buttons:
- LDassoc: `prostate_example_grch37.txt` and `prostate_example_grch38.txt`;
- LDscore: `22.bed`, `22.bim`, `22.fam`;
- LDherit / LDcorrelation: `BBJ_HDLC22.txt`, `BBJ_LDLC22.txt`.

With panels and examples inside, the app needs one data URL and nothing else. `$LD_PANELS` stays the pipeline's
source for the genotype build and the kernel verifier.
