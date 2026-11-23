-- Key and sort LDlink's annotation collections into intermediates (DuckDB; final files are written by
-- pipeline/annotations/build.mjs with hyparquet-writer). Run after pipeline/annotations/fetch.sh:
--   python -I tools/run_sql.py dbsnp_dump/work.duckdb pipeline/annotations/annotations.sql
-- Parity: rows keep dump order (ord) within equal keys, because LDlink's find/find_one return index order, which
-- for equal keys is insertion (dump) order. Composite key = chrom code (1-22, X=23, Y=24, MT=25) * 1e9 + position.
SET memory_limit = '32GB';
SET threads = 8;
SET temp_directory = 'dbsnp_dump/duckdb_tmp';
SET preserve_insertion_order = true;

CREATE OR REPLACE MACRO chrom_code(c) AS
  CASE upper(regexp_replace(CAST(c AS VARCHAR), '^(chr|CHR)', '')) WHEN 'X' THEN 23 WHEN 'Y' THEN 24 WHEN 'MT' THEN 25
       ELSE try_cast(regexp_replace(CAST(c AS VARCHAR), '^(chr|CHR)', '') AS INTEGER) END;
CREATE OR REPLACE MACRO mkkey(c, p) AS
  CASE WHEN chrom_code(c) BETWEEN 1 AND 25 AND try_cast(p AS BIGINT) BETWEEN 1 AND 999999999
       THEN CAST(chrom_code(c) AS BIGINT) * 1000000000 + try_cast(p AS BIGINT) END;
CREATE OR REPLACE MACRO tsv(path, cols) AS TABLE
  SELECT * FROM read_csv(path, delim = '\t', header = false, quote = '', escape = '', all_varchar = true, names = cols, null_padding = true);

-- ---------- recombination: rate (cM/Mb) by position, range queries ----------
CREATE OR REPLACE TABLE a_recomb AS SELECT * FROM tsv('annot_src/recomb.tsv', ['ord', 'rate', 'c37', 'p37', 'c38', 'p38']);
COPY (SELECT mkkey(c37, p37) AS key, CAST(rate AS DOUBLE) AS rate, CAST(ord AS BIGINT) AS ord FROM a_recomb WHERE mkkey(c37, p37) IS NOT NULL ORDER BY key, ord)
  TO 'annot_src/k_recomb37.parquet' (FORMAT parquet);
COPY (SELECT mkkey(c38, p38) AS key, CAST(rate AS DOUBLE) AS rate, CAST(ord AS BIGINT) AS ord FROM a_recomb WHERE mkkey(c38, p38) IS NOT NULL ORDER BY key, ord)
  TO 'annot_src/k_recomb38.parquet' (FORMAT parquet);

-- ---------- RegulomeDB: score by position, find_one (first doc in dump order) ----------
CREATE OR REPLACE TABLE a_regulome AS SELECT * FROM tsv('annot_src/regulome.tsv', ['ord', 'score', 'c37', 'p37', 'c38', 'p38']);
COPY (SELECT mkkey(c37, p37) AS key, score, CAST(ord AS BIGINT) AS ord FROM a_regulome WHERE mkkey(c37, p37) IS NOT NULL ORDER BY key, ord)
  TO 'annot_src/k_regulome37.parquet' (FORMAT parquet);
COPY (SELECT mkkey(c38, p38) AS key, score, CAST(ord AS BIGINT) AS ord FROM a_regulome WHERE mkkey(c38, p38) IS NOT NULL ORDER BY key, ord)
  TO 'annot_src/k_regulome38.parquet' (FORMAT parquet);

-- ---------- FORGEdb: score by rsID (find_one on snp_id) + position-keyed copies via the rsID index ----------
CREATE OR REPLACE TABLE a_forge AS SELECT CAST(ord AS BIGINT) AS ord, snp_id, score FROM tsv('annot_src/forge_score.tsv', ['ord', 'snp_id', 'score']);
COPY (SELECT try_cast(regexp_replace(snp_id, '^rs', '') AS BIGINT) AS rsid, snp_id, try_cast(score AS INTEGER) AS score, ord
      FROM a_forge WHERE try_cast(regexp_replace(snp_id, '^rs', '') AS BIGINT) IS NOT NULL ORDER BY rsid, ord)
  TO 'annot_src/k_forge_rsid.parquet' (FORMAT parquet);
-- position copies: the first FORGEdb doc per rsID, placed at the rsID's LDlink position (first dbSNP doc in dump order)
CREATE OR REPLACE TABLE a_forge_first AS SELECT * FROM (SELECT try_cast(regexp_replace(snp_id, '^rs', '') AS BIGINT) AS rsid, try_cast(score AS INTEGER) AS score, ord,
  row_number() OVER (PARTITION BY snp_id ORDER BY ord) AS rn FROM a_forge) WHERE rn = 1 AND rsid IS NOT NULL;
CREATE OR REPLACE TABLE a_dbsnp_first AS SELECT rsid, key37, key38 FROM (SELECT rsid, key37, key38, row_number() OVER (PARTITION BY rsid ORDER BY ord) AS rn
  FROM krs WHERE rsid IN (SELECT rsid FROM a_forge_first)) WHERE rn = 1;
COPY (SELECT d.key37 AS key, f.rsid, f.score FROM a_forge_first f JOIN a_dbsnp_first d USING (rsid) WHERE d.key37 IS NOT NULL ORDER BY key, f.rsid)
  TO 'annot_src/k_forge37.parquet' (FORMAT parquet);
COPY (SELECT d.key38 AS key, f.rsid, f.score FROM a_forge_first f JOIN a_dbsnp_first d USING (rsid) WHERE d.key38 IS NOT NULL ORDER BY key, f.rsid)
  TO 'annot_src/k_forge38.parquet' (FORMAT parquet);

-- ---------- refGene: transcripts, overlap queries on (chrom, txStart, txEnd) ----------
CREATE OR REPLACE TABLE a_refgene37 AS SELECT * FROM tsv('annot_src/refGene_grch37.tsv', ['ord', 'bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']);
CREATE OR REPLACE TABLE a_refgene38 AS SELECT * FROM tsv('annot_src/refGene_grch38.tsv', ['ord', 'bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']);
COPY (SELECT mkkey(chrom, CAST(txStart AS BIGINT) + 1) - 1 AS key, chrom, CAST(txStart AS INTEGER) txStart, CAST(txEnd AS INTEGER) txEnd, name, name2, strand,
      CAST(cdsStart AS INTEGER) cdsStart, CAST(cdsEnd AS INTEGER) cdsEnd, CAST(exonCount AS INTEGER) exonCount, exonStarts, exonEnds, exonFrames, cdsStartStat, cdsEndStat,
      CAST(bin AS INTEGER) bin, CAST(score AS INTEGER) score, CAST(ord AS BIGINT) ord
      FROM a_refgene37 WHERE chrom_code(chrom) BETWEEN 1 AND 25 ORDER BY chrom_code(chrom), txStart::INTEGER, txEnd::INTEGER, ord)
  TO 'annot_src/k_genes37.parquet' (FORMAT parquet);
COPY (SELECT mkkey(chrom, CAST(txStart AS BIGINT) + 1) - 1 AS key, chrom, CAST(txStart AS INTEGER) txStart, CAST(txEnd AS INTEGER) txEnd, name, name2, strand,
      CAST(cdsStart AS INTEGER) cdsStart, CAST(cdsEnd AS INTEGER) cdsEnd, CAST(exonCount AS INTEGER) exonCount, exonStarts, exonEnds, exonFrames, cdsStartStat, cdsEndStat,
      CAST(bin AS INTEGER) bin, CAST(score AS INTEGER) score, CAST(ord AS BIGINT) ord
      FROM a_refgene38 WHERE chrom_code(chrom) BETWEEN 1 AND 25 ORDER BY chrom_code(chrom), txStart::INTEGER, txEnd::INTEGER, ord)
  TO 'annot_src/k_genes38.parquet' (FORMAT parquet);

-- ---------- gene name → coordinates (LDassoc gene mode; find_one on name) ----------
CREATE OR REPLACE TABLE a_gnc AS SELECT * FROM tsv('annot_src/genes_name_coords.tsv', ['ord', 'name', 'c37', 'b37', 'e37', 'c38', 'b38', 'e38']);
COPY (SELECT name, c37 AS chromosome_grch37, try_cast(b37 AS INTEGER) AS begin_grch37, try_cast(e37 AS INTEGER) AS end_grch37,
             c38 AS chromosome_grch38, try_cast(b38 AS INTEGER) AS begin_grch38, try_cast(e38 AS INTEGER) AS end_grch38, CAST(ord AS BIGINT) ord
      FROM a_gnc ORDER BY name, ord) TO 'annot_src/k_gene_names.parquet' (FORMAT parquet);

-- ---------- GWAS Catalog (2025 backup): window queries on (chromosome, position) ----------
CREATE OR REPLACE TABLE a_gwas AS SELECT * FROM tsv('annot_src/gwas_catalog.tsv',
  ['ord', 'date_added', 'pubmedid', 'first_author', 'date', 'journal', 'study', 'trait', 'mapped_trait', 'mapped_gene', 'strongest_risk_allele', 'snps', 'merged',
   'snp_id_current', 'context', 'risk_allele_frequency', 'p_value', 'pvalue_mlog', 'p_value_text', 'or_beta', 'ci_95_text', 'study_accession', 'chromosome', 'p37', 'p38']);
COPY (SELECT mkkey(chromosome, p37) AS key, * EXCLUDE (p37, p38, ord), CAST(ord AS BIGINT) ord FROM a_gwas WHERE mkkey(chromosome, p37) IS NOT NULL ORDER BY key, ord)
  TO 'annot_src/k_gwas37.parquet' (FORMAT parquet);
COPY (SELECT mkkey(chromosome, p38) AS key, * EXCLUDE (p37, p38, ord), CAST(ord AS BIGINT) ord FROM a_gwas WHERE mkkey(chromosome, p38) IS NOT NULL ORDER BY key, ord)
  TO 'annot_src/k_gwas38.parquet' (FORMAT parquet);

-- ---------- GTEx v8 eQTL: match on tissue + position, inner-joined to tissues and genes (LDexpress $unwind) ----------
CREATE OR REPLACE TABLE a_gtex_t AS SELECT * FROM tsv('annot_src/gtex_tissues.tsv', ['ord', 'tissueSiteDetailId', 'tissueSiteDetail', 'tissueSiteDetailAbbr', 'tissueSite', 'colorHex', 'colorRgb', 'datasetId', 'hasEGenes', 'eGeneCount', 'rnaSeqAndGenotypeSampleCount', 'uberonId']);
CREATE OR REPLACE TABLE a_gtex_g AS SELECT * FROM tsv('annot_src/gtex_genes.tsv', ['ord', 'gene_id_v26', 'gene_name_v26', 'ncbi_id']);
COPY (SELECT * EXCLUDE (ord), CAST(ord AS BIGINT) ord FROM a_gtex_t ORDER BY tissueSiteDetailId, ord) TO 'annot_src/k_gtex_tissues.parquet' (FORMAT parquet);
COPY (SELECT * EXCLUDE (ord), CAST(ord AS BIGINT) ord FROM a_gtex_g ORDER BY gene_id_v26, ord) TO 'annot_src/k_gtex_genes.parquet' (FORMAT parquet);
CREATE OR REPLACE TABLE a_gtex AS SELECT * FROM tsv('annot_src/gtex_tissue_eqtl.tsv',
  ['ord', 'variant_id', 'gene_id', 'tss_distance', 'ma_samples', 'ma_count', 'maf', 'pval_nominal', 'slope', 'slope_se', 'pval_nominal_threshold', 'min_pval_nominal', 'pval_beta', 'tissue', 'c37', 'p37', 'c38', 'p38']);
COPY (SELECT mkkey(c37, p37) AS key, variant_id, gene_id, tissue, CAST(maf AS DOUBLE) maf, CAST(pval_nominal AS DOUBLE) pval_nominal, CAST(slope AS DOUBLE) slope,
             CAST(slope_se AS DOUBLE) slope_se, CAST(tss_distance AS INTEGER) tss_distance, CAST(ord AS BIGINT) ord
      FROM a_gtex WHERE mkkey(c37, p37) IS NOT NULL ORDER BY key, ord) TO 'annot_src/k_gtex37.parquet' (FORMAT parquet);
COPY (SELECT mkkey(c38, p38) AS key, variant_id, gene_id, tissue, CAST(maf AS DOUBLE) maf, CAST(pval_nominal AS DOUBLE) pval_nominal, CAST(slope AS DOUBLE) slope,
             CAST(slope_se AS DOUBLE) slope_se, CAST(tss_distance AS INTEGER) tss_distance, CAST(ord AS BIGINT) ord
      FROM a_gtex WHERE mkkey(c38, p38) IS NOT NULL ORDER BY key, ord) TO 'annot_src/k_gtex38.parquet' (FORMAT parquet);

-- ---------- genotyping arrays (SNPchip): platforms by position ----------
CREATE OR REPLACE TABLE a_plat AS SELECT * FROM tsv('annot_src/platforms.tsv', ['ord', 'platform', 'code']);
COPY (SELECT platform, code, CAST(ord AS BIGINT) ord FROM a_plat ORDER BY ord) TO 'annot_src/k_platforms.parquet' (FORMAT parquet);
CREATE OR REPLACE TABLE a_snpcol AS SELECT * FROM tsv('annot_src/snp_col.tsv', ['ord', 'platforms', 'c37', 'p37', 'c38', 'p38']);
COPY (SELECT mkkey(c37, p37) AS key, platforms, CAST(ord AS BIGINT) ord FROM a_snpcol WHERE mkkey(c37, p37) IS NOT NULL ORDER BY key, ord) TO 'annot_src/k_arrays37.parquet' (FORMAT parquet);
COPY (SELECT mkkey(c38, p38) AS key, platforms, CAST(ord AS BIGINT) ord FROM a_snpcol WHERE mkkey(c38, p38) IS NOT NULL ORDER BY key, ord) TO 'annot_src/k_arrays38.parquet' (FORMAT parquet);

-- ---------- counts: kept vs. not keyable (NA / unknown chromosome / impossible position) ----------
COPY (
  SELECT 'recomb' c, count(*) docs, count(mkkey(c37, p37)) keyed37, count(mkkey(c38, p38)) keyed38 FROM a_recomb
  UNION ALL SELECT 'regulome', count(*), count(mkkey(c37, p37)), count(mkkey(c38, p38)) FROM a_regulome
  UNION ALL SELECT 'gwas_catalog_2025', count(*), count(mkkey(chromosome, p37)), count(mkkey(chromosome, p38)) FROM a_gwas
  UNION ALL SELECT 'gtex_tissue_eqtl', count(*), count(mkkey(c37, p37)), count(mkkey(c38, p38)) FROM a_gtex
  UNION ALL SELECT 'snp_col', count(*), count(mkkey(c37, p37)), count(mkkey(c38, p38)) FROM a_snpcol
  UNION ALL SELECT 'forge_score', count(*), (SELECT count(*) FROM 'annot_src/k_forge37.parquet'), (SELECT count(*) FROM 'annot_src/k_forge38.parquet') FROM a_forge
  UNION ALL SELECT 'refGene_grch37', count(*), count(*) FILTER (WHERE chrom_code(chrom) BETWEEN 1 AND 25), NULL FROM a_refgene37
  UNION ALL SELECT 'refGene_grch38', count(*), NULL, count(*) FILTER (WHERE chrom_code(chrom) BETWEEN 1 AND 25) FROM a_refgene38
) TO 'annot_src/annotation_counts.csv' (HEADER);
