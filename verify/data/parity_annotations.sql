-- (Concatenations use coalesce: in SQL 'a' || NULL is NULL, which would hide records with empty fields.)
-- Expected LDlink answers for the annotation files, computed from the raw dump tables (a_*: stored strings) with
-- Mongo/LDlink query semantics — independent of the key-building macros. Run after annotations.sql:
--   python -I tools/run_sql.py dbsnp_dump/work.duckdb verify/data/parity_annotations.sql
SET memory_limit = '32GB';
SET threads = 8;
SET temp_directory = 'dbsnp_dump/duckdb_tmp';

-- RegulomeDB: get_regDB = find_one({chromosome_grchXX: "chr"+chr, position_grchXX: int}) → first doc in dump order
CREATE OR REPLACE MACRO reg_cases(cc, pp) AS TABLE
WITH docs AS (SELECT cc AS c, pp AS p, ord::BIGINT ord, score FROM a_regulome WHERE pp <> 'NA' AND cc LIKE 'chr%'),
multi AS (SELECT c, p FROM docs GROUP BY c, p HAVING count(*) > 1),
picks AS (
  (SELECT c, p FROM docs USING SAMPLE reservoir(1500 ROWS) REPEATABLE (7))
  UNION ALL (SELECT c, p FROM multi USING SAMPLE reservoir(300 ROWS) REPEATABLE (7))
  UNION ALL (SELECT c, CAST(CAST(p AS BIGINT) + 1 AS VARCHAR) FROM docs USING SAMPLE reservoir(100 ROWS) REPEATABLE (8))
)
SELECT DISTINCT q.c, q.p, (SELECT score FROM docs d WHERE d.c = q.c AND d.p = q.p ORDER BY ord LIMIT 1) AS expected FROM picks q;
COPY (SELECT * FROM reg_cases(c37, p37)) TO 'annot_src/parity_regulome37.json' (FORMAT json, ARRAY true);
COPY (SELECT * FROM reg_cases(c38, p38)) TO 'annot_src/parity_regulome38.json' (FORMAT json, ARRAY true);

-- FORGEdb: get_forgeDB = find_one({"snp_id": "rs..."}) → first doc in dump order
COPY (
  WITH multi AS (SELECT snp_id FROM a_forge GROUP BY snp_id HAVING count(*) > 1),
  picks AS (
    (SELECT snp_id FROM a_forge USING SAMPLE reservoir(2000 ROWS) REPEATABLE (7))
    UNION ALL (SELECT snp_id FROM multi USING SAMPLE reservoir(200 ROWS) REPEATABLE (7))
    UNION ALL (SELECT 'rs' || (2100000000 + i) FROM range(0, 50) t(i))
  )
  SELECT DISTINCT q.snp_id, (SELECT score FROM a_forge f WHERE f.snp_id = q.snp_id ORDER BY ord LIMIT 1) AS expected FROM picks q
) TO 'annot_src/parity_forge.json' (FORMAT json, ARRAY true);

-- recombination: getRecomb = find({chromosome_grchXX: str(chr), position_grchXX: {$gte, $lte}})
CREATE OR REPLACE MACRO recomb_cases(cc, pp) AS TABLE
WITH docs AS (SELECT cc AS c, try_cast(pp AS BIGINT) AS p, ord::BIGINT ord, rate FROM a_recomb WHERE try_cast(pp AS BIGINT) IS NOT NULL),
w AS (SELECT c, p - 250000 AS lo, p + 250000 AS hi FROM docs USING SAMPLE reservoir(100 ROWS) REPEATABLE (9))
SELECT w.c, w.lo, w.hi, (SELECT count(*) FROM docs d WHERE d.c = w.c AND d.p BETWEEN w.lo AND w.hi) AS n,
  (SELECT sum(CAST(rate AS DOUBLE)) FROM docs d WHERE d.c = w.c AND d.p BETWEEN w.lo AND w.hi) AS rate_sum,
  (SELECT list(p ORDER BY p, ord) FROM docs d WHERE d.c = w.c AND d.p BETWEEN w.lo AND w.hi) AS positions FROM w;
COPY (SELECT * FROM recomb_cases(c37, p37)) TO 'annot_src/parity_recomb37.json' (FORMAT json, ARRAY true);
COPY (SELECT * FROM recomb_cases(c38, p38)) TO 'annot_src/parity_recomb38.json' (FORMAT json, ARRAY true);

-- refGene: getRefGene = find({chrom: "chr"+chr, $or of 4 cases}) = inclusive overlap; order not guaranteed (no sort)
CREATE OR REPLACE MACRO gene_cases(tbl) AS TABLE
WITH docs AS (SELECT chrom, txStart::BIGINT s, txEnd::BIGINT e, name, name2, ord::BIGINT ord FROM query_table(tbl)),
w AS (SELECT chrom, s - 100000 AS lo,
  s + 400000 AS hi FROM docs WHERE chrom ~ '^chr([0-9]+|X|Y)$' USING SAMPLE reservoir(100 ROWS) REPEATABLE (11))
SELECT w.chrom, w.lo, w.hi,
  (SELECT list(coalesce(name, '') || '@' || s || '-' || e ORDER BY name, s, e,
    ord) FROM docs d WHERE d.chrom = w.chrom AND d.s <= w.hi AND d.e >= w.lo) AS transcripts FROM w;
COPY (SELECT * FROM gene_cases('a_refgene37')) TO 'annot_src/parity_genes37.json' (FORMAT json, ARRAY true);
COPY (SELECT * FROM gene_cases('a_refgene38')) TO 'annot_src/parity_genes38.json' (FORMAT json, ARRAY true);

-- GWAS Catalog: get_window_variants = find({chromosome, position_grchXX: {$gte: max(0, pos-w), $lte: pos+w}}) in index order
CREATE OR REPLACE MACRO gwas_cases(pp) AS TABLE
WITH docs AS (SELECT chromosome AS c, try_cast(pp AS BIGINT) AS p, ord::BIGINT ord, pubmedid, snp_id_current,
  trait FROM a_gwas WHERE try_cast(pp AS BIGINT) IS NOT NULL),
w AS (SELECT c, greatest(0, p - 500000) AS lo,
  p + 500000 AS hi FROM docs WHERE c ~ '^([0-9]+|X|Y)$' USING SAMPLE reservoir(100 ROWS) REPEATABLE (13))
SELECT w.c, w.lo, w.hi, (SELECT list(coalesce(pubmedid, '') || '|' || coalesce(snp_id_current, '') || '|' || coalesce(trait, '') ORDER BY p,
  ord) FROM docs d WHERE d.c = w.c AND d.p BETWEEN w.lo AND w.hi) AS records FROM w;
COPY (SELECT * FROM gwas_cases(p37)) TO 'annot_src/parity_gwas37.json' (FORMAT json, ARRAY true);
COPY (SELECT * FROM gwas_cases(p38)) TO 'annot_src/parity_gwas38.json' (FORMAT json, ARRAY true);

-- GTEx (LDexpress): match {tissue ∈ T, chromosome_grchXX: str(chr), position_grchXX: int}, inner-joined to tissues and genes
CREATE OR REPLACE MACRO gtex_cases(cc, pp) AS TABLE
WITH docs AS (SELECT cc AS c, try_cast(pp AS BIGINT) AS p, ord::BIGINT ord, gene_id, tissue, slope,
  pval_nominal FROM a_gtex WHERE try_cast(pp AS BIGINT) IS NOT NULL),
ts AS (SELECT list(tissueSiteDetailId ORDER BY tissueSiteDetailId) l FROM (SELECT tissueSiteDetailId FROM a_gtex_t USING SAMPLE reservoir(8 ROWS) REPEATABLE (17))),
q AS (SELECT DISTINCT c, p FROM docs USING SAMPLE reservoir(150 ROWS) REPEATABLE (19))
SELECT q.c, q.p, (SELECT l FROM ts) AS tissues,
  (SELECT list(coalesce(d.gene_id, '') || '|' || coalesce(d.tissue, '') || '|' || coalesce(d.slope, '') || '|' || coalesce(d.pval_nominal,
    '') || '|' || coalesce(g.gene_name_v26, '') || '|' || coalesce(t.tissueSiteDetail, '') ORDER BY d.ord)
   FROM docs d JOIN a_gtex_t t ON t.tissueSiteDetailId = d.tissue JOIN a_gtex_g g ON g.gene_id_v26 = d.gene_id
   WHERE d.c = q.c AND d.p = q.p AND list_contains((SELECT l FROM ts), d.tissue)) AS rows FROM q;
COPY (SELECT * FROM gtex_cases(c37, p37)) TO 'annot_src/parity_gtex37.json' (FORMAT json, ARRAY true);
COPY (SELECT * FROM gtex_cases(c38, p38)) TO 'annot_src/parity_gtex38.json' (FORMAT json, ARRAY true);

-- SNPchip: find({chromosome_grchXX: str(chr), position_grchXX: int(pos)}) → every doc's platform list, in dump order
CREATE OR REPLACE MACRO chip_cases(cc, pp) AS TABLE
WITH docs AS (SELECT cc AS c, try_cast(pp AS BIGINT) AS p, ord::BIGINT ord,
  platforms FROM a_snpcol WHERE try_cast(pp AS BIGINT) > 0 AND cc ~ '^([0-9]+|X|Y)$'),
picks AS ((SELECT c, p FROM docs USING SAMPLE reservoir(1000 ROWS) REPEATABLE (23)) UNION ALL (SELECT c,
  p + 1 FROM docs USING SAMPLE reservoir(100 ROWS) REPEATABLE (29)))
SELECT DISTINCT q.c, q.p, (SELECT list(coalesce(platforms,
  '') ORDER BY ord) FROM docs d WHERE d.c = q.c AND d.p = q.p) AS platforms FROM picks q;
COPY (SELECT * FROM chip_cases(c37, p37)) TO 'annot_src/parity_arrays37.json' (FORMAT json, ARRAY true);
COPY (SELECT * FROM chip_cases(c38, p38)) TO 'annot_src/parity_arrays38.json' (FORMAT json, ARRAY true);
