-- Compare LDlink's 2023 dbsnp dump (table d in dbsnp_dump/work.duckdb) with an NCBI dbSNP VCF rebuild.
-- Run with tools/run_sql.py (one statement per process; paths are relative to LD_DATA):
--   python -I tools/run_sql.py ncbi_dbsnp/<build>/compare.duckdb <compare.sql with {B} replaced>   (after sed s/{B}/<build>/)
-- Inputs: ncbi_dbsnp/{B}/vcf37.parquet, vcf38.parquet (pipeline/ncbi/vcf2tsv.js → DuckDB COPY), merged.tsv (pipeline/ncbi/merged2tsv.py).
-- LDlink's import (scripts/dbSNP_mongo_update/rsjson_parse_data_mongo.py) per refSNP JSON record:
--   one doc per rsID in [refsnp_id] + dbsnp1_merges, ref_id = refsnp_id; position = first number of the REF allele
--   HGVS on the GRCh37.p13 / GRCh38.p13 chromosome placement; type = variant_type; function = set of protein SO terms.
--   Then (rsjson_process_merges.py) merged rsIDs not yet present get a copy of their target's doc, ref_id = target.
SET memory_limit = '40GB';
SET temp_directory = 'ncbi_dbsnp/duckdb_tmp_cmp';
SET preserve_insertion_order = false;

ATTACH 'dbsnp_dump/work.duckdb' AS dump (READ_ONLY);

-- one row per rsID per assembly (first VCF row; count kept so multi-placement rsIDs are visible)
CREATE OR REPLACE TABLE n37 AS
SELECT rsid, count(*) AS n, arg_min(chrom, ord) AS chrom, arg_min(pos, ord) AS pos, arg_min(CASE WHEN vc = 'INS' THEN pos WHEN list_bool_and(list_transform(string_split(alt, ','), lambda a: a[1] = ref[1])) AND (length(ref) > 1 OR vc NOT IN ('SNV', 'MNV')) THEN pos + 1 ELSE pos END, ord) AS start,
  arg_min(ref, ord) AS ref, arg_min(alt, ord) AS alt, arg_min(vc, ord) AS vc, arg_min(fx, ord) AS fx
FROM 'ncbi_dbsnp/{B}/vcf37.parquet' GROUP BY rsid;

CREATE OR REPLACE TABLE n38 AS
SELECT rsid, count(*) AS n, arg_min(chrom, ord) AS chrom, arg_min(pos, ord) AS pos, arg_min(CASE WHEN vc = 'INS' THEN pos WHEN list_bool_and(list_transform(string_split(alt, ','), lambda a: a[1] = ref[1])) AND (length(ref) > 1 OR vc NOT IN ('SNV', 'MNV')) THEN pos + 1 ELSE pos END, ord) AS start,
  arg_min(ref, ord) AS ref, arg_min(alt, ord) AS alt, arg_min(vc, ord) AS vc, arg_min(fx, ord) AS fx
FROM 'ncbi_dbsnp/{B}/vcf38.parquet' GROUP BY rsid;

CREATE OR REPLACE TABLE ncbi AS
SELECT coalesce(a.rsid, b.rsid) AS rsid, coalesce(b.chrom, a.chrom) AS chrom, a.chrom AS chrom37, b.chrom AS chrom38,
  a.n AS n37, b.n AS n38, a.pos AS pos37, a.start AS start37, b.pos AS pos38, b.start AS start38,
  coalesce(b.vc, a.vc) AS vc, coalesce(b.fx, a.fx) AS fx, coalesce(b.ref, a.ref) AS ref38, coalesce(b.alt, a.alt) AS alt38,
  a.ref AS ref37, a.alt AS alt37
FROM n37 a FULL OUTER JOIN n38 b ON a.rsid = b.rsid;

CREATE OR REPLACE TABLE merged AS
SELECT * FROM read_csv('ncbi_dbsnp/{B}/merged.tsv', delim = '\t', header = false, quote = '', escape = '',
  columns = {'rsid': 'BIGINT', 'target': 'BIGINT', 'k': 'INTEGER', 'n_targets': 'INTEGER'});

-- primary dump docs (id == ref_id), first doc per rsID in dump order (what find_one returns)
CREATE OR REPLACE TABLE dp AS
SELECT rsid, count(*) AS ndocs, arg_min(chrom, ord) AS chrom, arg_min(p37, ord) AS p37, arg_min(p38, ord) AS p38,
  arg_min(type, ord) AS type, arg_min(function, ord) AS function
FROM dump.d WHERE rsid = ref_rsid GROUP BY rsid;

CREATE OR REPLACE TABLE cmp AS
SELECT coalesce(p.rsid, n.rsid) AS rsid, p.rsid IS NOT NULL AS in_dump, n.rsid IS NOT NULL AS in_ncbi,
  p.ndocs, p.chrom AS d_chrom, n.chrom AS n_chrom, n.chrom37, n.chrom38, n.n37, n.n38,
  p.p37, n.pos37, n.start37, p.p38, n.pos38, n.start38, p.type, n.vc, p.function, n.fx, n.ref38, n.alt38
FROM dp p FULL OUTER JOIN ncbi n ON p.rsid = n.rsid;

-- LDlink position rule (verified genome-wide on b155): HGVS start of the REF allele = POS for SNV/MNV and VC=INS,
-- POS+1 for any other anchored record (deletions, delins incl. insertions inside repeats).
-- Summary: every count is over rsIDs (primary docs only). Position classes per assembly:
--   eq_start = dump == LDlink-style start from VCF, eq_pos = dump == VCF POS (but not start), na_both, dump_na_only,
--   ncbi_na_only, differ.
CREATE OR REPLACE TABLE summary AS
WITH c AS (
  SELECT *,
    CASE WHEN p37 IS NULL AND pos37 IS NULL THEN 'na_both' WHEN p37 IS NULL THEN 'dump_na_only' WHEN pos37 IS NULL THEN 'ncbi_na_only'
         WHEN p37 = start37 THEN 'eq_start' WHEN p37 = pos37 THEN 'eq_pos' ELSE 'differ' END AS c37,
    CASE WHEN p38 IS NULL AND pos38 IS NULL THEN 'na_both' WHEN p38 IS NULL THEN 'dump_na_only' WHEN pos38 IS NULL THEN 'ncbi_na_only'
         WHEN p38 = start38 THEN 'eq_start' WHEN p38 = pos38 THEN 'eq_pos' ELSE 'differ' END AS c38,
    lower(CASE vc WHEN 'INDEL' THEN 'delins' ELSE vc END) AS ntype
  FROM cmp WHERE in_dump AND in_ncbi
)
SELECT 'rsids' AS metric, 'in_dump_and_ncbi' AS k, count(*) AS n FROM c
UNION ALL SELECT 'rsids', 'dump_only', count(*) FROM cmp WHERE in_dump AND NOT in_ncbi
UNION ALL SELECT 'rsids', 'ncbi_only', count(*) FROM cmp WHERE in_ncbi AND NOT in_dump
UNION ALL SELECT 'rsids', 'dump_multi_doc', count(*) FROM cmp WHERE ndocs > 1
UNION ALL SELECT 'rsids', 'ncbi_multi_row37', count(*) FROM cmp WHERE n37 > 1
UNION ALL SELECT 'rsids', 'ncbi_multi_row38', count(*) FROM cmp WHERE n38 > 1
UNION ALL SELECT 'chrom', CASE WHEN d_chrom = n_chrom THEN 'eq' ELSE 'differ' END, count(*) FROM c GROUP BY 2
UNION ALL SELECT 'pos37', c37, count(*) FROM c GROUP BY 2
UNION ALL SELECT 'pos38', c38, count(*) FROM c GROUP BY 2
UNION ALL SELECT 'pos37_by_type:' || coalesce(type, 'NULL'), c37, count(*) FROM c GROUP BY 1, 2
UNION ALL SELECT 'pos38_by_type:' || coalesce(type, 'NULL'), c38, count(*) FROM c GROUP BY 1, 2
UNION ALL SELECT 'type', CASE WHEN type = ntype THEN 'eq' ELSE coalesce(type, 'NULL') || '->' || coalesce(ntype, 'NULL') END, count(*) FROM c GROUP BY 2
UNION ALL SELECT 'dump_only_by_type', coalesce(type, 'NULL'), count(*) FROM cmp WHERE in_dump AND NOT in_ncbi GROUP BY 2
UNION ALL SELECT 'ncbi_only_by_vc', coalesce(vc, 'NULL'), count(*) FROM cmp WHERE in_ncbi AND NOT in_dump GROUP BY 2;

-- function vs. VCF flags, as sets, restricted to the four terms both sides can express
CREATE OR REPLACE TABLE summary_function AS
WITH c AS (
  SELECT list_sort(list_filter(string_split(coalesce(function, ''), ','), lambda t: t IN ('missense', 'synonymous', 'stop_gained', 'frameshift'))) AS df,
    list_sort(list_transform(list_filter(string_split(coalesce(fx, ''), ','), lambda f: f IN ('NSM', 'SYN', 'NSN', 'NSF')),
      lambda f: CASE f WHEN 'NSM' THEN 'missense' WHEN 'SYN' THEN 'synonymous' WHEN 'NSN' THEN 'stop_gained' ELSE 'frameshift' END)) AS nf
  FROM cmp WHERE in_dump AND in_ncbi
)
SELECT CASE WHEN df = nf THEN 'eq' WHEN len(df) = 0 THEN 'ncbi_flags_only' WHEN len(nf) = 0 THEN 'dump_only'
            WHEN list_has_all(nf, df) THEN 'ncbi_superset' WHEN list_has_all(df, nf) THEN 'dump_superset' ELSE 'differ' END AS k,
  count(*) AS n FROM c GROUP BY 1;

-- merged docs (id != ref_id): is (rsid -> ref_rsid) one of NCBI's merged_into pairs?
CREATE OR REPLACE TABLE summary_merged AS
WITH dm AS (SELECT DISTINCT rsid, ref_rsid FROM dump.d WHERE rsid <> ref_rsid),
j AS (SELECT dm.*, m.k, dm.rsid IN (SELECT rsid FROM merged) AS any_row FROM dm LEFT JOIN merged m ON m.rsid = dm.rsid AND m.target = dm.ref_rsid)
SELECT CASE WHEN k IS NOT NULL THEN 'pair_in_merged' WHEN any_row THEN 'rsid_merged_elsewhere' ELSE 'not_in_merged' END AS k, count(*) AS n FROM j GROUP BY 1
UNION ALL SELECT 'merged_rsids_ncbi', count(DISTINCT rsid) FROM merged
UNION ALL SELECT 'merged_rsids_target_in_ncbi', count(DISTINCT m.rsid) FROM merged m SEMI JOIN ncbi n ON n.rsid = m.target
UNION ALL SELECT 'merged_rsids_ncbi_missing_from_dump', count(*) FROM (SELECT DISTINCT m.rsid FROM merged m SEMI JOIN ncbi n ON n.rsid = m.target) x ANTI JOIN dm ON dm.rsid = x.rsid;

COPY (SELECT * FROM summary ORDER BY metric, n DESC) TO 'ncbi_dbsnp/{B}/summary.csv' (HEADER);
COPY (SELECT * FROM summary_function ORDER BY n DESC) TO 'ncbi_dbsnp/{B}/summary_function.csv' (HEADER);
COPY (SELECT * FROM summary_merged ORDER BY n DESC) TO 'ncbi_dbsnp/{B}/summary_merged.csv' (HEADER);
-- samples of each disagreement class for inspection
COPY (
  SELECT * FROM (SELECT *, row_number() OVER (PARTITION BY cls ORDER BY hash(rsid)) AS rn FROM (
    SELECT *, CASE WHEN NOT in_ncbi THEN 'dump_only' WHEN NOT in_dump THEN 'ncbi_only'
      WHEN d_chrom <> n_chrom THEN 'chrom' WHEN p37 IS DISTINCT FROM start37 AND p37 IS DISTINCT FROM pos37 THEN 'pos37'
      WHEN p38 IS DISTINCT FROM start38 AND p38 IS DISTINCT FROM pos38 THEN 'pos38' END AS cls FROM cmp) WHERE cls IS NOT NULL)
  WHERE rn <= 200
) TO 'ncbi_dbsnp/{B}/mismatch_samples.parquet' (FORMAT parquet);
