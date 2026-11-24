-- dbsnp_src/dup_ids.parquet: the rsIDs with more than one doc in LDlink's dump (296,308), sorted. Input of
-- pipeline/dbsnp/dup_rsids.js (the release's dbsnp/duplicated_rsids.parquet). Run after pipeline/dbsnp/sort.sql:
-- usage: python -I tools/run_sql.py dbsnp_dump/work.duckdb pipeline/dbsnp/dup_ids.sql
SET memory_limit = '32GB';
SET temp_directory = 'dbsnp_dump/duckdb_tmp';
COPY (SELECT rsid FROM krs GROUP BY rsid HAVING count(*) > 1 ORDER BY rsid) TO 'dbsnp_src/dup_ids.parquet' (FORMAT parquet);
