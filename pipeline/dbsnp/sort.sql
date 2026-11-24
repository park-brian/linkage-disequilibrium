-- External sort of the dbSNP dump TSV into intermediates (DuckDB; intermediates only — final files are written
-- by hyparquet-writer in pipeline/dbsnp/build.js). Run (the @for blocks need tools/run_sql.py, not the duckdb CLI):
--   python -I tools/run_sql.py dbsnp_dump/work.duckdb pipeline/dbsnp/sort.sql
-- Order matters for parity, so: keep insertion order on, materialize sorted tables, then write each partition with
-- its own explicit ORDER BY (PARTITION_BY writes do not guarantee row order).
-- Parity notes (server/LDcommon.py): rsID lookups are find_one({"id"}) → first doc in dump order;
-- chr:pos lookups are find({"chromosome", "position_grchXX": str(pos)}) → all docs in dump order.
SET memory_limit = '48GB';
SET temp_directory = 'dbsnp_dump/duckdb_tmp';
SET preserve_insertion_order = true;

CREATE OR REPLACE VIEW raw AS
SELECT * FROM read_csv('dbsnp_dump/dbsnp_dump.tsv', delim = '\t', header = false, quote = '', escape = '',
  columns = {'ord': 'BIGINT', 'id': 'VARCHAR', 'ref_id': 'VARCHAR', 'chromosome': 'VARCHAR',
             'position_grch37': 'VARCHAR', 'position_grch38': 'VARCHAR', 'type': 'VARCHAR', 'function': 'VARCHAR'});

CREATE OR REPLACE TABLE d AS
SELECT ord,
  try_cast(id AS BIGINT) AS rsid,
  try_cast(ref_id AS BIGINT) AS ref_rsid,
  id, ref_id, chromosome, position_grch37, position_grch38,
  CASE WHEN chromosome IN ('X', 'x') THEN 23 WHEN chromosome IN ('Y', 'y') THEN 24 WHEN chromosome = 'MT' THEN 25
       ELSE try_cast(chromosome AS INTEGER) END AS chrom,
  try_cast(position_grch37 AS BIGINT) AS p37,
  try_cast(position_grch38 AS BIGINT) AS p38,
  nullif(type, '') AS type,
  nullif(nullif(function, 'NA'), '') AS function
FROM raw;

-- Anything that cannot round-trip exactly (non-canonical id strings, positions that are not plain integers or
-- that would overflow the composite key, unknown chromosomes) is kept verbatim for inspection.
COPY (
  SELECT * FROM d WHERE rsid IS NULL OR id <> CAST(rsid AS VARCHAR)
     OR (ref_rsid IS NULL AND ref_id <> '') OR (ref_rsid IS NOT NULL AND ref_id <> CAST(ref_rsid AS VARCHAR))
     OR chrom IS NULL OR chrom NOT BETWEEN 1 AND 25
     OR (position_grch37 <> 'NA' AND (p37 IS NULL OR p37 <= 0 OR p37 >= 1000000000 OR position_grch37 <> CAST(p37 AS VARCHAR)))
     OR (position_grch38 <> 'NA' AND (p38 IS NULL OR p38 <= 0 OR p38 >= 1000000000 OR position_grch38 <> CAST(p38 AS VARCHAR)))
) TO 'dbsnp_dump/anomalies.parquet' (FORMAT parquet);

CREATE OR REPLACE VIEW k AS
SELECT ord, rsid, CASE WHEN ref_rsid = rsid THEN NULL ELSE ref_rsid END AS ref_rsid, chrom,
  CASE WHEN p37 > 0 AND p37 < 1000000000 AND chrom BETWEEN 1 AND 25 THEN CAST(chrom AS BIGINT) * 1000000000 + p37 END AS key37,
  CASE WHEN p38 > 0 AND p38 < 1000000000 AND chrom BETWEEN 1 AND 25 THEN CAST(chrom AS BIGINT) * 1000000000 + p38 END AS key38,
  (ref_rsid = rsid) AS is_ref, type, function
FROM d WHERE rsid IS NOT NULL;

-- sorted tables (rsID index: partition = rsid >> 26, 67,108,864 ids per file; position index: per chromosome)
CREATE OR REPLACE TABLE krs AS SELECT rsid >> 26 AS part, rsid, ref_rsid, chrom, key37, key38, ord FROM k ORDER BY rsid, ord;
CREATE OR REPLACE TABLE kp37 AS SELECT chrom, key37 AS key, rsid, is_ref, type, function, ord FROM k WHERE key37 IS NOT NULL ORDER BY key37, ord;
CREATE OR REPLACE TABLE kp38 AS SELECT chrom, key38 AS key, rsid, is_ref, type, function, ord FROM k WHERE key38 IS NOT NULL ORDER BY key38, ord;

@for p in 0..31
COPY (SELECT rsid, ref_rsid, chrom, key37, key38, ord FROM krs WHERE part = {p} ORDER BY rsid, ord) TO 'dbsnp_dump/sorted_rsid/part={p}.parquet' (FORMAT parquet, COMPRESSION 'zstd');
@end
@for c in 1..25
COPY (SELECT key, rsid, is_ref, type, function, ord FROM kp37 WHERE chrom = {c} ORDER BY key, ord) TO 'dbsnp_dump/sorted_pos37/chrom={c}.parquet' (FORMAT parquet, COMPRESSION 'zstd');
@end
@for c in 1..25
COPY (SELECT key, rsid, is_ref, type, function, ord FROM kp38 WHERE chrom = {c} ORDER BY key, ord) TO 'dbsnp_dump/sorted_pos38/chrom={c}.parquet' (FORMAT parquet, COMPRESSION 'zstd');
@end

-- summary for the docs
COPY (
  SELECT count(*) AS docs, count(DISTINCT rsid) AS distinct_rsid,
    count(*) - count(DISTINCT rsid) AS duplicate_rsid_docs,
    count(*) FILTER (WHERE NOT is_ref) AS merged_docs,
    count(*) FILTER (WHERE key37 IS NULL) AS no_key37, count(*) FILTER (WHERE key38 IS NULL) AS no_key38,
    (SELECT count(*) FROM 'dbsnp_dump/anomalies.parquet') AS anomalies
  FROM k
) TO 'dbsnp_dump/sort_summary.csv' (HEADER);
COPY (SELECT chromosome, count(*) n FROM d GROUP BY 1 ORDER BY 2 DESC) TO 'dbsnp_dump/chromosomes.csv' (HEADER);
