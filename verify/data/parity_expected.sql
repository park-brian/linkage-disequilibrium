-- Expected LDlink answers computed directly from the dump (Mongo semantics), for verify_parity.js.
-- Run after pipeline/dbsnp/sort.sql:  python -I tools/run_sql.py dbsnp_dump/work.duckdb verify/data/parity_expected.sql
-- (one DuckDB process per statement: DuckDB 2.0.0 segfaulted running these in one session)
-- Note: DuckDB applies USING SAMPLE to the FROM source before WHERE, so every filtered case is sampled from an
-- explicitly filtered subquery.
SET memory_limit = '32GB';
SET threads = 8;
SET temp_directory = 'dbsnp_dump/duckdb_tmp';

-- rsID inputs: random docs, merged docs, duplicated ids, ids with NA positions, and ids that do not exist
CREATE OR REPLACE TABLE dup_ids AS SELECT rsid FROM krs GROUP BY rsid HAVING count(*) > 1;
CREATE OR REPLACE TABLE s_rsid AS
SELECT DISTINCT rsid FROM (
  (SELECT rsid FROM krs USING SAMPLE reservoir(2000 ROWS) REPEATABLE (42))
  UNION ALL (SELECT rsid FROM (SELECT rsid FROM krs WHERE ref_rsid IS NOT NULL) USING SAMPLE reservoir(500 ROWS) REPEATABLE (42))
  UNION ALL (SELECT rsid FROM dup_ids USING SAMPLE reservoir(500 ROWS) REPEATABLE (42))
  UNION ALL (SELECT rsid FROM (SELECT rsid FROM krs WHERE key37 IS NULL) USING SAMPLE reservoir(300 ROWS) REPEATABLE (42))
  UNION ALL (SELECT rsid FROM (SELECT rsid FROM krs WHERE key38 IS NULL) USING SAMPLE reservoir(300 ROWS) REPEATABLE (42))
  UNION ALL (SELECT 2100208106 + i FROM range(1, 51) t(i))       -- beyond max id
  UNION ALL (SELECT i * 7 + 1 FROM range(0, 50) t(i))             -- small ids, some absent
);
CREATE OR REPLACE TABLE s_rsid_docs AS SELECT k.* FROM krs k JOIN s_rsid s ON k.rsid = s.rsid;
COPY (
  SELECT s.rsid, f.ref_rsid, f.chrom, f.key37, f.key38, f.rsid IS NOT NULL AS found, d.rsid IS NOT NULL AS duplicated
  FROM s_rsid s
  LEFT JOIN (SELECT * FROM s_rsid_docs QUALIFY row_number() OVER (PARTITION BY rsid ORDER BY ord) = 1) f ON f.rsid = s.rsid
  LEFT JOIN dup_ids d ON d.rsid = s.rsid
) TO 'dbsnp_dump/parity_rsid.json' (FORMAT json, ARRAY true);

-- chr:pos inputs per build: random positions, positions with several docs (with/without ref docs), absent positions
CREATE OR REPLACE MACRO expected_coords(g) AS TABLE
WITH picks AS (
  (SELECT key FROM query_table(g) USING SAMPLE reservoir(1500 ROWS) REPEATABLE (42))
  UNION ALL (SELECT key FROM (SELECT key FROM query_table(g) WHERE len(ids) > 1) USING SAMPLE reservoir(1000 ROWS) REPEATABLE (42))
  UNION ALL (SELECT key FROM (SELECT key FROM query_table(g) WHERE len(ids) > 1 AND coalesce(len(refs),
    0) = 0) USING SAMPLE reservoir(300 ROWS) REPEATABLE (42))
  UNION ALL (SELECT key FROM (SELECT key FROM query_table(g) WHERE coalesce(len(refs),
    0) > 1) USING SAMPLE reservoir(300 ROWS) REPEATABLE (42))
  UNION ALL (SELECT key + 1 FROM query_table(g) USING SAMPLE reservoir(200 ROWS) REPEATABLE (42))
)
SELECT p.key, x.ids, x.refs,
  CASE WHEN x.ids IS NULL THEN NULL
       WHEN len(x.ids) = 1 THEN x.ids[1]
       WHEN coalesce(len(x.refs), 0) > 1 THEN x.refs[1]
       WHEN coalesce(len(x.refs), 0) = 0 THEN x.ids[1]
       ELSE x.refs[1] END AS expected_rsid,
  (x.ids IS NOT NULL AND len(x.ids) > 1 AND coalesce(len(x.refs), 0) <> 1) AS expect_warning
FROM (SELECT DISTINCT key FROM picks) p LEFT JOIN query_table(g) x USING (key);

CREATE OR REPLACE TABLE g37 AS SELECT key, list(rsid ORDER BY ord) AS ids,
  list(rsid ORDER BY ord) FILTER (WHERE is_ref) AS refs FROM kp37 GROUP BY key;
COPY (SELECT * FROM expected_coords('g37')) TO 'dbsnp_dump/parity_pos37.json' (FORMAT json, ARRAY true);
CREATE OR REPLACE TABLE g38 AS SELECT key, list(rsid ORDER BY ord) AS ids,
  list(rsid ORDER BY ord) FILTER (WHERE is_ref) AS refs FROM kp38 GROUP BY key;
COPY (SELECT * FROM expected_coords('g38')) TO 'dbsnp_dump/parity_pos38.json' (FORMAT json, ARRAY true);
