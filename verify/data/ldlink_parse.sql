-- LDlink's own parser (rsjson_parse_data_mongo.py) re-run on NCBI b155 refsnp-chr22.json.bz2 vs. the dump's chr22 docs.
-- Run from $LD_DATA: duckdb -no-agent ncbi_dbsnp/b155/parse22.duckdb < <repo>/verify/data/ldlink_parse.sql
SET memory_limit = '24GB';
SET temp_directory = 'ncbi_dbsnp/duckdb_tmp_p22';
-- one thread + insertion order: row_number() OVER () then follows the parser's file order (pord is used for the order check)
SET threads = 1;
SET preserve_insertion_order = true;
ATTACH 'dbsnp_dump/work.duckdb' AS dump (READ_ONLY);

CREATE OR REPLACE TABLE p AS
SELECT row_number() OVER () - 1 AS pord, id, ref_id, chromosome, position_grch37, position_grch38, type,
  list_sort(string_split(coalesce(nullif(function, ''), 'NA'), ',')) AS fset -- same normalization as the dump side
FROM read_json('ncbi_dbsnp/b155/ldlink_parse/chr_22_filtered.json', format = 'newline_delimited',
  columns = {'id': 'VARCHAR', 'chromosome': 'VARCHAR', 'position_grch37': 'VARCHAR', 'position_grch38': 'VARCHAR',
             'function': 'VARCHAR', 'type': 'VARCHAR', 'ref_id': 'VARCHAR'});

CREATE OR REPLACE TABLE q AS
SELECT ord, id, ref_id, chromosome, position_grch37, position_grch38, type,
  list_sort(string_split(coalesce(function, 'NA'), ',')) AS fset
FROM dump.d WHERE chromosome = '22';

-- exact multiset comparison of the 7 fields (function as a set; Python set order is hash-randomized)
CREATE OR REPLACE TABLE x AS
WITH a AS (SELECT id, ref_id, position_grch37, position_grch38, type, fset, count(*) AS n FROM p GROUP BY ALL),
     b AS (SELECT id, ref_id, position_grch37, position_grch38, type, fset, count(*) AS n FROM q GROUP BY ALL)
SELECT coalesce(a.id, b.id) AS id, coalesce(a.ref_id, b.ref_id) AS ref_id, a.n AS n_parse, b.n AS n_dump,
  coalesce(a.position_grch37, b.position_grch37) AS p37, coalesce(a.position_grch38, b.position_grch38) AS p38,
  coalesce(a.type, b.type) AS type, coalesce(a.fset, b.fset) AS fset
FROM a FULL OUTER JOIN b ON a.id IS NOT DISTINCT FROM b.id AND a.ref_id IS NOT DISTINCT FROM b.ref_id AND a.position_grch37 IS NOT DISTINCT FROM b.position_grch37 AND a.position_grch38 IS NOT DISTINCT FROM b.position_grch38 AND a.type IS NOT DISTINCT FROM b.type AND a.fset IS NOT DISTINCT FROM b.fset; -- NULL-safe: USING never matches NULL to NULL

SELECT 'parse_rows' AS k, count(*) AS n FROM p
UNION ALL SELECT 'dump_rows', count(*) FROM q
UNION ALL SELECT 'identical_records', sum(least(n_parse, n_dump)) FROM x WHERE n_parse IS NOT NULL AND n_dump IS NOT NULL
UNION ALL SELECT 'parse_only_records', sum(n_parse - coalesce(n_dump, 0)) FROM x WHERE n_parse > coalesce(n_dump, 0)
UNION ALL SELECT 'dump_only_records', sum(n_dump - coalesce(n_parse, 0)) FROM x WHERE n_dump > coalesce(n_parse, 0)
UNION ALL SELECT 'dump_only_merged_copies(id<>ref_id)', sum(n_dump - coalesce(n_parse, 0)) FROM x WHERE n_dump > coalesce(n_parse,
  0) AND id <> ref_id;

-- order: does the dump keep the parser's order for the records both have?
SELECT count(*) AS pairs, sum((pord_rank = ord_rank)::INT) AS same_rank FROM (
  SELECT row_number() OVER (ORDER BY pord) AS pord_rank, row_number() OVER (ORDER BY ord) AS ord_rank
  FROM (SELECT p.pord,
    min(q.ord) AS ord FROM p JOIN q ON p.id IS NOT DISTINCT FROM q.id AND p.ref_id IS NOT DISTINCT FROM q.ref_id AND p.position_grch37 IS NOT DISTINCT FROM q.position_grch37 AND p.position_grch38 IS NOT DISTINCT FROM q.position_grch38 AND p.type IS NOT DISTINCT FROM q.type AND p.fset IS NOT DISTINCT FROM q.fset GROUP BY p.pord));

COPY (SELECT * FROM x WHERE n_parse IS DISTINCT FROM n_dump ORDER BY hash(id) LIMIT 300) TO 'ncbi_dbsnp/b155/parse22_diff_samples.csv' (HEADER);
