-- Phase B inputs for pipeline/genotypes/build.js, per chromosome (one statement per process: tools/run_sql.py).
--   phaseB/dbsnp_src/pos{37,38}_chr{N}.parquet  key{37,38}, rsid, ref_rsid, type, function — LDlink's dump docs with
--       a position 1..999,999,999 (the dbSNP index rule; drops the corrupt HGVS-concatenation positions), sorted
--       (key, rsid, dump order): duplicates of one rsID at one position keep LDlink's find_one order (DuckDB sorts are unstable)
--   phaseB/ncbi/{37,38}_chr{N}.parquet          rsid, pos, ref, alt (+ p_dump) — NCBI b155 VCF records in VCF order
-- usage: python -I tools/run_sql.py phaseB/slices.duckdb pipeline/genotypes/slices.sql
SET memory_limit = '64GB';
SET temp_directory = 'phaseB/duckdb_tmp';
SET preserve_insertion_order = true;
ATTACH 'dbsnp_dump/work.duckdb' AS dump (READ_ONLY);
CREATE OR REPLACE TABLE s37 AS SELECT CAST(chrom AS BIGINT) * 1000000000 + p37 AS key37, rsid, ref_rsid, type, function, ord FROM dump.d WHERE chrom BETWEEN 1 AND 25 AND p37 BETWEEN 1 AND 999999999 ORDER BY key37, rsid, ord;
CREATE OR REPLACE TABLE s38 AS SELECT CAST(chrom AS BIGINT) * 1000000000 + p38 AS key38, rsid, ref_rsid, type, function, ord FROM dump.d WHERE chrom BETWEEN 1 AND 25 AND p38 BETWEEN 1 AND 999999999 ORDER BY key38, rsid, ord;
CREATE OR REPLACE TABLE dp AS SELECT rsid, arg_min(p37, ord) AS p37, arg_min(p38, ord) AS p38 FROM dump.d WHERE rsid = ref_rsid GROUP BY rsid;
CREATE OR REPLACE TABLE a37 AS SELECT v.chrom, v.ord, v.rsid, v.pos, v.ref, v.alt, p.p37 AS p_dump FROM 'ncbi_dbsnp/b155/vcf37.parquet' v LEFT JOIN dp p USING (rsid) ORDER BY v.chrom, v.ord;
CREATE OR REPLACE TABLE a38 AS SELECT v.chrom, v.ord, v.rsid, v.pos, v.ref, v.alt, p.p38 AS p_dump FROM 'ncbi_dbsnp/b155/vcf38.parquet' v LEFT JOIN dp p USING (rsid) ORDER BY v.chrom, v.ord;
@for c in 1..23
COPY (SELECT key37, rsid, ref_rsid, type, function FROM s37 WHERE key37 BETWEEN {c}000000000 AND {c}999999999 ORDER BY key37, rsid, ord) TO 'phaseB/dbsnp_src/pos37_chr{c.chr}.parquet' (FORMAT parquet, COMPRESSION zstd);
COPY (SELECT rsid, pos, ref, alt, p_dump FROM a37 WHERE chrom = {c} ORDER BY ord) TO 'phaseB/ncbi/37_chr{c.chr}.parquet' (FORMAT parquet, COMPRESSION zstd);
COPY (SELECT key38, rsid, ref_rsid, type, function FROM s38 WHERE key38 BETWEEN {c}000000000 AND {c}999999999 ORDER BY key38, rsid, ord) TO 'phaseB/dbsnp_src/pos38_chr{c.chr}.parquet' (FORMAT parquet, COMPRESSION zstd);
COPY (SELECT rsid, pos, ref, alt, p_dump FROM a38 WHERE chrom = {c} ORDER BY ord) TO 'phaseB/ncbi/38_chr{c.chr}.parquet' (FORMAT parquet, COMPRESSION zstd);
@end
COPY (SELECT key37, rsid, ref_rsid, type, function FROM s37 WHERE key37 BETWEEN 24000000000 AND 24999999999 ORDER BY key37, rsid, ord) TO 'phaseB/dbsnp_src/pos37_chrY.parquet' (FORMAT parquet, COMPRESSION zstd);
COPY (SELECT rsid, pos, ref, alt, p_dump FROM a37 WHERE chrom = 24 ORDER BY ord) TO 'phaseB/ncbi/37_chrY.parquet' (FORMAT parquet, COMPRESSION zstd);
