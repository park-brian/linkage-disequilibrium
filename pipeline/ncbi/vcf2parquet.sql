-- NCBI b155 VCF records as Parquet, VCF order kept (ord): ncbi_dbsnp/b155/vcf{37,38}.tsv (pipeline/ncbi/vcf2tsv.js on
-- GCF_000001405.25.gz / .39.gz) → vcf{37,38}.parquet, the input of pipeline/genotypes/slices.sql (allele labels).
-- usage: python -I tools/run_sql.py ncbi_dbsnp/b155/convert.duckdb pipeline/ncbi/vcf2parquet.sql
SET memory_limit = '32GB';
SET temp_directory = 'ncbi_dbsnp/duckdb_tmp';
SET preserve_insertion_order = true;
@for a in 37..38
COPY (SELECT * FROM read_csv('ncbi_dbsnp/b155/vcf{a}.tsv', delim = '\t', header = false, quote = '', escape = '',
  columns = {'ord': 'BIGINT', 'chrom': 'UTINYINT', 'pos': 'INTEGER', 'rsid': 'BIGINT', 'ref': 'VARCHAR', 'alt': 'VARCHAR',
             'n_alt': 'UTINYINT', 'vc': 'VARCHAR', 'fx': 'VARCHAR', 'start_hgvs': 'INTEGER'}))
  TO 'ncbi_dbsnp/b155/vcf{a}.parquet' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 1000000);
@end
