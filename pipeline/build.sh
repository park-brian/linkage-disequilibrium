#!/usr/bin/env bash
# The runbook: a release built from its sources, stage by stage, in order. Each stage is skipped when its output exists
# (delete the output to redo it), so a failed run resumes where it stopped. Paths are data-relative ($LD_DATA).
#   sources: LDlink's 2023 mongodump (S3, needs $LD_AWS_ENV credentials), the 1000 Genomes VCFs ($LD_VCF), NCBI dbSNP
#   b155 (pipeline/ncbi/fetch.sh), LDlink's ldsc panels (ldscore_src/ldscore) and population lists ($LD_PANELS).
# stage        output (skip test)                        what
#   dump         dbsnp_dump/dbsnp_dump.tsv                 the dbsnp collection, streamed from S3 (pipeline/dbsnp/bson2tsv.js)
#   dbsnp_sort   dbsnp_dump/sorted_pos38/chrom=25.parquet  DuckDB external sort into intermediates (pipeline/dbsnp/sort.sql)
#   dup_ids      dbsnp_src/dup_ids.parquet                 rsIDs with several docs (pipeline/dbsnp/dup_ids.sql)
#   dbsnp        <release>/dbsnp/duplicated_rsids.parquet  the dbSNP index files (pipeline/dbsnp/build.js, dup_rsids.js)
#   ncbi         ncbi_dbsnp/b155/vcf38.parquet             NCBI VCFs → TSV → Parquet (pipeline/ncbi/*)
#   slices       phaseB/ncbi/37_chrY.parquet               per-chromosome inputs of the genotype build (slices.sql)
#   genotypes    <release>/genotypes/grch37/chr22.common.parquet (or with CHROMS, its first chromosome)
#                                                          build_all.sh into $OUT (default phaseB/geno), check.js, then
#                                                          to_release.sh
#   annotations  <release>/annotations/grch37/recomb.parquet  fetch.sh → annotations.sql → build.js
#   ldscore      <release>/ldscore/w_hm3.snplist.parquet   LDlink's ldsc panels (pipeline/ldscore/convert_ref.js)
#   extras       <release>/annotations/panels.parquet      panels, examples, ldsc counts (pipeline/release/extras.js)
#   manifest     <release>/manifest.json                   always rebuilt when run (pipeline/release/manifest.js)
# Environment: CHROMS="22 X" and BUILDS="grch37" limit the genotype stage (a scratch release, e.g. the gate run of
# docs/history: one chromosome rebuilt and compared byte for byte); P (default 3) genotype jobs at a time.
# usage: pipeline/build.sh <release (data-relative)> [stage ...]       no stages: all of them, in order
set -euo pipefail
. "$(dirname "$0")/../tools/config.sh"; cd "$LD_DATA"
REL=${1:?usage: pipeline/build.sh <release> [stage ...]}; shift
STAGES=("$@"); ALL=(dump dbsnp_sort dup_ids dbsnp ncbi slices genotypes annotations ldscore extras manifest)
[ ${#STAGES[@]} -eq 0 ] && STAGES=("${ALL[@]}")
P=${P:-3} OUT=${OUT:-phaseB/geno}
read -r -a GBUILDS <<< "${BUILDS:-grch37 grch38 grch38_high_coverage}"
S3=s3://nci-cbiit-dceg-data-nonprod/ldlink
sql() { python -I "$REPO/tools/run_sql.py" "$1" "$REPO/$2"; }
done_if() { [ -e "$1" ] && { echo "== $stage: done ($1 exists)"; return 0; }; echo "== $stage"; return 1; }
aws_env() { set -a; source "${LD_AWS_ENV:-$REPO/../.env}"; set +a; }
mkdir -p "$REL"

for stage in "${STAGES[@]}"; do
  case $stage in
    dump) done_if dbsnp_dump/dbsnp_dump.tsv || { mkdir -p dbsnp_dump; aws_env
      aws s3 cp "$S3/deploy/mongodb/dbsnp.bson.gz" - | node "$REPO/pipeline/dbsnp/bson2tsv.js" > dbsnp_dump/dbsnp_dump.tsv.part
      mv dbsnp_dump/dbsnp_dump.tsv.part dbsnp_dump/dbsnp_dump.tsv; } ;;
    dbsnp_sort) done_if dbsnp_dump/sorted_pos38/chrom=25.parquet || sql dbsnp_dump/work.duckdb pipeline/dbsnp/sort.sql ;;
    dup_ids) done_if dbsnp_src/dup_ids.parquet || { mkdir -p dbsnp_src; sql dbsnp_dump/work.duckdb pipeline/dbsnp/dup_ids.sql; } ;;
    dbsnp) done_if "$REL/dbsnp/duplicated_rsids.parquet" || {
      for k in rsid pos37 pos38; do node --max-old-space-size=16000 "$REPO/pipeline/dbsnp/build.js" $k "$REL"; done
      node "$REPO/pipeline/dbsnp/dup_rsids.js" "$REL"; } ;;
    ncbi) done_if ncbi_dbsnp/b155/vcf38.parquet || { bash "$REPO/pipeline/ncbi/fetch.sh"
      for a in 37:25 38:39; do # assembly : GCF_000001405 version
        tsv=ncbi_dbsnp/b155/vcf${a%:*}.tsv
        [ -e "$tsv" ] || node "$REPO/pipeline/ncbi/vcf2tsv.js" "ncbi_dbsnp/b155/GCF_000001405.${a#*:}.gz" "$tsv"
      done
      sql ncbi_dbsnp/b155/convert.duckdb pipeline/ncbi/vcf2parquet.sql; } ;;
    slices) done_if phaseB/ncbi/37_chrY.parquet || {
      mkdir -p phaseB/dbsnp_src phaseB/ncbi; sql phaseB/slices.duckdb pipeline/genotypes/slices.sql; } ;;
    genotypes) first=${CHROMS:-22}; first=${first%% *}
      done_if "$REL/genotypes/${GBUILDS[0]}/chr$first.common.parquet" || {
      OUT=$OUT CHROMS=${CHROMS:-} bash "$REPO/pipeline/genotypes/build_all.sh" "$P" "${GBUILDS[@]}"
      node "$REPO/pipeline/release/check.js" "$OUT" "$OUT/check_release.json"
      bash "$REPO/pipeline/genotypes/to_release.sh" "$OUT" "$REL" "${GBUILDS[@]}"; } ;;
    annotations) done_if "$REL/annotations/grch37/recomb.parquet" || { bash "$REPO/pipeline/annotations/fetch.sh"
      sql dbsnp_dump/work.duckdb pipeline/annotations/annotations.sql; node "$REPO/pipeline/annotations/build.js" "$REL"; } ;;
    ldscore) done_if "$REL/ldscore/w_hm3.snplist.parquet" || node "$REPO/pipeline/ldscore/convert_ref.js" ldscore_src/ldscore "$REL" ;;
    extras) done_if "$REL/annotations/panels.parquet" || node "$REPO/pipeline/release/extras.js" "$REL" ;;
    manifest) echo "== manifest"; node "$REPO/pipeline/release/manifest.js" "$REL" ;;
    *) echo "unknown stage $stage (stages: ${ALL[*]})"; exit 2 ;;
  esac
done
