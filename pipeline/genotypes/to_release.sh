#!/usr/bin/env bash
# Copy a genotype build tree (pipeline/genotypes/build_all.sh: <geno>/<build>/chr<N>/{common_bits_p64k,low_bits_p64k,
# rare_sparse_p16k}.parquet) into a release as <release>/genotypes/<build>/chr<N>.{common,low,rare}.parquet, byte for
# byte. Only chromosomes whose build_report.json exists (finished builds) are copied.
# release_v4 took grch37 from phaseB/geno_vcfid (the vcf_id rebuild, docs/history/CORE_REVIEW.md) and the GRCh38 builds
# from phaseB/geno.
# usage: pipeline/genotypes/to_release.sh <geno dir> <release> [build ...]    (data-relative; default: all three builds)
set -euo pipefail
. "$(dirname "$0")/../../tools/config.sh"; cd "$LD_DATA"
SRC=$1 REL=$2; shift 2
BUILDS=("$@"); [ ${#BUILDS[@]} -eq 0 ] && BUILDS=(grch37 grch38 grch38_high_coverage)
for b in "${BUILDS[@]}"; do
  mkdir -p "$REL/genotypes/$b"
  for dir in "$SRC/$b"/chr*/; do
    c=$(basename "$dir"); [ -f "$dir/build_report.json" ] || { echo "skip $b $c (unfinished)"; continue; }
    cp "$dir/common_bits_p64k.parquet" "$REL/genotypes/$b/$c.common.parquet"
    cp "$dir/low_bits_p64k.parquet" "$REL/genotypes/$b/$c.low.parquet"
    cp "$dir/rare_sparse_p16k.parquet" "$REL/genotypes/$b/$c.rare.parquet"
    echo "$b $c"
  done
done
