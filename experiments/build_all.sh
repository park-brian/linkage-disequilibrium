#!/usr/bin/env bash
# Phase B (ROADMAP §3): build the genotype files for every chromosome of every build with build_genotypes_v2.mjs
# (FINAL layouts, NCBI b155 allele labels). Resumable: a chromosome whose build_report.json exists is skipped.
# Jobs run largest-first, P at a time; each job's log and exit status go to out/phaseB/geno/<build>/chr<N>.log.
# usage (from experiments/): ./build_all.sh [P=3] [build ...]      builds: grch37 grch38 grch38_high_coverage
set -uo pipefail
P=${1:-3}; shift || true
BUILDS=("$@"); [ ${#BUILDS[@]} -eq 0 ] && BUILDS=(grch37 grch38 grch38_high_coverage)
G=/d/Projects/LDlink/db/1000G
PANELS=/d/Projects/LDlink/db/data/1000G_population_samples
OUT=${OUT:-out/phaseB/geno} # OUT=<dir> builds elsewhere (e.g. a rebuild next to the current files)
code() { case $1 in X) echo 23;; Y) echo 24;; *) echo $1;; esac; }
vcf() { # build chrom → VCF path
  case $1 in
    grch37) echo "$G/GRCh37/ALL.chr$2.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz";;
    grch38) echo "$G/GRCh38/ALL.chr$2.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz";;
    grch38_high_coverage) echo "$G/GRCh38_High_Coverage/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$2.filtered.shapeit2-duohmm-phased.vcf.gz";;
  esac
}
jobs=()
for b in "${BUILDS[@]}"; do
  chroms="$(seq 1 22) X"; [ "$b" = grch37 ] && chroms="$chroms Y"
  a=$([ "$b" = grch37 ] && echo 37 || echo 38)
  for c in $chroms; do
    v=$(vcf $b $c); [ -f "$v" ] || { echo "missing VCF $v"; exit 1; }
    for f in out/phaseB/dbsnp_src/pos${a}_chr$c.parquet out/phaseB/ncbi/${a}_chr$c.parquet; do [ -f "$f" ] || { echo "missing slice $f"; exit 1; }; done
    jobs+=("$(stat -c %s "$v") $b $c $a $v")
  done
done
run() { # size build chrom assembly vcf
  local b=$2 c=$3 a=$4 v=$5 dir=$OUT/$2/chr$3
  [ -f "$dir/build_report.json" ] && { echo "skip $b chr$c"; return 0; }
  mkdir -p "$dir"; rm -f "$dir"/*.parquet
  local t0=$(date +%s)
  FINAL=1 ALLELES=out/phaseB/ncbi/${a}_chr$c.parquet node --max-old-space-size=48000 build_genotypes_v2.mjs \
    "$v" "$dir" $PANELS out/phaseB/dbsnp_src/pos${a}_chr$c.parquet $(code $c) > "$OUT/$b/chr$c.log" 2>&1
  local s=$?
  echo "$b chr$c exit=$s $(( $(date +%s) - t0 ))s" | tee -a "$OUT/progress.log"
  return $s
}
export -f run code; export OUT PANELS
mkdir -p "$OUT"; for b in "${BUILDS[@]}"; do mkdir -p "$OUT/$b"; done
printf '%s\n' "${jobs[@]}" | sort -rn | xargs -P "$P" -L 1 bash -c 'run "$@"' _
echo "failed: $(grep -c 'exit=[1-9]' "$OUT/progress.log" 2>/dev/null || echo 0)"
