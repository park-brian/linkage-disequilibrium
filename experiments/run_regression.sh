#!/usr/bin/env bash
# The safety net for changes to src/: unit tests, the data verifiers, and every comparison against the local LDlink
# (../ldlink-local must be up for the compare_* steps). Each step logs to out/regression/<step>.log; the last line of
# each log (its totals) is collected into out/regression/summary.txt.
# usage (from experiments/): ./run_regression.sh [ROOT=out/release_v4]
set -uo pipefail
export ROOT=${1:-out/release_v4}
G=/d/Projects/LDlink/db/1000G/GRCh37 C=out/ldlink_compare L=out/regression
mkdir -p "$L"; : > "$L/summary.txt"
step() { # name cmd...
  local name=$1; shift; local t0=$(date +%s)
  "$@" > "$L/$name.log" 2>&1; local s=$?
  printf '%-18s exit=%s %5ss  %s\n' "$name" "$s" "$(( $(date +%s) - t0 ))" "$(tail -1 "$L/$name.log" | cut -c1-200)" | tee -a "$L/summary.txt"
}
step tests          node --test test/phaseA.test.mjs
step verify_planner node verify_planner.mjs "$ROOT"
step verify_parity  node verify_parity.mjs "$ROOT"
step verify_annot   node verify_annotations.mjs "$ROOT"
step kernel_check   node kernel_check.mjs "$ROOT/genotypes/grch37/chr1" "$G/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz" \
                      "$G/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz.tbi" 1 207000000 "$L/kernel_check.json"
step verify_select  node verify_select.mjs "$G/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz" \
                      out/phaseB/geno_vcfid/grch37/chr22 out/phaseB/ncbi/37_chr22.parquet "$L/verify_select.json" 20000
step ldsc_munge     node ldsc/verify_munge.mjs out/ldscore_src/ldscore out/ldsc_runs/h2rg "$ROOT/ldscore"
step ldsc_regress   node ldsc/verify_regress.mjs out/ldsc_runs/h2rg "$ROOT/ldscore"
step ldmatrix       node compare_ldmatrix.mjs $C/ldmatrix_cases.json $C/ldmatrix_report.json
step ldpair         node compare_ldpair.mjs $C/ldmatrix_cases.json $C/ldpair_report.json ldpair
step ldpop          node compare_ldpair.mjs $C/ldmatrix_cases.json $C/ldpop_report.json ldpop
step ldhap          node compare_ldhap.mjs $C/ldmatrix_cases.json $C/ldhap_report.json
step ldproxy        node compare_ldproxy.mjs $C/ldmatrix_cases.json $C/ldproxy_report.json
step snpclip        node compare_snpclip.mjs $C/ldmatrix_cases.json $C/snpclip_report.json
step snpchip        node compare_snpchip.mjs $C/ldmatrix_cases.json $C/snpchip_report.json
step ldtrait        node compare_ldtrait.mjs $C/ldmatrix_cases.json $C/ldtrait_report.json
step ldexpress      node compare_ldexpress.mjs $C/ldmatrix_cases.json $C/ldexpress_report.json
step ldassoc        node compare_ldassoc.mjs $C/ldassoc_report.json
echo "done: $L/summary.txt"
