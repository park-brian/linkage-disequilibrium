#!/usr/bin/env bash
# The regression suite: unit tests, data verifiers, and (with `ldlink`) every module against a local LDlink
# ($LDLINK_URL; tools/ldlink-local/README.md). Each step logs to $LD_DATA/verify/<step>.log; the last line of each log
# is collected in $LD_DATA/verify/summary.txt.
# usage: verify/run.sh [data|ldlink|all]     (npm run verify = data, npm run verify:ldlink = ldlink)
set -uo pipefail
. "$(dirname "$0")/../tools/config.sh"
WHAT=${1:-data}
V="$REPO/verify" L="$LD_DATA/verify" C="$REPO/verify/ldlink/cases/ldmatrix_cases.json" R="$LD_DATA/verify/ldlink"
G="$LD_VCF/GRCh37" REL=${LD_RELEASE#$LD_DATA/}
mkdir -p "$L" "$R"; cd "$LD_DATA"; [ "$WHAT" = ldlink ] || : > "$L/summary.txt"
step() { # name cmd...
  local name=$1; shift; local t0=$(date +%s)
  "$@" > "$L/$name.log" 2>&1; local s=$?
  printf '%-18s exit=%s %5ss  %s\n' "$name" "$s" "$(( $(date +%s) - t0 ))" "$(tail -1 "$L/$name.log" | cut -c1-200)" | tee -a "$L/summary.txt"
}
vcf() { echo "$G/ALL.chr$1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz"; }
if [ "$WHAT" != ldlink ]; then
  # from the repo (node --test globs are cwd-relative); zero tests run is a failure, not a pass
  step tests          bash -c 'cd "$REPO" && node --test "test/*.test.mjs" | tee /dev/stderr | grep -qE "^ℹ tests [1-9]"'
  step verify_planner node "$V/data/planner.mjs" "$REL"
  step verify_parity  node "$V/data/parity.mjs" "$REL"
  step verify_annot   node "$V/data/annotations.mjs" "$REL"
  step kernel_chr1    node "$V/data/kernel.mjs" "$REL/genotypes/grch37/chr1" "$(vcf 1)" "$(vcf 1).tbi" 1 207000000 "verify/kernel_chr1.json"
  step kernel_chrX    node "$V/data/kernel.mjs" "$REL/genotypes/grch37/chrX" "$(vcf X)" "$(vcf X).tbi" X 100100000 "verify/kernel_chrX.json"
  step verify_select  node "$V/data/select.mjs" "$(vcf 22)" phaseB/geno_vcfid/grch37/chr22 phaseB/ncbi/37_chr22.parquet "verify/select.json" 20000
  step ldsc_munge     node "$V/ldsc/munge.mjs" ldscore_src/ldscore ldsc_runs/h2rg "$REL/ldscore"
  step ldsc_regress   node "$V/ldsc/regress.mjs" ldsc_runs/h2rg "$REL/ldscore"
fi
if [ "$WHAT" != data ]; then
  step ldmatrix       node "$V/ldlink/compare_ldmatrix.mjs" "$C" "$R/ldmatrix.json"
  step ldpair         node "$V/ldlink/compare_ldpair.mjs" "$C" "$R/ldpair.json" ldpair
  step ldpop          node "$V/ldlink/compare_ldpair.mjs" "$C" "$R/ldpop.json" ldpop
  step ldhap          node "$V/ldlink/compare_ldhap.mjs" "$C" "$R/ldhap.json"
  step ldproxy        node "$V/ldlink/compare_ldproxy.mjs" "$C" "$R/ldproxy.json"
  step snpclip        node "$V/ldlink/compare_snpclip.mjs" "$C" "$R/snpclip.json"
  step snpchip        node "$V/ldlink/compare_snpchip.mjs" "$C" "$R/snpchip.json"
  step ldtrait        node "$V/ldlink/compare_ldtrait.mjs" "$C" "$R/ldtrait.json"
  step ldexpress      node "$V/ldlink/compare_ldexpress.mjs" "$C" "$R/ldexpress.json"
  step ldassoc        node "$V/ldlink/compare_ldassoc.mjs" "$R/ldassoc.json"
fi
echo "done: $L/summary.txt"
