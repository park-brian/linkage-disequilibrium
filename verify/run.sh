#!/usr/bin/env bash
# The regression suite: unit tests, data verifiers, and (with `ldlink`) every module against a local LDlink
# ($LDLINK_URL; tools/ldlink-local/README.md). Each step logs to $LD_DATA/verify/<step>.log; the last line of each log
# is collected in $LD_DATA/verify/summary.txt.
# usage: verify/run.sh [data|ldlink|all]     (npm run verify = data, npm run verify:ldlink = ldlink)
set -uo pipefail
. "$(dirname "$0")/../tools/config.sh"
WHAT=${1:-data}
V="$REPO/verify" L="$LD_DATA/verify"
G="$LD_VCF/GRCh37" REL=${LD_RELEASE#$LD_DATA/}
FAILED=""
mkdir -p "$L"; cd "$LD_DATA"; [ "$WHAT" = ldlink ] || : > "$L/summary.txt"
step() { # name cmd...
  local name=$1; shift; local t0=$(date +%s)
  "$@" > "$L/$name.log" 2>&1; local s=$?; [ $s -eq 0 ] || FAILED="$FAILED $name"
  printf '%-18s exit=%s %5ss  %s\n' "$name" "$s" "$(( $(date +%s) - t0 ))" "$(tail -1 "$L/$name.log" | cut -c1-200)" |
    tee -a "$L/summary.txt"
}
vcf() { echo "$G/ALL.chr$1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz"; }
if [ "$WHAT" != ldlink ]; then
  # from the repo (node --test globs are cwd-relative); zero tests run is a failure, not a pass
  step tests          bash -c 'cd "$REPO" && node --test "test/*.test.js" | tee /dev/stderr | grep -qE "^ℹ tests [1-9]"'
  step verify_planner node "$V/data/planner.js" "$REL"
  step verify_parity  node "$V/data/parity.js" "$REL"
  step verify_annot   node "$V/data/annotations.js" "$REL"
  step kernel_chr1    node "$V/data/kernel.js" "$REL/genotypes/grch37/chr1" "$(vcf 1)" "$(vcf 1).tbi" 1 207000000 "verify/kernel_chr1.json"
  step kernel_chrX    node "$V/data/kernel.js" "$REL/genotypes/grch37/chrX" "$(vcf X)" "$(vcf X).tbi" X 100100000 "verify/kernel_chrX.json"
  step verify_select  node "$V/data/select.js" "$(vcf 22)" "$REL/genotypes/grch37/chr22" phaseB/ncbi/37_chr22.parquet \
    "verify/select.json" 20000
  step ldsc_munge     node "$V/ldsc/munge.js" ldscore_src/ldscore ldsc_runs/h2rg "$REL/ldscore"
  step ldsc_regress   node "$V/ldsc/regress.js" ldsc_runs/h2rg "$REL/ldscore"
fi
if [ "$WHAT" != data ]; then
  for m in ldmatrix ldpair ldpop ldhap ldproxy snpclip snpchip ldtrait ldexpress ldassoc sumstats ldherit ldcorrelation ldscore; do
    step "$m" node "$V/ldlink/compare.js" "$m"
  done
fi
echo "done: $L/summary.txt"
if [ -n "${FAILED:-}" ]; then echo "FAILED:$FAILED"; exit 1; fi
