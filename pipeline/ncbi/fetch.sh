#!/bin/bash
# Download NCBI dbSNP VCFs (+tbi) and refsnp-merged JSON for b155 and b157, resumable, md5-verified.
# usage: pipeline/ncbi/fetch.sh [outDir=ncbi_dbsnp (data-relative)]   → <outDir>/b155, <outDir>/b157
set -u
. "$(dirname "$0")/../../tools/config.sh"; cd "$LD_DATA"
out=${1:-ncbi_dbsnp}; mkdir -p "$out/b155" "$out/b157"
get() { # url dir
  local f=$(basename "$1")
  curl -sS -f --retry 20 --retry-delay 10 -C - -o "$2/$f" "$1" || { echo "FAIL $1"; return 1; }
  curl -sS -f -o "$2/$f.md5" "$1.md5" 2>/dev/null && (cd "$2" && want=$(awk '{print $1}' "$f.md5") && got=$(md5sum "$f" | awk '{print $1}') && [ "$want" = "$got" ] && echo "OK $2/$f" || echo "MD5MISMATCH $2/$f")
}
B=https://ftp.ncbi.nih.gov/snp
for f in VCF/GCF_000001405.25.gz VCF/GCF_000001405.39.gz JSON/refsnp-merged.json.bz2 release_notes.txt; do get $B/archive/b155/$f $out/b155; done
for f in VCF/GCF_000001405.25.gz VCF/GCF_000001405.25.gz.tbi VCF/GCF_000001405.40.gz VCF/GCF_000001405.40.gz.tbi JSON/refsnp-merged.json.bz2 release_notes.txt; do get $B/latest_release/$f $out/b157; done
echo DONE
