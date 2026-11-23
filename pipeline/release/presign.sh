#!/usr/bin/env bash
# Presign every release file in the private bucket (GET, 12 h) → presigned.json for bench/net.mjs.
# usage: pipeline/release/presign.sh <localReleaseRoot> <s3 prefix> <profile> > $LD_DATA/presigned.json
set -euo pipefail
root=$1; prefix=$2; profile=$3
echo "{"
first=1
while IFS= read -r f; do
  rel=${f#"$root"/}
  url=$(aws --profile "$profile" --region us-east-1 s3 presign "$prefix/$rel" --expires-in 43200)
  [ $first -eq 1 ] && first=0 || echo ","
  printf '  "%s": "%s"' "$rel" "$url"
done < <(find "$root" -type f \( -name '*.parquet' -o -name 'manifest.json' \) | sort)
echo; echo "}"
