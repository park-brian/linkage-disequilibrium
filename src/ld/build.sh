#!/usr/bin/env bash
# Build the standalone LD kernel (no libc, no JS glue) with the emsdk clang. npm run build:kernel (any directory).
# Scalar i64.popcnt only: V8's SIMD popcount is not faster (docs/CORE_REVIEW.md §2).
#   CLANG=<emsdk>/upstream/bin/clang npm run build:kernel
set -euo pipefail
cd "$(dirname "$0")"
CLANG=${CLANG:-clang}
$CLANG --target=wasm32 -O3 -nostdlib -fno-builtin -mno-simd128 -Wl,--no-entry -Wl,--export=__heap_base \
  -Wl,--growable-table -Wl,--initial-memory=1048576 -o kernel.wasm kernel.c
ls -l kernel.wasm
