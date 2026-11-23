#!/usr/bin/env bash
# Build the standalone LD kernel (no libc, no JS glue) with the emsdk clang. Run from src/ld/.
# Scalar i64.popcnt only: V8's SIMD popcount is not faster (docs/CORE_REVIEW.md §2).
#   CLANG=/e/Projects/emsdk/upstream/bin/clang ./build.sh
set -euo pipefail
CLANG=${CLANG:-clang}
$CLANG --target=wasm32 -O3 -nostdlib -fno-builtin -mno-simd128 -Wl,--no-entry -Wl,--export=__heap_base \
  -Wl,--growable-table -Wl,--initial-memory=1048576 -o kernel.wasm kernel.c
ls -l kernel.wasm
