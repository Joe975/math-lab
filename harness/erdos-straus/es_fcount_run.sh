#!/bin/bash
# Chunked, checkpointed driver for es_fcount: f(p) for all primes
# p ≡ 1 (mod 24), p < 10^6.  Each completed chunk lands atomically in
# $MATHLAB_OUT (default ./out) as f1mod24_<lo>_<hi>.txt; rerunning skips
# chunks whose files already exist, so a restart never loses finished work.
set -euo pipefail
TOOLS="$(cd "$(dirname "$0")" && pwd)"
DATA="${MATHLAB_OUT:-$PWD/out}"
mkdir -p "$DATA"
BOUNDS=(2 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000)
for ((i = 0; i < ${#BOUNDS[@]} - 1; i++)); do
    lo=${BOUNDS[i]}; hi=${BOUNDS[i+1]}
    out="$DATA/f1mod24_${lo}_${hi}.txt"
    if [ -s "$out" ]; then echo "skip $out (exists)"; continue; fi
    echo "chunk [$lo, $hi) ..."
    t0=$(date +%s)
    "$TOOLS/es_fcount" "$lo" "$hi" 24 1 > "$out.tmp"
    mv "$out.tmp" "$out"
    echo "  done: $(wc -l < "$out") primes in $(( $(date +%s) - t0 ))s"
done
echo "ALL CHUNKS DONE"
