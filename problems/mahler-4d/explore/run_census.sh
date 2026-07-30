#!/bin/bash
# Restartable eval driver for the {0,+-1}^4 symmetric-polytope census.
#
# Usage: run_census.sh KMIN KMAX [SLICE]
#   For each level k in KMIN..KMAX, evaluates data/level{k}.masks in slices
#   of SLICE masks (default 200000) with `census eval`, writing each slice
#   atomically to data/eval{k}_{start}.txt (via .partial + rename).  A slice
#   file that already exists with the expected line count is skipped, so the
#   run can be killed and restarted at any moment.
#
# Run from problems/mahler-4d/.  Deterministic: slice boundaries depend only
# on SLICE, and census eval is single-threaded exact enumeration.
set -u
KMIN=$1
KMAX=$2
SLICE=${3:-200000}
cd "$(dirname "$0")/.."   # problems/mahler-4d
CEN=explore/census

for k in $(seq "$KMIN" "$KMAX"); do
    kk=$(printf %02d "$k")
    IN=data/level$kk.masks
    [ -f "$IN" ] || { echo "missing $IN" >&2; exit 1; }
    N=$(wc -l < "$IN")
    start=0
    while [ "$start" -lt "$N" ]; do
        want=$SLICE
        [ $((start + want)) -gt "$N" ] && want=$((N - start))
        OUT=data/eval${kk}_$(printf %08d "$start").txt
        if [ -f "$OUT" ] && [ "$(wc -l < "$OUT")" -eq "$want" ]; then
            echo "level $k slice $start: complete, skipping"
        else
            rm -f "$OUT.partial"
            "$CEN" eval "$IN" "$OUT.partial" "$start" "$SLICE" || exit 1
            [ "$(wc -l < "$OUT.partial")" -eq "$want" ] || { echo "slice $start bad count" >&2; exit 1; }
            mv "$OUT.partial" "$OUT"
            echo "level $k slice $start: done ($want lines)"
        fi
        start=$((start + SLICE))
    done
done
echo "ALL EVAL SLICES COMPLETE k=$KMIN..$KMAX"
