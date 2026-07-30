#!/bin/bash
# Collision-safe eval worker for one level: processes slices FIRST..LAST
# (inclusive, slice indices, slice size fixed at 200000), skipping complete
# files.  Partial output carries $$ so concurrent workers never share a
# temp file; the final mv is an atomic rename and slice content is
# deterministic, so a duplicate computation is only wasted work, never
# corruption.
#
# Usage: eval_worker.sh LEVEL FIRST LAST
set -u
LEVEL=$1; FIRST=$2; LAST=$3
SLICE=200000
cd "$(dirname "$0")/.."
kk=$(printf %02d "$LEVEL")
IN=data/level$kk.masks
N=$(wc -l < "$IN")
for s in $(seq "$FIRST" "$LAST"); do
    start=$((s * SLICE))
    [ "$start" -ge "$N" ] && break
    want=$SLICE
    [ $((start + want)) -gt "$N" ] && want=$((N - start))
    OUT=data/eval${kk}_$(printf %08d "$start").txt
    if [ -f "$OUT" ] && [ "$(wc -l < "$OUT")" -eq "$want" ]; then
        continue
    fi
    TMP="$OUT.partial.$$"
    explore/census eval "$IN" "$TMP" "$start" "$SLICE" 2>/dev/null || exit 1
    [ "$(wc -l < "$TMP")" -eq "$want" ] || { echo "bad count slice $s" >&2; exit 1; }
    mv "$TMP" "$OUT"
    echo "worker level $LEVEL slice $s done"
done
echo "worker level $LEVEL range $FIRST-$LAST complete"
