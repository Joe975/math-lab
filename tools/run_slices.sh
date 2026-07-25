#!/bin/bash
# Restartable slice runner for the girth>=5 Erdos-Gyarfas search.
#
# Usage: run_slices.sh N MOD [JOBS]
#   Generates connected cubic girth>=5 graphs on N vertices in MOD geng
#   slices (res/mod), filters each through tools/cycle_filter, and writes
#   each slice's full output atomically to
#     attempts/erdos-gyarfas/data/g5n{N}_s{R}of{MOD}.txt
#   A slice file is only moved into place after cycle_filter emits its
#   final SUMMARY line, so any file that exists and contains SUMMARY is
#   complete and is skipped on re-run.  geng's stderr (which reports its
#   own generated-graph count, ">Z ... graphs generated") is kept beside
#   it as .gengerr for an independent generation-count cross-check.
#
# Safe to re-run after a crash or container restart: completed slices are
# never re-run; partial files (*.partial) are discarded and redone.

set -u
N=$1
MOD=$2
JOBS=${3:-4}
TOOLS="$(cd "$(dirname "$0")" && pwd)"
DATA="$TOOLS/../attempts/erdos-gyarfas/data"
mkdir -p "$DATA"

run_slice() {
    local r=$1
    local out="$DATA/g5n${N}_s${r}of${MOD}.txt"
    local err="$DATA/g5n${N}_s${r}of${MOD}.gengerr"
    if [ -f "$out" ] && grep -q '^SUMMARY' "$out"; then
        echo "slice $r/$MOD: already complete, skipping"
        return 0
    fi
    local tmp="$out.partial"
    rm -f "$tmp"
    nauty-geng -c -d3 -D3 -tf "$N" "$r/$MOD" 2>"$err" \
        | "$TOOLS/cycle_filter" > "$tmp"
    if grep -q '^SUMMARY' "$tmp"; then
        mv "$tmp" "$out"
        echo "slice $r/$MOD: done ($(grep -c '' "$out") lines)"
    else
        echo "slice $r/$MOD: FAILED (no SUMMARY; partial kept as $tmp)" >&2
        return 1
    fi
}

fail=0
active=0
for r in $(seq 0 $((MOD - 1))); do
    run_slice "$r" &
    active=$((active + 1))
    if [ "$active" -ge "$JOBS" ]; then
        wait -n || fail=1
        active=$((active - 1))
    fi
done
while [ "$active" -gt 0 ]; do
    wait -n || fail=1
    active=$((active - 1))
done

if [ "$fail" -eq 0 ]; then
    echo "ALL SLICES COMPLETE n=$N mod=$MOD"
else
    echo "SOME SLICES FAILED n=$N mod=$MOD" >&2
fi
exit $fail
