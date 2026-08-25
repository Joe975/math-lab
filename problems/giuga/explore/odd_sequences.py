"""Enumerate ODD Giuga sequences, and check whether any is all-prime.

This is the widest handle on the open odd-Giuga-number question. An odd Giuga
number is exactly an all-prime odd Giuga sequence, and odd sequences are far
more plentiful than odd prime sets because composite odd entries (9, 15, 21, ...)
help the reciprocal sum reach 1 with fewer terms. So enumerating odd sequences
and filtering for all-prime is a strictly better-covered search than enumerating
odd prime sets directly -- if an odd Giuga number of this length exists, it
appears here.

Checkpoints after every length.

Usage:
    python problems/giuga/explore/odd_sequences.py --min-len 7 --max-len 10 \
        --out problems/giuga/data/odd-sequences.json
"""

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "harness"))
from giuga.sequences import (  # noqa: E402
    all_prime, enumerate_sequences, k_upper_bound, verify_sequence,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-len", type=int, default=7)
    ap.add_argument("--max-len", type=int, default=10)
    ap.add_argument("--min-entry", type=int, default=3)
    ap.add_argument("--odd", action="store_true", default=True)
    ap.add_argument("--all-integers", dest="odd", action="store_false")
    ap.add_argument("--rho-budget", type=int, default=400_000)
    ap.add_argument("--window-cap", type=int, default=100_000)
    ap.add_argument("--out", type=str, default=None)
    a = ap.parse_args()

    rows = []
    for m in range(a.min_len, a.max_len + 1):
        nodes, unresolved = [0], []
        t0 = time.time()
        seqs = sorted(enumerate_sequences(
            m, min_entry=a.min_entry, odd_only=a.odd, unresolved=unresolved,
            rho_budget=a.rho_budget, window_cap=a.window_cap,
            on_node=lambda: nodes.__setitem__(0, nodes[0] + 1)))
        dt = time.time() - t0
        for s in seqs:
            assert verify_sequence(s), s
        primes = [s for s in seqs if all_prime(s)]
        status = ("EXHAUSTIVE" if not unresolved
                  else f"COMPLETE EXCEPT {len(unresolved)} UNFACTORED")
        print(f"len={m:>3}  k<={k_upper_bound(m, a.min_entry, a.odd)}  "
              f"nodes={nodes[0]:>10}  {dt:9.2f}s  sequences={len(seqs)}  "
              f"all-prime={len(primes)}  [{status}]", flush=True)
        for s in primes:
            n = 1
            for x in s:
                n *= x
            print(f"        *** ODD GIUGA NUMBER {n} = {' * '.join(map(str, s))}",
                  flush=True)
        for s in seqs[:8]:
            print(f"        {s}", flush=True)
        rows.append({"len": m, "nodes": nodes[0], "seconds": round(dt, 2),
                     "count": len(seqs),
                     "sequences": [list(s) for s in seqs],
                     "all_prime": [list(s) for s in primes],
                     "unresolved_count": len(unresolved)})
        if a.out:
            with open(a.out, "w") as f:
                json.dump({"odd_only": a.odd, "min_entry": a.min_entry,
                           "rows": rows}, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
