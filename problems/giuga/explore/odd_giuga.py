"""Search for an odd Giuga number by exhausting small factor counts.

An odd Giuga number is an open *object*: nobody has exhibited one and nobody has
ruled one out.  Unlike a Giuga counterexample it is not known to be enormous, so
it is worth attacking directly rather than bounding.

The reciprocal sum forces the factor count: sum 1/p - 1/n = k >= 1 needs
sum 1/p > 1, and the eight smallest odd primes give only 0.998956, so an odd
Giuga number has at least 9 prime factors.

Checkpoints after every factor count, so a long run can be resumed or read
while it is still going.

Usage:
    python problems/giuga/explore/odd_giuga.py --min-factors 9 --max-factors 14 --fast \
        --out problems/giuga/data/odd-giuga-search.json
"""

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "harness"))
from giuga.enumerate_giuga import (  # noqa: E402
    SieveTooLarge, enumerate_by_prime_sets, k_upper_bound, product, verify_set,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-factors", type=int, default=9)
    ap.add_argument("--max-factors", type=int, default=14)
    ap.add_argument("--min-prime", type=int, default=3)
    ap.add_argument("--fast", action="store_true",
                    help="use the P8 factorisation closed form for the last two primes")
    ap.add_argument("--rho-budget", type=int, default=400_000)
    ap.add_argument("--budget-seconds", type=float, default=1e9)
    ap.add_argument("--out", type=str, default=None)
    a = ap.parse_args()

    rows = []
    t_start = time.time()
    for m in range(a.min_factors, a.max_factors + 1):
        if time.time() - t_start > a.budget_seconds:
            print(f"# budget spent, stopping before m={m}", flush=True)
            break
        nodes, unresolved = [0], []
        kw = {"use_p8": True, "unresolved": unresolved,
              "rho_budget": a.rho_budget} if a.fast else {}
        t0 = time.time()
        try:
            sets = sorted(enumerate_by_prime_sets(
                m, min_prime=a.min_prime,
                on_node=lambda: nodes.__setitem__(0, nodes[0] + 1), **kw))
        except SieveTooLarge as e:
            print(f"m={m:>3}  ABORTED after {time.time() - t0:.1f}s, "
                  f"{nodes[0]} nodes: {e}", flush=True)
            rows.append({"m": m, "aborted": str(e), "nodes": nodes[0]})
            continue
        dt = time.time() - t0
        for s in sets:
            assert verify_set(s), s
        status = ("EXHAUSTIVE" if not unresolved
                  else f"COMPLETE EXCEPT {len(unresolved)} UNFACTORED BRANCHES")
        print(f"m={m:>3}  k<={k_upper_bound(m, a.min_prime)}  nodes={nodes[0]:>10}  "
              f"{dt:9.2f}s  found={len(sets)}  [{status}]", flush=True)
        for s in sets:
            print(f"        FOUND {product(s)} = {' * '.join(map(str, s))}", flush=True)
        for x in unresolved[:3]:
            print(f"        unresolved: chosen={x['chosen']} "
                  f"cofactor has {len(str(x['cofactor']))} digits", flush=True)
        rows.append({
            "m": m, "nodes": nodes[0], "seconds": round(dt, 2),
            "found": [list(s) for s in sets],
            "unresolved_count": len(unresolved),
            "unresolved": [{"chosen": x["chosen"], "a": x["a"],
                            "cofactor": str(x["cofactor"])} for x in unresolved],
        })
        if a.out:
            with open(a.out, "w") as f:
                json.dump({"min_prime": a.min_prime, "fast": a.fast,
                           "rho_budget": a.rho_budget, "rows": rows}, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
