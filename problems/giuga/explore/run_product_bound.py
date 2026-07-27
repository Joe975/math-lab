"""Lower bound on the *size* of a composite solution, not just its factor count.

Same split search, different objective: minimise the product of the admissible
set rather than its cardinality.  Reports decimal digits.

Usage:
    python run_product_bound.py --pool 3000000 --cap 8000 --xs 3,13,31,53,73,101,151,199
"""

import argparse
import json
import os
import sys
import time

sys.set_int_max_str_digits(200000)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "harness"))
from giuga.admissible import Pool, search_product  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", type=int, default=3_000_000)
    ap.add_argument("--cap", type=int, default=8000)
    ap.add_argument("--xs", type=str, default="3,13,31,53,73,101,151,199")
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--budget-seconds", type=float, default=1200.0)
    a = ap.parse_args()

    rows = []
    t_start = time.time()
    for X in [int(x) for x in a.xs.split(",")]:
        if time.time() - t_start > a.budget_seconds:
            print(f"# time budget spent, stopping before X={X}")
            break
        pool = Pool(a.pool, X)
        t0 = time.time()
        best, witness, nodes = search_product(pool, cap=a.cap)
        dt = time.time() - t0
        digits = len(str(best))
        T, t = witness
        rows.append({"X": X, "digits": digits, "nodes": nodes,
                     "seconds": round(dt, 2), "cap_hits": pool.cap_hits,
                     "witness_T_len": len(T), "witness_tail": t})
        print(f"X={X:>5}  n >= 10^{digits - 1} ({digits} digits)  nodes={nodes:>8}  "
              f"{dt:8.2f}s  cap_hits={pool.cap_hits}  |T|={len(T)} tail={t}")
    if a.out:
        with open(a.out, "w") as f:
            json.dump({"pool_limit": a.pool, "cap": a.cap, "rows": rows}, f, indent=1)
        print(f"wrote {a.out}")


if __name__ == "__main__":
    sys.exit(main())
