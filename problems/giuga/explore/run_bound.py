"""Run the admissible-set factor-count bound at a series of branch limits X.

Usage:
    python run_bound.py --pool 5000000 --cap 5000 --xs 3,5,7,11,13,17,19,23,29,31
"""

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "harness"))
from giuga.admissible import Pool, PoolExhausted, search, leaf_bound_exact  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", type=int, default=5_000_000)
    ap.add_argument("--cap", type=int, default=5000)
    ap.add_argument("--xs", type=str, default="3,5,7,11,13,17,19,23,29,31,37,41")
    ap.add_argument("--node-cap", type=int, default=None)
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--budget-seconds", type=float, default=1200.0)
    a = ap.parse_args()

    rows = []
    t_start = time.time()
    for X in [int(x) for x in a.xs.split(",")]:
        if time.time() - t_start > a.budget_seconds:
            print(f"# time budget spent, stopping before X={X}")
            break
        t0 = time.time()
        pool = Pool(a.pool, X)
        t_build = time.time() - t0
        t0 = time.time()
        try:
            bound, witness, nodes = search(pool, best=a.cap, node_cap=a.node_cap)
        except (PoolExhausted, TimeoutError) as e:
            print(f"X={X:>5}  ABORTED: {type(e).__name__}: {e}")
            rows.append({"X": X, "aborted": f"{type(e).__name__}: {e}"})
            continue
        dt = time.time() - t0
        T, tleaf = witness if witness else ([], None)
        row = {"X": X, "bound": bound, "nodes": nodes, "seconds": round(dt, 2),
               "build_seconds": round(t_build, 2), "witness_T": T,
               "witness_tail": tleaf, "branchable": pool.nbranch}
        rows.append(row)
        print(f"X={X:>5}  L>={bound:>6}  nodes={nodes:>9}  {dt:8.2f}s "
              f"(+{t_build:.1f}s build)  witness T={T} + {tleaf} primes above X")
        if bound >= a.cap:
            print(f"          (hit the cap {a.cap}; rerun with a larger --cap)")

    if a.out:
        with open(a.out, "w") as f:
            json.dump({"pool_limit": a.pool, "cap": a.cap, "rows": rows}, f, indent=1)
        print(f"wrote {a.out}")


if __name__ == "__main__":
    sys.exit(main())
