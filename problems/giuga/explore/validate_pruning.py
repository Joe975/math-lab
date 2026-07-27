"""Check the pruned search against an unpruned run of the same code.

R2/R3 are the only place a wrong number can come from without leaving a trace:
an unsound prune silently discards the branch that would have produced a
*smaller* bound, and the search reports a larger one.  So the test is that
turning the pruning off does not change the answer.

Usage:
    python validate_pruning.py --xs 3,5,7,11,13,17,19,23,29,31,37,43,53 --pool 300000
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "harness"))
from giuga.admissible import Pool, search, leaf_bound_exact  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xs", type=str, default="3,5,7,11,13,17,19,23,29,31,37,43,53")
    ap.add_argument("--pool", type=int, default=300000)
    ap.add_argument("--cap", type=int, default=3000)
    a = ap.parse_args()

    ok = True
    for X in [int(x) for x in a.xs.split(",")]:
        pool = Pool(a.pool, X)
        t0 = time.time()
        b1, w1, n1 = search(pool, best=a.cap, prune=True)
        t1 = time.time() - t0
        t0 = time.time()
        b2, w2, n2 = search(pool, best=a.cap, prune=False)
        t2 = time.time() - t0
        # third opinion: recompute the winning leaf with exact Fractions
        exact = leaf_bound_exact(w1[0], X, a.pool) if w1 else None
        agree = (b1 == b2 == exact)
        ok &= agree
        print(f"X={X:>4}  pruned={b1:>6} ({n1:>7} nodes, {t1:6.2f}s)   "
              f"unpruned={b2:>6} ({n2:>7} nodes, {t2:6.2f}s)   "
              f"exact-Fraction recheck of witness={exact:>6}   "
              f"{'OK' if agree else '*** DISAGREE ***'}")
    print("all agree" if ok else "MISMATCH")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
