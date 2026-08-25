"""Second BvLS control: complete the real orbit matrix from a deep seed.

Seeds orbit7.py with rows 0..N-1 of BvLS(243)'s real order-11 orbit matrix
(in the graph's own labelling, so symmetry breaking is OFF -- the lex rule
is validated separately on small profiles) and checks that the search
from there finds the real matrix among its completions. This exercises the
row-by-row search, the (S) transposition, the inner-product constraints
and the diagonal sum-of-squares on a 23-orbit fixed-point profile, which
the order-7 run at 99 also relies on.

    python problems/conway-99/explore/bvls_deep_seed.py --rows 12 [--node-budget N]
"""

import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import orbit7  # noqa: E402
import bvls_control  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=int, default=12)
    ap.add_argument("--node-budget", type=int, default=None)
    args = ap.parse_args()
    adj, perm, orbs, R = bvls_control.build()
    sizes = [len(o) for o in orbs]
    t = len(sizes)
    seed = {(i, j): R[i][j] for i in range(args.rows) for j in range(t)}
    t0 = time.time()
    found, stats = orbit7.enumerate_orbit_matrices(
        sizes, 22, 1, 2, seed=seed, symmetry=False, node_budget=args.node_budget)
    hit = any(M == R for M in found)
    print(f"seeded rows 0..{args.rows - 1}: {len(found)} completions, nodes "
          f"{stats['nodes']}, exhausted {stats['exhausted']}, rows "
          f"{stats.get('row_completions')}, {time.time() - t0:.0f}s")
    print(f"REAL MATRIX AMONG COMPLETIONS: {hit}")
    return 0 if hit else 1


if __name__ == "__main__":
    sys.exit(main())
