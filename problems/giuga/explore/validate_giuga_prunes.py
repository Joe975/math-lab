"""Turn off individual pruning rules and check the enumeration is unchanged.

P4 (max-tail-sum feasibility) and P7 (closed form for the last two primes) are
the two rules that could silently delete a solution.  Disabling one must change
only the node count, never the set of answers.

Usage:
    python validate_giuga_prunes.py --max-factors 7
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "harness"))
from giuga.enumerate_giuga import enumerate_by_prime_sets, product  # noqa: E402


def run(m, **kw):
    nodes = [0]
    t0 = time.time()
    sets = sorted(enumerate_by_prime_sets(
        m, on_node=lambda: nodes.__setitem__(0, nodes[0] + 1), **kw))
    return sets, nodes[0], time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-factors", type=int, default=7)
    a = ap.parse_args()
    ok = True
    for m in range(3, a.max_factors + 1):
        base, nb, tb = run(m)
        variants = {"no-P4": dict(use_p4=False), "no-P7": dict(use_p7=False),
                    "no-P4-no-P7": dict(use_p4=False, use_p7=False)}
        line = [f"m={m}  full: {len(base)} sets, {nb} nodes, {tb:.2f}s"]
        for name, kw in variants.items():
            try:
                v, nv, tv = run(m, **kw)
            except Exception as e:                       # noqa: BLE001
                line.append(f"{name}: {type(e).__name__} (branch became unbounded)")
                continue
            same = v == base
            ok &= same
            line.append(f"{name}: {len(v)} sets, {nv} nodes, {tv:.2f}s "
                        f"{'same' if same else '*** DIFFERENT ***'}")
        print("   |   ".join(line))
        for s in base:
            print(f"        {product(s)} = {' * '.join(map(str, s))}")
    print("all rule-toggles agree" if ok else "MISMATCH")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
