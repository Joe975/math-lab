"""Enumerate Giuga numbers by prime-factor count, with a two-pass factoriser.

The last-two-primes step (P8) needs the full divisor list of N^2 - a.  Most of
those factor instantly; a few are hard semiprimes.  Rather than let one hard
integer set the cost of the whole run, pass 1 uses a small Pollard-Brent budget
and records what it could not split, then pass 2 retries only those with a much
larger budget.  Whatever still resists is reported as an unresolved branch, so
the claim is "exhaustive except for these listed N^2 - a", never a bare "none".

Checkpoints after every factor count.

Usage:
    python problems/giuga/explore/giuga_by_factor_count.py --m 8 \
        --out problems/giuga/data/giuga-m8.json
"""

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "harness"))
from giuga.enumerate_giuga import enumerate_by_prime_sets, product, verify_set  # noqa: E402
from giuga.factorint import divisors, factorize_full  # noqa: E402
from giuga.conditions import is_prime  # noqa: E402


def resolve_branch(x, budget):
    """Retry one unresolved branch with a bigger budget; return the sets it yields."""
    N, a, M = x["N"], x["a"], x["M"]
    last = max(x["chosen"])
    fac, cofactor = factorize_full(M, budget=budget)
    if cofactor != 1:
        return None
    out = []
    seen = set()
    for d in divisors(fac):
        for sd in (d, -d):
            num = sd + N
            if num <= 0 or num % a:
                continue
            p = num // a
            if p <= last:
                continue
            num2 = (M // d if sd > 0 else -(M // d)) + N
            if num2 <= 0 or num2 % a:
                continue
            q = num2 // a
            if q <= p or (p, q) in seen:
                continue
            if is_prime(p) and is_prime(q):
                seen.add((p, q))
                out.append(tuple(sorted(x["chosen"] + [p, q])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, nargs="+", default=[8])
    ap.add_argument("--min-prime", type=int, default=2)
    ap.add_argument("--pass1-budget", type=int, default=25_000)
    ap.add_argument("--window-cap", type=int, default=100_000)
    ap.add_argument("--pass2-budget", type=int, default=3_000_000)
    ap.add_argument("--out", type=str, default=None)
    a = ap.parse_args()

    rows = []
    for m in a.m:
        nodes, unresolved = [0], []
        t0 = time.time()
        sets = set(enumerate_by_prime_sets(
            m, min_prime=a.min_prime, use_p8=True, unresolved=unresolved,
            rho_budget=a.pass1_budget, window_cap=a.window_cap,
            on_node=lambda: nodes.__setitem__(0, nodes[0] + 1)))
        t1 = time.time()
        print(f"m={m}  pass 1: {nodes[0]} nodes, {t1 - t0:.1f}s, "
              f"{len(sets)} found, {len(unresolved)} branches unresolved", flush=True)

        still = []
        for i, x in enumerate(unresolved):
            got = resolve_branch(x, a.pass2_budget)
            if got is None:
                still.append(x)
            else:
                sets.update(got)
            if (i + 1) % 25 == 0:
                print(f"    pass 2: {i + 1}/{len(unresolved)} retried, "
                      f"{len(still)} still unresolved", flush=True)
        t2 = time.time()
        sets = sorted(sets)
        for s in sets:
            assert verify_set(s), s
        status = "EXHAUSTIVE" if not still else f"COMPLETE EXCEPT {len(still)} BRANCHES"
        print(f"m={m}  pass 2: {t2 - t1:.1f}s, {len(still)} still unresolved  "
              f"-> {len(sets)} Giuga numbers with {m} prime factors  [{status}]",
              flush=True)
        for s in sets:
            print(f"        {product(s)} = {' * '.join(map(str, s))}", flush=True)
        for x in still[:10]:
            print(f"        UNRESOLVED chosen={x['chosen']} a={x['a']} "
                  f"cofactor={x['cofactor']}", flush=True)

        rows.append({"m": m, "nodes": nodes[0],
                     "pass1_seconds": round(t1 - t0, 1),
                     "pass2_seconds": round(t2 - t1, 1),
                     "count": len(sets),
                     "sets": [list(s) for s in sets],
                     "products": [str(product(s)) for s in sets],
                     "unresolved_after_pass1": len(unresolved),
                     "unresolved_final": [
                         {"chosen": x["chosen"], "a": x["a"],
                          "M": str(x["M"]), "cofactor": str(x["cofactor"])}
                         for x in still],
                     "status": status})
        if a.out:
            with open(a.out, "w") as f:
                json.dump({"min_prime": a.min_prime,
                           "pass1_budget": a.pass1_budget,
                           "pass2_budget": a.pass2_budget, "rows": rows}, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
