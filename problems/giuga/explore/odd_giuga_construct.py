"""Try to CONSTRUCT an odd Giuga number, rather than bound one away.

The exhaustive search over odd prime sets stalls around 12 factors.  This goes
after longer prefixes it will never reach, by giving up completeness *between*
prefixes while staying complete *within* each one.

For a prefix p_1 < ... < p_{m-2} of odd primes with product N and
rem = 1 - sum 1/p_i = a/N (the denominator is always exactly N -- see
harness/giuga/sequences.py), the last two primes satisfy

    a*p*q - N*p - N*q + 1 = 0,   so   q = (N*p - 1) / (a*p - N),

and 1/p + 1/q = a/N + 1/(N p q) with q > p forces

    N/a  <  p  <=  2N/a.

So scanning every prime p in that window and testing whether q comes out an
integer, odd, prime and greater than p **settles the prefix completely** -- with
no factorisation anywhere.  The window has width N/a = 1/rem, so a prefix is
usable exactly when its reciprocal sum is not too close to 1.

Each prefix that gets scanned is therefore a finished piece of exhaustion; what
is incomplete is only which prefixes were sampled.  The script reports both, so
a negative result has an exact scope.

Deterministic: seeded PRNG, no wall-clock or hash dependence.

Usage:
    python problems/giuga/explore/odd_giuga_construct.py --seed 1 --trials 20000 \
        --window-cap 3000000 --out problems/giuga/data/odd-construct.json
"""

import argparse
import json
import os
import random
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "harness"))
from giuga.admissible import primes_upto  # noqa: E402
from giuga.conditions import is_prime, is_giuga_number  # noqa: E402


def scan_prefix(prefix, N, a, window_cap):
    """Every completion of ``prefix`` by two more primes, or None if too wide.

    Complete for this prefix when it returns a list.
    """
    lo = N // a
    hi = (2 * N) // a + 1
    last = prefix[-1]
    lo = max(lo, last)
    if hi - lo > window_cap:
        return None
    out = []
    for p in range(lo + 1, hi + 1):
        d = a * p - N
        if d <= 0:
            continue
        num = N * p - 1
        if num % d:
            continue
        q = num // d
        if q <= p or q % 2 == 0:
            continue
        if is_prime(p) and is_prime(q):
            out.append((p, q))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--trials", type=int, default=20000)
    ap.add_argument("--pool", type=int, default=200_000)
    ap.add_argument("--window-cap", type=int, default=3_000_000)
    ap.add_argument("--max-len", type=int, default=60)
    ap.add_argument("--skip-prob", type=float, default=0.25,
                    help="chance of skipping a candidate during greedy extension; "
                         "without this only the 2^14 small-prime subsets vary")
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--report-every", type=float, default=60.0)
    a = ap.parse_args()

    pool = [p for p in primes_upto(a.pool) if p > 2]
    rng = random.Random(a.seed)
    small = [p for p in pool if p < 60]

    seen_prefixes = set()  # (N, a) pairs already settled
    scanned = 0            # prefixes settled completely
    window_total = 0       # integers examined inside those windows
    by_len = {}
    found = []
    t0 = time.time()
    t_report = t0

    for trial in range(a.trials):
        # random subset of the small odd primes, then greedily extend upward
        base = [p for p in small if rng.random() < 0.55]
        s = sum((Fraction(1, p) for p in base), Fraction(0))
        if s >= 1 or not base:
            continue
        prefix = list(base)
        N = 1
        for p in prefix:
            N *= p
        idx = 0
        while len(prefix) < a.max_len:
            # every prefix along the way is a candidate worth settling
            rem = 1 - s
            if rem <= 0:
                break
            aa = rem.numerator
            assert rem.denominator == N, (prefix, rem)
            if (N, aa) in seen_prefixes:        # same prefix reached twice
                hits = None
            else:
                seen_prefixes.add((N, aa))
                hits = scan_prefix(prefix, N, aa, a.window_cap)
            if hits is not None:
                scanned += 1
                window_total += (2 * N) // aa - max(N // aa, prefix[-1])
                by_len[len(prefix) + 2] = by_len.get(len(prefix) + 2, 0) + 1
                for (p, q) in hits:
                    seq = tuple(sorted(prefix + [p, q]))
                    n = 1
                    for x in seq:
                        n *= x
                    if n % 2 == 1 and is_giuga_number(n):
                        found.append({"n": str(n), "primes": list(seq)})
                        print(f"*** ODD GIUGA NUMBER {n} = "
                              f"{' * '.join(map(str, seq))}", flush=True)
            # extend: smallest unused prime that keeps the sum below 1
            nxt = None
            while idx < len(pool):
                p = pool[idx]
                idx += 1
                if p <= prefix[-1] or p in prefix:
                    continue
                if s + Fraction(1, p) < 1:
                    if rng.random() < a.skip_prob:
                        continue
                    nxt = p
                    break
            if nxt is None:
                break
            prefix.append(nxt)
            s += Fraction(1, nxt)
            N *= nxt

        if time.time() - t_report > a.report_every:
            t_report = time.time()
            print(f"  trial {trial + 1}/{a.trials}: {scanned} prefixes settled, "
                  f"{window_total:,} integers scanned, {len(found)} found "
                  f"({time.time() - t0:.0f}s)", flush=True)

    dt = time.time() - t0
    print(f"seed={a.seed} trials={a.trials}: {scanned} prefixes settled "
          f"completely, {window_total:,} integers scanned inside their windows, "
          f"{dt:.0f}s")
    print(f"prefix lengths settled (as full factor counts): "
          f"{dict(sorted(by_len.items()))}")
    print(f"odd Giuga numbers found: {len(found)}")
    if a.out:
        with open(a.out, "w") as f:
            json.dump({"seed": a.seed, "trials": a.trials, "pool": a.pool,
                       "window_cap": a.window_cap,
                       "prefixes_settled": scanned,
                       "integers_scanned": window_total,
                       "by_factor_count": {str(k): v for k, v in sorted(by_len.items())},
                       "found": found, "seconds": round(dt, 1)}, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
