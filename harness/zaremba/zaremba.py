#!/usr/bin/env python3
"""Reference implementation for Zaremba's conjecture (harness, tier 0).

Canonical continued fractions of reduced fractions m/n, 0 < m < n:
m/n = [0; a1, ..., ak] with ak >= 2, which is exactly what the Euclidean
algorithm produces and is unique. K(m/n) = max ai; z(n) = min K over
coprime m. See problems/zaremba/PROBLEM.md for the verification contract.

Subcommands:
    cf <n> <m>              print the canonical partial quotients of m/n
    verify <n> <m> [--max A]  check gcd(m,n)=1 and K(m/n) <= A (default 5)
    z <n>                   exact z(n) by exhaustive scan over m
    scan <N> [--max A]      first witness m with K <= A for every n in 2..N
    zscan <N>               exact z(n) for every n in 2..N (reference; slow)
    verify-file <path> [--max A]  verify "n m" witness lines; report failures

Standard library only. All arithmetic exact. Exit status 0 iff every check
requested passed.
"""

from __future__ import annotations

import argparse
import sys


def cf(n: int, m: int) -> list[int]:
    """Canonical partial quotients of m/n for 0 < m < n, gcd(m, n) = 1.

    Raises ValueError if the inputs are out of range or not coprime, so a
    caller cannot mistake a garbage pair for a verified witness.
    """
    if not 0 < m < n:
        raise ValueError(f"need 0 < m < n, got m={m} n={n}")
    quots = []
    p, q = n, m
    while q:
        a, r = divmod(p, q)
        quots.append(a)
        p, q = q, r
    if p != 1:
        raise ValueError(f"gcd({m}, {n}) = {p} != 1")
    return quots


def k_of(n: int, m: int) -> int:
    """K(m/n) = max canonical partial quotient."""
    return max(cf(n, m))


def z_exact(n: int) -> int:
    """z(n) by exhaustive scan over all coprime m, with early aborts.

    The last canonical quotient is always >= 2, so z(n) >= 2 for n >= 2 and
    the scan stops if 2 is reached.
    """
    if n < 2:
        return 1
    best = n  # K(1/n) = n, so this is attained
    for m in range(1, n):
        p, q, worst = n, m, 0
        while q:
            a, r = divmod(p, q)
            if a >= best:
                worst = best  # can't improve; abandon this m
                break
            if a > worst:
                worst = a
            p, q = q, r
        else:
            if p != 1:
                continue  # not coprime
            if worst < best:
                best = worst
                if best == 2:
                    break
    return best


def first_witness(n: int, bound: int) -> tuple[int, int] | None:
    """Smallest m with gcd(m,n)=1 and K(m/n) <= bound, plus tries used."""
    for tries, m in enumerate(range(1, n), start=1):
        p, q = n, m
        ok = True
        while q:
            a, r = divmod(p, q)
            if a > bound:
                ok = False
                break
            p, q = q, r
        if ok and p == 1:
            return m, tries
    return None


def main() -> int:
    top = argparse.ArgumentParser(description=__doc__)
    sub = top.add_subparsers(dest="cmd", required=True)

    p_cf = sub.add_parser("cf")
    p_cf.add_argument("n", type=int)
    p_cf.add_argument("m", type=int)

    p_verify = sub.add_parser("verify")
    p_verify.add_argument("n", type=int)
    p_verify.add_argument("m", type=int)
    p_verify.add_argument("--max", type=int, default=5, dest="bound")

    p_z = sub.add_parser("z")
    p_z.add_argument("n", type=int)

    p_scan = sub.add_parser("scan")
    p_scan.add_argument("N", type=int)
    p_scan.add_argument("--max", type=int, default=5, dest="bound")

    p_zscan = sub.add_parser("zscan")
    p_zscan.add_argument("N", type=int)

    p_vf = sub.add_parser("verify-file")
    p_vf.add_argument("path")
    p_vf.add_argument("--max", type=int, default=5, dest="bound")

    args = top.parse_args()

    if args.cmd == "cf":
        print(" ".join(map(str, cf(args.n, args.m))))
        return 0

    if args.cmd == "verify":
        quots = cf(args.n, args.m)
        worst = max(quots)
        ok = worst <= args.bound
        print(f"{args.m}/{args.n} = [0; {', '.join(map(str, quots))}]  "
              f"K = {worst}  {'OK' if ok else f'EXCEEDS {args.bound}'}")
        return 0 if ok else 1

    if args.cmd == "z":
        print(z_exact(args.n))
        return 0

    if args.cmd == "scan":
        failures = 0
        for n in range(2, args.N + 1):
            hit = first_witness(n, args.bound)
            if hit is None:
                print(f"FAIL n={n}: no witness with K <= {args.bound}")
                failures += 1
            else:
                m, tries = hit
                print(f"{n} {m} {tries}")
        return 0 if failures == 0 else 1

    if args.cmd == "zscan":
        for n in range(2, args.N + 1):
            print(f"{n} {z_exact(n)}")
        return 0

    if args.cmd == "verify-file":
        bad = total = 0
        with open(args.path, encoding="utf-8") as fh:
            for line in fh:
                parts = line.split()
                if not parts or parts[0].startswith("#"):
                    continue
                n, m = int(parts[0]), int(parts[1])
                if (n, m) == (1, 0):
                    continue  # conventional z(1) = 1 entry; no fraction exists
                total += 1
                try:
                    worst = k_of(n, m)
                except ValueError as exc:
                    print(f"BAD n={n} m={m}: {exc}")
                    bad += 1
                    continue
                if worst > args.bound:
                    print(f"BAD n={n} m={m}: K = {worst} > {args.bound}")
                    bad += 1
        print(f"checked {total} witnesses, {bad} bad")
        return 0 if bad == 0 else 1

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
