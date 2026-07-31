#!/usr/bin/env python3
"""Independent count of B4-orbits of k-subsets of the 40 antipodal pairs of
nonzero points of {0,+-1}^4, by Burnside's lemma.

This deliberately shares no code with census.c: the group is built from
itertools, the action is on frozensets of points, and the count comes from
cycle-index generating functions.  If census.c's expand levels report the
same per-k orbit counts, both the group action and the canonicalization are
corroborated.

Usage: burnside_check.py [KMAX]
"""
import itertools
import sys
from math import comb


def pairs():
    out = []
    for p in itertools.product((-1, 0, 1), repeat=4):
        if p == (0, 0, 0, 0):
            continue
        fz = next(x for x in p if x)
        if fz == 1:
            out.append(p)
    assert len(out) == 40
    return out


def normalize(p):
    fz = next((x for x in p if x), 0)
    return tuple(-x for x in p) if fz < 0 else p


def group_perms(reps):
    idx = {p: i for i, p in enumerate(reps)}
    perms = []
    for perm in itertools.permutations(range(4)):
        for signs in itertools.product((1, -1), repeat=4):
            img = []
            for p in reps:
                q = tuple(signs[c] * p[perm[c]] for c in range(4))
                img.append(idx[normalize(q)])
            perms.append(tuple(img))
    # The action on antipodal PAIRS is not faithful: -identity acts
    # trivially, so exactly 192 distinct permutations appear, each twice.
    # Burnside is still a sum over all 384 group elements.
    assert len(perms) == 384 and len(set(perms)) == 192
    return perms


def cycle_lengths(perm):
    seen = [False] * len(perm)
    out = []
    for i in range(len(perm)):
        if seen[i]:
            continue
        n, j = 0, i
        while not seen[j]:
            seen[j] = True
            j = perm[j]
            n += 1
        out.append(n)
    return out


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    reps = pairs()
    perms = group_perms(reps)
    # sum over g of the generating function prod (1 + x^len)
    total = [0] * 41
    for perm in perms:
        poly = [1]
        for ln in cycle_lengths(perm):
            new = poly + [0] * ln
            for i, c in enumerate(poly):
                new[i + ln] += c
            poly = new
        for k, c in enumerate(poly):
            total[k] += c
    print(" k   subsets C(40,k)   B4-orbits (Burnside)")
    for k in range(kmax + 1):
        assert total[k] % 384 == 0
        print(f"{k:2d} {comb(40, k):16d} {total[k] // 384:12d}")


if __name__ == "__main__":
    main()
