#!/usr/bin/env python3
"""uc_search.py -- search union-closed families minimizing max element frequency.

Union-closed sets conjecture (Frankl): every finite union-closed family
F != {emptyset} has an element in >= 1/2 of its sets.  This tool searches for
families minimizing (max element frequency) / |F|:

  * n <= 4 : EXHAUSTIVE over all 2^(2^n) subfamilies of 2^[n] (filtered for
             union-closure).  Ground truth for small n.
  * n = 5,6: heuristic -- (i) sweep of all generator sets of size <= 3
             closed under union, (ii) random generator sets closed under
             union, (iii) hill-climbing local search (moves: remove a set
             whose removal preserves closure / add a random set and re-close).

Any family found with max frequency < 1/2 would be a counterexample to the
conjecture; such a find is re-verified from scratch (independent closure
check, exact Fraction arithmetic, closure fixed-point check) and loudly
flagged.  Expected result if the conjecture holds: minimum = 1/2 exactly
(achieved e.g. by {emptyset, {1}} or any power set).

Sets are bitmasks in [0, 2^n).  Family = python set of masks.

Usage: python3 uc_search.py [--nmax 6] [--seed 0] [--restarts 60] [--iters 400]
Re-runnable and deterministic for a fixed seed.
"""

import argparse
import random
from fractions import Fraction
from itertools import combinations


# ----------------------------------------------------------------- primitives

def is_union_closed(fam):
    """Independent O(|F|^2) verification."""
    ms = sorted(fam)
    s = set(fam)
    for i, a in enumerate(ms):
        for b in ms[i:]:
            if (a | b) not in s:
                return False
    return True


def closure(gen):
    """Union-closure of a set of masks (fixed point of pairwise OR)."""
    s = set(gen)
    while True:
        new = set()
        ls = list(s)
        for i, a in enumerate(ls):
            for b in ls[i:]:
                u = a | b
                if u not in s:
                    new.add(u)
        if not new:
            return s
        s |= new


def max_freq(fam, n):
    """(max_i #{S in fam : i in S}) / |fam| as an exact Fraction."""
    counts = [0] * n
    for mk in fam:
        for i in range(n):
            if (mk >> i) & 1:
                counts[i] += 1
    return Fraction(max(counts), len(fam))


def degenerate(fam):
    """Families the conjecture does not talk about: empty, or {emptyset}."""
    return len(fam) == 0 or set(fam) == {0}


def triple_check(fam, n):
    """Independent re-verification of a candidate family.  Returns the
    verified Fraction objective, or raises AssertionError if the family is
    not actually union-closed (i.e. the 'candidate' was bogus)."""
    fam = set(fam)
    assert not degenerate(fam), "degenerate family"
    # check 1: pairwise closure
    assert is_union_closed(fam), "NOT union-closed (pairwise check failed)"
    # check 2: closure fixed point
    assert closure(fam) == fam, "NOT union-closed (closure grew the family)"
    # check 3: exact rational objective, recomputed per element
    m = len(fam)
    fr = max(Fraction(sum(1 for S in fam if (S >> i) & 1), m)
             for i in range(n))
    assert fr == max_freq(fam, n)
    return fr


# ----------------------------------------------------------------- exhaustive

def exhaustive(n):
    """Minimize over ALL union-closed families on ground set [n] (n <= 4).
    Returns (best_fraction, best_family, num_uc_families, num_achieving_half).
    """
    N = 1 << n
    best_fr, best_fam = None, None
    n_uc = 0
    n_half = 0
    half = Fraction(1, 2)
    for code in range(1, 1 << N):
        masks = [m for m in range(N) if (code >> m) & 1]
        if masks == [0]:
            continue  # {emptyset}
        s = set(masks)
        ok = True
        for i, a in enumerate(masks):
            for b in masks[i:]:
                if (a | b) not in s:
                    ok = False
                    break
            if not ok:
                break
        if not ok:
            continue
        n_uc += 1
        fr = max_freq(masks, n)
        if fr == half:
            n_half += 1
        if best_fr is None or fr < best_fr:
            best_fr, best_fam = fr, set(masks)
    return best_fr, best_fam, n_uc, n_half


# ----------------------------------------------------------------- heuristics

def generator_sweep(n, max_gens=3):
    """Closure of every generator set of size <= max_gens (deterministic)."""
    N = 1 << n
    best_fr, best_fam = None, None
    for k in range(1, max_gens + 1):
        for gens in combinations(range(N), k):
            fam = closure(gens)
            if degenerate(fam):
                continue
            fr = max_freq(fam, n)
            if best_fr is None or fr < best_fr:
                best_fr, best_fam = fr, fam
    return best_fr, best_fam


def removable(fam, S):
    """Can S be deleted from union-closed fam leaving it union-closed?
    True iff S is not the union of two other members."""
    rest = [x for x in fam if x != S]
    for i, a in enumerate(rest):
        for b in rest[i:]:
            if (a | b) == S:
                return False
    return True


def local_search(n, rng, iters=400, cap=400, start=None):
    """Hill-climb (with small sideways/uphill tolerance) on max_freq."""
    N = 1 << n
    if start is None:
        k = rng.randint(1, 7)
        cur = closure(rng.sample(range(N), k))
    else:
        cur = set(start)
    if degenerate(cur):
        cur = {0, 1}
    cur_fr = max_freq(cur, n)
    best_fr, best_fam = cur_fr, set(cur)
    for _ in range(iters):
        cand = None
        if rng.random() < 0.6 and len(cur) > 2:
            S = rng.choice(sorted(cur))
            if removable(cur, S):
                cand = set(x for x in cur if x != S)
                if degenerate(cand):
                    cand = None
        elif len(cur) < cap:
            cand = closure(cur | {rng.randrange(N)})
            if len(cand) > cap or degenerate(cand):
                cand = None
        if cand is None:
            continue
        fr = max_freq(cand, n)
        if fr <= cur_fr or rng.random() < 0.03:
            cur, cur_fr = cand, fr
            if fr < best_fr:
                best_fr, best_fam = fr, set(cand)
    return best_fr, best_fam


def heuristic(n, rng, restarts=60, iters=400):
    best_fr, best_fam = generator_sweep(n)
    # structured seeds: power set of first k coords; {0, singleton}; a slice's closure
    seeds = []
    for k in range(1, min(n, 4) + 1):
        seeds.append(set(range(1 << k)))                       # power set 2^[k]
    seeds.append({0, 1})
    slice_gens = [sum(1 << i for i in c) for c in combinations(range(n), 2)]
    seeds.append(closure(slice_gens))                          # closure of 2-slice
    for sd in seeds:
        fr = max_freq(sd, n)
        if fr < best_fr:
            best_fr, best_fam = fr, set(sd)
        fr2, fam2 = local_search(n, rng, iters=iters, start=sd)
        if fr2 < best_fr:
            best_fr, best_fam = fr2, fam2
    for _ in range(restarts):
        fr, fam = local_search(n, rng, iters=iters)
        if fr < best_fr:
            best_fr, best_fam = fr, fam
    return best_fr, best_fam


# ----------------------------------------------------------------------- main

def fmt_family(fam, n, limit=12):
    def one(mk):
        return "{" + ",".join(str(i + 1) for i in range(n) if (mk >> i) & 1) + "}"
    ms = sorted(fam)
    body = " ".join(one(m) for m in ms[:limit])
    if len(ms) > limit:
        body += f" ... ({len(ms)} sets total)"
    return body


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--nmax", type=int, default=6)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--restarts", type=int, default=60)
    ap.add_argument("--iters", type=int, default=400)
    args = ap.parse_args()

    half = Fraction(1, 2)
    overall_ok = True
    print(f"# uc_search  seed={args.seed} restarts={args.restarts} iters={args.iters}")
    print("# objective: min over union-closed F (F != {0}, F nonempty) of "
          "max_i freq(i)/|F|")
    for n in range(1, args.nmax + 1):
        if n <= 4:
            fr, fam, n_uc, n_half = exhaustive(n)
            mode = f"EXHAUSTIVE ({n_uc} union-closed families; "\
                   f"{n_half} achieve exactly 1/2)"
        else:
            rng = random.Random(args.seed + n)
            fr, fam = heuristic(n, rng, restarts=args.restarts, iters=args.iters)
            mode = "heuristic (generator sweep <=3 + local search)"
        vfr = triple_check(fam, n)  # independent re-verification, always
        assert vfr == fr
        tag = "OK (>= 1/2, conjecture-consistent)"
        if vfr < half:
            overall_ok = False
            tag = ("*** BELOW 1/2 -- POTENTIAL COUNTEREXAMPLE -- "
                   "family passed independent union-closure re-checks; "
                   "DO NOT TRUST WITHOUT FURTHER INDEPENDENT VERIFICATION ***")
        print(f"n={n}: min max-frequency = {vfr} = {float(vfr):.6f}  [{mode}]")
        print(f"      witness (|F|={len(fam)}): {fmt_family(fam, n)}   {tag}")
    print("RESULT:", "no union-closed family with max frequency < 1/2 found"
          if overall_ok else "!!! candidate below 1/2 found -- see above !!!")


if __name__ == "__main__":
    main()
