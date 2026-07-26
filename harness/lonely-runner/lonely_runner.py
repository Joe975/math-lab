#!/usr/bin/env python3
"""Exact maximum-loneliness (gap) tool for the Lonely Runner Conjecture,
plus a k=8 near-tight scan over bounded speed tuples.

Reduction (standard).  For k runners with pairwise distinct constant speeds
on the unit circle, subtracting one runner's speed from all (loneliness of
that runner only depends on speed differences) reduces to: one runner is
stationary at 0, the other k-1 runners have pairwise distinct nonzero real
speeds.  A standard approximation/limit argument (Bohl; see Bienia et al.)
further reduces to *integer* speeds, and negating a speed does not change
||v t||, so we may take distinct positive integers v_1 < ... < v_{k-1}.
The conjecture (for the stationary runner) then reads

    ML(v) := sup_{t in R} min_i ||v_i t|| >= 1/k,

where ||x|| is the distance from x to the nearest integer.  For k = 8 this
means: for every 7-tuple of distinct positive integers, ML(v) >= 1/8.
ML is invariant under scaling v -> c*v (substitute t -> t/c), so it is
enough to consider tuples with gcd(v_1,...,v_7) = 1.

Exact algorithm.  Let g_i(t) = ||v_i t|| and f(t) = min_i g_i(t).  Since the
v_i are integers, f has period 1; it is continuous and piecewise linear, and
every linear piece of every g_i has the form (v_i t - a) [slope +v_i] or
(b - v_i t) [slope -v_i] with integers a, b.  f attains its max on [0,1).
At a global maximizer t*, the left derivative of f is >= 0 and the right
derivative is <= 0, so there are active pieces p(t) = v_i t - a (achieving
f just left of t*) and q(t) = b - v_j t (achieving f just right of t*) with
p(t*) = q(t*) = f(t*); the case i = j (a kink of a single g_i at its peak)
is included.  Solving v_i t* - a = b - v_j t* gives

    t* = (a+b)/(v_i+v_j),  i.e.  t* in T := { m/(v_i+v_j) : i <= j, 0 <= m < v_i+v_j }.

Conversely every t in T gives a valid value f(t) <= ML(v).  Hence

    ML(v) = max_{t in T} f(t),   a finite exact computation.

At t = m/d, ||v t|| = min(v*m mod d, d - v*m mod d) / d, so everything is
integer arithmetic; results are returned as exact Fractions.

Usage:
  lonely_runner.py --validate                  run the validation suite
  lonely_runner.py --ml 1,2,3,4,5,6,7          exact ML of one tuple (+witness)
  lonely_runner.py --scan V [--jobs N] [--out FILE.jsonl]
        scan all gcd=1 7-subsets of {1..V} and record every tuple with
        ML < 1/8 + 1/200 = 13/100 ("near-tight"), with exact ML + witness.
"""

import argparse
import itertools
import json
import multiprocessing as mp
import random
import sys
import time
from fractions import Fraction
from math import gcd

# near-tight threshold 1/8 + 1/200 = 13/100 (strict <)
THR_NUM, THR_DEN = 13, 100


# ----------------------------------------------------------------------
# exact ML
# ----------------------------------------------------------------------

def candidate_denominators(speeds):
    """Denominators v_i + v_j, i <= j (i = j covers the kink peaks m/(2 v_i))."""
    return sorted({a + b for i, a in enumerate(speeds) for b in speeds[i:]})


def ml_exact(speeds):
    """Exact ML(speeds) as (Fraction value, Fraction witness_time).

    Max of f over the (superset of critical points) T described in the
    module docstring.  Pure integer arithmetic.
    """
    best_n, best_d, best_t = 0, 1, Fraction(0)
    for d in candidate_denominators(speeds):
        for m in range(1, d):
            e = d  # min_i min(r, d-r) for this t = m/d
            for v in speeds:
                r = (v * m) % d
                if r > d - r:
                    r = d - r
                if r < e:
                    e = r
                    if r == 0:
                        break
            if e * best_d > best_n * d:  # e/d > best
                best_n, best_d, best_t = e, d, Fraction(m, d)
    return Fraction(best_n, best_d), best_t


def ml_exact_paranoid(speeds):
    """Independent recomputation over a strictly larger candidate set
    (sums AND absolute differences of pairs, plus 2*v_i) -- used only for
    validation cross-checks of ml_exact's candidate-set argument."""
    dens = set(candidate_denominators(speeds))
    dens |= {abs(a - b) for a in speeds for b in speeds if a != b}
    dens |= {2 * v for v in speeds}
    dens.discard(0)
    best = Fraction(0)
    for d in sorted(dens):
        for m in range(1, d):
            e = min(min((v * m) % d, d - (v * m) % d) for v in speeds)
            if Fraction(e, d) > best:
                best = Fraction(e, d)
    return best


def ml_grid_numeric(speeds, n_grid=200_001):
    """Dense-grid numeric lower bound on ML (floating point).  f is
    Lipschitz with constant max(speeds), so ML <= grid_max + vmax/(2*n)."""
    vmax = max(speeds)
    best = 0.0
    for j in range(n_grid):
        t = j / n_grid
        e = min(abs(v * t - round(v * t)) for v in speeds)
        if e > best:
            best = e
    return best, vmax / (2 * n_grid)


# ----------------------------------------------------------------------
# scan with early bail-out
# ----------------------------------------------------------------------

# Probe denominators d with an attainable value e/d >= 13/100; checking a few
# arbitrary times t = m/d is sound for bailing because ANY f(t) is a lower
# bound on ML.  Small d kill most tuples in O(1).
PROBE_DENOMS = (7, 13, 14, 15, 9, 11, 20, 22, 23, 26, 27, 29, 31)


def is_near_tight(speeds, tn=THR_NUM, td=THR_DEN):
    """Return None if ML(speeds) >= tn/td (witnessed, early bail), else
    (ml_num, ml_den, wit_m, wit_d) with ML = ml_num/ml_den exactly < tn/td
    witnessed at t = wit_m/wit_d (reduced fractions not required here)."""
    # probe phase: cheap certified lower bounds (any f(t) lower-bounds ML)
    for d in PROBE_DENOMS:
        for m in range(1, (d // 2) + 1):
            e = d
            for v in speeds:
                r = (v * m) % d
                if r > d - r:
                    r = d - r
                if r < e:
                    e = r
                    if td * e < tn * d:
                        break
            if td * e >= tn * d:
                return None
    # full exact phase with early bail
    best_n, best_d, wit = 0, 1, (0, 1)
    for d in candidate_denominators(speeds):
        for m in range(1, d):
            e = d
            for v in speeds:
                r = (v * m) % d
                if r > d - r:
                    r = d - r
                if r < e:
                    e = r
                    if e == 0:
                        break
            if e * best_d > best_n * d:
                best_n, best_d, wit = e, d, (m, d)
                if td * e >= tn * d:
                    return None  # ML >= tn/td: not near-tight
    return best_n, best_d, wit[0], wit[1]


def _scan_chunk(args):
    """Worker: scan all 7-subsets of {1..V} whose leading pair (v1,v2) is
    assigned to this chunk (round-robin over pairs, for load balance and to
    avoid redundant enumeration).  Scaled (gcd>1) copies of primitive tuples
    are scanned too and deduplicated in post.  Returns (n_scanned, hits)."""
    V, k, chunk, jobs, tn, td = args
    hits = []
    n = 0
    check = is_near_tight
    combos = itertools.combinations
    pair_idx = 0
    for v1 in range(1, V - k + 3):
        for v2 in range(v1 + 1, V - k + 4):
            pair_idx += 1
            if pair_idx % jobs != chunk:
                continue
            for rest in combos(range(v2 + 1, V + 1), k - 3):
                n += 1
                tup = (v1, v2) + rest
                res = check(tup, tn, td)
                if res is not None:
                    bn, bd, wm, wd = res
                    ml = Fraction(bn, bd)
                    hits.append({
                        "speeds": list(tup),
                        "ml": f"{ml.numerator}/{ml.denominator}",
                        "ml_float": float(ml),
                        "witness_t": f"{wm}/{wd}",
                    })
    return n, hits


def scan(V, k=8, jobs=4, out_path=None, tn=THR_NUM, td=THR_DEN):
    t0 = time.time()
    with mp.Pool(jobs) as pool:
        results = pool.map(_scan_chunk,
                           [(V, k, c, jobs, tn, td) for c in range(jobs)])
    n_scanned = sum(r[0] for r in results)
    hits = [h for r in results for h in r[1]]
    # dedupe scaled copies: keep primitive (gcd=1) representatives, but note
    # every raw hit was independently computed (a consistency check in itself)
    prim = {}
    for h in hits:
        v = h["speeds"]
        g = 0
        for x in v:
            g = gcd(g, x)
        key = tuple(x // g for x in v)
        h["primitive"] = list(key)
        h["gcd"] = g
        if g == 1:
            prim[key] = h
        else:
            # scaled copy must have identical ML to its primitive if present
            if key in prim and prim[key]["ml"] != h["ml"]:
                print(f"WARNING: scale inconsistency {v} vs {list(key)}")
    hits.sort(key=lambda h: (h["ml_float"], h["speeds"]))
    dt = time.time() - t0
    print(f"V={V} k={k}: scanned {n_scanned} tuples in {dt:.1f}s; "
          f"{len(hits)} near-tight hits, {len(prim)} primitive (ML < {tn}/{td})")
    for key in sorted(prim, key=lambda x: (prim[x]['ml_float'], x)):
        h = prim[key]
        print(f"  {tuple(h['speeds'])}  ML = {h['ml']} = {h['ml_float']:.6f}"
              f"  at t = {h['witness_t']}")
    if out_path:
        with open(out_path, "w") as f:
            json.dump({"V": V, "k": k, "threshold": f"{tn}/{td}",
                       "n_scanned": n_scanned, "seconds": round(dt, 1),
                       "n_hits_raw": len(hits), "n_primitive": len(prim),
                       "near_tight_raw": hits}, f, indent=1)
        print(f"wrote {out_path}")
    return hits


# ----------------------------------------------------------------------
# validation
# ----------------------------------------------------------------------

def validate():
    ok = True
    # 1. Known tight cases: (1,...,k-1) must give exactly 1/k for k = 3..8.
    for k in range(3, 9):
        v = tuple(range(1, k))
        ml, t = ml_exact(v)
        good = (ml == Fraction(1, k))
        ok &= good
        print(f"[tight] k={k} v={v}: ML = {ml} (expect 1/{k}) "
              f"witness t={t}  {'OK' if good else 'FAIL'}")
    # scaled copy: ML invariant under scaling
    ml, _ = ml_exact(tuple(3 * i for i in range(1, 8)))
    good = ml == Fraction(1, 8)
    ok &= good
    print(f"[scale] 3*(1..7): ML = {ml} (expect 1/8)  {'OK' if good else 'FAIL'}")

    rng = random.Random(20260725)
    # 2. Candidate-set cross-check vs strictly larger candidate set.
    for trial in range(150):
        kk = rng.randint(3, 7)
        v = tuple(sorted(rng.sample(range(1, 30), kk)))
        a, _ = ml_exact(v)
        b = ml_exact_paranoid(v)
        if a != b:
            ok = False
            print(f"[paranoid] MISMATCH v={v}: {a} vs {b}")
    print("[paranoid] 150 random tuples: sums-only candidate set == "
          "sums+diffs+2v set" if ok else "[paranoid] FAILURES above")

    # 3. Dense-grid numeric cross-check (Lipschitz bound both ways).
    worst = 0.0
    for trial in range(25):
        v = tuple(sorted(rng.sample(range(1, 40), 7)))
        a, _ = ml_exact(v)
        g, half_step = ml_grid_numeric(v, 100_001)
        # grid value is a genuine f-value so g <= ML (+ float eps);
        # conversely ML <= g + vmax/n_grid... use Lipschitz slack.
        lo_ok = g <= float(a) + 1e-9
        hi_ok = float(a) <= g + 2 * half_step + 1e-9
        worst = max(worst, abs(float(a) - g))
        if not (lo_ok and hi_ok):
            ok = False
            print(f"[grid] MISMATCH v={v}: exact={float(a):.8f} grid={g:.8f} "
                  f"slack={2*half_step:.2e}")
    print(f"[grid] 25 random 7-tuples vs 1e5-point grid: max |diff| = "
          f"{worst:.2e} (within Lipschitz slack)")

    # 4. is_near_tight agrees with ml_exact on random tuples.
    for trial in range(400):
        v = tuple(sorted(rng.sample(range(1, 26), 7)))
        a, _ = ml_exact(v)
        r = is_near_tight(v)
        near = a < Fraction(THR_NUM, THR_DEN)
        if near != (r is not None) or (r and Fraction(r[0], r[1]) != a):
            ok = False
            print(f"[bail] MISMATCH v={v}: exact={a} bailres={r}")
    print("[bail] 400 random 7-tuples: is_near_tight consistent with ml_exact")

    print("VALIDATION", "PASSED" if ok else "FAILED")
    return ok


# ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--ml", type=str, help="comma-separated speeds")
    ap.add_argument("--scan", type=int, metavar="V")
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--out", type=str)
    ap.add_argument("--thr", type=str, default="13/100",
                    metavar="NUM/DEN", help="near-tight threshold (strict <)")
    args = ap.parse_args()
    if args.validate:
        sys.exit(0 if validate() else 1)
    if args.ml:
        v = tuple(sorted(int(x) for x in args.ml.split(",")))
        assert len(set(v)) == len(v) and min(v) >= 1
        ml, t = ml_exact(v)
        print(f"v={v}  ML = {ml} = {float(ml):.8f}  witness t = {t}")
        return
    if args.scan:
        tn, td = (int(x) for x in args.thr.split("/"))
        scan(args.scan, k=args.k, jobs=args.jobs, out_path=args.out,
             tn=tn, td=td)
        return
    ap.print_help()


if __name__ == "__main__":
    main()
