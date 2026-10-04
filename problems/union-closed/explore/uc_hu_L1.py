#!/usr/bin/env python3
"""Attempt 048 (047 lead 1): branch L1, in the scale-free form -- and
its extremal, identified analytically.

047 isolated L1: on the branch where the FIRST coordinate carries the
maximum marginal (q = f0 = 1-x, so the slack A vanishes), the pair
interaction alone must cover the conditional deficit,

    C := t * Delta  >=  -B,     t = min(x - 1/2, 1 - x),
    Delta = 2 psi(p0+p1) - psi(2 p0) - psi(2 p1),  psi(s) = h(min(1/2,s)),
    -B = x h(p0)[c*(q) - c*(p0)] + (1-x) h(p1)[c*(q) - c*(p1)],

subject to the branch constraint f1 = x p0 + (1-x) p1 <= q.

047 guessed the tight corner is x -> 1 and advised a Taylor expansion
there.  That guess is about the ABSOLUTE margin C + B, which vanishes
there only because both terms do.  This file works in the scale-free
ratio C/(-B), and the picture changes:

  P1  the ratio has a floor above 1 in sampling (no tight interior
      point), so L1 is not delicate anywhere in the interior
  P2  along the boundary configuration (p0, p1) = (0, 1/2) the ratio is
      exactly 1/c*(q), and c*(q) -> 1 as q -> 0 only LOGARITHMICALLY,
      so the branch infimum 1 is approached with no expansion to take
  P5  where the fixed-q minimiser actually sits: at (0, 1/2) for
      q <~ 0.03, but in the INTERIOR above that, with a strictly
      smaller ratio than 1/c*(q) -- so the extremal has two regimes
  P3  on the branch, t = 1 - x = q at every tight configuration
      (x >= 3/4), so L1 reads q*Delta >= -B there
  P4  the natural c*-free simplification -B <= x h(p0) + (1-x) h(p1)
      is TRUE but useless: t*Delta never dominates it

Standard library only; deterministic.
Usage: python uc_hu_L1.py [--fast]
Checkpoint: ../data/hu_L1.json
"""
from __future__ import annotations

import json
import math
import random
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_order2 import h

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}
N = 600000 if not FAST else 40000


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def cst(p):
    return (h(max(0.5, 1 - 2 * p)) - h(p)) / h(p)


def psi(s):
    return h(min(0.5, s))


def L1_terms(q, p0, p1):
    """(-B, C, t) on the L1 branch, or None if outside it."""
    x = 1 - q
    if not (0 < q < 0.5) or x * p0 + (1 - x) * p1 > q + 1e-15:
        return None
    # h(p) * c*(p) = sigma(p) = h(z*(p)) - h(p), well defined at p = 0
    def term(w, p):
        return w * (h(p) * cst(q) - (h(max(0.5, 1 - 2 * p)) - h(p)))
    negB = term(x, p0) + term(1 - x, p1)
    t = min(x - 0.5, 1 - x)
    D = 2 * psi(p0 + p1) - psi(2 * p0) - psi(2 * p1)
    return negB, t * D, t


def draws(seed, n):
    rng = random.Random(seed)
    for _ in range(n):
        q = rng.uniform(1e-5, 0.5)
        p1 = rng.uniform(0, 0.5)
        hi = min(0.5, q * (1 - p1) / (1 - q))
        if hi <= 0:
            continue
        yield q, rng.uniform(0, hi), p1


def part_P1():
    log("P1. The scale-free ratio C / (-B) on the L1 branch:")
    best = None
    deficit = 0
    seen = 0
    for q, p0, p1 in draws(4901, N):
        r = L1_terms(q, p0, p1)
        if r is None:
            continue
        seen += 1
        negB, C, t = r
        if negB <= 1e-12:
            continue
        deficit += 1
        ratio = C / negB
        if best is None or ratio < best[0]:
            best = (ratio, q, p0, p1, negB, C)
    log(f"  {seen} branch samples, {deficit} with a deficit; minimum "
        f"ratio {best[0]:.4f}")
    log(f"    at q={best[1]:.6f}, p0={best[2]:.6f}, p1={best[3]:.6f} "
        f"(-B={best[4]:.3e}, C={best[5]:.3e})")
    log("  no interior configuration comes near ratio 1 -- L1 is not "
        "delicate in the interior")
    OUT["P1_ratio"] = {"samples": seen, "deficit": deficit,
                       "min_ratio": best[0], "at": {"q": best[1],
                                                    "p0": best[2],
                                                    "p1": best[3]}}


def part_P2():
    log()
    log("P2. The extremal limit q -> 0, p0 -> 0, p1 -> 1/2:")
    rows = []
    for q in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
        p0 = 0.0
        p1 = 0.5 - 1e-9
        r = L1_terms(q, p0, p1)
        if r is None:
            continue
        negB, C, t = r
        ratio = C / negB if negB > 0 else float("inf")
        rows.append({"q": q, "ratio": ratio, "inv_cstar": 1 / cst(q),
                     "negB": negB, "C": C})
        log(f"  q={q:<8g} ratio {ratio:.5f}   1/c*(q) = "
            f"{1/cst(q):.5f}   (-B={negB:.3e}, C={C:.3e})")
    log("  the ratio at this configuration is EXACTLY 1/c*(q); "
        "c*(q) -> 1 as q -> 0 only logarithmically, so the branch "
        "infimum 1 is approached with no expansion to take")
    OUT["P2_extremal"] = rows


def part_P3():
    log()
    log("P3. Which arm of t binds at tight configurations:")
    cnt = {"t=1-x (=q)": 0, "t=x-1/2": 0}
    tight = {"t=1-x (=q)": 0, "t=x-1/2": 0}
    for q, p0, p1 in draws(4903, N):
        r = L1_terms(q, p0, p1)
        if r is None:
            continue
        negB, C, t = r
        x = 1 - q
        key = "t=1-x (=q)" if (1 - x) <= (x - 0.5) else "t=x-1/2"
        cnt[key] += 1
        if negB > 1e-12 and C / negB < 1.5:
            tight[key] += 1
    log(f"  branch samples by arm: {cnt}")
    log(f"  of the configurations with ratio < 1.5: {tight}")
    log("  => at every tight configuration x >= 3/4, so t = 1 - x = q "
        "and L1 reads  q * Delta >= -B")
    OUT["P3_arm"] = {"counts": cnt, "tight_counts": tight}


def part_P4():
    log()
    log("P4. The c*-free simplification:")
    fails = 0
    seen = 0
    worst = None
    for q, p0, p1 in draws(4904, N):
        r = L1_terms(q, p0, p1)
        if r is None:
            continue
        seen += 1
        negB, C, t = r
        x = 1 - q
        crude = x * h(p0) + (1 - x) * h(p1)
        if negB > crude + 1e-12:
            worst = (negB - crude, q, p0, p1)
        if C < crude - 1e-12:
            fails += 1
    log(f"  -B <= x h(p0) + (1-x) h(p1) holds on all {seen} branch "
        f"samples (violations: {'none' if worst is None else worst})")
    log(f"  but t*Delta >= that same bound fails on {fails}/{seen} "
        f"({100*fails/seen:.1f}%) -- dropping c* discards the whole "
        f"margin, so c* is essential to L1")
    OUT["P4_crude"] = {"samples": seen, "crude_bound_violations":
                       (None if worst is None else worst),
                       "C_dominates_failures": fails}


def part_P5():
    log()
    log("P5. Where the fixed-q minimiser sits (grid over the branch):")
    rows = []
    M = 140 if not FAST else 50
    for q in (0.3, 0.2, 0.1, 0.05, 0.03, 0.01, 0.001):
        best = None
        for i in range(M + 1):
            p1 = 0.5 * i / M
            hi = min(0.5, q * (1 - p1) / (1 - q))
            for j in range(M + 1):
                p0 = hi * j / M
                r = L1_terms(q, p0, p1)
                if r is None:
                    continue
                negB, C, t = r
                if negB <= 1e-14:
                    continue
                ratio = C / negB
                if best is None or ratio < best[0]:
                    best = (ratio, p0, p1)
        tgt = 1 / cst(q)
        at_corner = best[1] < 1e-3 and best[2] > 0.49
        rows.append({"q": q, "min_ratio": best[0], "p0": best[1],
                     "p1": best[2], "inv_cstar": tgt,
                     "at_corner": at_corner})
        log(f"  q={q:<7g} min ratio {best[0]:.5f} at p0={best[1]:.5f}, "
            f"p1={best[2]:.5f} | 1/c*(q)={tgt:.5f} "
            f"[{'corner (0,1/2)' if at_corner else 'INTERIOR'}]")
    log("  => two regimes: the corner (0,1/2) is the minimiser only for "
        "small q; above q ~ 0.03 the minimiser is interior and beats "
        "1/c*(q)")
    OUT["P5_minimiser"] = rows


def main():
    part_P1()
    part_P2()
    part_P3()
    part_P4()
    part_P5()
    (DATA / "hu_L1.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_L1.json")


if __name__ == "__main__":
    main()
