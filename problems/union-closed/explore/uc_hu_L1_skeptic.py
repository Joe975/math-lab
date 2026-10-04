#!/usr/bin/env python3
"""Skeptic pass on uc_hu_L1.py (attempt 048).  Stance: refute.

Independence: c*, psi and the L1 terms are rebuilt here in NATS from
the record's prose, and the branch constraint is re-derived rather than
copied.  The checks that matter are the two the record leans on: that
the ratio never drops below 1 anywhere on the branch (attacked by
descent, not just sampling), and that the boundary configuration
(0, 1/2) really evaluates to 1/c*(q) exactly.

  T1  the ratio floor: descent that actively minimises C/(-B)
  T2  the identity ratio(q, 0, 1/2) = 1/c*(q), independently
  T3  P5's two-regime claim: re-locate the fixed-q minimiser
  T4  P4's negative: the c*-free bound really does fail everywhere

Usage: python uc_hu_L1_skeptic.py     (exit 0 iff nothing refuted)
"""
from __future__ import annotations

import json
import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
LN2 = math.log(2.0)
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def hb(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * math.log(p) + (1 - p) * math.log(1 - p)) / LN2


def cstar(p):
    return (hb(max(0.5, 1 - 2 * p)) - hb(p)) / hb(p) if 0 < p < 1 else None


def psi(s):
    return hb(min(0.5, s))


def ratio(q, p0, p1):
    x = 1 - q
    if not (0 < q < 0.5):
        return None
    if x * p0 + (1 - x) * p1 > q + 1e-15:
        return None
    cq = cstar(q)
    if cq is None:
        return None
    def sig(p):
        return hb(max(0.5, 1 - 2 * p)) - hb(p)
    negB = x * (hb(p0) * cq - sig(p0)) + (1 - x) * (hb(p1) * cq - sig(p1))
    if negB <= 1e-14:
        return None
    t = min(x - 0.5, 1 - x)
    C = t * (2 * psi(p0 + p1) - psi(2 * p0) - psi(2 * p1))
    return C / negB


def main():
    out = json.loads((DATA / "hu_L1.json").read_text())
    rng = random.Random(5150)

    print("T1. The ratio floor, attacked by descent:")
    best = None
    for _ in range(200000):
        q = rng.uniform(1e-6, 0.5)
        p1 = rng.uniform(0, 0.5)
        hi = min(0.5, q * (1 - p1) / (1 - q))
        if hi <= 0:
            continue
        r = ratio(q, rng.uniform(0, hi), p1)
        if r is not None and (best is None or r < best[0]):
            best = (r, q, p1)
    q, p1 = best[1], best[2]
    p0 = min(0.5, q * (1 - p1) / (1 - q)) / 2
    cur = ratio(q, p0, p1) or best[0]
    step = 0.05
    while step > 1e-8:
        improved = False
        for d in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0),
                  (0, 0, 1), (0, 0, -1)):
            c = (min(0.4999, max(1e-7, q + step * d[0] * q)),
                 max(0.0, min(0.5, p0 + step * d[1])),
                 max(0.0, min(0.5, p1 + step * d[2])))
            v = ratio(*c)
            if v is not None and v < cur - 1e-15:
                q, p0, p1, cur = c[0], c[1], c[2], v
                improved = True
        if not improved:
            step /= 2
    check("C / (-B) >= 1 survives a descent that minimises it",
          cur >= 1 - 1e-9,
          f"min ratio {cur:.6f} at q={q:.3e}, p0={p0:.3e}, p1={p1:.5f}")

    print("T2. The boundary identity ratio(q, 0, 1/2) = 1/c*(q):")
    worst = 0.0
    for q in (0.2, 0.1, 0.01, 1e-3, 1e-4, 1e-5):
        r = ratio(q, 0.0, 0.5 - 1e-12)
        if r is None:
            continue
        worst = max(worst, abs(r - 1 / cstar(q)))
    check("identity holds to 1e-6", worst < 1e-6, f"max |diff| {worst:.2e}")

    print("T3. P5's two-regime claim, re-located independently:")
    mism = []
    for row in out["P5_minimiser"]:
        q = row["q"]
        best2 = None
        M = 90
        for i in range(M + 1):
            p1 = 0.5 * i / M
            hi = min(0.5, q * (1 - p1) / (1 - q))
            for j in range(M + 1):
                p0 = hi * j / M
                r = ratio(q, p0, p1)
                if r is not None and (best2 is None or r < best2[0]):
                    best2 = (r, p0, p1)
        corner = best2[1] < 5e-3 and best2[2] > 0.487
        if corner != row["at_corner"]:
            mism.append((q, best2, row["at_corner"]))
    check("the corner/interior classification reproduces at every q",
          not mism, f"mismatches: {mism[:2]}")

    print("T4. P4's negative (the c*-free bound is useless):")
    fails = 0
    seen = 0
    for _ in range(100000):
        q = rng.uniform(1e-5, 0.5)
        p1 = rng.uniform(0, 0.5)
        hi = min(0.5, q * (1 - p1) / (1 - q))
        if hi <= 0:
            continue
        p0 = rng.uniform(0, hi)
        x = 1 - q
        if x * p0 + (1 - x) * p1 > q + 1e-15:
            continue
        seen += 1
        t = min(x - 0.5, 1 - x)
        C = t * (2 * psi(p0 + p1) - psi(2 * p0) - psi(2 * p1))
        if C < x * hb(p0) + (1 - x) * hb(p1) - 1e-12:
            fails += 1
    check("t*Delta fails to dominate the c*-free bound on ~all of the "
          "branch", seen and fails / seen > 0.99,
          f"{fails}/{seen} = {100*fails/max(seen,1):.1f}%")

    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("No refutation found.")


if __name__ == "__main__":
    main()
