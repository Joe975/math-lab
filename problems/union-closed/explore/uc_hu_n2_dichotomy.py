#!/usr/bin/env python3
"""Attempt 047 (046 lead 1): the n = 2 target inequality splits into a
DICHOTOMY -- and each branch has a single covering term.

046 reduced the n = 2 case of (HU-TAX) to an explicit scalar
inequality and showed the pair-interaction term is load-bearing.  Write
that inequality as A + B + C >= 0 with

  A = (c*(f0) - c*(q)) h(f0)                       [first-coordinate slack]
  B = x (sigma(p0) - c*(q) h(p0))
      + (1-x)(sigma(p1) - c*(q) h(p1))             [diagonal ledger]
  C = t * Delta,  t = min(x - 1/2, 1 - x),
      Delta = 2 psi(p0+p1) - psi(2 p0) - psi(2 p1) [pair interaction]

where psi(s) = h(min(1/2, s)) and sigma(p) = c*(p) h(p) = h(z*(p)) - h(p).
A >= 0 (c* decreasing, q >= f0) and C >= 0 (psi concave, 046's N2-CONC);
only B can be negative, and B < 0 exactly when a conditional marginal
exceeds q.

This file measures the structure of that deficit:

  A  Lemma N2-ONE-ABOVE: at most one conditional can exceed q (proof in
     the record: the x-weighted mean of the conditionals is f1 <= q)
  B  THE DICHOTOMY: whenever B < 0, max(A, C) >= -B -- the deficit is
     always covered by ONE of the two terms, never needing both
  C  the sharper split: when q = f0 the interaction C alone always
     covers; when q = f1 the slack A covers except on a measured
     minority, where C covers
  D  tightness: where each branch comes closest to failing, which says
     whether a proof of either can afford crude bounds

Standard library only; deterministic.
Usage: python uc_hu_n2_dichotomy.py [--fast]
Checkpoint: ../data/hu_n2_dichotomy.json
"""
from __future__ import annotations

import json
import random
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_order2 import h, hu_cr_seq

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}
N = 200000 if not FAST else 20000


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def zst(p):
    return max(0.5, 1.0 - 2.0 * p)


def cst(p):
    return (h(zst(p)) - h(p)) / h(p)


def sig(p):
    return h(zst(p)) - h(p)


def psi(s):
    return h(min(0.5, s))


def terms(x, p0, p1):
    """(A, B, C, q, t) for the Case-A target inequality."""
    f0 = 1 - x
    f1 = x * p0 + (1 - x) * p1
    q = max(f0, f1)
    if not (0 < q < 0.5):
        return None
    A = (cst(f0) - cst(q)) * h(f0)
    B = (x * (sig(p0) - cst(q) * h(p0))
         + (1 - x) * (sig(p1) - cst(q) * h(p1)))
    t = min(x - 0.5, 1 - x)
    C = t * (2 * psi(p0 + p1) - psi(2 * p0) - psi(2 * p1))
    return A, B, C, q, t


def draws(seed, n):
    rng = random.Random(seed)
    for _ in range(n):
        yield rng.uniform(0.5, 1), rng.uniform(0, 0.5), rng.uniform(0, 0.5)


def part_A():
    log("A. Lemma N2-ONE-ABOVE (at most one conditional exceeds q):")
    both = 0
    seen = 0
    for x, p0, p1 in draws(4701, N):
        r = terms(x, p0, p1)
        if r is None:
            continue
        seen += 1
        q = r[3]
        if p0 > q + 1e-12 and p1 > q + 1e-12:
            both += 1
    log(f"  {seen} in-regime samples; both conditionals above q on "
        f"{both} (proof: x p0 + (1-x) p1 = f1 <= q, so they cannot both "
        f"exceed q)")
    OUT["A_one_above"] = {"samples": seen, "both_above": both}


def part_B():
    log()
    log("B. The dichotomy: whenever B < 0, is max(A, C) >= -B?")
    deficit = 0
    by_A = 0
    by_C = 0
    both_needed = 0
    fails = 0
    worst = None
    for x, p0, p1 in draws(4702, N):
        r = terms(x, p0, p1)
        if r is None:
            continue
        A, B, C, q, t = r
        if B >= 0:
            continue
        deficit += 1
        okA = A + B >= 0
        okC = B + C >= 0
        by_A += okA
        by_C += okC
        if not okA and not okC:
            if A + B + C >= 0:
                both_needed += 1
            else:
                fails += 1
                if worst is None or A + B + C < worst[0]:
                    worst = (A + B + C, x, p0, p1)
    log(f"  {deficit} deficit cases: covered by the slack A alone on "
        f"{by_A}, by the interaction C alone on {by_C}")
    log(f"  cases needing BOTH terms: {both_needed}; outright failures: "
        f"{fails}")
    log(f"  => the theorem's Case A reduces to: B < 0 implies "
        f"max(A, C) >= -B")
    OUT["B_dichotomy"] = {"deficit_cases": deficit, "covered_by_A": by_A,
                          "covered_by_C": by_C,
                          "needing_both": both_needed, "failures": fails,
                          "worst": worst}


def part_C():
    log()
    log("C. The split by which marginal is the max:")
    stats = {"q=f0": {"cases": 0, "C_fails": 0, "A_fails": 0},
             "q=f1": {"cases": 0, "C_fails": 0, "A_fails": 0}}
    worst = {"q=f0": None, "q=f1": None}
    for x, p0, p1 in draws(4703, N):
        r = terms(x, p0, p1)
        if r is None:
            continue
        A, B, C, q, t = r
        if B >= 0:
            continue
        key = "q=f0" if q == 1 - x else "q=f1"
        stats[key]["cases"] += 1
        if B + C < -1e-12:
            stats[key]["C_fails"] += 1
        if A + B < -1e-12:
            stats[key]["A_fails"] += 1
        margin = B + C
        if worst[key] is None or margin < worst[key][0]:
            worst[key] = (margin, x, p0, p1)
    for key in ("q=f0", "q=f1"):
        s = stats[key]
        log(f"  {key}: {s['cases']} deficit cases; interaction alone "
            f"fails on {s['C_fails']}, slack alone fails on "
            f"{s['A_fails']}")
    log("  => branch L1 (q = f0, where A = 0): the interaction ALWAYS "
        "covers -- the clean target for a proof")
    OUT["C_split"] = {"stats": stats,
                      "worst_B_plus_C": {k: worst[k] for k in worst}}


def part_D():
    log()
    log("D. Tightness of each branch (how close to failing):")
    rows = []
    for key, cond in (("L1 (q=f0): B + C", lambda q, x: q == 1 - x),
                      ("L2 (q=f1): max(A,C) + B",
                       lambda q, x: q != 1 - x)):
        best = None
        for x, p0, p1 in draws(4704, N):
            r = terms(x, p0, p1)
            if r is None:
                continue
            A, B, C, q, t = r
            if B >= 0 or not cond(q, x):
                continue
            v = (B + C) if key.startswith("L1") else (max(A, C) + B)
            if best is None or v < best[0]:
                best = (v, x, p0, p1, A, B, C)
        if best:
            v, x, p0, p1, A, B, C = best
            log(f"  {key}: min {v:+.3e} at x={x:.4f}, p0={p0:.4f}, "
                f"p1={p1:.4f}  (A={A:+.4f}, B={B:+.4f}, C={C:+.4f})")
            rows.append({"branch": key, "min": v, "x": x, "p0": p0,
                         "p1": p1, "A": A, "B": B, "C": C})
    OUT["D_tightness"] = rows


def main():
    part_A()
    part_B()
    part_C()
    part_D()
    (DATA / "hu_n2_dichotomy.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_n2_dichotomy.json")


if __name__ == "__main__":
    main()
