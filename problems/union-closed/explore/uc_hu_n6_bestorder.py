#!/usr/bin/env python3
"""Attempt 063 (059 lead 1): best-order at n = 6, made affordable.

n = 6 is the size where the promoted conjecture had essentially no
adversarial evidence: the best-order objective costs 720 CR
evaluations per candidate, so 054 could only afford a single-start
spot check and 059 stopped at n = 5.

The cheap route 059 suggested: rollout is a LOWER BOUND on best-order
(the best order is at least as good as the one rollout picks), so if
the ROLLOUT margin is already >= 0 the best-order bound holds with no
enumeration at all.  Enumerate the 720 orders only on instances where
rollout leaves the margin negative.

Parts:
  A  census at n = 6, cap 0.49: how often does the lower bound
     suffice, and are there violations among the rest
  B  the same at cap 0.497, where the margin is tighter
  C  a control at n = 5, where full enumeration is affordable, to
     confirm the shortcut agrees with the enumerated answer

Standard library only; deterministic (seeds 9100+).
Usage: python uc_hu_n6_bestorder.py [--fast]
Checkpoint: ../data/hu_n6_bestorder.json
"""
from __future__ import annotations

import itertools
import json
import random
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_order2 import hu_cr_seq, seq_roll, h
from uc_hu_rollcensus import gen_instance
from uc_hu_anneal_essential import margins

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}
MIN_MARG = 0.03


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def cstar(p):
    return (h(max(0.5, 1.0 - 2.0 * p)) - h(p)) / h(p)


def census(n, cap, count, seed, enumerate_all=False):
    rng = random.Random(seed)
    checked = enum_needed = viol = 0
    worst = None
    tried = 0
    agree = 0
    while checked < count and tried < 40 * count:
        tried += 1
        mu = gen_instance(rng, n, cap)
        tot = sum(mu.values())
        m = {a: w / tot for a, w in mu.items() if w / tot > 1e-14}
        mg = margins(n, m)
        if min(mg) < MIN_MARG or max(mg) >= 0.5:
            continue
        cr, H = hu_cr_seq(n, m, seq_roll(n, m))
        if H < 0.2:
            continue
        checked += 1
        c = cstar(max(mg))
        lb = cr / H - c                      # lower bound on best-order
        if lb >= 0 and not enumerate_all:
            if worst is None or lb < worst[0]:
                worst = (lb, "rollout-suffices")
            continue
        best = max(hu_cr_seq(n, m, s)[0]
                   for s in itertools.permutations(range(n)))
        v = best / H - c
        if lb < 0:
            enum_needed += 1
        if enumerate_all and lb >= 0 and v >= 0:
            agree += 1
        if v < -1e-12:
            viol += 1
        if worst is None or v < worst[0]:
            worst = (v, "enumerated")
    pct = 100 * (checked - enum_needed) / max(checked, 1)
    log(f"  n={n} cap {cap}: {checked} essential in-regime instances; "
        f"rollout alone certified the bound on {checked - enum_needed} "
        f"({pct:.1f}%); enumeration needed on {enum_needed}; "
        f"{viol} violations; min margin {worst[0]:+.6f} ({worst[1]})")
    return {"checked": checked, "lower_bound_sufficed": checked - enum_needed,
            "enumeration_needed": enum_needed, "violations": viol,
            "min_margin": worst[0], "min_source": worst[1],
            "agreements": agree if enumerate_all else None}


def main():
    N = 3000 if not FAST else 200
    log("A. n = 6, cap 0.49 (rollout as a lower bound on best-order):")
    OUT["A_n6_cap_0.49"] = census(6, 0.49, N, 9100)
    log()
    log("B. n = 6, cap 0.497:")
    OUT["B_n6_cap_0.497"] = census(6, 0.497, N, 9101)
    log()
    log("C. Control at n = 5 -- enumerate everything and confirm the "
        "shortcut never disagrees:")
    OUT["C_n5_control"] = census(5, 0.49, N // 3, 9102, enumerate_all=True)
    log(f"    (instances where the lower bound sufficed AND the "
        f"enumerated answer agreed: "
        f"{OUT['C_n5_control']['agreements']})")
    (DATA / "hu_n6_bestorder.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_n6_bestorder.json")


if __name__ == "__main__":
    main()
