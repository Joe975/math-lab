#!/usr/bin/env python3
"""Attempt 038 (037 leads 3 + 1): how close is ROLLOUT to the best
order, and does rollout-HU survive at n = 6 and 7?

037 established rollout-order HU at n <= 5: best order found on all
three record witnesses, descents survive to cap 0.499.  Two open edges:

  P1  CENSUS (037 lead 3): random in-regime measures at n = 4, 5, full
      order enumeration alongside -- how often is roll exactly best?
      what is its worst rank / worst gap?  canonical ranked for
      contrast.  Falsifiable reading: a census family where roll ranks
      LOW is the counterexample seed 037 asked for.
  P2  n = 6 descents (037 lead 1): the three starts 037 skipped for
      budget, plus fresh n = 6 hostile seeds, caps 0.49 / 0.497.
  P3  n = 7 spot descents (reduced rounds), cap 0.49.

Standard library only; deterministic (fixed seeds 938000+/941000+).
Usage: python uc_hu_rollcensus.py [--fast]
Checkpoint: ../data/hu_rollcensus.json
"""
from __future__ import annotations

import itertools
import json
import math
import random
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_order2 import (hu_cr_seq, seq_roll, seq_can, rule_ratio, descend,
                          h)
from uc_hu_attack import starts

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}

def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def gen_instance(rng, n, cap):
    """Random in-regime measure, 035-style delta_0 projection."""
    k = rng.randint(3, min(14, 1 << n))
    supp = rng.sample(range(1 << n), k)
    mu = {a: math.exp(rng.uniform(-2.5, 2.5)) for a in supp}
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    fm = max(sum(w for a, w in mu.items() if (a >> i) & 1)
             for i in range(n))
    if fm >= cap - 0.005:
        al = 1.0 - (cap - 0.005) / fm
        mu = {a: w * (1 - al) for a, w in mu.items()}
        mu[0] = mu.get(0, 0.0) + al
    return mu


def part_P1():
    log("P1. Census: roll rank vs full order enumeration:")
    N = 300 if not FAST else 30
    for n, cap in ((4, 0.45), (4, 0.49), (5, 0.45), (5, 0.49)):
        rows = []
        rng = random.Random(938000 + n * 1000 + int(cap * 1000))
        tried = 0
        while len(rows) < N and tried < 20 * N:
            tried += 1
            mu = gen_instance(rng, n, cap)
            H = -sum(w * math.log2(w) for w in mu.values() if w > 0)
            if H < 0.2:
                continue
            allcr = sorted(hu_cr_seq(n, mu, s)[0]
                           for s in itertools.permutations(range(n)))
            rollcr, _ = hu_cr_seq(n, mu, seq_roll(n, mu))
            cancr, _ = hu_cr_seq(n, mu, seq_can(n, mu))
            rank = 1 + sum(1 for c in allcr if c < rollcr - 1e-12)
            crank = 1 + sum(1 for c in allcr if c < cancr - 1e-12)
            ties_at_top = sum(1 for c in allcr if c > allcr[-1] - 1e-12)
            rows.append({
                "roll_CR": rollcr, "best_CR": allcr[-1],
                "worst_CR": allcr[0], "H": H,
                "roll_rank": rank, "canon_rank": crank,
                "roll_is_best": bool(allcr[-1] - rollcr <= 1e-12),
                "gap_over_H": (allcr[-1] - rollcr) / H,
                "any_neg": bool(allcr[0] < 0),
                "roll_neg": bool(rollcr < 0),
                "top_ties": ties_at_top})
        nb = sum(1 for r in rows if r["roll_is_best"])
        worst_rank = min(r["roll_rank"] for r in rows)
        wgap = max(r["gap_over_H"] for r in rows)
        rollneg = sum(1 for r in rows if r["roll_neg"])
        anyneg = sum(1 for r in rows if r["any_neg"])
        canon_below = sum(1 for r in rows
                          if r["canon_rank"] < r["roll_rank"])
        log(f"  n={n} cap {cap}: {len(rows)} instances; roll=best on "
            f"{nb} ({100*nb/len(rows):.1f}%), worst roll rank "
            f"{worst_rank}/{math.factorial(n)}, worst (best-roll)/H "
            f"{wgap:.2e}; roll negative on {rollneg} "
            f"(some order negative on {anyneg}); canonical strictly "
            f"below roll on {canon_below}")
        OUT[f"P1_n{n}_cap{cap}"] = {
            "instances": len(rows), "roll_best": nb,
            "worst_roll_rank": worst_rank, "worst_gap_over_H": wgap,
            "roll_negative": rollneg, "any_order_negative": anyneg,
            "canon_below_roll": canon_below,
            "rows": rows}


def part_P2():
    log()
    log("P2. n = 6 rollout descents (the starts 037 skipped + hostile "
        "seeds):")
    caps = (0.49, 0.497) if not FAST else (0.49,)
    for cap in caps:
        seeds = [(name, n, mu) for name, n, mu in starts(cap) if n == 6]
        rng = random.Random(941000 + int(cap * 1000))
        for j in range(4 if not FAST else 1):
            seeds.append((f"n6seed{j}", 6, gen_instance(rng, 6, cap)))
        rows, gmin = [], 1e9
        for name, n, mu in seeds:
            r = descend(name, n, mu, cap,
                        lambda n_, m_, c_: rule_ratio(n_, m_, c_, seq_roll),
                        rounds_max=120)
            if r is None:
                continue
            rows.append(r)
            gmin = min(gmin, r["floor"])
            log(f"    {name:12s}: floor {r['floor']:+.6f}")
        viol = [r for r in rows if r["floor"] < 0]
        log(f"  roll n=6 cap {cap}: {len(rows)} starts, global floor "
            f"{gmin:+.6f}, {len(viol)} violations "
            f"[product extremal {(1 - h(cap)) / h(cap):+.6f}]")
        OUT[f"P2_cap_{cap}"] = {"rows": rows, "global_floor": gmin,
                                "violations": viol}


def part_P3():
    log()
    log("P3. n = 7 spot descents, cap 0.49 (reduced rounds):")
    cap = 0.49
    rng = random.Random(941700)
    seeds = [(f"n7seed{j}", 7, gen_instance(rng, 7, cap))
             for j in range(3 if not FAST else 1)]
    rows, gmin = [], 1e9
    for name, n, mu in seeds:
        r = descend(name, n, mu, cap,
                    lambda n_, m_, c_: rule_ratio(n_, m_, c_, seq_roll),
                    rounds_max=50, step_floor=0.05)
        if r is None:
            continue
        rows.append(r)
        gmin = min(gmin, r["floor"])
        log(f"    {name:12s}: floor {r['floor']:+.6f}")
    viol = [r for r in rows if r["floor"] < 0]
    log(f"  roll n=7 cap {cap}: {len(rows)} starts, global floor "
        f"{gmin:+.6f}, {len(viol)} violations "
        f"[product extremal {(1 - h(cap)) / h(cap):+.6f}]")
    OUT["P3_cap_0.49"] = {"rows": rows, "global_floor": gmin,
                          "violations": viol}


def main():
    part_P1()
    part_P2()
    part_P3()
    (DATA / "hu_rollcensus.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_rollcensus.json")


if __name__ == "__main__":
    main()
