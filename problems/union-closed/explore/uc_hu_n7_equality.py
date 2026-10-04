#!/usr/bin/env python3
"""Attempt 057 (056 lead 2): can a constrained adversary reach the
equality family at n = 7, or only at n = 6?

055's constrained anneal landed exactly on the equality family at
n = 6 (certified in 056: margin enclosed within 4.9e-32 of 0) but
stopped at +2.484e-04 at n = 7.  Two readings, and they differ in
consequence:

  (i)  the anneal simply did not converge at n = 7 -- more budget or a
       better seed reaches the family, and "a constrained adversary
       can always reach equality" survives;
  (ii) the family is genuinely harder to reach at n = 7, which would
       mean the essential floor there is strictly positive and the
       route has a size-dependent gap it has never seen.

The falsifiable separation: SEED the anneal at the family itself
(a diagonal at the cap) and see whether it stays.  If the family is a
local minimum of the constrained objective at n = 7, the anneal will
sit there and reading (i) is right; if it drifts away and settles
higher, something is wrong with the framing rather than the budget.
Then give the unseeded search a longer budget and see if it arrives.

Parts:
  A  seeded-at-the-family runs at n = 7 (and n = 6 as a control)
  B  unseeded long runs at n = 7 (3x the 055 budget)
  C  the diagonal's exact margin at n = 7 by the 042 identity, as the
     target the search is trying to hit

Standard library only; deterministic (seeds 8300+).
Usage: python uc_hu_n7_equality.py [--fast]
Checkpoint: ../data/hu_n7_equality.json
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
from uc_hu_order2 import hu_cr_seq, seq_roll, h
from uc_hu_rollcensus import gen_instance
from uc_hu_anneal_essential import anneal, essential_margin, margins

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def cstar(p):
    return (h(max(0.5, 1.0 - 2.0 * p)) - h(p)) / h(p)


def diagonal(n, p):
    return {0: 1.0 - p, (1 << n) - 1: p}


def part_C():
    log("C. The diagonal's margin at n = 6 and 7 (042 identity):")
    rows = []
    for n in (6, 7):
        for p in (0.3826, 0.45, 0.487):
            mu = diagonal(n, p)
            cr, H = hu_cr_seq(n, mu, tuple(range(n)))
            mg = margins(n, {a: w for a, w in mu.items()})
            m = cr / H - cstar(max(mg))
            rows.append({"n": n, "p": p, "margin": m, "CR": cr, "H": H})
            log(f"  n={n} p={p}: CR={cr:+.6f} H={H:.6f} margin "
                f"{m:+.3e}  (identity predicts exactly 0)")
    OUT["C_diagonal_identity"] = rows


def part_A():
    log()
    log("A. Anneal SEEDED at the family -- does it stay?")
    steps = 900 if not FAST else 150
    rows = []
    for n in (6, 7):
        for p in (0.3826, 0.487):
            rng = random.Random(8300 + n * 10 + int(p * 1000))
            mu = diagonal(n, p)
            if essential_margin(n, mu) is None:
                log(f"  n={n} p={p}: seed not admissible, skipped")
                continue
            r = anneal(f"seeded_p{p}", n, mu, steps, rng)
            if r is None:
                continue
            stayed = r["floor"] < 1e-9
            rows.append({"n": n, "p": p, "floor": r["floor"],
                         "effective_n": r["effective_n"],
                         "stayed_on_family": stayed})
            log(f"  n={n} p={p}: floor after {steps} steps "
                f"{r['floor']:+.3e}, effective n={r['effective_n']} -> "
                f"{'STAYED on the family' if stayed else 'drifted off'}")
    OUT["A_seeded"] = rows


def part_B():
    log()
    log("B. Unseeded long runs at n = 7 (3x the 055 budget):")
    steps = 2700 if not FAST else 200
    rng = random.Random(8400)
    rows, gmin = [], 1e9
    tries = 0
    while len(rows) < (3 if not FAST else 1) and tries < 60:
        tries += 1
        mu = gen_instance(rng, 7, 0.49)
        if essential_margin(7, mu) is None:
            continue
        r = anneal(f"long{len(rows)}", 7, mu, steps, rng)
        if r is None:
            continue
        rows.append(r)
        gmin = min(gmin, r["floor"])
        log(f"    {r['start']}: floor {r['floor']:+.3e}, effective "
            f"n={r['effective_n']}")
    reached = gmin < 1e-9
    log(f"  n=7 cap 0.49, {len(rows)} runs x {steps} steps: floor "
        f"{gmin:+.3e} -> "
        f"{'REACHED the family' if reached else 'did NOT reach it'} "
        f"(055 at 900 steps: +2.484e-04)")
    OUT["B_unseeded_long"] = {"rows": rows, "floor": gmin,
                              "reached_family": reached, "steps": steps}


def main():
    part_C()
    part_A()
    part_B()
    (DATA / "hu_n7_equality.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_n7_equality.json")


if __name__ == "__main__":
    main()
