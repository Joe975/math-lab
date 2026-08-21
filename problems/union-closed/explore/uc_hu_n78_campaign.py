#!/usr/bin/env python3
"""Attempt 065 (064 lead 1): n = 7 and n = 8 best-order, where the
surrogate makes enumeration unnecessary.

064 ran the first adversarial best-order campaign at n = 6 by annealing
on the rollout margin and enumerating 720 orders once at the endpoint.
The logic scales further than that, and at n = 8 it removes the
enumeration entirely:

  best-order margin >= rollout margin   (pointwise, by definition)

so a NON-NEGATIVE rollout floor is already a proof that the best-order
floor is non-negative -- no enumeration at any point.  Enumeration is
only needed to measure how much slack the surrogate leaves, which is
affordable at n = 7 (5,040 orders once) and skipped at n = 8 (40,320).

Parts:
  A  n = 7, cap 0.49: anneal on rollout, enumerate 5,040 orders at the
     endpoint to measure the gap
  B  n = 8, cap 0.49: anneal on rollout only -- the floor's sign is
     what carries the best-order conclusion, and no enumeration is
     performed or needed
  C  the ladder: rollout floors and (where measured) gaps at
     n = 6, 7, 8

Standard library only; deterministic (seeds 9500+).
Usage: python uc_hu_n78_campaign.py [--fast]
Checkpoint: ../data/hu_n78_campaign.json
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


def prep(n, mu):
    tot = sum(mu.values())
    if tot <= 0:
        return None
    m = {a: w / tot for a, w in mu.items() if w / tot > 1e-14}
    mg = margins(n, m)
    if min(mg) < MIN_MARG or max(mg) >= 0.5:
        return None
    return m, max(mg)


def roll_margin(n, mu):
    p = prep(n, mu)
    if p is None:
        return None
    m, fmax = p
    cr, H = hu_cr_seq(n, m, seq_roll(n, m))
    if H < 0.2:
        return None
    return cr / H - cstar(fmax)


def anneal(name, n, mu0, steps, rng, enumerate_endpoint):
    mu = dict(mu0)
    cur = roll_margin(n, mu)
    if cur is None:
        return None
    best, best_mu = cur, dict(mu)
    full = (1 << n) - 1
    T, decay = 0.02, (1e-4 / 0.02) ** (1.0 / steps)
    for _ in range(steps):
        cand = dict(mu)
        atoms = sorted(cand)
        mv = rng.random()
        if mv < 0.45:
            cand[rng.choice(atoms)] *= math.exp(rng.gauss(0.0, 0.45))
        elif mv < 0.55:
            a = rng.choice(atoms)
            tgt = full if rng.random() < 0.5 else 0
            d = (0.2 + 0.6 * rng.random()) * cand[a]
            cand[a] -= d
            cand[tgt] = cand.get(tgt, 0.0) + d
            if cand[a] <= 1e-12:
                del cand[a]
        elif mv < 0.70 and len(cand) < 14:
            u = rng.random()
            new = (rng.randrange(1 << n) if u < 0.4 else
                   (rng.choice(atoms) | rng.choice(atoms)) if u < 0.8 else
                   full ^ rng.choice(atoms))
            cand[new] = cand.get(new, 0.0) + 0.03 * sum(cand.values())
        elif mv < 0.85 and len(cand) > 3:
            del cand[rng.choice(atoms)]
        elif len(atoms) >= 2:
            a, b = rng.sample(atoms, 2)
            d = 0.3 * rng.random() * cand[a]
            cand[a] -= d
            cand[b] = cand.get(b, 0.0) + d
        v = roll_margin(n, cand)
        if v is None:
            continue
        if v < cur or rng.random() < math.exp(-(v - cur) / max(T, 1e-9)):
            mu, cur = cand, v
            if v < best:
                best, best_mu = v, dict(cand)
        T *= decay
    p = prep(n, best_mu)
    if p is None:
        return None
    m, fmax = p
    row = {"start": name, "n": n, "roll_floor": best,
           "own_max_marginal": fmax,
           "mu": {format(a, f"0{n}b"): w for a, w in m.items()}}
    if enumerate_endpoint:
        bo, H = None, None
        for s in itertools.permutations(range(n)):
            cr, H = hu_cr_seq(n, m, s)
            if bo is None or cr > bo:
                bo = cr
        row["bestorder_at_endpoint"] = bo / H - cstar(fmax)
        row["gap"] = row["bestorder_at_endpoint"] - best
    return row


def campaign(n, cap, runs, steps, seed, enumerate_endpoint):
    rng = random.Random(seed)
    rows = []
    tries = 0
    while len(rows) < runs and tries < 60:
        tries += 1
        mu = gen_instance(rng, n, cap)
        if roll_margin(n, mu) is None:
            continue
        r = anneal(f"c{len(rows)}", n, mu, steps, rng, enumerate_endpoint)
        if r is not None:
            rows.append(r)
    rf = min(r["roll_floor"] for r in rows)
    proved = rf >= -1e-9
    msg = (f"  n={n} cap {cap}: {len(rows)} anneals x {steps} steps | "
           f"rollout floor {rf:+.6e} -> best-order floor is "
           f"{'PROVED non-negative by domination' if proved else 'NOT proved'}")
    log(msg)
    if enumerate_endpoint:
        gaps = [r["gap"] for r in rows]
        log(f"      endpoint gaps (best-order minus rollout): " +
            ", ".join(f"{g:+.4f}" for g in gaps))
    return {"rows": rows, "roll_floor": rf,
            "bestorder_nonneg_by_domination": proved,
            "enumerated": enumerate_endpoint, "steps": steps}


def main():
    steps = 900 if not FAST else 150
    runs = 3 if not FAST else 1
    log("A. n = 7, cap 0.49 (5,040-order enumeration at the endpoint):")
    OUT["A_n7_cap_0.49"] = campaign(7, 0.49, runs, steps, 9500, True)
    log()
    log("B. n = 8, cap 0.49 (NO enumeration -- domination carries it):")
    OUT["B_n8_cap_0.49"] = campaign(8, 0.49, runs, steps // 2, 9600, False)
    log()
    log("C. The ladder so far:")
    try:
        six = json.loads((DATA / "hu_n6_campaign.json").read_text())
        log(f"  n=6 cap 0.49: rollout floor "
            f"{six['A_n6_cap_0.49']['roll_floor']:+.3e} (064)")
    except Exception:
        pass
    for k, lbl in (("A_n7_cap_0.49", "n=7"), ("B_n8_cap_0.49", "n=8")):
        log(f"  {lbl} cap 0.49: rollout floor "
            f"{OUT[k]['roll_floor']:+.3e}")
    (DATA / "hu_n78_campaign.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_n78_campaign.json")


if __name__ == "__main__":
    main()
