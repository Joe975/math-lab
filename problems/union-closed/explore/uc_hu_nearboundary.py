#!/usr/bin/env python3
"""Attempt 060 (correction to 059): the caps constrained only the START,
so 059 never tested the boundary it claimed to.

059 reported "at the hardest cap tested -- 0.499 -- the anneal lands
exactly on the equality family".  Inspecting the endpoints shows the
cap constrains only the instances the campaign STARTS from: the
admissibility test is max marginal < 1/2, so the anneal is free to
wander to small marginals, and it does.  The sharpest endpoint of every
block has drifted far below its cap:

    cap 0.499 (n=4): sharpest endpoint own max marginal 0.0360
    cap 0.497 (n=4): 0.1217        cap 0.495 (n=5): 0.0930
    cap 0.499 (n=5): 0.0390 and 0.0347

while the endpoints that stay near the cap have floors three orders of
magnitude LARGER (+1.2e-02 at 0.4998).  So 059's zero-violation result
stands, but it is a statement about the whole in-regime space, not
about the boundary.

This file adds the missing constraint -- a marginal FLOOR, max marginal
>= cap - 0.005, so the adversary must stay near the boundary -- and
reports the real near-boundary floors.

Parts:
  A  near-boundary anneals at n = 4, caps 0.495 / 0.497 / 0.499
  B  the same at n = 5, cap 0.499
  C  the drift table from 059's checkpoint, as the evidence for the
     correction

Standard library only; deterministic (seeds 8900+).
Usage: python uc_hu_nearboundary.py [--fast]
Checkpoint: ../data/hu_nearboundary.json
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
from uc_hu_order2 import hu_cr_seq, h
from uc_hu_rollcensus import gen_instance
from uc_hu_anneal_essential import margins

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}
MIN_MARG = 0.03
BAND = 0.005
ORDERS = {n: list(itertools.permutations(range(n))) for n in (4, 5)}


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def cstar(p):
    return (h(max(0.5, 1.0 - 2.0 * p)) - h(p)) / h(p)


def near_margin(n, mu, cap):
    """Best-order own-constant margin, admissible only if the instance is
    essential AND its max marginal stays within BAND of the cap."""
    tot = sum(mu.values())
    if tot <= 0:
        return None
    m = {a: w / tot for a, w in mu.items() if w / tot > 1e-14}
    mg = margins(n, m)
    if min(mg) < MIN_MARG or max(mg) >= 0.5 or max(mg) < cap - BAND:
        return None
    best, H = None, None
    for s in ORDERS[n]:
        cr, H = hu_cr_seq(n, m, s)
        if best is None or cr > best:
            best = cr
    if H is None or H < 0.2:
        return None
    return best / H - cstar(max(mg))


def anneal(name, n, mu0, cap, steps, rng):
    mu = dict(mu0)
    cur = near_margin(n, mu, cap)
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
        v = near_margin(n, cand, cap)
        if v is None:
            continue
        if v < cur or rng.random() < math.exp(-(v - cur) / max(T, 1e-9)):
            mu, cur = cand, v
            if v < best:
                best, best_mu = v, dict(cand)
        T *= decay
    tot = sum(best_mu.values())
    m = {a: w / tot for a, w in best_mu.items()}
    mg = margins(n, m)
    return {"start": name, "n": n, "cap": cap, "floor": best,
            "own_max_marginal": max(mg), "effective_n":
                sum(1 for x in mg if x > 1e-12),
            "mu": {format(a, f"0{n}b"): w for a, w in m.items()}}


def campaign(n, cap, runs, steps, seed):
    rng = random.Random(seed)
    rows, gmin = [], 1e9
    tries = 0
    while len(rows) < runs and tries < 120:
        tries += 1
        mu = gen_instance(rng, n, cap)
        if near_margin(n, mu, cap) is None:
            continue
        r = anneal(f"nb{len(rows)}", n, mu, cap, steps, rng)
        if r is None:
            continue
        rows.append(r)
        gmin = min(gmin, r["floor"])
    viol = [r for r in rows if r["floor"] < -1e-9]
    mm = [r["own_max_marginal"] for r in rows] or [float("nan")]
    log(f"  n={n} cap {cap}: {len(rows)} anneals x {steps} steps, floor "
        f"{gmin:+.3e}, {len(viol)} VIOLATIONS; endpoint marginals stay in "
        f"[{min(mm):.4f}, {max(mm):.4f}] (band floor {cap - BAND:.4f})")
    return {"rows": rows, "floor": gmin, "violations": viol,
            "steps": steps, "band_floor": cap - BAND}


def part_C():
    log()
    log("C. The drift in 059's endpoints (evidence for the correction):")
    try:
        old = json.loads((DATA / "hu_bestorder_anneal.json").read_text())
    except Exception:
        return
    rows = []
    for key in sorted(k for k in old if k.startswith(("A_", "B_"))):
        cap = float(key.split("_")[-1])
        r = min(old[key]["rows"], key=lambda z: z["floor"])
        rows.append({"block": key, "cap": cap,
                     "sharpest_floor": r["floor"],
                     "sharpest_own_marginal": r["max_marginal"],
                     "drifted": r["max_marginal"] < cap - BAND})
        log(f"  {key}: sharpest floor {r['floor']:+.3e} at own marginal "
            f"{r['max_marginal']:.4f} -> "
            f"{'DRIFTED off the cap' if r['max_marginal'] < cap - BAND else 'stayed near it'}")
    OUT["C_drift"] = rows


def main():
    steps = 1500 if not FAST else 200
    runs = 3 if not FAST else 1
    log("A. Near-boundary anneals at n = 4 (marginal held within "
        f"{BAND} of the cap):")
    for cap in (0.495, 0.497, 0.499):
        OUT[f"A_n4_cap_{cap}"] = campaign(4, cap, runs, steps,
                                          8900 + int(cap * 1000))
    log()
    log("B. The same at n = 5, cap 0.499:")
    OUT["B_n5_cap_0.499"] = campaign(5, 0.499, runs, steps // 2, 8999)
    part_C()
    (DATA / "hu_nearboundary.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_nearboundary.json")


if __name__ == "__main__":
    main()
