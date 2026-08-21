#!/usr/bin/env python3
"""Attempt 059 (058 lead 1): the strongest adversary yet, pointed at the
caps where the promoted conjecture is closest to its boundary.

045 attacked best-order HU at caps 0.495/0.497/0.499 with pattern-search
DESCENTS only.  Since then the adversary has improved three times:
  051  essentiality (no escaping to a lower-dimensional face)
  055  annealing instead of descent (five orders of magnitude sharper
       at n = 6)
  058  the collapse-toward-a-diagonal move (reaches the equality family
       at n = 7 where tripling the budget did not)

None of those existed when 045 ran, and 052 showed 4 of its 21-23
endpoints per cap had dropped a coordinate.  This re-runs 045's caps
with all three improvements, against the BEST-ORDER objective -- the
route's last unrefuted for-all-mu positivity statement.

Parts:
  A  n = 4 (24 orders per candidate), caps 0.495 / 0.497 / 0.499
  B  n = 5 (120 orders), caps 0.495 / 0.499
  C  comparison against 045's recorded floors, with each recorded
     endpoint labelled admissible or excluded by essentiality

A violation here would refute the promoted conjecture outright; the
expected outcome is the equality family, i.e. floors at 0.

Standard library only; deterministic (seeds 8700+).
Usage: python uc_hu_bestorder_anneal.py [--fast]
Checkpoint: ../data/hu_bestorder_anneal.json
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
ORDERS = {n: list(itertools.permutations(range(n))) for n in (4, 5)}


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def cstar(p):
    return (h(max(0.5, 1.0 - 2.0 * p)) - h(p)) / h(p)


def best_margin(n, mu):
    """Own-constant margin of the BEST order, or None if not essential."""
    tot = sum(mu.values())
    if tot <= 0:
        return None
    m = {a: w / tot for a, w in mu.items() if w / tot > 1e-14}
    mg = margins(n, m)
    if min(mg) < MIN_MARG or max(mg) >= 0.5:
        return None
    best, H = None, None
    for s in ORDERS[n]:
        cr, H = hu_cr_seq(n, m, s)
        if best is None or cr > best:
            best = cr
    if H is None or H < 0.2:
        return None
    return best / H - cstar(max(mg))


def anneal(name, n, mu0, steps, rng):
    """055's constrained anneal plus 058's collapse move."""
    mu = dict(mu0)
    cur = best_margin(n, mu)
    if cur is None:
        return None
    best, best_mu = cur, dict(mu)
    full = (1 << n) - 1
    T = 0.02
    decay = (1e-4 / T) ** (1.0 / steps)
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
        else:
            if len(atoms) >= 2:
                a, b = rng.sample(atoms, 2)
                d = 0.3 * rng.random() * cand[a]
                cand[a] -= d
                cand[b] = cand.get(b, 0.0) + d
        v = best_margin(n, cand)
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
    live = {a: w for a, w in m.items() if w > 1e-6}
    return {"start": name, "n": n, "floor": best,
            "effective_n": sum(1 for x in mg if x > 1e-12),
            "max_marginal": max(mg),
            "is_diagonal": set(live) <= {0, full} and len(live) == 2,
            "mu": {format(a, f"0{n}b"): w for a, w in m.items()}}


def campaign(n, cap, runs, steps, seed):
    rng = random.Random(seed)
    rows, gmin = [], 1e9
    tries = 0
    while len(rows) < runs and tries < 60:
        tries += 1
        mu = gen_instance(rng, n, cap)
        if best_margin(n, mu) is None:
            continue
        r = anneal(f"a{len(rows)}", n, mu, steps, rng)
        if r is None:
            continue
        rows.append(r)
        gmin = min(gmin, r["floor"])
    viol = [r for r in rows if r["floor"] < -1e-9]
    ndiag = sum(1 for r in rows if r["is_diagonal"])
    log(f"  n={n} cap {cap}: {len(rows)} anneals x {steps} steps, floor "
        f"{gmin:+.3e}, {len(viol)} VIOLATIONS, {ndiag} endpoints on the "
        f"diagonal family")
    return {"rows": rows, "floor": gmin, "violations": viol,
            "diagonal_endpoints": ndiag, "steps": steps}


def part_C():
    log()
    log("C. Against the 045 floors (descent-only, pre-improvements):")
    try:
        old = json.loads((DATA / "hu_bestorder.json").read_text())
    except Exception:
        return
    rows = []
    for cap in (0.495, 0.497, 0.499):
        key = f"cap_{cap}"
        newkey = f"A_n4_cap_{cap}"
        if key not in old or newkey not in OUT:
            continue
        r = min(old[key]["rows"], key=lambda z: z["floor"])
        n = r["n"]
        m = {int(s, 2): w for s, w in r["mu"].items()}
        t = sum(m.values())
        m = {a: w / t for a, w in m.items()}
        mg = margins(n, m)
        rows.append({"cap": cap, "new_floor": OUT[newkey]["floor"],
                     "recorded_floor": r["floor"],
                     "recorded_min_marginal": min(mg),
                     "recorded_admissible": min(mg) >= MIN_MARG})
        log(f"  cap {cap}: anneal {OUT[newkey]['floor']:+.3e} vs 045 "
            f"descent {r['floor']:+.3e} | recorded endpoint min marginal "
            f"{min(mg):.4f} -> "
            f"{'admissible' if min(mg) >= MIN_MARG else 'EXCLUDED here'}")
    OUT["C_comparison"] = rows


def main():
    steps = 1500 if not FAST else 200
    runs = 3 if not FAST else 1
    log("A. n = 4, all 24 orders per candidate:")
    for cap in (0.495, 0.497, 0.499):
        OUT[f"A_n4_cap_{cap}"] = campaign(4, cap, runs, steps,
                                          8700 + int(cap * 1000))
    log()
    log("B. n = 5, all 120 orders per candidate:")
    for cap in (0.495, 0.499):
        OUT[f"B_n5_cap_{cap}"] = campaign(5, cap, runs, steps // 2,
                                          8800 + int(cap * 1000))
    part_C()
    (DATA / "hu_bestorder_anneal.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_bestorder_anneal.json")


if __name__ == "__main__":
    main()
