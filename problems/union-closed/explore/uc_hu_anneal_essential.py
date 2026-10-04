#!/usr/bin/env python3
"""Attempt 055 (054 lead 2): the 040 anneal, re-run with no degenerate
escape -- the last campaign still lacking a constrained counterpart.

052 found the route's worst dimension overstatement in 040: its
"n = 7, floor +0.000488" endpoint is a TWO-dimensional instance, and
the n = 6 and n = 7 runs reported the same floor because they had
collapsed to the same instance.  054 re-ran 038's descents under 051's
essentiality constraint and every floor rose; the anneal was left
undone because it is a different search (temperature, atom add/drop,
pairwise transfer) and could behave differently -- an anneal accepts
uphill moves, so it can wander into a degenerate face and stay there
even when a descent would not.

This file re-runs it with the constraint enforced INSIDE the move
acceptance, so no accepted state ever loses a coordinate.

Parts:
  A  constrained anneal at n = 6, caps 0.49 and 0.497
  B  constrained anneal at n = 7, cap 0.49
  C  the comparison: constrained floors against 040's recorded ones,
     with the effective dimension of each endpoint reported

Standard library only; deterministic (fixed seeds 8100+).
Usage: python uc_hu_anneal_essential.py [--fast]
Checkpoint: ../data/hu_anneal_essential.json
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

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}
MIN_MARG = 0.03


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def cstar(p):
    return (h(max(0.5, 1.0 - 2.0 * p)) - h(p)) / h(p)


def margins(n, m):
    return [sum(w for a, w in m.items() if (a >> i) & 1) for i in range(n)]


def essential_margin(n, mu):
    """Own-constant margin, or None if the state is not essential.
    The constraint lives here, so an accepted anneal state can never
    have lost a coordinate."""
    tot = sum(mu.values())
    if tot <= 0:
        return None
    m = {a: w / tot for a, w in mu.items() if w / tot > 1e-14}
    mg = margins(n, m)
    if min(mg) < MIN_MARG or max(mg) >= 0.5:
        return None
    cr, H = hu_cr_seq(n, m, seq_roll(n, m))
    if H < 0.2:
        return None
    return cr / H - cstar(max(mg))


def anneal(name, n, mu0, steps, rng):
    """040's move set, with essentiality enforced in acceptance."""
    mu = dict(mu0)
    cur = essential_margin(n, mu)
    if cur is None:
        return None
    best, best_mu = cur, dict(mu)
    T = 0.02
    decay = (1e-4 / T) ** (1.0 / steps)
    rejected_ess = 0
    for _ in range(steps):
        cand = dict(mu)
        atoms = sorted(cand)
        mv = rng.random()
        if mv < 0.55:
            cand[rng.choice(atoms)] *= math.exp(rng.gauss(0.0, 0.45))
        elif mv < 0.70 and len(cand) < 14:
            u = rng.random()
            if u < 0.4:
                new = rng.randrange(1 << n)
            elif u < 0.8:
                new = rng.choice(atoms) | rng.choice(atoms)
            else:
                new = ((1 << n) - 1) ^ rng.choice(atoms)
            cand[new] = cand.get(new, 0.0) + 0.03 * sum(cand.values())
        elif mv < 0.85 and len(cand) > 3:
            del cand[rng.choice(atoms)]
        else:
            if len(atoms) >= 2:
                a, b = rng.sample(atoms, 2)
                d = 0.3 * rng.random() * cand[a]
                cand[a] -= d
                cand[b] = cand.get(b, 0.0) + d
        v = essential_margin(n, cand)
        if v is None:
            rejected_ess += 1
            continue
        if v < cur or rng.random() < math.exp(-(v - cur) / max(T, 1e-9)):
            mu, cur = cand, v
            if v < best:
                best, best_mu = v, dict(cand)
        T *= decay
    tot = sum(best_mu.values())
    m = {a: w / tot for a, w in best_mu.items()}
    mg = margins(n, m)
    return {"start": name, "n": n, "floor": best,
            "effective_n": sum(1 for x in mg if x > 1e-12),
            "min_marginal": min(mg), "max_marginal": max(mg),
            "rejected_for_essentiality": rejected_ess,
            "mu": {format(a, f"0{n}b"): w for a, w in m.items()}}


def campaign(n, cap, runs, steps, seed):
    rng = random.Random(seed)
    rows, gmin = [], 1e9
    tries = 0
    while len(rows) < runs and tries < 40 * runs:
        tries += 1
        mu = gen_instance(rng, n, cap)
        if essential_margin(n, mu) is None:
            continue
        r = anneal(f"ess{len(rows)}", n, mu, steps, rng)
        if r is None:
            continue
        rows.append(r)
        gmin = min(gmin, r["floor"])
    viol = [r for r in rows if r["floor"] < -1e-12]
    lost = [r for r in rows if r["effective_n"] < n]
    log(f"  n={n} cap {cap}: {len(rows)} anneals x {steps} steps, floor "
        f"{gmin:+.6f}, {len(viol)} violations, {len(lost)} endpoints that "
        f"lost a coordinate (should be 0); moves rejected for "
        f"essentiality: {sum(r['rejected_for_essentiality'] for r in rows)}")
    return {"rows": rows, "floor": gmin, "violations": viol,
            "lost_dimension": len(lost)}


def part_C():
    log()
    log("C. Constrained floors against 040's recorded ones:")
    try:
        old = json.loads((DATA / "hu_roll_anneal.json").read_text())
    except Exception:
        log("  (040 checkpoint unavailable)")
        return
    rows = []
    pairs = [("A_n6_cap_0.49", "A_n6_cap_0.49"),
             ("A_n6_cap_0.497", "A_n6_cap_0.497"),
             ("B_n7_cap_0.49", "A_n7_cap_0.49")]
    for newkey, oldkey in pairs:
        if newkey not in OUT or oldkey not in old:
            continue
        o = old[oldkey]
        oldfloor = o.get("global_floor", o.get("floor"))
        r = min(o["rows"], key=lambda z: z.get("floor", 1e9))
        n = r["n"]
        m = {int(s, 2): w for s, w in r["mu"].items()}
        t = sum(m.values())
        m = {a: w / t for a, w in m.items()}
        mg_old = margins(n, m)
        eff_old = sum(1 for x in mg_old if x > 1e-12)
        # would the recorded endpoint even be admissible here?  If its
        # smallest marginal is below MIN_MARG the constraint excludes it,
        # and the floors are then not comparable like for like -- say so.
        admissible = min(mg_old) >= MIN_MARG
        rows.append({"block": newkey, "constrained_floor": OUT[newkey]["floor"],
                     "recorded_floor": oldfloor, "recorded_effective_n": eff_old,
                     "claimed_n": n, "recorded_min_marginal": min(mg_old),
                     "recorded_endpoint_admissible": admissible})
        log(f"  {newkey}: constrained {OUT[newkey]['floor']:+.6f} vs "
            f"recorded {oldfloor:+.6f} | recorded endpoint: effectively "
            f"n={eff_old} of {n}, min marginal {min(mg_old):.4f} -> "
            f"{'admissible here' if admissible else 'EXCLUDED by the constraint'}")
    OUT["C_comparison"] = rows


def main():
    steps = 900 if not FAST else 120
    log("A. Constrained anneal at n = 6 (040 rerun):")
    for cap in (0.49, 0.497):
        OUT[f"A_n6_cap_{cap}"] = campaign(6, cap, 3 if not FAST else 1,
                                          steps, 8100 + int(cap * 1000))
    log()
    log("B. Constrained anneal at n = 7 (040 rerun -- the worst "
        "overstatement of record):")
    OUT["B_n7_cap_0.49"] = campaign(7, 0.49, 3 if not FAST else 1,
                                    steps, 8200)
    part_C()
    (DATA / "hu_anneal_essential.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_anneal_essential.json")


if __name__ == "__main__":
    main()
