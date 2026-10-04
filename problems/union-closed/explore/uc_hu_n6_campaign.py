#!/usr/bin/env python3
"""Attempt 064 (063 lead 1): an ADVERSARIAL best-order campaign at n = 6,
at roughly rollout cost.

063 censused n = 6 best-order cheaply by noting that rollout
lower-bounds it, and found the bound sufficed on 100% of 6,000 random
instances.  A census is not a campaign, though: random measures are
easy, and every floor this route has recorded came from a search.

The same shortcut makes the campaign affordable.  Anneal on the
ROLLOUT margin -- one order per candidate instead of 720 -- and
enumerate all 720 orders only at the ENDPOINT, to convert the rollout
floor into a best-order floor.  Since best-order >= rollout pointwise,
a non-negative rollout floor already proves the best-order floor is
non-negative; the endpoint enumeration then reports how much slack the
shortcut left.

Parts:
  A  constrained anneals at n = 6, caps 0.49 and 0.497 (essentiality
     of 051 plus the collapse move of 058), objective = rollout margin
  B  endpoint enumeration: the best-order margin at each anneal
     endpoint, and the gap between it and the rollout margin there
  C  the comparison against 063's census minima, which is what the
     campaign has to beat to be worth running

Standard library only; deterministic (seeds 9300+).
Usage: python uc_hu_n6_campaign.py [--fast]
Checkpoint: ../data/hu_n6_campaign.json
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


def best_margin_at(n, m, fmax):
    best, H = None, None
    for s in itertools.permutations(range(n)):
        cr, H = hu_cr_seq(n, m, s)
        if best is None or cr > best:
            best = cr
    return best / H - cstar(fmax)


def anneal(name, n, mu0, steps, rng):
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
    bo = best_margin_at(n, m, fmax)          # the 720-order evaluation
    mg = margins(n, m)
    return {"start": name, "n": n, "roll_floor": best,
            "bestorder_at_endpoint": bo, "gap": bo - best,
            "own_max_marginal": fmax, "min_marginal": min(mg),
            "mu": {format(a, f"0{n}b"): w for a, w in m.items()}}


def campaign(n, cap, runs, steps, seed):
    rng = random.Random(seed)
    rows = []
    tries = 0
    while len(rows) < runs and tries < 60:
        tries += 1
        mu = gen_instance(rng, n, cap)
        if roll_margin(n, mu) is None:
            continue
        r = anneal(f"c{len(rows)}", n, mu, steps, rng)
        if r is not None:
            rows.append(r)
    rf = min(r["roll_floor"] for r in rows)
    bf = min(r["bestorder_at_endpoint"] for r in rows)
    viol = [r for r in rows if r["bestorder_at_endpoint"] < -1e-9]
    log(f"  n={n} cap {cap}: {len(rows)} anneals x {steps} steps | "
        f"rollout floor {rf:+.6f}, best-order at those endpoints "
        f"{bf:+.6f}, {len(viol)} violations")
    log(f"      endpoint gaps (best-order minus rollout): " +
        ", ".join(f"{r['gap']:+.4f}" for r in rows))
    return {"rows": rows, "roll_floor": rf, "bestorder_floor": bf,
            "violations": viol, "steps": steps}


def part_C():
    log()
    log("C. Against the 063 census minima (what the campaign must beat):")
    try:
        cen = json.loads((DATA / "hu_n6_bestorder.json").read_text())
    except Exception:
        return
    rows = []
    for cap, key in ((0.49, "A_n6_cap_0.49"), (0.497, "B_n6_cap_0.497")):
        ck = f"A_n6_cap_{cap}"
        if ck not in OUT or key not in cen:
            continue
        rows.append({"cap": cap, "campaign": OUT[ck]["bestorder_floor"],
                     "census": cen[key]["min_margin"]})
        log(f"  cap {cap}: campaign {OUT[ck]['bestorder_floor']:+.6f} vs "
            f"census {cen[key]['min_margin']:+.6f} "
            f"({'campaign is sharper' if OUT[ck]['bestorder_floor'] < cen[key]['min_margin'] else 'census was already sharper'})")
    OUT["C_vs_census"] = rows


def main():
    steps = 1200 if not FAST else 200
    runs = 4 if not FAST else 1
    log("A/B. Adversarial best-order campaign at n = 6 (anneal on the "
        "rollout margin, enumerate 720 orders at the endpoint):")
    for cap in (0.49, 0.497):
        OUT[f"A_n6_cap_{cap}"] = campaign(6, cap, runs, steps,
                                          9300 + int(cap * 1000))
    part_C()
    (DATA / "hu_n6_campaign.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_n6_campaign.json")


if __name__ == "__main__":
    main()
