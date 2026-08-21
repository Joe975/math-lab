#!/usr/bin/env python3
"""Attempt 053 (052 lead 1): the high-n floors, re-run with no
degenerate escape allowed.

052 showed that the route's n = 6 and n = 7 claims rest on endpoints
that had dropped coordinates -- most severely an "n = 7" anneal floor
that is a two-dimensional instance.  This file re-runs those campaigns
under 051's essentiality constraint (every coordinate must keep
marginal >= MIN_MARG) so the adversary cannot leave the dimension it
was launched in, and asks the falsifiable question: do the floors stay
positive?

At n = 3 the same constraint moved the floor from +2.6e-09 (degenerate)
to +1.8e-05 (essential), i.e. UP -- the degenerate escape was making
the problem look harder than it is.  If that pattern holds at n = 6/7
the floors rise and stay positive; if the constraint instead exposes a
genuinely n-dimensional violation, the promoted conjecture is in
trouble at high n, which would be the most consequential outcome this
route has produced.

Parts:
  A  rollout descents at n = 6, caps 0.49 / 0.497 (038 P2 rerun)
  B  rollout descents at n = 7, cap 0.49 (038 P3 rerun)
  C  a best-order SPOT CHECK at n = 6 (the promoted conjecture
     itself), cap 0.49, all 720 orders enumerated per candidate: one
     start on a short budget, since 720 CR evaluations per candidate
     make a full campaign hours long

All margins are own-constant (044 standard) and every endpoint is
checked essential before it is recorded.

Standard library only; deterministic.
Usage: python uc_hu_essential_highn.py [--fast]
Checkpoint: ../data/hu_essential_highn.json
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
    margs = [sum(w for a, w in m.items() if (a >> i) & 1) for i in range(n)]
    if min(margs) < MIN_MARG or max(margs) >= 0.5:
        return None
    return m, max(margs)


def roll_margin(n, mu):
    p = prep(n, mu)
    if p is None:
        return None
    m, fmax = p
    cr, H = hu_cr_seq(n, m, seq_roll(n, m))
    if H < 0.2:
        return None
    return cr / H - cstar(fmax)


def best_margin_sampled(n, mu, orders):
    p = prep(n, mu)
    if p is None:
        return None
    m, fmax = p
    best, H = None, None
    for s in orders:
        cr, H = hu_cr_seq(n, m, s)
        if best is None or cr > best:
            best = cr
    if H is None or H < 0.2:
        return None
    return best / H - cstar(fmax)


def descend(name, n, mu, obj, rounds=60):
    cur = obj(n, mu)
    if cur is None:
        return None
    step = 0.4
    r = 0
    while step > 0.02 and r < (rounds if not FAST else 12):
        improved = False
        for a in sorted(mu):
            for f in (1 + step, 1 / (1 + step)):
                cand = dict(mu)
                cand[a] = cand[a] * f
                v = obj(n, cand)
                if v is not None and v < cur - 1e-12:
                    mu, cur, improved = cand, v, True
        tot = sum(mu.values())
        for extra in (0, (1 << n) - 1):
            cand = dict(mu)
            cand[extra] = cand.get(extra, 0.0) + 0.08 * tot
            v = obj(n, cand)
            if v is not None and v < cur - 1e-12:
                mu, cur, improved = cand, v, True
        if not improved:
            step *= 0.5
        r += 1
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items()}
    margs = [sum(w for a, w in m.items() if (a >> i) & 1) for i in range(n)]
    return {"start": name, "n": n, "floor": cur,
            "min_marginal": min(margs), "max_marginal": max(margs),
            "effective_n": sum(1 for x in margs if x > 1e-12),
            "mu": {format(a, f"0{n}b"): w for a, w in m.items()}}


def campaign(tag, n, cap, obj, starts, seed, rounds=60):
    rng = random.Random(seed)
    rows, gmin = [], 1e9
    tries = 0
    while len(rows) < starts and tries < 40 * starts:
        tries += 1
        mu = gen_instance(rng, n, cap)
        if prep(n, mu) is None:
            continue
        r = descend(f"{tag}{len(rows)}", n, mu, obj, rounds=rounds)
        if r is None:
            continue
        rows.append(r)
        gmin = min(gmin, r["floor"])
    viol = [r for r in rows if r["floor"] < -1e-12]
    bad_dim = [r for r in rows if r["effective_n"] < n]
    log(f"  {tag} n={n} cap {cap}: {len(rows)} descents, floor "
        f"{gmin:+.6f}, {len(viol)} violations, {len(bad_dim)} endpoints "
        f"that lost a coordinate (should be 0)")
    return {"rows": rows, "floor": gmin, "violations": viol,
            "lost_dimension": len(bad_dim)}


def main():
    log("A. Rollout descents at n = 6, essential support (038 P2 rerun):")
    for cap in (0.49, 0.497):
        OUT[f"A_n6_cap_{cap}"] = campaign("ess", 6, cap, roll_margin,
                                          6 if not FAST else 2,
                                          7100 + int(cap * 1000))
    log()
    log("B. Rollout descents at n = 7, essential support (038 P3 rerun):")
    OUT["B_n7_cap_0.49"] = campaign("ess", 7, 0.49, roll_margin,
                                    4 if not FAST else 2, 7200)
    log()
    log("C. BEST-ORDER descents at n = 6 (the promoted conjecture), "
        "all 720 orders:")
    orders = list(itertools.permutations(range(6)))
    # one start, short budget: at n = 6 the best-order objective costs
    # 720 CR evaluations per candidate, so a full campaign is hours.
    # Reported as a spot check, not a campaign.
    OUT["C_best_n6_cap_0.49"] = campaign(
        "bo", 6, 0.49, lambda n, mu: best_margin_sampled(n, mu, orders),
        1, 7300, rounds=12)
    (DATA / "hu_essential_highn.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_essential_highn.json")


if __name__ == "__main__":
    main()
