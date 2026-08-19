#!/usr/bin/env python3
"""Attempt 023 (022 lead 2): adversarial search for CR < 0 — attack Gap 2's
(TAX at p) candidate where Gap 1 died.

The candidate: for every mu with H(mu) > 0 and element marginals < 0.38271,
the recipe coupling has CR(mu, pi(mu)) > 0 (009/011).  The recipe takes the
best lambda, so a kill needs sup_lambda CR < 0.  Seeds: the two 020 kill
geometries (which annihilate every per-history OR condition but leave CR
at +0.40/+0.80), plus random supports.  Objective: minimize
sup_{lambda in grid} CR + marginal penalty, free-support anneal with
atom-count mutations.

Standard library only; deterministic.
Usage: python uc_cr_attack.py [--fast]
Checkpoint: ../data/cr_attack.json
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
import uc_or_avg as ORA
from uc_mitax import cr_eval

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
MARG = 0.38271
T0 = time.monotonic()
LAM_GRID = (0.05, 0.15, 0.35, 0.8, 1.8, 3.5)


def log(m):
    print(f"[{time.monotonic() - T0:6.0f}s] {m}", flush=True)


def sup_cr(mu: dict[int, float], n: int):
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    best = None
    for lam in LAM_GRID:
        pi, resid = ORA.sinkhorn(mu, lam, tol=1e-12, iters=30000)
        if resid > 1e-8:
            continue
        ev = cr_eval(n, pi)
        if best is None or ev["CR"] > best["CR"]:
            best = ev
            best["lambda"] = lam
    return best


def seeds():
    out = []
    d = json.loads((DATA / "gap1deep_partB.json").read_text())
    rB = [x for x in d["descents"] if x["a2_end"] < 0][0]
    out.append(("windowkill", 4,
                {int(s, 2): w for s, w in rB["mu"].items()}))
    d2 = json.loads((DATA / "gap1deep_partC.json").read_text())
    allv = sorted(d2["verified"],
                  key=lambda v: min(p["mm_abs"] for p in v["points"]))
    src = [q for q in d2["designs"]
           if abs(q["model_c1"] - allv[0]["model_c1"]) < 1e-12][0]
    out.append(("mmabskill", 5,
                {int(s, 2): w for s, w in src["mu"].items()}))
    return out


def main():
    rows = []
    steps = 1200 if not FAST else 40
    jobs = []
    for name, n, mu in seeds():
        jobs.append((name, n, mu))
        jobs.append((name + "_n6", 6, {a: w for a, w in mu.items()}))
    for k in range(4 if not FAST else 1):
        rng = random.Random(707000 + k)
        n = (4, 5, 6, 6)[k]
        jobs.append((f"random{k}", n,
                     {a: math.exp(rng.uniform(-1.5, 1.5))
                      for a in rng.sample(range(1 << n), 10)}))

    global_best = None
    for name, n, mu0 in jobs:
        rng = random.Random(808000 + hash(name) % 9999)
        mu = dict(mu0)
        ev = sup_cr(mu, n)
        if ev is None:
            continue

        def obj(e):
            return e["CR"] + 50.0 * max(0.0, e["max_freq"] - MARG)

        cur = obj(ev)
        best_local = (dict(mu), ev) if ev["max_freq"] <= MARG else (None, None)
        for step in range(steps):
            sigma = 0.6 * (1 - step / steps) + 0.05
            cand = dict(mu)
            r = rng.random()
            if r < 0.65 or len(cand) <= 4:
                a = rng.choice(list(cand))
                cand[a] = max(cand[a] * math.exp(rng.gauss(0, sigma)), 1e-6)
            elif r < 0.85 and len(cand) < 14:
                cand[rng.getrandbits(n)] = math.exp(rng.uniform(-2, 0))
            else:
                del cand[rng.choice(list(cand))]
            ev2 = sup_cr(cand, n)
            if ev2 is None:
                continue
            v2 = obj(ev2)
            if v2 < cur:
                mu, cur, ev = cand, v2, ev2
                if (ev["max_freq"] <= MARG
                        and (best_local[1] is None
                             or ev["CR"] < best_local[1]["CR"])):
                    best_local = (dict(mu), dict(ev))
        u_best, e_best = best_local
        if e_best is None:
            log(f"{name}: no in-regime endpoint")
            continue
        row = {"seed": name, "n": n, "sup_cr": e_best["CR"],
               "at_lambda": e_best["lambda"], "gain": e_best["gain"],
               "second_tax": e_best["second_tax"],
               "max_freq": e_best["max_freq"],
               "H_mu": e_best["H_mu"], "atoms": len(u_best),
               "mu": {format(a, f"0{n}b"): w for a, w in u_best.items()}}
        rows.append(row)
        log(f"{name:12s} n={n}: best sup_lam CR {e_best['CR']:+.6e} "
            f"(gain {e_best['gain']:+.4e}, H {e_best['H_mu']:.3f}, "
            f"mf {e_best['max_freq']:.3f}, {len(u_best)} atoms)")
        if global_best is None or e_best["CR"] < global_best["sup_cr"]:
            global_best = row
    out = {"rows": rows, "best": global_best, "lam_grid": LAM_GRID,
           "steps": steps}
    (DATA / "cr_attack.json").write_text(json.dumps(out, indent=1) + "\n")
    if global_best:
        log(f"GLOBAL BEST sup_lam CR {global_best['sup_cr']:+.6e} "
            f"({global_best['seed']}, n={global_best['n']})")


if __name__ == "__main__":
    main()
