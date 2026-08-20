#!/usr/bin/env python3
"""Attempt 033 part A (queue 1(b1)): adversarial campaign against the
half-union coupling as primary target -- minimize CR_HU/H(mu) over
(mu, coordinate order) jointly, in-regime, H >= 0.2.

(HU-TAX) is posed for every order, so the adversary minimizes over
orders too: the objective at a candidate mu is min over an order pool
(identity + orders found by greedy adjacent swaps).  Weight moves are
deterministic pattern search; starts are the recorded hard instances
plus fixed-seed hostile shapes.  No anneal, no unseeded RNG.

Result (033): no violation at the working cap 0.38271; a certified
order-dependent violation at cap 0.45 (witness saved in the checkpoint,
certified in uc_hu_certify.py).

Usage: python uc_hu_attack.py [--fast]
Checkpoint: ../data/hu_attack.json
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
from uc_mitax import cr_eval, crash_family, h
from uc_hu_probe import half_union_pairs

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}
HMIN = 0.2


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def permute_mu(n, mu, perm):
    out = {}
    for a, w in mu.items():
        b = 0
        for i in range(n):
            if (a >> i) & 1:
                b |= 1 << perm[i]
        out[b] = out.get(b, 0.0) + w
    return out


def hu_cr_h(n, mu, cap):
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items() if w / tot > 1e-12}
    fmax = max(sum(w for a, w in m.items() if (a >> i) & 1)
               for i in range(n))
    if fmax >= cap:
        return None
    r = cr_eval(n, half_union_pairs(n, m))
    if r["H_mu"] < HMIN:
        return None
    return r["CR"] / r["H_mu"]


def min_over_orders(n, mu, cap, pool):
    vals = []
    for perm in pool:
        v = hu_cr_h(n, permute_mu(n, mu, perm), cap)
        if v is None:
            return None
        vals.append(v)
    return min(vals)


def improve_order_pool(n, mu, cap, pool):
    """Greedy adjacent swaps from the pool's current best order."""
    best_perm = min(pool, key=lambda p: hu_cr_h(
        n, permute_mu(n, mu, p), cap) or 1e9)
    perm = list(best_perm)
    changed = True
    while changed:
        changed = False
        base = hu_cr_h(n, permute_mu(n, mu, tuple(perm)), cap)
        if base is None:
            break
        for i in range(n - 1):
            cand = perm[:]
            cand[i], cand[i + 1] = cand[i + 1], cand[i]
            v = hu_cr_h(n, permute_mu(n, mu, tuple(cand)), cap)
            if v is not None and v < base - 1e-10:
                perm, changed = cand, True
                break
    t = tuple(perm)
    if t not in pool:
        pool.append(t)
    return pool


def descend(name, n, mu0, cap):
    pool = [tuple(range(n))]
    mu = dict(mu0)
    pool = improve_order_pool(n, mu, cap, pool)
    best = min_over_orders(n, mu, cap, pool)
    if best is None:
        return None
    step = 0.5
    rounds = 0
    while step >= 0.03 and rounds < (200 if not FAST else 20):
        improved = False
        for a in sorted(mu):
            for f in (1 + step, 1 / (1 + step)):
                cand = dict(mu)
                cand[a] = cand[a] * f
                v = min_over_orders(n, cand, cap, pool)
                if v is not None and v < best - 1e-9:
                    mu, best, improved = cand, v, True
        # atom moves: add delta_0 / full-set / drop lightest
        tot = sum(mu.values())
        for extra in (0, (1 << n) - 1):
            cand = dict(mu)
            cand[extra] = cand.get(extra, 0.0) + 0.1 * tot
            v = min_over_orders(n, cand, cap, pool)
            if v is not None and v < best - 1e-9:
                mu, best, improved = cand, v, True
        if len(mu) > 3:
            a0 = min(mu, key=mu.get)
            cand = {a: w for a, w in mu.items() if a != a0}
            v = min_over_orders(n, cand, cap, pool)
            if v is not None and v < best - 1e-9:
                mu, best, improved = cand, v, True
        pool = improve_order_pool(n, mu, cap, pool)
        v = min_over_orders(n, mu, cap, pool)
        if v is not None:
            best = min(best, v)
        if not improved:
            step *= 0.5
        rounds += 1
    tot = sum(mu.values())
    worst_perm = min(pool, key=lambda pp: hu_cr_h(
        n, permute_mu(n, mu, pp), cap) if hu_cr_h(
        n, permute_mu(n, mu, pp), cap) is not None else 1e9)
    return {"start": name, "n": n, "cap": cap, "floor": best,
            "atoms": len(mu), "orders_tried": len(pool),
            "mu": {format(a, f"0{n}b"): w / tot for a, w in mu.items()},
            "worst_order": list(worst_perm),
            "identity_value": hu_cr_h(n, mu, cap)}


def starts(cap):
    out = []
    d = json.loads((DATA / "cr_attack.json").read_text())
    for row in d["rows"]:
        if row["H_mu"] < 0.5:
            continue
        out.append((f"floor:{row['seed']}", row["n"],
                    {int(s, 2): w for s, w in row["mu"].items()}))
    for n in (5, 8):
        atoms, logm = crash_family(n)
        out.append((f"crash{n}", n,
                    {a: 2.0 ** lm for a, lm in zip(atoms, logm)}))
    # near-product at the cap (the conjectured extremal), perturbed
    for n in (4, 5):
        p = cap - 0.01
        mu = {a: p ** bin(a).count("1") * (1 - p) ** (n - bin(a).count("1"))
              for a in range(1 << n)}
        out.append((f"product{n}", n, mu))
    # fixed-seed hostile shapes
    for seed in range(8 if not FAST else 2):
        rng = random.Random(917000 + seed)
        n = rng.choice([4, 5, 6])
        k = rng.randint(3, min(12, 1 << n))
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
        out.append((f"seed{seed}", n, mu))
    return out


def main():
    for cap in ((0.38271, 0.45) if not FAST else (0.38271,)):
        log(f"cap = {cap}: joint (mu, order) descent, floors:")
        rows = []
        gmin = 1e9
        for name, n, mu in starts(cap):
            r = descend(name, n, mu, cap)
            if r is None:
                log(f"  {name:20s}: start left the feasible region")
                continue
            rows.append(r)
            gmin = min(gmin, r["floor"])
            log(f"  {name:20s} n={r['n']}: floor CR_HU/H = "
                f"{r['floor']:+.5f} ({r['orders_tried']} orders)")
        ext = (1 - h(cap)) / h(cap)
        log(f"  GLOBAL floor at cap {cap}: {gmin:+.5f}   "
            f"[product extremal {ext:+.5f}]")
        viol = [r for r in rows if r["floor"] < 0]
        OUT[f"cap_{cap}"] = {"rows": rows, "global_floor": gmin,
                             "product_extremal": ext,
                             "violations": viol}
        if viol:
            log(f"  VIOLATION at cap {cap}: {len(viol)} start(s); "
                f"witness in the checkpoint (see 033 record)")
    (DATA / "hu_attack.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log("checkpoint: data/hu_attack.json")


if __name__ == "__main__":
    main()
