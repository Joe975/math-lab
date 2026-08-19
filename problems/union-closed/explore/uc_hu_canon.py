#!/usr/bin/env python3
"""Attempt 035 (033 lead 4 / queue 1(b1) follow-up): CANONICAL-ORDER HU.

033 killed the for-all-orders form of (HU-TAX) at cap 9/20 with a witness
whose bad order reveals its deterministic coordinate LAST (a
surplus-sterile slot).  The dissection's own heuristic -- reveal
(near-)deterministic coordinates FIRST -- becomes a definition here:

  CANONICAL ORDER pi_can(mu): greedily pick, among the unrevealed
  coordinates, one minimising the conditional entropy H(A_i | A_S)
  where S is the already-chosen set (ties -> lowest index).

It depends on mu and the chosen SET only, never on realised values, so
both copies of the coupling use the same order: canonical-HU is a
genuine coupling and still total + closed-form.

Question this file answers: does canonical-HU survive where the
for-all-orders form died?  Same joint adversary as 033 but with the
order fixed by the rule (the adversary keeps the weights).

Standard library only; deterministic (fixed seeds for start shapes only).
Usage: python uc_hu_canon.py [--fast]
Checkpoint: ../data/hu_canon.json
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
from uc_mitax import cr_eval, crash_family, h
from uc_hu_probe import half_union_pairs
from uc_hu_attack import permute_mu, starts

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}
HMIN = 0.2


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def cond_entropy(n, mu, i, S):
    """H(A_i | A_S) in bits, exact over the atom measure."""
    groups = {}
    for a, w in mu.items():
        key = tuple((a >> j) & 1 for j in S)
        g = groups.setdefault(key, [0.0, 0.0])
        g[0] += w
        if not (a >> i) & 1:
            g[1] += w
    tot = 0.0
    for m_, z_ in groups.values():
        if m_ > 0:
            tot += m_ * h(z_ / m_)
    return tot


TIE_TOL = 1e-12


def canonical_order(n, mu):
    """Greedy most-predictable-first; returns perm with perm[i] = slot of
    coordinate i (the convention permute_mu uses).

    TIES ARE REAL (symmetric measures give exactly equal conditional
    entropies) and float noise breaks them differently in different
    implementations -- the 035 skeptic caught exactly that.  So ties are
    detected with an explicit tolerance and then broken by LOWEST INDEX,
    which makes the rule reproducible across implementations."""
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    S = []
    remaining = set(range(n))
    while remaining:
        scored = sorted((cond_entropy(n, mu, i, S), i)
                        for i in sorted(remaining))
        lo = scored[0][0]
        best = min(i for v, i in scored if v <= lo + TIE_TOL)
        S.append(best)
        remaining.discard(best)
    # S[k] is the coordinate revealed at slot k -> perm[coord] = slot
    perm = [0] * n
    for slot, coord in enumerate(S):
        perm[coord] = slot
    return tuple(perm)


def cr_ratio(n, mu, cap, perm=None):
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items() if w / tot > 1e-12}
    fmax = max(sum(w for a, w in m.items() if (a >> i) & 1)
               for i in range(n))
    if fmax >= cap:
        return None
    if perm is not None:
        m = permute_mu(n, m, perm)
    r = cr_eval(n, half_union_pairs(n, m))
    if r["H_mu"] < HMIN:
        return None
    return r["CR"] / r["H_mu"]


def canon_ratio(n, mu, cap):
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items() if w / tot > 1e-12}
    return cr_ratio(n, m, cap, canonical_order(n, m))


def descend_canon(name, n, mu0, cap, rounds_max=200):
    """Weight-only descent; the order is re-derived by the rule at every
    candidate (so the adversary cannot exploit a stale order)."""
    mu = dict(mu0)
    best = canon_ratio(n, mu, cap)
    if best is None:
        return None
    step = 0.5
    rounds = 0
    while step >= 0.03 and rounds < (rounds_max if not FAST else 20):
        improved = False
        for a in sorted(mu):
            for f in (1 + step, 1 / (1 + step)):
                cand = dict(mu)
                cand[a] = cand[a] * f
                v = canon_ratio(n, cand, cap)
                if v is not None and v < best - 1e-9:
                    mu, best, improved = cand, v, True
        tot = sum(mu.values())
        for extra in (0, (1 << n) - 1):
            cand = dict(mu)
            cand[extra] = cand.get(extra, 0.0) + 0.1 * tot
            v = canon_ratio(n, cand, cap)
            if v is not None and v < best - 1e-9:
                mu, best, improved = cand, v, True
        if len(mu) > 3:
            a0 = min(mu, key=mu.get)
            cand = {a: w for a, w in mu.items() if a != a0}
            v = canon_ratio(n, cand, cap)
            if v is not None and v < best - 1e-9:
                mu, best, improved = cand, v, True
        if not improved:
            step *= 0.5
        rounds += 1
    tot = sum(mu.values())
    return {"start": name, "n": n, "cap": cap, "floor": best,
            "atoms": len(mu),
            "mu": {format(a, f"0{n}b"): w / tot for a, w in mu.items()},
            "canon_order": list(canonical_order(n, mu))}


def part_A():
    """The 033 witness under the canonical rule."""
    log("A. The 033 order-kill witness under the canonical rule:")
    d = json.loads((DATA / "hu_attack.json").read_text())
    viol = d["cap_0.45"]["violations"]
    rows = []
    for v in viol:
        n = v["n"]
        mu = {int(s, 2): w for s, w in v["mu"].items()}
        allv = []
        for perm in itertools.permutations(range(n)):
            r = cr_ratio(n, mu, v["cap"], perm)
            if r is not None:
                allv.append((r, perm))
        allv.sort()
        cperm = canonical_order(n, mu)
        cval = canon_ratio(n, mu, v["cap"])
        rank = sum(1 for r, _ in allv if r < cval - 1e-12)
        log(f"  witness ({v['start']}, cap {v['cap']}): "
            f"worst order {allv[0][0]:+.5f} {allv[0][1]}, "
            f"canonical {cval:+.5f} {cperm} "
            f"(rank {rank + 1}/{len(allv)}, best {allv[-1][0]:+.5f})")
        rows.append({"start": v["start"], "cap": v["cap"],
                     "worst_ratio": allv[0][0],
                     "worst_order": list(allv[0][1]),
                     "canon_ratio": cval, "canon_order": list(cperm),
                     "canon_rank": rank + 1, "orders": len(allv),
                     "best_ratio": allv[-1][0]})
    OUT["A_witness"] = rows


def part_B():
    """Adversarial descent against canonical-HU at several caps."""
    log()
    log("B. Weight-descent against CANONICAL-HU (order re-derived each step):")
    caps = (0.38271, 0.45, 0.49) if not FAST else (0.45,)
    for cap in caps:
        rows = []
        gmin = 1e9
        for name, n, mu in starts(cap):
            r = descend_canon(name, n, mu, cap)
            if r is None:
                continue
            rows.append(r)
            gmin = min(gmin, r["floor"])
        ext = (1 - h(cap)) / h(cap)
        viol = [r for r in rows if r["floor"] < 0]
        log(f"  cap {cap}: {len(rows)} starts, global floor "
            f"{gmin:+.6f} [product extremal {ext:+.5f}], "
            f"{len(viol)} violations")
        OUT[f"B_cap_{cap}"] = {"rows": rows, "global_floor": gmin,
                               "product_extremal": ext,
                               "violations": viol}


def part_C():
    """Random in-regime sweep: canonical vs worst vs identity order."""
    log()
    log("C. Fixed-seed sweep, canonical vs worst-order vs identity:")
    for cap in ((0.45, 0.49) if not FAST else (0.45,)):
        nviol_can = nviol_any = 0
        worst_can = (1e9, None)
        gaps = []
        tested = 0
        for seed in range(120 if not FAST else 20):
            rng = random.Random(921000 + seed)
            n = rng.choice([3, 4])
            k = rng.randint(3, min(10, 1 << n))
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
            cval = canon_ratio(n, mu, cap)
            if cval is None:
                continue
            tested += 1
            allv = [cr_ratio(n, mu, cap, p)
                    for p in itertools.permutations(range(n))]
            allv = [v for v in allv if v is not None]
            if not allv:
                continue
            if cval < 0:
                nviol_can += 1
            if min(allv) < 0:
                nviol_any += 1
            gaps.append(cval - min(allv))
            if cval < worst_can[0]:
                worst_can = (cval, seed)
        log(f"  cap {cap}: {tested} instances; canonical negative on "
            f"{nviol_can}; SOME order negative on {nviol_any}; "
            f"worst canonical {worst_can[0]:+.6f} (seed {worst_can[1]}); "
            f"mean canonical-minus-worst gap {sum(gaps)/len(gaps):+.5f}")
        OUT[f"C_cap_{cap}"] = {"tested": tested, "canon_violations": nviol_can,
                               "any_order_violations": nviol_any,
                               "worst_canon": worst_can[0],
                               "worst_seed": worst_can[1],
                               "mean_gap": sum(gaps) / len(gaps)}


def main():
    part_A()
    part_B()
    part_C()
    (DATA / "hu_canon.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_canon.json")


if __name__ == "__main__":
    main()
