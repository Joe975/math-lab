#!/usr/bin/env python3
"""Attempt 040 (038 lead 3): a STRONGER adversary against rollout-HU at
n = 6 and 7.

038's n >= 6 descents used 035's pattern-descent move set and stalled
~40x above the n <= 5 floors, so "zero violations" there was weak
evidence.  This file attacks the rollout rule with simulated annealing
over measures -- multiplicative weight kicks, atom add (random mask /
union of two existing / complement), atom drop, pairwise transfer --
with delta_0 projection into the regime, followed by a 035-style
pattern-descent polish from the annealed best.  Hostile seeds include
the 035 cap-0.49 kill measures tensored with near-cap Bernoulli
coordinates up to n = 6 (the known-hard geometry, embedded).

Any negative endpoint here is a kill of the rollout rule (subject to
skeptic re-evaluation + certification); zero negatives upgrade 038's
n = 6/7 story from "descent stalled" to "descent AND anneal stalled".

Standard library only; deterministic (seeds 947000+).
Usage: python uc_hu_roll_anneal.py [--fast]
Checkpoint: ../data/hu_roll_anneal.json
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
from uc_hu_order2 import seq_roll, rule_ratio, descend, h

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def roll_ratio(n, mu, cap):
    return rule_ratio(n, mu, cap, seq_roll)


def project(n, mu, cap):
    """delta_0-mix into the regime, normalized."""
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items() if w > 0}
    fm = max(sum(w for a, w in mu.items() if (a >> i) & 1)
             for i in range(n))
    if fm >= cap - 0.003:
        al = 1.0 - (cap - 0.003) / fm
        mu = {a: w * (1 - al) for a, w in mu.items()}
        mu[0] = mu.get(0, 0.0) + al
    return mu


def anneal(name, n, mu0, cap, steps, rng):
    mu = project(n, dict(mu0), cap)
    cur = roll_ratio(n, mu, cap)
    if cur is None:
        return None
    best, best_mu = cur, dict(mu)
    T = 0.02
    decay = (1e-4 / T) ** (1.0 / steps)
    for t in range(steps):
        cand = dict(mu)
        mv = rng.random()
        atoms = sorted(cand)
        if mv < 0.55:                                # weight kick
            a = rng.choice(atoms)
            cand[a] *= math.exp(rng.gauss(0.0, 0.45))
        elif mv < 0.70 and len(cand) < 14:           # add atom
            u = rng.random()
            if u < 0.4:
                new = rng.randrange(1 << n)
            elif u < 0.8:
                new = rng.choice(atoms) | rng.choice(atoms)
            else:
                new = ((1 << n) - 1) ^ rng.choice(atoms)
            cand[new] = cand.get(new, 0.0) + 0.03 * sum(cand.values())
        elif mv < 0.85 and len(cand) > 3:            # drop atom
            del cand[rng.choice(atoms)]
        else:                                        # pairwise transfer
            a, b = rng.sample(atoms, 2) if len(atoms) >= 2 else (atoms[0],) * 2
            d = 0.3 * rng.random() * cand[a]
            cand[a] -= d
            cand[b] = cand.get(b, 0.0) + d
        cand = project(n, cand, cap)
        v = roll_ratio(n, cand, cap)
        if v is None:
            continue
        if v < cur or rng.random() < math.exp(-(v - cur) / max(T, 1e-9)):
            mu, cur = cand, v
            if v < best:
                best, best_mu = v, dict(cand)
        T *= decay
    return best, best_mu


def seeds_for(n, cap, rng):
    """Hostile + random anneal seeds."""
    can = json.loads((DATA / "hu_canon.json").read_text())
    out = []
    q = cap - 0.003
    for r in can["B_cap_0.49"]["violations"]:
        base = {int(s, 2): w for s, w in r["mu"].items()}
        nb = r["n"]
        if nb >= n:
            continue
        ext = {}
        for a, w in base.items():
            for mask in range(1 << (n - nb)):
                pm = 1.0
                for j in range(n - nb):
                    pm *= q if (mask >> j) & 1 else (1.0 - q)
                ext[a | (mask << nb)] = ext.get(a | (mask << nb), 0.0) + w * pm
        # keep support manageable: top 14 atoms, renormalized
        top = sorted(ext.items(), key=lambda kv: -kv[1])[:14]
        tot = sum(w for _, w in top)
        out.append((f"embed:{r['start']}-n{r['n']}",
                    {a: w / tot for a, w in top}))
    for j in range(4 if not FAST else 1):
        k = rng.randint(4, 12)
        supp = rng.sample(range(1 << n), k)
        out.append((f"rand{j}",
                    {a: math.exp(rng.uniform(-2.5, 2.5)) for a in supp}))
    return out


def campaign(n, cap, steps, polish_rounds):
    rng = random.Random(947000 + n * 100 + int(cap * 1000))
    rows, gmin = [], 1e9
    for name, mu0 in seeds_for(n, cap, rng):
        r = anneal(name, n, mu0, cap, steps, rng)
        if r is None:
            continue
        av, amu = r
        d = descend(f"{name}+polish", n, amu, cap,
                    lambda n_, m_, c_: rule_ratio(n_, m_, c_, seq_roll),
                    rounds_max=polish_rounds)
        floor = d["floor"] if d is not None and d["floor"] < av else av
        endmu = d["mu"] if d is not None and d["floor"] < av else {
            format(a, f"0{n}b"): w / sum(amu.values()) for a, w in amu.items()}
        rows.append({"start": name, "n": n, "cap": cap,
                     "anneal_best": av, "floor": floor, "mu": endmu,
                     "atoms": len(endmu)})
        gmin = min(gmin, floor)
        log(f"    {name:22s}: anneal {av:+.6f}, polished {floor:+.6f}")
    viol = [r for r in rows if r["floor"] < 0]
    log(f"  n={n} cap {cap}: {len(rows)} runs, global floor {gmin:+.6f}, "
        f"{len(viol)} violations [product extremal "
        f"{(1 - h(cap)) / h(cap):+.6f}]")
    OUT[f"A_n{n}_cap_{cap}"] = {"rows": rows, "global_floor": gmin,
                                "violations": viol}


def main():
    log("A. Anneal + polish against rollout-HU:")
    if FAST:
        campaign(6, 0.49, 300, 15)
    else:
        campaign(6, 0.49, 4000, 80)
        campaign(6, 0.497, 4000, 80)
        campaign(7, 0.49, 1200, 40)
    (DATA / "hu_roll_anneal.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_roll_anneal.json")


if __name__ == "__main__":
    main()
