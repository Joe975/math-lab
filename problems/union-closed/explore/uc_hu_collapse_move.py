#!/usr/bin/env python3
"""Attempt 058 (057 lead 2): is the n = 6 -> n = 7 difficulty jump pure
volume, or a missing move?

057 showed the equality family is a stable attractor at n = 7 but that
the anneal does not FIND it unseeded, even at 3x budget.  Its move set
(weight kick, atom add, atom drop, pairwise transfer) never proposes a
diagonal directly: reaching {empty, full} from a 10-atom support needs
eight coordinated drops, each of which the temperature must accept.
So the failure may be a missing move rather than a hard problem.

Falsifiable test: add ONE move -- "collapse toward a diagonal", which
transfers a slice of mass from a random atom onto the empty set or the
full set -- and re-run the same unseeded campaign at n = 7.  If the
anneal now reaches the family, the difficulty was the move set; if it
still does not, 057's volume reading stands.

Parts:
  A  the augmented anneal at n = 7, cap 0.49, same budget as 057
  B  a control: the SAME code with the collapse move disabled, so the
     comparison isolates the move and not the reimplementation
  C  the augmented anneal at n = 6, to confirm the move does not break
     what already worked

Standard library only; deterministic (seeds 8500+).
Usage: python uc_hu_collapse_move.py [--fast]
Checkpoint: ../data/hu_collapse_move.json
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
from uc_hu_anneal_essential import essential_margin, margins
from uc_hu_rollcensus import gen_instance

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def anneal(name, n, mu0, steps, rng, collapse=True):
    """040's move set plus, optionally, a collapse-toward-diagonal move.
    Essentiality is enforced in acceptance exactly as in 055."""
    mu = dict(mu0)
    cur = essential_margin(n, mu)
    if cur is None:
        return None
    best, best_mu = cur, dict(mu)
    full = (1 << n) - 1
    T = 0.02
    decay = (1e-4 / T) ** (1.0 / steps)
    used = 0
    for _ in range(steps):
        cand = dict(mu)
        atoms = sorted(cand)
        mv = rng.random()
        thresh = 0.45 if collapse else 0.55
        if mv < thresh:
            cand[rng.choice(atoms)] *= math.exp(rng.gauss(0.0, 0.45))
        elif collapse and mv < 0.55:
            # THE NEW MOVE: pull a slice of one atom onto empty or full
            a = rng.choice(atoms)
            tgt = full if rng.random() < 0.5 else 0
            d = (0.2 + 0.6 * rng.random()) * cand[a]
            cand[a] -= d
            cand[tgt] = cand.get(tgt, 0.0) + d
            if cand[a] <= 1e-12:
                del cand[a]
            used += 1
        elif mv < 0.70 and len(cand) < 14:
            u = rng.random()
            if u < 0.4:
                new = rng.randrange(1 << n)
            elif u < 0.8:
                new = rng.choice(atoms) | rng.choice(atoms)
            else:
                new = full ^ rng.choice(atoms)
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
            "collapse_moves_used": used,
            "is_diagonal": set(live) <= {0, full} and len(live) == 2,
            "mu": {format(a, f"0{n}b"): w for a, w in m.items()}}


def campaign(tag, n, cap, runs, steps, seed, collapse):
    rng = random.Random(seed)
    rows, gmin = [], 1e9
    tries = 0
    while len(rows) < runs and tries < 60:
        tries += 1
        mu = gen_instance(rng, n, cap)
        if essential_margin(n, mu) is None:
            continue
        r = anneal(f"{tag}{len(rows)}", n, mu, steps, rng, collapse)
        if r is None:
            continue
        rows.append(r)
        gmin = min(gmin, r["floor"])
    reached = gmin < 1e-9
    ndiag = sum(1 for r in rows if r["is_diagonal"])
    log(f"  {tag} n={n} cap {cap} ({'WITH' if collapse else 'without'} "
        f"the collapse move): {len(rows)} runs x {steps} steps, floor "
        f"{gmin:+.3e}, {ndiag} endpoints are diagonals -> "
        f"{'REACHED the family' if reached else 'did not reach it'}")
    return {"rows": rows, "floor": gmin, "reached_family": reached,
            "diagonal_endpoints": ndiag, "collapse": collapse,
            "steps": steps}


def main():
    steps = 2700 if not FAST else 250
    runs = 3 if not FAST else 1
    log("A. n = 7 WITH the collapse move (057's campaign, one move added):")
    OUT["A_n7_collapse"] = campaign("coll", 7, 0.49, runs, steps, 8500, True)
    log()
    log("B. Control -- same code, collapse move disabled:")
    OUT["B_n7_control"] = campaign("ctrl", 7, 0.49, runs, steps, 8500, False)
    log()
    log("C. n = 6 with the move, to check it does not break what worked:")
    OUT["C_n6_collapse"] = campaign("coll6", 6, 0.49, runs, steps, 8600, True)
    (DATA / "hu_collapse_move.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_collapse_move.json")


if __name__ == "__main__":
    main()
