#!/usr/bin/env python3
"""Attempt 028 (027 lead 3 / queue 1(b')): recipe v2 and its battery.

Recipe v2 (branch map; full statement in the record).  On the genres of
record, mu is given WITH its structure (genre detection stays Gap 3):

  G-mix  (product mixtures, incl. delta_0 components and the
          {0}u[1/2,1) genre as the 2-block special case):
          LATENT-MIXING, q optimized -- couple the component labels
          (k_A, k_B) with the correct marginals and free
          P(k_A = k_B = light) = q; equal labels use the adaptive
          per-coordinate joint (009 part D), differing labels take the
          two products independent.  Specializations: q = P_light is
          exactly the adaptive-block coupling (cross cells vanish);
          a delta_0 light component makes it 025's empty-mixing (EM:
          ST = 0, CR closed-form, coverage lemma).  Currently stated
          for two components (+ an optional delta_0 atom); >= 3
          nondegenerate components remain Gap 3 residue.
  G-gen  (generic / slice measures, no declared mixture structure):
          lambda-swept Sinkhorn tilt, grid including lambda = 0.

(TAX at p) v2: for every mu in a stated genre with H(mu) > 0 and all
element marginals < p, CR(mu, pi_v2(mu)) > 0.

This script runs every hard instance of record through the v2
assignment:
  - the two 020 kill geometries and 023's committed floor endpoints
    (generic genre -> Sinkhorn sweep);
  - 025's three certified sliver instances ({0}u[1/2,1) -> latent-mixing
    with delta_0 light component = empty-mixing);
  - 025/026's block kill mixtures (G-mix -> latent-mixing);
  - 027's weakest lo/hi probe points (G-mix -> latent-mixing).
Failure of any instance (CR <= 0) refutes v2 where stated.

Standard library only; deterministic; ~90 s (Sinkhorn rows dominate).
Usage: python uc_recipe_v2.py [--fast]
Checkpoint: ../data/recipe_v2_battery.json
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import uc_or_avg as ORA
from uc_mitax import cr_eval
from uc_branch_attack import cr_mix, F_empty_mix
from uc_adaptive_cell_attack import latent_mixing_pairs, bern
import uc_cr_attack as CA

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
MARG = 0.38271
T0 = time.monotonic()
ROWS = []


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def record(name, genre, branch, cr, detail=""):
    ok = cr > 0
    ROWS.append({"instance": name, "genre": genre, "branch": branch,
                 "CR": cr, "positive": ok})
    log(f"  {name:26s} {genre:5s} -> {branch:14s} CR = {cr:+.5f} "
        f"{'ok' if ok else 'VIOLATION'} {detail}")
    return ok


def sinkhorn_sup(n, mu):
    mu = {a: w for a, w in mu.items() if w > 1e-13}
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    best = None
    for lam in (0.0, 0.05, 0.15, 0.35, 0.8, 1.8, 3.5):
        if lam == 0.0:
            pairs, resid = ({(a, b): wa * wb for a, wa in mu.items()
                             for b, wb in mu.items()}, 0.0)
        else:
            pairs, resid = ORA.sinkhorn(mu, lam, tol=1e-12, iters=30000)
        if resid > 1e-8:
            continue
        r = cr_eval(n, pairs)
        if best is None or r["CR"] > best["CR"]:
            best = r
    return best["CR"]


def latent_sup(n, plo, phi, Plo, qsteps=60):
    qlo = max(0.0, 2 * Plo - 1.0)
    best = -1e9
    for j in range(qsteps + 1):
        q = qlo + (Plo - qlo) * j / qsteps
        r = cr_eval(n, latent_mixing_pairs(n, plo, phi, Plo, q))
        best = max(best, r["CR"])
    return best


def main():
    ok = True
    log("Recipe v2 battery: every hard instance of record, v2 assignment")
    log("=" * 74)

    # -- generic genre: 020 kill geometries (022's check, re-run) -------
    for name, n, mu in CA.seeds():
        tot = sum(mu.values())
        mu = {a: w / tot for a, w in mu.items()}
        ok &= record(name, "G-gen", "sinkhorn-sweep",
                     sinkhorn_sup(n, mu))

    # -- generic genre: 023's committed floor endpoints -----------------
    d = json.loads((DATA / "cr_attack.json").read_text())
    for row in d["rows"]:
        if row["H_mu"] < 0.5:
            continue           # the degenerate corner, excluded by H > 0
        mu = {int(s, 2): w for s, w in row["mu"].items()}
        ok &= record(f"floor:{row['seed']}", "G-gen", "sinkhorn-sweep",
                     sinkhorn_sup(row["n"], mu))

    # -- {0}u[1/2,1) genre: 025/026 sliver instances --------------------
    for n, ph in [(2, 0.9935), (6, 0.9985), (16, 0.99965)]:
        P1 = 0.64
        qlo = max(0.0, 2 * P1 - 1.0)
        best = max(cr_mix(n, ph, P1, qlo + (P1 - qlo) * j / 400)
                   for j in range(401))
        ok &= record(f"sliver n={n}", "G-mix", "latent(EM)", best,
                     f"(closed form; certified at rational neighbors "
                     f"in 026)" if n != 16 else "(closed form)")

    # -- product mixtures: 025/026 block kills --------------------------
    # 2block {0: .62, .95: .38}: delta_0 light component -> EM closed
    # form with nu = Bern(.95)^n
    for n in (4, 8):
        P1 = 0.62
        qlo = max(0.0, 2 * P1 - 1.0)
        best = max(cr_mix(n, 0.95, P1, qlo + (P1 - qlo) * j / 400)
                   for j in range(401))
        ok &= record(f"2block n={n}", "G-mix", "latent(EM)", best)
    # 3block {0:.60, .90:.30, .99:.10}: delta_0 light vs hi-mixture --
    # EM with mixture s-law (025 B3 value at n=6, q swept here)
    from uc_branch_attack import profile_H_mixture_with_empty
    scomps = [(0.75, 0.90), (0.25, 0.99)]
    for n in (4, 6, 8):
        P1 = 0.60
        qlo = max(0.0, 2 * P1 - 1.0)
        best = max(profile_H_mixture_with_empty(
            n, scomps, qlo + (P1 - qlo) * j / 400)
            - profile_H_mixture_with_empty(n, scomps, P1)
            for j in range(401))
        ok &= record(f"3block n={n}", "G-mix", "latent(EM)", best)

    # -- lo/hi cell: 027's weakest probe points -------------------------
    for plo, n in ([(1e-5, 2), (1e-5, 6), (1e-3, 4)] if not FAST
                   else [(1e-3, 4)]):
        ok &= record(f"lo/hi p_lo={plo:.0e} n={n}", "G-mix",
                     "latent-mixing", latent_sup(n, plo, 0.85, 0.58))
    # and the tightest scan instance from 027
    ok &= record("lo/hi tightest(027)", "G-mix", "latent-mixing",
                 latent_sup(2, 0.001, 0.85, 0.58))

    log()
    npos = sum(1 for r in ROWS if r["positive"])
    log(f"{npos}/{len(ROWS)} battery instances positive under the v2 "
        f"assignment")
    out = {"rows": ROWS, "all_positive": ok}
    (DATA / "recipe_v2_battery.json").write_text(
        json.dumps(out, indent=1, default=float) + "\n")
    log("checkpoint: data/recipe_v2_battery.json")
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
