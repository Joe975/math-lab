#!/usr/bin/env python3
"""Attempt 027 skeptic pass.  Default stance: refute.

Independent of the engine (uc_adaptive_cell_attack.py):
  - couplings rebuilt from 009 part D's verbal rule by a direct
    (A, B)-mask double loop with per-coordinate cell probabilities
    multiplied out (no product-DP, no shared builder);
  - CR evaluated with uc_branch_skeptic.cr_chain (the 025 skeptic's
    natural-log, B-major evaluator -- shares nothing with uc_mitax);
  - marginal symmetry of the latent-mixing build checked explicitly;
  - EM limits recomputed by brute-force entropy difference.

Reads data/adaptive_cell.json only to compare claims.
Usage: python uc_adaptive_cell_skeptic.py
Output: ../data/adaptive_cell_skeptic_out.txt; exit 1 on refutation.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_branch_skeptic import cr_chain, bern_measure, mix, H_of, LN2

DATA = HERE.parent / "data"
FAILS = []
LINES = []


def say(m=""):
    print(m, flush=True)
    LINES.append(m)


def check(name, ok, detail=""):
    say(f"  [{'ok' if ok else 'REFUTED'}] {name}"
        + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def cell(p):
    """Adaptive per-coordinate joint, from 009 part D's verbal rule."""
    if p > 0.5:
        return {(0, 0): 1.0 - p, (1, 1): p}
    m = 2.0 * p if p <= 0.25 else 0.5     # adaptive union marginal
    z = 1.0 - m
    off = 1.0 - p - z
    return {(0, 0): z, (0, 1): off, (1, 0): off, (1, 1): p - off}


def pair_mass(a, b, n, cj):
    w = 1.0
    for i in range(n):
        w *= cj.get(((a >> i) & 1, (b >> i) & 1), 0.0)
        if w == 0.0:
            return 0.0
    return w


def adaptive_pairs(n, comps):
    pairs = {}
    for P, p in comps:
        cj = cell(p)
        for a in range(1 << n):
            for b in range(1 << n):
                w = P * pair_mass(a, b, n, cj)
                if w > 0.0:
                    pairs[(a, b)] = pairs.get((a, b), 0.0) + w
    return pairs


def latent_pairs(n, plo, phi, Plo, q):
    pairs = {}
    blo, bhi = bern_measure(n, plo), bern_measure(n, phi)
    cells = [(q, cell(plo), None), (1.0 - 2.0 * Plo + q, cell(phi), None),
             (Plo - q, None, (blo, bhi)), (Plo - q, None, (bhi, blo))]
    for wt, cj, cross in cells:
        if wt <= 0.0:
            continue
        if cj is not None:
            for a in range(1 << n):
                for b in range(1 << n):
                    w = wt * pair_mass(a, b, n, cj)
                    if w > 0.0:
                        pairs[(a, b)] = pairs.get((a, b), 0.0) + w
        else:
            ma, mb = cross
            for a, wa in ma.items():
                for b, wb in mb.items():
                    w = wt * wa * wb
                    if w > 0.0:
                        pairs[(a, b)] = pairs.get((a, b), 0.0) + w
    return pairs


def main():
    say("Skeptic pass on uc_adaptive_cell_attack.py")
    say("=" * 70)
    ref = json.loads((DATA / "adaptive_cell.json").read_text())

    say("\nS1. The ten tightest scan instances (claimed CR > 0):")
    for row in ref["tightest"]:
        comps = [(1.0 - row["P_hi"], row["p_lo"]),
                 (row["P_hi"], row["p_hi"])]
        r = cr_chain(row["n"], adaptive_pairs(row["n"], comps))
        check(f"p_lo={row['p_lo']}, p_hi={row['p_hi']}, "
              f"P_hi={row['P_hi']}, n={row['n']}",
              r["CR"] > 0 and abs(r["CR"] - row["CR"]) < 1e-10,
              f"CR={r['CR']:+.6f} vs engine {row['CR']:+.6f}")

    say("\nS2. Scaling probe rows (adaptive CR -> 0, latent-mix holds):")
    for row in [r for r in ref["scaling_probe"]
                if r["p_lo"] in (1e-5, 1e-3) and r["n"] in (2, 6)]:
        comps = [(0.58, row["p_lo"]), (0.42, 0.85)]
        r = cr_chain(row["n"], adaptive_pairs(row["n"], comps))
        check(f"adaptive CR at p_lo={row['p_lo']:.0e}, n={row['n']}",
              abs(r["CR"] - row["CR_adaptive"]) < 1e-9
              and r["CR"] < 0.01 * row["H_mu"] + 1e-9,
              f"CR={r['CR']:+.7f}, H={r['H_mu']:.3f}")
        lp = latent_pairs(row["n"], row["p_lo"], 0.85, 0.58,
                          row["q_star"])
        # marginal symmetry of the independent latent build
        ma, mb = {}, {}
        for (a, b), w in lp.items():
            ma[a] = ma.get(a, 0.0) + w
            mb[b] = mb.get(b, 0.0) + w
        mdev = max(abs(ma.get(k, 0.0) - mb.get(k, 0.0))
                   for k in set(ma) | set(mb))
        r2 = cr_chain(row["n"], lp)
        check(f"latent-mix CR at p_lo={row['p_lo']:.0e}, n={row['n']}, "
              f"q={row['q_star']:.2f}",
              abs(r2["CR"] - row["CR_latmix"]) < 1e-9 and mdev < 1e-12
              and r2["CR"] > 0.2,
              f"CR={r2['CR']:+.5f} vs engine {row['CR_latmix']:+.5f}, "
              f"marg_dev={mdev:.1e}")

    say("\nS3. EM limits by brute-force entropy difference:")
    for n, lim in sorted(ref["em_limit"].items(), key=lambda kv: int(kv[0])):
        n = int(n)
        s = bern_measure(n, 0.85)
        qlo = max(0.0, 2 * 0.58 - 1.0)
        best = max((H_of(mix([({0: 1.0}, q), (s, 1.0 - q)]))
                    - H_of(mix([({0: 1.0}, 0.58), (s, 0.42)]))) / LN2
                   for q in [qlo + j * (0.58 - qlo) / 400 for j in
                             range(401)])
        check(f"EM limit at n={n}", abs(best - lim) < 2e-4,
              f"{best:+.5f} vs engine {lim:+.5f}")
        # convergence: latent-mix at p_lo=1e-5 within 2e-3 of the limit
        pr = [r for r in ref["scaling_probe"]
              if r["p_lo"] == 1e-5 and r["n"] == n]
        if pr:
            check(f"latent-mix(p_lo=1e-5) -> EM limit, n={n}",
                  abs(pr[0]["CR_latmix"] - lim) < 2e-3,
                  f"delta={pr[0]['CR_latmix'] - lim:+.1e}")

    say("\n" + "=" * 70)
    if FAILS:
        say(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            say(f"  - {f}")
    else:
        say("No refutation found.")
    (DATA / "adaptive_cell_skeptic_out.txt").write_text(
        "\n".join(LINES) + "\n")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
