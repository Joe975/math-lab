#!/usr/bin/env python3
"""Attempt 027 (025 lead 4 / queue 1(e)): adversarial scan of the last
branch region with neither proof nor kill -- the ADAPTIVE block coupling
on two-component lo/hi product mixtures

    mu = P_lo Bern(p_lo)^n + P_hi Bern(p_hi)^n,
    p_lo in (0, psi], p_hi in (1/2, 1), in-regime marginal
    P_lo p_lo + P_hi p_hi < 0.38271.

(The {0}u[1/2,1) genre and delta_0-plus-anything mixtures are covered by
025's EM lemma; all-above-half mixtures are the certified 025/026 kills;
this lo/hi cell is what remains.)

Mechanism under test (the 025 sliver logic transplanted): the adaptive
coupling's per-coordinate gain slope P_lo [h(m_lo) - h(p_lo)] collapses
as p_lo -> 0 while the mixture-boundary cost stays O(1), and at
p_lo = 0 exactly the coupling degenerates to comonotone (CR = 0) -- so
the question is the sign of CR at small p_lo, small n.

Also swept on any violating instance:
  - the lambda-grid Sinkhorn branch (incl. lambda = 0, the iid coupling);
  - a new LATENT-MIXING family generalizing 025's empty-mixing to
    nondegenerate light components: couple (k_A, k_B) with the correct
    marginals and P(k_A = k_B = lo) = q free; given equal ks use the
    adaptive per-coordinate joint, given differing ks the two products
    are independent.  (No ST = 0 lemma here -- ks differing breaks the
    EM proof -- so its CR is measured, not closed-form.)

Standard library only; deterministic (grid scan, no RNG).
Usage: python uc_adaptive_cell_attack.py [--fast]
Checkpoint: ../data/adaptive_cell.json
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
from uc_mitax import cr_eval, h
from uc_branch_attack import adaptive_block_pairs

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
MARG = 0.38271
PSI = (3.0 - math.sqrt(5.0)) / 2.0
T0 = time.monotonic()
OUT = {}


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def bern(n, p):
    return {a: p ** bin(a).count("1") * (1 - p) ** (n - bin(a).count("1"))
            for a in range(1 << n)}


def cell_joint(p):
    """Adaptive per-coordinate joint for a p-component (009 part D rule)."""
    if p > 0.5:
        return {(1, 1): p, (0, 0): 1.0 - p}
    z = 1.0 - 2.0 * p if p <= 0.25 else 0.5
    off = 1.0 - p - z
    return {(0, 0): z, (0, 1): off, (1, 0): off, (1, 1): p - off}


def product_pairs(n, cell, weight):
    """Coordinate-product of a 2x2 joint, as a pairs dict scaled by
    weight."""
    cur = {(0, 0): weight}
    for i in range(n):
        nxt = {}
        for (a, b), w in cur.items():
            for (ai, bi), cw in cell.items():
                if cw <= 0.0:
                    continue
                k = (a | (ai << i), b | (bi << i))
                nxt[k] = nxt.get(k, 0.0) + w * cw
        cur = nxt
    return cur


def latent_mixing_pairs(n, plo, phi, Plo, q):
    """(k_A, k_B) joint: P(lo,lo)=q, P(lo,hi)=P(hi,lo)=Plo-q,
    P(hi,hi)=1-2Plo+q; equal ks -> adaptive cell product; differing ks
    -> independent products."""
    assert max(0.0, 2 * Plo - 1.0) - 1e-12 <= q <= Plo + 1e-12
    pairs = {}

    def acc(d):
        for k, w in d.items():
            if w > 1e-17:
                pairs[k] = pairs.get(k, 0.0) + w

    acc(product_pairs(n, cell_joint(plo), q))
    acc(product_pairs(n, cell_joint(phi), 1.0 - 2.0 * Plo + q))
    if Plo - q > 1e-15:
        blo, bhi = bern(n, plo), bern(n, phi)
        for a, wa in blo.items():
            for b, wb in bhi.items():
                w = (Plo - q) * wa * wb
                if w > 1e-17:
                    pairs[(a, b)] = pairs.get(pairs_key := (a, b), 0.0) + w
                    pairs[(b, a)] = pairs.get((b, a), 0.0) + w
    return pairs


def sup_sinkhorn_cr(n, mu, lams=(0.0, 0.05, 0.15, 0.35, 0.8, 1.8, 3.5)):
    mu = {a: w for a, w in mu.items() if w > 1e-13}
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    best = None
    for lam in lams:
        if lam == 0.0:
            pairs, resid = ({(a, b): wa * wb for a, wa in mu.items()
                             for b, wb in mu.items()}, 0.0)
        else:
            pairs, resid = ORA.sinkhorn(mu, lam, tol=1e-12, iters=30000)
        if resid > 1e-8:
            continue
        r = cr_eval(n, pairs)
        r["lambda"] = lam
        if best is None or r["CR"] > best["CR"]:
            best = r
    return best


def main():
    log("=" * 74)
    log("Scan: adaptive block coupling on lo/hi mixtures, in-regime")
    log("=" * 74)
    plos = [0.001, 0.003, 0.01, 0.03, 0.1, 0.2, 0.3, 0.38]
    phis = [0.55, 0.7, 0.85, 0.95, 0.99, 0.999]
    Phis = [0.05, 0.1, 0.2, 0.3, 0.42]
    ns = [2, 4, 6] if not FAST else [2, 4]
    rows, neg = [], []
    tested = 0
    for plo in plos:
        for phi in phis:
            for Phi in Phis:
                Plo = 1.0 - Phi
                marg = Plo * plo + Phi * phi
                if marg >= MARG:
                    continue
                comps = [(Plo, plo), (Phi, phi)]
                for n in ns:
                    tested += 1
                    r = cr_eval(n, adaptive_block_pairs(n, comps))
                    row = {"p_lo": plo, "p_hi": phi, "P_hi": Phi, "n": n,
                           "marg": marg, "CR": r["CR"], "gain": r["gain"],
                           "H_mu": r["H_mu"]}
                    rows.append(row)
                    if r["CR"] < -1e-12:
                        neg.append(row)
    log(f"{tested} instances evaluated; {len(neg)} with CR < 0")
    OUT["scan_rows"] = rows
    OUT["negatives"] = sorted(neg, key=lambda r: r["CR"])[:40]
    if neg:
        worst = OUT["negatives"][0]
        log(f"WORST: CR = {worst['CR']:+.5f} at p_lo={worst['p_lo']}, "
            f"p_hi={worst['p_hi']}, P_hi={worst['P_hi']}, n={worst['n']} "
            f"(marg {worst['marg']:.4f}, H {worst['H_mu']:.3f})")
        # refine the worst family: n = 8 and denser p_lo
        log()
        log("Refinement around the worst instance:")
        ref = []
        for n in [2, 3, 4, 6, 8]:
            comps = [(1.0 - worst["P_hi"], worst["p_lo"]),
                     (worst["P_hi"], worst["p_hi"])]
            r = cr_eval(n, adaptive_block_pairs(n, comps))
            ref.append({"n": n, "CR": r["CR"], "gain": r["gain"],
                        "H_mu": r["H_mu"]})
            log(f"    n={n}: CR = {r['CR']:+.6f} (gain {r['gain']:+.5f}, "
                f"H {r['H_mu']:.3f})")
        OUT["refine_n"] = ref

        # rescues on the worst instance at its scan n
        n = worst["n"]
        comps = [(1.0 - worst["P_hi"], worst["p_lo"]),
                 (worst["P_hi"], worst["p_hi"])]
        mu = {}
        for P, p in comps:
            for a, w in bern(n, p).items():
                mu[a] = mu.get(a, 0.0) + P * w
        best = sup_sinkhorn_cr(n, mu)
        log(f"    Sinkhorn sweep rescue: sup_lam CR = {best['CR']:+.5f} "
            f"at lam={best['lambda']}")
        Plo = 1.0 - worst["P_hi"]
        qlo = max(0.0, 2 * Plo - 1.0)
        bq = None
        for j in range(101):
            q = qlo + (Plo - qlo) * j / 100.0
            r = cr_eval(n, latent_mixing_pairs(
                n, worst["p_lo"], worst["p_hi"], Plo, q))
            if bq is None or r["CR"] > bq[0]:
                bq = (r["CR"], q, r["marg_dev"])
        log(f"    latent-mixing rescue: sup_q CR = {bq[0]:+.5f} at "
            f"q={bq[1]:.3f} (marg_dev {bq[2]:.1e})")
        OUT["rescue"] = {"sinkhorn_CR": best["CR"],
                         "sinkhorn_lam": best["lambda"],
                         "latent_mixing_CR": bq[0],
                         "latent_mixing_q": bq[1]}
    else:
        # survival: report the tightest margins
        tight = sorted(rows, key=lambda r: r["CR"])[:10]
        log("No violation; tightest instances:")
        for r in tight:
            log(f"    CR = {r['CR']:+.5f} at p_lo={r['p_lo']}, "
                f"p_hi={r['p_hi']}, P_hi={r['P_hi']}, n={r['n']}")
        OUT["tightest"] = tight
    # scaling probe: p_lo -> 0 at fixed (p_hi, P_hi).  The adaptive CR
    # must -> 0 (comonotone limit) while H(mu) stays Theta(1) -- so the
    # branch's CR/H collapses.  Does the latent-mixing family hold the
    # 0-limit's empty-mixing value instead?
    log()
    log("Scaling probe at (p_hi=0.85, P_hi=0.42):")
    log(f"    {'p_lo':>8} {'n':>2} {'CR_adaptive':>12} {'CR/H':>9} "
        f"{'CR_latmix':>10} {'q*':>5} {'H_mu':>6}")
    phi, Phi = 0.85, 0.42
    Plo = 1.0 - Phi
    probe = []
    for plo in [1e-5, 1e-4, 1e-3, 1e-2, 3e-2]:
        for n in ([2, 4] if FAST else [2, 4, 6]):
            comps = [(Plo, plo), (Phi, phi)]
            r = cr_eval(n, adaptive_block_pairs(n, comps))
            qlo = max(0.0, 2 * Plo - 1.0)
            bq = None
            for j in range(51):
                q = qlo + (Plo - qlo) * j / 50.0
                rr = cr_eval(n, latent_mixing_pairs(n, plo, phi, Plo, q))
                if bq is None or rr["CR"] > bq[0]:
                    bq = (rr["CR"], q)
            probe.append({"p_lo": plo, "n": n, "CR_adaptive": r["CR"],
                          "cr_over_h": r["CR"] / r["H_mu"],
                          "CR_latmix": bq[0], "q_star": bq[1],
                          "H_mu": r["H_mu"]})
            log(f"    {plo:8.0e} {n:2d} {r['CR']:+12.6f} "
                f"{r['CR']/r['H_mu']:+9.5f} {bq[0]:+10.5f} {bq[1]:5.2f} "
                f"{r['H_mu']:6.3f}")
    OUT["scaling_probe"] = probe
    # the p_lo = 0 limit is the empty-mixing genre: closed-form value
    from uc_branch_attack import cr_mix
    qgrid = [max(0.0, 2 * Plo - 1.0) + j / 200.0 *
             (Plo - max(0.0, 2 * Plo - 1.0)) for j in range(201)]
    for n in [2, 4, 6]:
        lim = max(cr_mix(n, phi, Plo, q) for q in qgrid)
        log(f"    p_lo=0 limit (EM closed form), n={n}: sup_q CR = "
            f"{lim:+.5f}")
        OUT.setdefault("em_limit", {})[str(n)] = lim

    (DATA / "adaptive_cell.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/adaptive_cell.json")


if __name__ == "__main__":
    main()
