#!/usr/bin/env python3
"""Attempt 030 (028 lead 2): K-component latent-mixing -- definition and
first adversarial battery.  Lives in explore/; checkpoint in ../data/.

Definition (K components (P_k, p_k)): couple the labels (k_A, k_B) by a
KxK matrix Q with both marginals P; given (j, j) use the adaptive
per-coordinate joint for p_j (009 part D rule); given (j, k), j != k,
take Bern(p_j)^n x Bern(p_k)^n independent.  Q = diag(P) is exactly the
adaptive-block coupling; K = 2 with a delta_0 light component is 025's
empty-mixing.  Parametrize Q by symmetric transfers t_jk >= 0:
q_jk = t_jk (j != k), q_jj = P_j - sum_{k != j} t_jk >= 0.

Optimizer: deterministic coarse grid over the feasible t-simplex,
then one local refinement around the best (no RNG).

Standard library only; deterministic.
"""
from __future__ import annotations

import itertools
import json
import math
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_mitax import cr_eval, h
from uc_branch_attack import cr_mix

T0 = time.monotonic()
MARG = 0.38271
OUT = {}


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def bern(n, p):
    if p <= 0.0:
        return {0: 1.0}
    if p >= 1.0:
        return {(1 << n) - 1: 1.0}
    return {a: p ** bin(a).count("1") * (1 - p) ** (n - bin(a).count("1"))
            for a in range(1 << n)}


def cell_joint(p):
    if p > 0.5:
        return {(1, 1): p, (0, 0): 1.0 - p}
    if p <= 0.0:
        return {(0, 0): 1.0}
    z = 1.0 - 2.0 * p if p <= 0.25 else 0.5
    off = 1.0 - p - z
    return {(0, 0): z, (0, 1): off, (1, 0): off, (1, 1): p - off}


def product_pairs(n, cell, weight):
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


def latentK_pairs(n, comps, Q):
    """comps: [(P_k, p_k)]; Q: KxK matrix (list of lists), both marginals
    must equal P."""
    K = len(comps)
    pairs = {}

    def acc(d):
        for k, w in d.items():
            if w > 1e-17:
                pairs[k] = pairs.get(k, 0.0) + w

    margs = [bern(n, p) for _, p in comps]
    for j in range(K):
        for k in range(K):
            w = Q[j][k]
            if w <= 1e-15:
                continue
            if j == k:
                acc(product_pairs(n, cell_joint(comps[j][1]), w))
            else:
                for a, wa in margs[j].items():
                    for b, wb in margs[k].items():
                        ww = w * wa * wb
                        if ww > 1e-17:
                            pairs[(a, b)] = pairs.get((a, b), 0.0) + ww
    return pairs


def q_from_t(P, t):
    """t: dict (j,k) j<k -> transfer mass.  Returns Q or None if
    infeasible."""
    K = len(P)
    Q = [[0.0] * K for _ in range(K)]
    for (j, k), v in t.items():
        Q[j][k] = Q[k][j] = v
    for j in range(K):
        d = P[j] - sum(Q[j][k] for k in range(K) if k != j)
        if d < -1e-12:
            return None
        Q[j][j] = max(d, 0.0)
    return Q


def optimize_t(n, comps, steps=7, refine=2):
    """Deterministic grid + refinement over the transfer simplex.
    Returns (best_CR, best_t, diag_CR)."""
    P = [c[0] for c in comps]
    K = len(comps)
    pairs_idx = [(j, k) for j in range(K) for k in range(K) if j < k]
    caps = {jk: min(P[jk[0]], P[jk[1]]) for jk in pairs_idx}

    diag = cr_eval(n, latentK_pairs(n, comps, q_from_t(P, {})))["CR"]

    best = (diag, {jk: 0.0 for jk in pairs_idx})
    grids = [[caps[jk] * i / steps for i in range(steps + 1)]
             for jk in pairs_idx]
    for _ in range(refine + 1):
        for combo in itertools.product(*grids):
            t = dict(zip(pairs_idx, combo))
            Q = q_from_t(P, t)
            if Q is None:
                continue
            cr = cr_eval(n, latentK_pairs(n, comps, Q))["CR"]
            if cr > best[0]:
                best = (cr, t)
        # refine around best
        newgrids = []
        for gi, jk in enumerate(pairs_idx):
            c = best[1][jk]
            w = (grids[gi][1] - grids[gi][0]) if len(grids[gi]) > 1 else \
                caps[jk] / steps
            lo, hi = max(0.0, c - w), min(caps[jk], c + w)
            newgrids.append([lo + (hi - lo) * i / steps
                             for i in range(steps + 1)])
        grids = newgrids
    return best[0], best[1], diag


def main():
    log("K-component latent-mixing: battery (K = 3), in-regime")
    log("=" * 72)

    # K=2 anchor: must reproduce 027's latent-mixing values
    log("Anchor vs 027 (K=2, p_lo/p_hi cell):")
    for plo, n, ref in [(1e-5, 2, 0.22098), (1e-5, 6, 1.18930),
                        (1e-3, 6, 1.16759)]:
        comps = [(0.58, plo), (0.42, 0.85)]
        cr, t, diag = optimize_t(n, comps)
        log(f"  p_lo={plo:.0e} n={n}: CR* = {cr:+.5f} (027: {ref:+.5f}), "
            f"diag(adaptive) = {diag:+.6f}")
        OUT.setdefault("anchor", []).append(
            {"p_lo": plo, "n": n, "CR": cr, "ref_027": ref, "diag": diag})

    insts = [
        ("lo+mid+hi", [(0.5, 0.01), (0.3, 0.55), (0.2, 0.9)]),
        ("spread", [(0.6, 0.05), (0.25, 0.7), (0.15, 0.99)]),
        ("near0+2hi", [(0.62, 0.001), (0.22, 0.9), (0.16, 0.99)]),
        ("all-low+1mid", [(0.4, 0.2), (0.35, 0.3), (0.25, 0.55)]),
        ("near-psi", [(0.4, 0.36), (0.3, 0.39), (0.3, 0.40)]),
    ]
    log()
    log("Battery (CR* = optimized transfers; diag = adaptive-block):")
    rows = []
    for tag, comps in insts:
        marg = sum(P * p for P, p in comps)
        assert marg < MARG, (tag, marg)
        for n in (2, 4, 6):
            cr, t, diag = optimize_t(n, comps)
            r = cr_eval(n, latentK_pairs(
                n, comps, q_from_t([c[0] for c in comps], t)))
            tstr = ", ".join(f"t{j}{k}={v:.3f}"
                             for (j, k), v in t.items() if v > 1e-4)
            log(f"  {tag:12s} n={n}: CR* = {cr:+.5f}  diag = {diag:+.5f}  "
                f"H = {r['H_mu']:.3f}  CR*/H = {cr / r['H_mu']:+.4f}  "
                f"[{tstr or 'diagonal'}]")
            rows.append({"tag": tag, "n": n, "marg": marg, "CR_opt": cr,
                         "CR_diag": diag, "H_mu": r["H_mu"],
                         "cr_over_h": cr / r["H_mu"],
                         "t": {f"{j}{k}": v for (j, k), v in t.items()}})
    OUT["battery"] = rows
    neg = [r for r in rows if r["CR_opt"] <= 0]
    log()
    log(f"{len(rows)} battery rows; {len(neg)} with optimized CR <= 0")

    # degeneration probe: near0+2hi with the light p -> 0; compare the
    # K=3 optimum against the 2-block EM limit where the two hi
    # components fuse into the conditional mixture (delta_0 light).
    log()
    log("Degeneration probe (near0+2hi, light p -> 0), n = 4:")
    from uc_branch_attack import profile_H_mixture_with_empty as PH
    scomps = [(0.22 / 0.38, 0.9), (0.16 / 0.38, 0.99)]
    P1 = 0.62
    qlo = max(0.0, 2 * P1 - 1.0)
    em_lim = max(PH(4, scomps, qlo + (P1 - qlo) * j / 400)
                 - PH(4, scomps, P1) for j in range(401))
    for plight in (1e-2, 1e-3, 1e-4):
        comps = [(0.62, plight), (0.22, 0.9), (0.16, 0.99)]
        cr, t, diag = optimize_t(4, comps)
        log(f"  p_light={plight:.0e}: CR* = {cr:+.5f} (EM limit "
            f"{em_lim:+.5f}), diag = {diag:+.6f}, "
            f"t = {[f'{v:.3f}' for v in t.values()]}")
        OUT.setdefault("degeneration", []).append(
            {"p_light": plight, "CR_opt": cr, "em_limit": em_lim,
             "diag": diag})

    out_path = Path(__file__).resolve().parent.parent / "data" / "latentK_probe.json"
    out_path.write_text(json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log(f"checkpoint: {out_path}")


if __name__ == "__main__":
    main()
