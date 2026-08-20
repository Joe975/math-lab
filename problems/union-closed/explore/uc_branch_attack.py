#!/usr/bin/env python3
"""Attempt 025 (023 lead 2 / queue 1b pre-step): adversarial pass on the
009 recipe's two non-Sinkhorn branches -- half-mixing and block-adaptive --
which no attack had exercised, plus the CR/H floor ratios (lead 1 pre-work).

Half-mixing genre closed form (ST = 0 by 009's HM lemma, VERIFIED in 011,
extended here to arbitrary (eps_A, eps_B) joints and arbitrary s-laws):
for mu = P1 delta_empty + (1-P1) law(s) and the coupling A = eps_A s,
B = eps_B s with P(eps_A = eps_B = 0) = q,

    CR(q) = F(q) - F(P1),   F(w) = H(w delta_empty + (1-w) law(s)).

The 004/009 branch fixes q = P1/2.  F is strictly concave in w, so
sup_q CR > 0 iff F'(P1) < 0 iff the empty set's surprisal under mu is
below the nu-average surprisal (nu = law(s)); the sliver where the stated
branch fails is p_h -> 1 in-regime.

Block branch: per-coordinate iid-given-k gain is h((1-p)^2) - h(p), which
changes sign exactly at p = psi = (3-sqrt5)/2; mixtures with above-psi
components go Theta(n) negative at profile level (Gain >= CR kills CR).

Standard library only; deterministic (no RNG anywhere).
Usage: python uc_branch_attack.py [--fast]
Checkpoint: ../data/branch_attack.json
"""
from __future__ import annotations

import json
import math
import sys
import time
from math import log2
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import uc_or_avg as ORA
from uc_mitax import cr_eval, h

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
MARG = 0.38271
PSI = (3.0 - math.sqrt(5.0)) / 2.0
T0 = time.monotonic()
OUT = {}


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def xlog2x(x):
    return 0.0 if x <= 0.0 else x * log2(x)


# ---------------------------------------------------------------- closed forms
def F_empty_mix(n, ph, w):
    """H( w*delta_empty + (1-w)*Bern(ph)^n ), exactly (float)."""
    beta = (1.0 - ph) ** n
    m0 = w + (1.0 - w) * beta
    return (-xlog2x(m0)
            - (1.0 - w) * (1.0 - beta) * log2(1.0 - w)
            + (1.0 - w) * (n * h(ph) + xlog2x(beta)))


def cr_mix(n, ph, P1, q):
    """CR of the generalized empty-mixing coupling (ST = 0 lemma)."""
    return F_empty_mix(n, ph, q) - F_empty_mix(n, ph, P1)


def gen_mixing_pairs(n, ph, P1, q):
    """A = eps_A s, B = eps_B s; P(0,0)=q, P(0,1)=P(1,0)=P1-q,
    P(1,1)=1-2P1+q; s ~ Bern(ph)^n.  q = P1/2 is 004's half-mixing."""
    assert max(0.0, 2 * P1 - 1.0) - 1e-12 <= q <= P1 + 1e-12
    pairs = {}

    def add(k, w):
        pairs[k] = pairs.get(k, 0.0) + w

    for s in range(1 << n):
        w1 = bin(s).count("1")
        ps = ph ** w1 * (1.0 - ph) ** (n - w1)
        add((0, 0), ps * q)
        add((0, s), ps * (P1 - q))
        add((s, 0), ps * (P1 - q))
        add((s, s), ps * (1.0 - 2.0 * P1 + q))
    return pairs


def solve_ph(target_h):
    """ph in [1/2, 1) with h(ph) = target_h (bisection, 200 steps)."""
    lo, hi = 0.5, 1.0 - 1e-15
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if h(mid) > target_h:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ------------------------------------------------------------------- part A
def part_A():
    log("=" * 74)
    log("A. Half-mixing branch (q = P1/2): in-regime failure map + repair")
    log("=" * 74)

    # A1: cross-validate the closed form against the full-history evaluator
    # (cr_eval: independent implementation, validated in 009/011/024), and
    # check ST = 0 numerically for generalized q (extended HM lemma).
    battery, worst_cr, worst_st, worst_md = [], 0.0, 0.0, 0.0
    for n, ph, P1, q in [(4, 0.60, 0.50, 0.25), (6, 0.60, 0.50, 0.30),
                         (6, 0.99, 0.64, 0.32), (6, 0.99, 0.64, 0.50),
                         (2, 0.9935, 0.64, 0.32), (2, 0.9935, 0.64, 0.50),
                         (8, 0.75, 0.40, 0.10), (5, 0.51, 0.163, 0.0815),
                         (6, 0.9985, 0.64, 0.32), (6, 0.9985, 0.64, 0.50),
                         (7, 0.55, 0.6666, 0.3333), (7, 0.55, 0.6666, 0.50)]:
        r = cr_eval(n, gen_mixing_pairs(n, ph, P1, q))
        cf = cr_mix(n, ph, P1, q)
        dev = abs(r["CR"] - cf)
        worst_cr = max(worst_cr, dev)
        worst_st = max(worst_st, abs(r["second_tax"]))
        worst_md = max(worst_md, r["marg_dev"])
        battery.append({"n": n, "ph": ph, "P1": P1, "q": q,
                        "CR_eval": r["CR"], "CR_closed": cf,
                        "second_tax": r["second_tax"],
                        "tax_A": r["tax_A"], "marg_dev": r["marg_dev"]})
    log(f"A1 battery ({len(battery)} points): max |CR_eval - closed| = "
        f"{worst_cr:.3e}; max |ST| = {worst_st:.3e} (extended-HM lemma "
        f"prediction: 0); max marginal dev = {worst_md:.3e}")
    OUT["A1_battery"] = battery
    OUT["A1_max_dev"] = {"cr": worst_cr, "st": worst_st, "md": worst_md}

    # A2: the failure sliver, per n.  In-regime: (1-P1) ph < MARG,
    # genre ph in [1/2, 1); branch validity q = P1/2 needs P1 <= 2/3.
    log()
    log("A2: q = P1/2 failure instances, one per n (P1 = 0.64; ph chosen so")
    log("    n h(ph) < 2[h(P1)-h(P1/2)]/P1; margin shown is CR_hm < 0):")
    P1 = 0.64
    cost = h(P1) - h(P1 / 2.0)  # fixed mixing cost of the stated branch
    fam = []
    ns = [1, 2, 3, 4, 6, 8, 12, 16, 24, 40, 64, 100]
    for n in ns:
        # aim h(ph) at 60% of the failure threshold
        ph = solve_ph(0.6 * 2.0 * cost / (P1 * n))
        crv = cr_mix(n, ph, P1, P1 / 2.0)
        best_q = max((cr_mix(n, ph, P1, q / 1000.0)
                      for q in range(int(1000 * max(0.0, 2 * P1 - 1)),
                                     int(1000 * P1) + 1)))
        marg = (1.0 - P1) * ph
        Hmu = F_empty_mix(n, ph, P1)
        row = {"n": n, "ph": ph, "P1": P1, "marg": marg, "H_mu": Hmu,
               "CR_halfmix": crv, "CR_best_q": best_q,
               "CR_q_half": cr_mix(n, ph, P1, 0.5)}
        fam.append(row)
        log(f"    n={n:3d}: ph={ph:.6f} marg={marg:.4f} H(mu)={Hmu:.4f} "
            f"CR_hm={crv:+.5f}  CR(q=1/2)={row['CR_q_half']:+.5f} "
            f"CR(q*)={best_q:+.5f}")
    OUT["A2_family"] = fam
    assert all(r["CR_halfmix"] < 0 for r in fam), "sliver construction failed"
    assert all(r["marg"] < MARG for r in fam)
    assert all(r["CR_q_half"] > 0 for r in fam), "q=1/2 repair failed?"

    # A3: map the whole in-regime failure region at a few n (grid scan).
    log()
    log("A3: in-regime grid scan (ph x P1), fraction of genre grid where")
    log("    the stated branch fails vs where sup_q CR <= 0:")
    grid_rows = []
    for n in ([1, 2, 4, 8, 16] if not FAST else [1, 4]):
        tot = fail_hm = fail_all = 0
        worst = None
        for iph in range(500, 1000):
            ph = (iph + 0.5) / 1000.0
            for iP in range(1, 667):
                P1 = iP / 1000.0
                if (1.0 - P1) * ph >= MARG:
                    continue
                tot += 1
                crv = cr_mix(n, ph, P1, P1 / 2.0)
                if crv < 0:
                    fail_hm += 1
                    if worst is None or crv < worst["CR_halfmix"]:
                        worst = {"n": n, "ph": ph, "P1": P1,
                                 "CR_halfmix": crv}
                    # sup over q on the fine grid
                    qlo = max(0.0, 2 * P1 - 1.0)
                    best = max(cr_mix(n, ph, P1, qlo + (P1 - qlo) * j / 400)
                               for j in range(401))
                    if best <= 0:
                        fail_all += 1
        grid_rows.append({"n": n, "grid_points": tot, "halfmix_fail": fail_hm,
                          "allq_fail": fail_all, "worst": worst})
        log(f"    n={n:3d}: {tot} in-regime grid pts; stated branch fails on "
            f"{fail_hm} ({100.0*fail_hm/tot:.2f}%); sup_q CR <= 0 on "
            f"{fail_all}; worst CR_hm = "
            f"{worst['CR_halfmix'] if worst else 0:+.5f} at "
            f"(ph={worst['ph'] if worst else 0:.4f}, "
            f"P1={worst['P1'] if worst else 0:.3f})")
    OUT["A3_grid"] = grid_rows

    # A4: concavity criterion F'(P1) < 0 (sup_q CR > 0 iff, since q <= P1):
    # F'(w) = -[log2 m_w(empty) - E_nu log2 m_w].  Verify sign agreement
    # with the grid at n = 4.
    log()
    n = 4
    agree = tested = 0
    for iph in range(500, 1000, 7):
        ph = (iph + 0.5) / 1000.0
        for iP in range(10, 667, 13):
            P1 = iP / 1000.0
            if (1.0 - P1) * ph >= MARG:
                continue
            beta = (1.0 - ph) ** n
            m0 = P1 + (1.0 - P1) * beta
            # E_nu log2 m_w: nu = Bern(ph)^n; mass under mu of t != 0 is
            # (1-P1) nu(t); nu-avg log2-mass = beta log2 m0 (t=0 part)
            #   + sum_{t!=0} nu(t) log2((1-P1) nu(t))
            s_ne = -n * h(ph) - xlog2x(beta)   # sum_{t!=0} nu log2 nu
            enu = beta * log2(m0) + (1.0 - beta) * log2(1.0 - P1) + s_ne
            fprime = -(log2(m0) - enu)
            qlo = max(0.0, 2 * P1 - 1.0)
            best = max(cr_mix(n, ph, P1, qlo + (P1 - qlo) * j / 400)
                       for j in range(401))
            tested += 1
            if (fprime < -1e-9) == (best > 1e-9):
                agree += 1
    log(f"A4: concavity criterion F'(P1) < 0  <=>  sup_q CR > 0: "
        f"{agree}/{tested} grid points agree (n = 4)")
    OUT["A4_criterion"] = {"n": n, "agree": agree, "tested": tested}

    # A5: Sinkhorn rescue on violating instances (does the tilt branch
    # cover the sliver if the recipe routes it there?)
    log()
    log("A5: sup_lambda Sinkhorn CR on violating instances:")
    resc = []
    for n, ph, P1 in [(2, 0.9935, 0.64), (6, 0.9985, 0.64)]:
        mu = {}
        for s in range(1 << n):
            w1 = bin(s).count("1")
            mu[s] = mu.get(s, 0.0) + (1.0 - P1) * ph ** w1 \
                * (1.0 - ph) ** (n - w1)
        mu[0] = mu.get(0, 0.0) + P1
        mu = {a: w for a, w in mu.items() if w > 1e-13}
        tot = sum(mu.values())
        mu = {a: w / tot for a, w in mu.items()}
        best = None
        for lam in (0.0, 0.05, 0.15, 0.35, 0.8, 1.8, 3.5):
            if lam == 0.0:
                pairs = {(a, b): wa * wb for a, wa in mu.items()
                         for b, wb in mu.items()}
                resid = 0.0
            else:
                pairs, resid = ORA.sinkhorn(mu, lam, tol=1e-12, iters=30000)
            if resid > 1e-8:
                continue
            r = cr_eval(n, pairs)
            r["lambda"] = lam
            if best is None or r["CR"] > best["CR"]:
                best = r
        crhm = cr_mix(n, ph, P1, P1 / 2.0)
        log(f"    (n={n}, ph={ph}, P1={P1}): CR_halfmix={crhm:+.5f}; "
            f"sup_lam Sinkhorn CR = {best['CR']:+.5f} at lam="
            f"{best['lambda']}; CR(q=1/2) = "
            f"{cr_mix(n, ph, P1, 0.5):+.5f}")
        resc.append({"n": n, "ph": ph, "P1": P1, "CR_halfmix": crhm,
                     "sinkhorn_best_CR": best["CR"],
                     "sinkhorn_best_lam": best["lambda"],
                     "CR_q_half": cr_mix(n, ph, P1, 0.5)})
    OUT["A5_rescue"] = resc


def part_A6():
    """Coverage lemma for the generalized empty-mixing branch.

    Claim: for every mu = P1 delta_empty + (1-P1) Bern(ph)^n with
    ph in [1/2, 1), 0 < P1 < 1 and (1-P1) ph < MARG, sup_q CR > 0.

    Proof structure (record has the derivations):
      sup_q CR > 0  <=>  F'(P1) < 0  <=>  G > 0, where
      G = (1-b)[log2 m0 - log2(1-P1)] + n h(ph) + b log2 b,
      b = (1-ph)^n, m0 = P1 + (1-P1) b.
      (F strictly concave along the delta_empty--nu segment; q < P1 is
      always feasible.)
      (1) G is strictly increasing in P1 (both bracket terms are).
      (2) P1 >= 1/2: m0 >= P1 >= 1-P1 makes the bracket >= 0, and
          n h(ph) + b log2 b > 0 for ph in (0,1) -- so G > 0, every n.
      (3) P1 < 1/2 (in-regime forces P1 > 1-2*MARG = 0.23458 and
          ph < 2*MARG, so h(ph) > h(2*MARG) = 0.78591 and
          P1/(1-P1) > 0.30647): bracket >= -log2(1/0.30647) = -1.7062
          and n h(ph) + b log2 b >= 0.78591 n - 0.5307, positive for
          n >= 3.
      (4) n in {1, 2}, P1 < 1/2: by (1) it suffices at the in-regime
          infimum P1*(ph) = 1 - MARG/ph; this sweep checks that 1-D
          boundary on a dense grid (float-level).
    """
    log()
    log("A6: coverage-lemma boundary sweep, n in {1, 2} at "
        "P1 = 1 - MARG/ph + 1e-12:")

    def G(n, ph, P1):
        b = (1.0 - ph) ** n
        m0 = P1 + (1.0 - P1) * b
        return ((1.0 - b) * (log2(m0) - log2(1.0 - P1))
                + n * h(ph) + xlog2x(b))

    res = []
    for n in (1, 2):
        worst = None
        steps = 200001 if not FAST else 2001
        for j in range(steps):
            ph = 0.5 + (2 * MARG - 0.5) * j / (steps - 1)
            P1 = max(1e-12, 1.0 - MARG / ph + 1e-12)
            g = G(n, ph, P1)
            if worst is None or g < worst[0]:
                worst = (g, ph, P1)
        log(f"    n={n}: min G over boundary grid ({steps} pts) = "
            f"{worst[0]:+.6f} at ph={worst[1]:.6f} (P1*={worst[2]:.6f})")
        res.append({"n": n, "min_G": worst[0], "at_ph": worst[1],
                    "at_P1": worst[2], "grid": steps})
        assert worst[0] > 0.0
    # analytic-case constants used in (2)/(3), recorded for the skeptic
    OUT["A6_boundary"] = {
        "rows": res,
        "case3_consts": {"h_2marg": h(2 * MARG),
                         "log2_r_min": log2((1 - 2 * MARG) / (2 * MARG)),
                         "max_neg_xlog2x": math.log2(math.e) / math.e,
                         "n3_margin": 3 * h(2 * MARG)
                         - math.log2(math.e) / math.e
                         - (-log2((1 - 2 * MARG) / (2 * MARG)))}}
    log(f"    case-(3) n>=3 analytic margin: 3*h(2*MARG) - 0.5307 - 1.7062"
        f" = {OUT['A6_boundary']['case3_consts']['n3_margin']:+.4f} > 0")


# ------------------------------------------------------------------- part B
def profile_H_mixture(n, comps):
    """H of sum_k P_k Bern(p_k)^n by weight-profile arithmetic (exact in
    profile space; float).  comps: list of (P_k, p_k)."""
    # log-mass of a single mask at weight w
    H = tot = 0.0
    logmass = []
    for w in range(n + 1):
        terms = []
        for P, p in comps:
            if p <= 0.0:
                if w == 0:
                    terms.append(log2(P))
                continue
            if p >= 1.0:
                if w == n:
                    terms.append(log2(P))
                continue
            terms.append(log2(P) + w * log2(p) + (n - w) * log2(1.0 - p))
        if not terms:
            logmass.append(None)
            continue
        m = max(terms)
        la = m + log2(sum(2.0 ** (t - m) for t in terms))
        logmass.append(la)
        lc = math.lgamma(n + 1) - math.lgamma(w + 1) - math.lgamma(n - w + 1)
        lc /= math.log(2.0)
        mass = 2.0 ** (lc + la)
        tot += mass
        H -= mass * la
    assert abs(tot - 1.0) < 1e-9, tot
    return H


def block_pairs(n, comps):
    """Shared-k block coupling, iid given k: pairs (A,B) -> weight."""
    pairs = {}
    for P, p in comps:
        if p <= 0.0:
            pairs[(0, 0)] = pairs.get((0, 0), 0.0) + P
            continue
        pa = {}
        for a in range(1 << n):
            w1 = bin(a).count("1")
            pa[a] = p ** w1 * (1.0 - p) ** (n - w1)
        for a, wa in pa.items():
            if wa < 1e-16:
                continue
            for b, wb in pa.items():
                w = P * wa * wb
                if w >= 1e-16:
                    pairs[(a, b)] = pairs.get((a, b), 0.0) + w
    return pairs


def part_B():
    log()
    log("=" * 74)
    log("B. Block branch on product mixtures with above-psi components")
    log(f"   (psi = (3-sqrt5)/2 = {PSI:.6f}: per-coordinate iid-given-k")
    log("   gain h((1-p)^2) - h(p) changes sign exactly there)")
    log("=" * 74)

    # B0: the sign flip at psi (sanity of the mechanism statement)
    eps = 1e-9
    gL = h((1.0 - (PSI - 1e-4)) ** 2) - h(PSI - 1e-4)
    gR = h((1.0 - (PSI + 1e-4)) ** 2) - h(PSI + 1e-4)
    gAt = h((1.0 - PSI) ** 2) - h(PSI)
    log(f"B0: gain/coord at psi-1e-4: {gL:+.3e}; at psi: {gAt:+.3e}; "
        f"at psi+1e-4: {gR:+.3e}")
    OUT["B0_psi_flip"] = {"below": gL, "at": gAt, "above": gR}
    assert gL > 0 > gR and abs(gAt) < 1e-7

    # B1: in-regime mixtures, profile-level Gain (>= CR) vs n.
    insts = [("2block", [(0.62, 0.0), (0.38, 0.95)]),
             ("3block", [(0.60, 0.0), (0.30, 0.90), (0.10, 0.99)]),
             ("hi+lo", [(0.50, 0.10), (0.50, 0.60)])]
    rows = []
    log()
    log("B1: Gain_block = H(U-mixture) - H(mu) (profile arithmetic);")
    log("    CR <= Gain, so negative Gain kills CR at that n:")
    for tag, comps in insts:
        marg = sum(P * p for P, p in comps)
        assert marg < MARG, (tag, marg)
        gr = []
        for n in [2, 4, 6, 8, 12, 16, 24, 40] if not FAST else [2, 6]:
            Hm = profile_H_mixture(n, comps)
            Hu = profile_H_mixture(
                n, [(P, 1.0 - (1.0 - p) ** 2) for P, p in comps])
            gr.append({"n": n, "H_mu": Hm, "gain_block": Hu - Hm})
        slope = (gr[-1]["gain_block"] - gr[-2]["gain_block"]) / \
            (gr[-1]["n"] - gr[-2]["n"])
        pred = sum(P * (h(1.0 - (1.0 - p) ** 2) - h(p)) for P, p in comps)
        log(f"    {tag:7s} (marg {marg:.3f}): gain at n=2..: "
            + " ".join(f"{r['gain_block']:+.3f}" for r in gr)
            + f"  slope {slope:+.4f} (per-coord prediction {pred:+.4f})")
        rows.append({"tag": tag, "comps": comps, "marg": marg, "rows": gr,
                     "slope": slope, "slope_pred": pred})
    OUT["B1_gain"] = rows

    # B2: actual full-history CR of the iid-block coupling at n <= 8
    # (is CR even worse than Gain?), plus comonotone check CR = 0.
    log()
    log("B2: full-history CR (cr_eval) of the iid-block coupling:")
    b2 = []
    for tag, comps in insts:
        for n in ([4, 6, 8] if not FAST else [4]):
            r = cr_eval(n, block_pairs(n, comps))
            log(f"    {tag:7s} n={n}: CR={r['CR']:+.4f} gain={r['gain']:+.4f}"
                f" ST={r['second_tax']:+.4f} tax_A={r['tax_A']:+.4f}"
                f" (marg_dev {r['marg_dev']:.1e})")
            b2.append({"tag": tag, "n": n, "CR": r["CR"], "gain": r["gain"],
                       "second_tax": r["second_tax"], "tax_A": r["tax_A"]})
    OUT["B2_full_history"] = b2

    # comonotone (A=B) coupling: CR must be exactly 0 (identity check)
    n = 6
    comps = insts[1][1]
    pairs = {}
    for P, p in comps:
        for a in range(1 << n):
            w1 = bin(a).count("1")
            w = P * (p ** w1 * (1.0 - p) ** (n - w1) if 0.0 < p < 1.0
                     else (1.0 if (p == 0.0) == (a == 0) else 0.0))
            if w > 0:
                pairs[(a, a)] = pairs.get((a, a), 0.0) + w
    r = cr_eval(n, pairs)
    log(f"    comonotone (adaptive branch on all-above-psi blocks), 3block "
        f"n=6: CR = {r['CR']:+.2e} (exactly 0 predicted)")
    OUT["B2_comonotone_CR"] = r["CR"]

    # B2b: the ADAPTIVE block coupling (009 part D's m_k rule): per
    # coordinate, given k: p <= 1/4 -> antithetic to union marginal 2p;
    # 1/4 < p <= 1/2 -> union marginal 1/2; p > 1/2 -> comonotone (A=B).
    # Slope is >= 0 always, but the O(1) mixture cross-terms can push
    # small-n CR negative -- measure the n-threshold on hi+lo.
    log()
    log("B2b: ADAPTIVE block coupling on hi+lo (slope > 0, but small n?):")
    b2b = []
    for n in ([2, 3, 4, 5, 6, 8, 10] if not FAST else [2, 4]):
        pairs = adaptive_block_pairs(n, insts[2][1])
        r = cr_eval(n, pairs)
        log(f"    hi+lo n={n:2d}: CR={r['CR']:+.4f} gain={r['gain']:+.4f} "
            f"ST={r['second_tax']:+.4f} (marg_dev {r['marg_dev']:.1e})")
        b2b.append({"n": n, "CR": r["CR"], "gain": r["gain"],
                    "second_tax": r["second_tax"]})
    OUT["B2b_adaptive_hilo"] = b2b

    # B3: Sinkhorn rescue + empty-mixing rescue on the 3block instance, n=6
    log()
    log("B3: rescues on the 3block instance at n = 6:")
    comps = insts[1][1]
    n = 6
    mu = {}
    for P, p in comps:
        for a in range(1 << n):
            w1 = bin(a).count("1")
            w = P * (p ** w1 * (1.0 - p) ** (n - w1) if 0.0 < p < 1.0
                     else (1.0 if (p == 0.0) == (a == 0) else 0.0))
            mu[a] = mu.get(a, 0.0) + w
    mu = {a: w for a, w in mu.items() if w > 1e-13}
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    best = None
    for lam in (0.0, 0.05, 0.15, 0.35, 0.8, 1.8, 3.5):
        if lam == 0.0:
            pairs = {(a, b): wa * wb for a, wa in mu.items()
                     for b, wb in mu.items()}
            resid = 0.0
        else:
            pairs, resid = ORA.sinkhorn(mu, lam, tol=1e-12, iters=30000)
        if resid > 1e-8:
            continue
        r = cr_eval(n, pairs)
        r["lambda"] = lam
        if best is None or r["CR"] > best["CR"]:
            best = r
    log(f"    sup_lam Sinkhorn CR = {best['CR']:+.5f} at lam="
        f"{best['lambda']}")
    # empty-mixing rescue: A = eps_A s, B = eps_B s, s ~ mu|nonempty-part;
    # P1 = mass of the 0-component = 0.6; q swept.
    P1 = comps[0][0]
    scomps = [(P / (1.0 - P1), p) for P, p in comps[1:]]
    qgrid = [max(0.0, 2 * P1 - 1.0) + j / 400.0 *
             (P1 - max(0.0, 2 * P1 - 1.0)) for j in range(401)]
    crq = []
    for q in qgrid:
        Hq = profile_H_mixture_with_empty(n, scomps, q)
        HP = profile_H_mixture_with_empty(n, scomps, P1)
        crq.append((Hq - HP, q))
    bq = max(crq)
    log(f"    empty-mixing rescue: sup_q CR = {bq[0]:+.5f} at q={bq[1]:.4f}"
        f" (closed form; ST = 0 by the extended lemma)")
    # confirm with cr_eval at the optimal q
    pr = general_empty_mixing_pairs(n, scomps, P1, bq[1])
    r = cr_eval(n, pr)
    log(f"    cr_eval confirmation at q*: CR = {r['CR']:+.5f}, "
        f"ST = {r['second_tax']:+.2e}, marg_dev = {r['marg_dev']:.1e}")
    OUT["B3_rescue"] = {"sinkhorn_best_CR": best["CR"],
                        "sinkhorn_best_lam": best["lambda"],
                        "empty_mixing_CR": bq[0], "empty_mixing_q": bq[1],
                        "empty_mixing_cr_eval": r["CR"],
                        "empty_mixing_ST": r["second_tax"]}


def adaptive_block_pairs(n, comps):
    """Shared-k block coupling with 009 part D's adaptive per-coordinate
    joint: p <= 1/4 -> P(0,0) = 1-2p (union marginal 2p); 1/4 < p <= 1/2
    -> P(0,0) = 1/2; p > 1/2 -> comonotone."""
    pairs = {}
    for P, p in comps:
        if p <= 0.0:
            pairs[(0, 0)] = pairs.get((0, 0), 0.0) + P
            continue
        if p > 0.5:
            cell = {(1, 1): p, (0, 0): 1.0 - p}
        else:
            z = 1.0 - 2.0 * p if p <= 0.25 else 0.5
            # marginals p each, both-zero z: P(1,0)=P(0,1)=1-p-z,
            # P(1,1)=p-(1-p-z)
            off = 1.0 - p - z
            cell = {(0, 0): z, (0, 1): off, (1, 0): off, (1, 1): p - off}
            assert min(cell.values()) > -1e-12
        # product over coordinates of the cell joint
        cur = {(0, 0): P}
        for i in range(n):
            nxt = {}
            for (a, b), w in cur.items():
                for (ai, bi), cw in cell.items():
                    if cw <= 0.0:
                        continue
                    k = (a | (ai << i), b | (bi << i))
                    nxt[k] = nxt.get(k, 0.0) + w * cw
            cur = nxt
        for k, w in cur.items():
            if w > 1e-17:
                pairs[k] = pairs.get(k, 0.0) + w
    return pairs


def profile_H_mixture_with_empty(n, scomps, w):
    """H( w*delta_empty + (1-w) * sum_k S_k Bern(p_k)^n )."""
    H = tot = 0.0
    for wt in range(n + 1):
        terms = []
        for S, p in scomps:
            terms.append(math.log2(S * (1.0 - w))
                         + wt * log2(p) + (n - wt) * log2(1.0 - p))
        m = max(terms)
        la = m + log2(sum(2.0 ** (t - m) for t in terms))
        if wt == 0:
            la = log2(w + 2.0 ** la)
        lc = (math.lgamma(n + 1) - math.lgamma(wt + 1)
              - math.lgamma(n - wt + 1)) / math.log(2.0)
        mass = 2.0 ** (lc + la)
        tot += mass
        H -= mass * la
    assert abs(tot - 1.0) < 1e-9, tot
    return H


def general_empty_mixing_pairs(n, scomps, P1, q):
    """A = eps_A s, B = eps_B s with s ~ mixture scomps (shared), joint
    eps with P(0,0)=q."""
    pairs = {}

    def add(k, w):
        if w > 0:
            pairs[k] = pairs.get(k, 0.0) + w

    for S, p in scomps:
        for s in range(1 << n):
            w1 = bin(s).count("1")
            ps = S * p ** w1 * (1.0 - p) ** (n - w1)
            add((0, 0), ps * q)
            add((0, s), ps * (P1 - q))
            add((s, 0), ps * (P1 - q))
            add((s, s), ps * (1.0 - 2.0 * P1 + q))
    return pairs


# ------------------------------------------------------------------- part C
def part_C():
    log()
    log("=" * 74)
    log("C. CR/H ratios on the 023 floor endpoints (lead 1 pre-work)")
    log("=" * 74)
    path = DATA / "cr_attack.json"
    d = json.loads(path.read_text())
    rows = []
    # seed H values for comparison (was the anneal grinding CR or H?)
    import uc_cr_attack as CA
    seedH = {}
    for name, nn, mu in CA.seeds():
        tot = sum(mu.values())
        seedH[name] = -sum(w / tot * log2(w / tot) for w in mu.values())
    log(f"    {'seed':14s} {'n':>2} {'sup_CR':>10} {'H_end':>7} "
        f"{'CR/H':>8} {'H_seed':>7}")
    for r in d["rows"]:
        base = r["seed"].replace("_n6", "")
        hs = seedH.get(base)
        ratio = r["sup_cr"] / r["H_mu"]
        rows.append({"seed": r["seed"], "n": r["n"], "sup_cr": r["sup_cr"],
                     "H_end": r["H_mu"], "cr_over_h": ratio,
                     "H_seed": hs})
        log(f"    {r['seed']:14s} {r['n']:2d} {r['sup_cr']:+10.4f} "
            f"{r['H_mu']:7.3f} {ratio:+8.4f} "
            f"{hs if hs is not None else float('nan'):7.3f}")
    OUT["C_ratios"] = rows


def main():
    part_A()
    part_A6()
    part_B()
    part_C()
    (DATA / "branch_attack.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/branch_attack.json")


if __name__ == "__main__":
    main()
