#!/usr/bin/env python3
"""Attempt 029 (reviewer pass on 025-028): independent re-implementation.

Fresh-context reviewer code.  Shares NOTHING with the explore modules
under review (no imports from uc_mitax / uc_or_avg / uc_branch_* /
uc_adaptive_* / uc_recipe_v2): every coupling is rebuilt from the verbal
definitions in the records (004 R1 half-mixing, 025 Lemma EM, 009 part D
adaptive rule as quoted in 027, 027's latent-mixing description), and CR
is evaluated by this file's own chain-rule accounting:

  - couplings accumulated in exact Fractions where parameters are
    rational (the 026 instances), floats elsewhere;
  - per-coordinate conditional z computed as P(U_i = 1 | prefix) --
    the complementary probability to both prior evaluators' P(U_i = 0);
  - math.fsum accumulation; entropies in bits directly;
  - own IPF/Sinkhorn for the tilt kernel t^{|A cap B|} (G-gen rows);
  - own profile-sum entropy for the closed forms (no algebraic
    F_empty_mix reuse: mass list built and summed directly).

Reads committed checkpoints only to COMPARE against their claims.
Usage: python uc_reviewer029_reimpl.py
Output: ../data/uc_reviewer029_reimpl_out.txt; exit 1 on any refutation.
"""
from __future__ import annotations

import json
import math
import sys
from fractions import Fraction
from math import comb, fsum, log2
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
MARG = 0.38271
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


# ---------------------------------------------------------------- basics
def h(p):
    p = float(p)
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * log2(p) - (1.0 - p) * log2(1.0 - p)


def H_bits(measure):
    """Entropy (bits) of a dict mass function; accepts Fraction weights."""
    return -fsum(float(w) * log2(float(w)) for w in measure.values()
                 if float(w) > 0.0)


def bern(n, p):
    """Bern(p)^n as mask -> weight; keeps Fraction weights if p is one."""
    one = Fraction(1) if isinstance(p, Fraction) else 1.0
    out = {}
    for a in range(1 << n):
        k = bin(a).count("1")
        w = (p ** k) * ((one - p) ** (n - k))
        if float(w) > 0.0:
            out[a] = w
    return out


# ------------------------------------------------------- my CR evaluator
def cr_mine(n, pairs):
    """Chain-rule CR from the definitions of record (009/023 as quoted in
    025's notation block): CR = sum_i E[h(z_i)] - H(mu),
    z_i = P(U_i | A_<i, B_<i}), U = A|B; Gain = H(U) - H(mu);
    ST = Gain - CR.  z computed as P(U_i = 1 | prefix)."""
    ma, mb, ulaw = {}, {}, {}
    for (A, B), w in pairs.items():
        ma[A] = ma.get(A, 0) + w
        mb[B] = mb.get(B, 0) + w
        u = A | B
        ulaw[u] = ulaw.get(u, 0) + w
    marg_dev = max(abs(float(ma.get(k, 0)) - float(mb.get(k, 0)))
                   for k in set(ma) | set(mb))
    Hmu = H_bits(ma)
    HU = H_bits(ulaw)
    terms = []
    for i in range(n):
        mask = (1 << i) - 1
        cells = {}
        for (A, B), w in pairs.items():
            key = (A & mask, B & mask)
            c = cells.get(key)
            if c is None:
                c = cells[key] = [0, 0]
            c[0] += w
            if (A | B) >> i & 1:
                c[1] += w      # U_i = 1
        terms.extend(float(m) * h(float(o) / float(m))
                     for m, o in cells.values() if float(m) > 0.0)
    CR = fsum(terms) - Hmu
    gain = HU - Hmu
    return {"CR": CR, "gain": gain, "ST": gain - CR, "H_mu": Hmu,
            "marg_dev": marg_dev}


# ------------------------------------------------- builders (from records)
def eps_mixing(n, slaw, P1, q):
    """025 Lemma EM coupling: A = eps_A s, B = eps_B s, s ~ slaw,
    (eps_A, eps_B) with equal marginals P(eps=0) = P1, P(0,0) = q."""
    one = Fraction(1) if isinstance(P1, Fraction) else 1.0
    joint = {(0, 0): q, (0, 1): P1 - q, (1, 0): P1 - q,
             (1, 1): one - 2 * P1 + q}
    pairs = {}
    for (ea, eb), we in joint.items():
        if float(we) <= 0.0:
            continue
        for s, ws in slaw.items():
            # A = eps_A * s with eps in {0,1}, P(eps = 0) = P1 the
            # empty-branch probability: eps = 1 keeps s.
            k = (s if ea else 0, s if eb else 0)
            pairs[k] = pairs.get(k, 0) + we * ws
    return pairs


def block_iid(n, comps):
    """iid-given-k: shared component label, A and B iid Bern(p_k)^n."""
    pairs = {}
    for P, p in comps:
        if float(p) == 0.0:
            pairs[(0, 0)] = pairs.get((0, 0), 0) + P
            continue
        bm = bern(n, p)
        for a, wa in bm.items():
            for b, wb in bm.items():
                w = P * wa * wb
                if float(w) > 0.0:
                    pairs[(a, b)] = pairs.get((a, b), 0) + w
    return pairs


def adaptive_cell_joint(p):
    """009 part D rule as quoted in 025/027: p <= 1/4 -> antithetic,
    union marginal 2p; 1/4 < p <= 1/2 -> union marginal 1/2;
    p > 1/2 -> comonotone."""
    p = float(p)
    if p > 0.5:
        return {(0, 0): 1.0 - p, (1, 1): p}
    m = 2.0 * p if p <= 0.25 else 0.5      # union marginal
    z00 = 1.0 - m
    off = 1.0 - p - z00                    # P(0,1) = P(1,0)
    z11 = p - off
    cell = {(0, 0): z00, (0, 1): off, (1, 0): off, (1, 1): z11}
    assert min(cell.values()) > -1e-12
    return {k: v for k, v in cell.items() if v > 0.0}


def cell_product(n, cell, wt):
    cur = {(0, 0): wt}
    for i in range(n):
        nxt = {}
        for (a, b), w in cur.items():
            for (ai, bi), cw in cell.items():
                k = (a | (ai << i), b | (bi << i))
                nxt[k] = nxt.get(k, 0.0) + w * cw
        cur = nxt
    return cur


def adaptive_block(n, comps):
    pairs = {}
    for P, p in comps:
        for k, w in cell_product(n, adaptive_cell_joint(p), float(P)).items():
            pairs[k] = pairs.get(k, 0.0) + w
    return pairs


def latent_mixing(n, plo, phi, Plo, q):
    """027's description: labels (k_A, k_B) with P(lo) = Plo marginals,
    P(lo,lo) = q free; equal labels -> adaptive cell product; differing
    labels -> the two Bern products independent."""
    pairs = {}

    def acc(d):
        for k, w in d.items():
            if w > 0.0:
                pairs[k] = pairs.get(k, 0.0) + w

    acc(cell_product(n, adaptive_cell_joint(plo), q))
    acc(cell_product(n, adaptive_cell_joint(phi), 1.0 - 2.0 * Plo + q))
    cross = Plo - q
    if cross > 1e-15:
        blo, bhi = bern(n, plo), bern(n, phi)
        for a, wa in blo.items():
            for b, wb in bhi.items():
                w = cross * wa * wb
                if w > 0.0:
                    pairs[(a, b)] = pairs.get((a, b), 0.0) + w
                    pairs[(b, a)] = pairs.get((b, a), 0.0) + w
    return pairs


# ------------------------------------------- closed forms via profile sums
def F_profile(n, w, scomps):
    """H( w*delta_0 + (1-w)*sum_k S_k Bern(p_k)^n ) by direct profile
    summation (floats; own derivation, no algebraic shortcuts)."""
    masses = []
    for t in range(n + 1):
        a = fsum(float(S) * float(p) ** t * (1.0 - float(p)) ** (n - t)
                 for S, p in scomps)
        m = (1.0 - w) * a
        if t == 0:
            masses.append((1, m + w))
        else:
            masses.append((comb(n, t), m))
    tot = fsum(mult * m for mult, m in masses)
    assert abs(tot - 1.0) < 1e-9, tot
    return -fsum(mult * m * log2(m) for mult, m in masses if m > 0.0)


def sup_q_em(n, P1, scomps, steps=400):
    qlo = max(0.0, 2 * P1 - 1.0)
    FP1 = F_profile(n, P1, scomps)
    return max(F_profile(n, qlo + (P1 - qlo) * j / steps, scomps) - FP1
               for j in range(steps + 1))


# --------------------------------------------------------- my IPF Sinkhorn
def tilt_coupling(mu, lam, iters=40000, tol=1e-13):
    """pi(a,b) = r(a) c(b) t^{|a&b|} with both marginals mu; my own IPF."""
    t = 2.0 ** lam
    supp = sorted(mu)
    K = {(a, b): t ** bin(a & b).count("1") for a in supp for b in supp}
    r = {a: mu[a] for a in supp}
    c = {b: 1.0 for b in supp}
    for _ in range(iters):
        # column fit
        for b in supp:
            s = fsum(r[a] * K[(a, b)] for a in supp)
            c[b] = mu[b] / s
        # row fit
        dev = 0.0
        for a in supp:
            s = fsum(c[b] * K[(a, b)] for b in supp)
            dev = max(dev, abs(r[a] * s - mu[a]))
            r[a] = mu[a] / s
        if dev < tol:
            break
    pairs = {(a, b): r[a] * K[(a, b)] * c[b] for a in supp for b in supp}
    resid = max(abs(fsum(pairs[(a, b)] for a in supp) - mu[b])
                for b in supp)
    return pairs, resid


def sup_sinkhorn(n, mu, lams=(0.0, 0.05, 0.15, 0.35, 0.8, 1.8, 3.5)):
    mu = {a: w for a, w in mu.items() if w > 1e-13}
    tot = fsum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    best = None
    for lam in lams:
        if lam == 0.0:
            pairs, resid = ({(a, b): wa * wb for a, wa in mu.items()
                             for b, wb in mu.items()}, 0.0)
        else:
            pairs, resid = tilt_coupling(mu, lam)
        if resid > 1e-8:
            continue
        r = cr_mine(n, pairs)
        if best is None or r["CR"] > best["CR"]:
            best = r
    return best["CR"]


# ====================================================================
def main():
    say("Reviewer 029: independent re-implementation vs 025/026/027/028")
    say("=" * 72)
    ba = json.loads((DATA / "branch_attack.json").read_text())
    bc = json.loads((DATA / "branch_certify.json").read_text())
    ac = json.loads((DATA / "adaptive_cell.json").read_text())
    rv = json.loads((DATA / "recipe_v2_battery.json").read_text())
    cert = {c["name"].strip(): c for c in bc["certified"]}

    # ---------------------------------------------------------------- R1
    say("\nR1. Half-mixing CR at one violating instance per n in {2, 6}")
    say("    (exact-Fraction couplings at 026's rational instances; my")
    say("    own evaluator; vs 025 floats AND 026 certified enclosures):")
    for n, phr, ref_float in [
            (2, Fraction(1987, 2000), ba["A5_rescue"][0]["CR_halfmix"]),
            (6, Fraction(1997, 2000), ba["A5_rescue"][1]["CR_halfmix"])]:
        P1 = Fraction(16, 25)
        pairs = eps_mixing(n, bern(n, phr), P1, P1 / 2)
        r = cr_mine(n, pairs)
        marg = float((1 - P1) * phr)
        c = cert[f"CR_halfmix(n={n:2d}, ph={phr}, P1=16/25)"]
        mid = 0.5 * (c["lo"] + c["hi"])
        check(f"n={n}: CR < 0, in-regime, ST = 0, matches 025 + 026",
              r["CR"] < 0 and marg < MARG and abs(r["ST"]) < 1e-11
              and abs(r["CR"] - ref_float) < 1e-11
              and abs(r["CR"] - mid) < 1e-11 and r["marg_dev"] == 0,
              f"CR={r['CR']:+.12f} certmid={mid:+.12f} ST={r['ST']:+.1e}")
        # repair at q = 1/2 vs 026 C2
        r2 = cr_mine(n, eps_mixing(n, bern(n, phr), P1, Fraction(1, 2)))
        c2 = cert[f"CR(q=1/2) (n={n:2d}) repair"]
        mid2 = 0.5 * (c2["lo"] + c2["hi"])
        check(f"n={n}: q=1/2 repair positive, matches 026",
              r2["CR"] > 0 and abs(r2["CR"] - mid2) < 1e-11,
              f"CR={r2['CR']:+.12f} certmid={mid2:+.12f}")

    # A2-family float instances (025's table values)
    say("\n    025 A2 family spot checks (float instances, my evaluator):")
    for n in (1, 2, 6):
        row = next(r for r in ba["A2_family"] if r["n"] == n)
        pairs = eps_mixing(n, bern(n, row["ph"]), 0.64, 0.32)
        r = cr_mine(n, pairs)
        check(f"A2 n={n}: CR_hm reproduces and < 0",
              r["CR"] < 0 and abs(r["CR"] - row["CR_halfmix"]) < 1e-10,
              f"mine={r['CR']:+.6f} engine={row['CR_halfmix']:+.6f}")
    # n = 16, 100 via my own profile closed form (EM: CR = F(q)-F(P1))
    for n in (16, 100):
        row = next(r for r in ba["A2_family"] if r["n"] == n)
        sc = [(1.0, row["ph"])]
        crv = F_profile(n, 0.32, sc) - F_profile(n, 0.64, sc)
        check(f"A2 n={n}: closed-form CR_hm reproduces and < 0",
              crv < 0 and abs(crv - row["CR_halfmix"]) < 1e-9,
              f"mine={crv:+.6f} engine={row['CR_halfmix']:+.6f}")

    # ---------------------------------------------------------------- R2
    say("\nR2. Block Gain at n = 8, 2block (exact-Fraction mixture, my")
    say("    entropies by direct 2^n enumeration; vs 025 float + 026):")
    n = 8
    comps = [(Fraction(31, 50), Fraction(0)), (Fraction(19, 50),
                                               Fraction(19, 20))]
    mu, uu = {}, {}
    for P, p in comps:
        for a, w in bern(n, p).items():
            mu[a] = mu.get(a, 0) + P * w
        pu = 1 - (1 - p) ** 2
        for a, w in bern(n, pu).items():
            uu[a] = uu.get(a, 0) + P * w
    gain = H_bits(uu) - H_bits(mu)
    ref = next(r["gain_block"] for r in ba["B1_gain"][0]["rows"]
               if r["n"] == 8)
    c = cert["Gain_block(2block, n=8)"]
    check("2block Gain(n=8) < 0, matches 025 float and 026 enclosure",
          gain < 0 and abs(gain - ref) < 1e-10
          and abs(gain - 0.5 * (c["lo"] + c["hi"])) < 1e-10,
          f"mine={gain:+.10f} cert={0.5*(c['lo']+c['hi']):+.10f}")
    # and CR itself (full-history) is even below Gain: ST >= 0 check
    r = cr_mine(n, block_iid(n, comps))
    check("2block full-history CR <= Gain (ST >= 0) and CR < 0",
          r["CR"] < 0 and r["ST"] > -1e-12 and abs(r["gain"] - gain) < 1e-9,
          f"CR={r['CR']:+.5f} gain={r['gain']:+.5f} ST={r['ST']:+.2e}")

    # ---------------------------------------------------------------- R3
    say("\nR3. Empty-mixing rescue, 3block n = 6, q = 1/5 (exact-Fraction")
    say("    coupling from the EM definition; my evaluator):")
    n = 6
    slaw = {}
    for S, p in [(Fraction(3, 4), Fraction(9, 10)),
                 (Fraction(1, 4), Fraction(99, 100))]:
        for a, w in bern(n, p).items():
            slaw[a] = slaw.get(a, 0) + S * w
    pairs = eps_mixing(n, slaw, Fraction(3, 5), Fraction(1, 5))
    r = cr_mine(n, pairs)
    c = cert["CR_empty_mixing(3block, n=6, q=1/5)"]
    mid = 0.5 * (c["lo"] + c["hi"])
    check("rescue CR > 0, ST = 0 (Lemma EM), matches 025 B3 + 026 C4",
          r["CR"] > 0.5 and abs(r["ST"]) < 1e-11
          and abs(r["CR"] - ba["B3_rescue"]["empty_mixing_CR"]) < 1e-10
          and abs(r["CR"] - mid) < 1e-11 and r["marg_dev"] == 0,
          f"CR={r['CR']:+.12f} certmid={mid:+.12f} ST={r['ST']:+.1e}")

    # ---------------------------------------------------------------- R4
    say("\nR4. Lemma EM stress (my builders): ST = 0 on fresh arbitrary")
    say("    (eps-joint, s-law) instances never run by 025/026:")
    worst = 0.0
    mixed = {}
    for a, w in bern(5, 0.87).items():
        mixed[a] = mixed.get(a, 0.0) + 0.6 * w
    for a, w in bern(5, 0.52).items():
        mixed[a] = mixed.get(a, 0.0) + 0.4 * w
    for n, slaw, P1, q in [
            (5, mixed, 0.55, 0.13),
            (4, bern(4, 0.51), 0.71, 0.47),   # q > max(0, 2P1-1) = 0.42
            (3, {0b101: 0.5, 0b010: 0.3, 0b111: 0.2}, 0.44, 0.0)]:
        r = cr_mine(n, eps_mixing(n, slaw, P1, q))
        worst = max(worst, abs(r["ST"]))
    check("ST = 0 on 3 fresh EM instances (incl. non-product s-law)",
          worst < 1e-12, f"max |ST| = {worst:.2e}")

    # ---------------------------------------------------------------- R5
    say("\nR5. EM-coverage boundary sweep, n in {1, 2} (my own G):")

    def G(n, ph, P1):
        b = (1.0 - ph) ** n
        m0 = P1 + (1.0 - P1) * b
        return ((1.0 - b) * (log2(m0) - log2(1.0 - P1)) + n * h(ph)
                + (b * log2(b) if b > 0 else 0.0))

    for n, ref_min in [(1, 0.344848410901304), (2, 0.8657862643611467)]:
        vals = []
        steps = 20001
        for j in range(steps):
            ph = 0.5 + (2 * MARG - 0.5) * j / (steps - 1)
            P1 = max(1e-12, 1.0 - MARG / ph + 1e-12)
            vals.append(G(n, ph, P1))
        check(f"min G > 0 on n={n} boundary (my 20001-pt grid)",
              min(vals) > 0.3 and abs(min(vals) - ref_min) < 1e-3,
              f"min={min(vals):+.6f} engine={ref_min:+.6f}")
    # constants of case (3), computed from scratch
    h2m = h(2 * MARG)
    ratio = log2(2 * MARG / (1 - 2 * MARG))
    xmax = log2(math.e) / math.e
    margin3 = 3 * h2m - xmax - ratio
    check("case-(3) constants: h(2M)=0.78591, ratio=1.7062, "
          "max(-xlog2x)=0.53074, margin=+0.1208",
          abs(h2m - 0.78591) < 5e-6 and abs(ratio - 1.7062) < 5e-5
          and abs(xmax - 0.530738) < 1e-6 and abs(margin3 - 0.1208) < 5e-5,
          f"{h2m:.6f} {ratio:.6f} {xmax:.6f} {margin3:+.6f}")
    # psi threshold: h((1-p)^2) = h(p) at p = (3-sqrt5)/2 -- algebraic:
    # 2p - p^2 = 1 - p  <=>  p^2 - 3p + 1 = 0
    psi = (3 - math.sqrt(5)) / 2
    check("psi root identity 2p-p^2 = 1-p at (3-sqrt5)/2 (exact algebra)",
          abs((2 * psi - psi * psi) - (1 - psi)) < 1e-15
          and abs(h((1 - psi) ** 2) - h(psi)) < 1e-12)

    # ---------------------------------------------------------------- R6
    say("\nR6. 027 spot checks (my adaptive builder + evaluator):")
    for row in ac["tightest"][:4]:
        comps = [(1.0 - row["P_hi"], row["p_lo"]),
                 (row["P_hi"], row["p_hi"])]
        r = cr_mine(row["n"], adaptive_block(row["n"], comps))
        check(f"scan p_lo={row['p_lo']}, p_hi={row['p_hi']}, "
              f"P_hi={row['P_hi']}, n={row['n']}: CR > 0, reproduces",
              r["CR"] > 0 and abs(r["CR"] - row["CR"]) < 1e-10,
              f"mine={r['CR']:+.8f} engine={row['CR']:+.8f}")
    # the missing n = 8 refinement (027's record claims it; the code
    # path never ran) -- run it here on the tightest instance:
    tt = ac["tightest"][0]
    comps = [(1.0 - tt["P_hi"], tt["p_lo"]), (tt["P_hi"], tt["p_hi"])]
    for n in (3, 8):
        r = cr_mine(n, adaptive_block(n, comps))
        check(f"tightest instance at n={n}: adaptive CR > 0 "
              f"(029's own refinement -- absent from 027's checkpoint)",
              r["CR"] > 0, f"CR={r['CR']:+.6f}")
    # scaling probe rows + latent-mixing values
    for plo, n in [(1e-5, 6), (1e-3, 6)]:
        row = next(r for r in ac["scaling_probe"]
                   if r["p_lo"] == plo and r["n"] == n)
        comps = [(0.58, plo), (0.42, 0.85)]
        r = cr_mine(n, adaptive_block(n, comps))
        check(f"probe p_lo={plo:.0e} n={n}: adaptive CR reproduces, "
              f"CR/H tiny",
              abs(r["CR"] - row["CR_adaptive"]) < 1e-9
              and r["CR"] / r["H_mu"] < 0.003,
              f"CR={r['CR']:+.7f} CR/H={r['CR']/r['H_mu']:.1e}")
        lp = latent_mixing(n, plo, 0.85, 0.58, row["q_star"])
        r2 = cr_mine(n, lp)
        check(f"probe p_lo={plo:.0e} n={n}: latent-mixing CR reproduces",
              abs(r2["CR"] - row["CR_latmix"]) < 1e-9
              and r2["marg_dev"] < 1e-12 and r2["CR"] > 0.2,
              f"mine={r2['CR']:+.6f} engine={row['CR_latmix']:+.6f} "
              f"marg_dev={r2['marg_dev']:.1e}")
    # EM limit via my profile form
    for n in (2, 6):
        lim = sup_q_em(n, 0.58, [(1.0, 0.85)], steps=200)
        check(f"EM p_lo=0 limit at n={n} reproduces",
              abs(lim - ac["em_limit"][str(n)]) < 1e-9,
              f"mine={lim:+.6f} engine={ac['em_limit'][str(n)]:+.6f}")

    # ---------------------------------------------------------------- R7
    say("\nR7. 028 battery sample (6 rows, my implementations):")
    bat = {r["instance"]: r for r in rv["rows"]}
    # G-mix rows via my profile closed form / latent build
    mine = sup_q_em(2, 0.64, [(1.0, 0.9935)])
    check("sliver n=2 battery row",
          abs(mine - bat["sliver n=2"]["CR"]) < 1e-9 and mine > 0.068,
          f"mine={mine:+.6f} battery={bat['sliver n=2']['CR']:+.6f}")
    mine = sup_q_em(16, 0.64, [(1.0, 0.99965)])
    check("sliver n=16 battery row",
          abs(mine - bat["sliver n=16"]["CR"]) < 1e-9 and mine > 0.067,
          f"mine={mine:+.6f} battery={bat['sliver n=16']['CR']:+.6f}")
    mine = sup_q_em(4, 0.62, [(1.0, 0.95)])
    check("2block n=4 battery row",
          abs(mine - bat["2block n=4"]["CR"]) < 1e-9,
          f"mine={mine:+.6f} battery={bat['2block n=4']['CR']:+.6f}")
    mine = sup_q_em(8, 0.60, [(0.75, 0.90), (0.25, 0.99)])
    check("3block n=8 battery row",
          abs(mine - bat["3block n=8"]["CR"]) < 1e-9,
          f"mine={mine:+.6f} battery={bat['3block n=8']['CR']:+.6f}")
    # 2block n=4 full-history confirmation at the closed-form q*
    # (matches 028 section 3's transcript-only cr_chain cross-check)
    qs = [max(0.0, 2 * 0.62 - 1.0) + j / 400.0
          * (0.62 - max(0.0, 2 * 0.62 - 1.0)) for j in range(401)]
    Fs = [F_profile(4, q, [(1.0, 0.95)]) for q in qs]
    qstar = qs[Fs.index(max(Fs))]
    slaw4 = bern(4, 0.95)
    r = cr_mine(4, eps_mixing(4, slaw4, 0.62, qstar))
    check("2block n=4 full-history CR at q* matches closed form, ST = 0",
          abs(r["CR"] - bat["2block n=4"]["CR"]) < 1e-9
          and abs(r["ST"]) < 1e-12,
          f"CR={r['CR']:+.6f} ST={r['ST']:+.1e}")
    # G-gen rows via my own IPF Sinkhorn on the committed endpoint mus
    catt = json.loads((DATA / "cr_attack.json").read_text())
    for seed in ("mmabskill", "windowkill"):
        row = next(r for r in catt["rows"] if r["seed"] == seed)
        mu = {int(s, 2): w for s, w in row["mu"].items()}
        mine = sup_sinkhorn(row["n"], mu)
        ref = bat[f"floor:{seed}"]["CR"]
        check(f"floor:{seed} battery row (my own IPF)",
              abs(mine - ref) < 1e-6 and mine > 0.05,
              f"mine={mine:+.6f} battery={ref:+.6f}")

    say("\n" + "=" * 72)
    if FAILS:
        say(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            say(f"  - {f}")
    else:
        say("No numerical refutation: every re-implemented value agrees.")
    (DATA / "uc_reviewer029_reimpl_out.txt").write_text(
        "\n".join(LINES) + "\n")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
