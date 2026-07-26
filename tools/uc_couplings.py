#!/usr/bin/env python3
"""uc_couplings.py -- dependent / family-adaptive couplings (Idea C of
attempts/union-closed/001-entropy-barrier-map.md; evaluated in
attempts/union-closed/003-dependent-couplings.md).

Statement under test (the coupling analogue of the (S_c) ladder of 002):

  (S-coup at p, class C):  for every mu on 2^[n] with H(mu) > 0 and all
  marginals Pr[i in A] < p, there is a coupling pi in C(mu) -- a joint law
  of (A,B) with BOTH marginals equal to mu -- such that  H_pi(A u B) > H(mu).

Licensing lemma (family side; one line): if F is exactly union-closed and
mu = Unif(F), then EVERY coupling of mu with itself keeps A u B in F, so
H_pi(A u B) <= log|F| = H(mu).  Hence (S-coup at p) implies Frankl with
constant p.  Approximate closure licenses NONE of this (a coupling can put
Theta(1) mass on the o(1) fraction of escaping pairs), which is why the
Chase-Lovett barrier and the 002 no-go do not apply at this interface.

Gain functional (per class C containing the diagonal):
  Gain_C(mu) = sup_{pi in C(mu)} H_pi(A u B) - H(mu)   ( >= 0, diag gives 0 ).
mu is a counterexample to (S-coup, C) iff Gain_C(mu) = 0 (sup not strict).

Parts:
  A. per-coordinate overlap-tilt calculus on Bernoulli products (closed
     forms): Plackett coupling z_rho(x,y); a product Bern(p)^n stops
     obstructing the tilt exactly at odds ratio rho*(p) = p(3p-1)/(1-2p)^2;
     the fixed-rho ceiling is p*(rho) = (4rho-1-sqrt(4rho+1))/(2(4rho-3)),
     = psi at rho=1 (Gilmer/AHS endpoint), -> 1/2 as rho -> inf.
  B. the PURE Chase-Lovett slice under the pure overlap tilt
     pi_rho(A,B) prop rho^{|A n B|} on slice x slice (marginals are
     automatically uniform; no Sinkhorn needed): exact gain sweep in rho.
  B2. brute-force validation at n=8 of the exchangeable bookkeeping
     (U uniform on each slice; H(U) formula) against direct enumeration.
  C. Sinkhorn overlap-tilt pi_lam(A,B) = 2^(g|A| + g|B| + lam|A n B|) with
     profile potentials g fitted so both marginals equal mu, evaluated on
     (i) a product (engine validation against part A closed form),
     (ii) Sawin Prop-6 geometric-mixture gadgets (the exact family that
     killed idea B), (iii) the smoothed-slice family of 002 part D.
  D. block-adaptive couplings (conditionally-iid in the Sawin-gadget case),
     exact large-n evaluation on 002's certificate instances.
  E. adversary hunt for the *fixed global tilt*: minimize over 2-block
     product mixtures (marginal-capped) the best tilt gain, in the
     asymptotic block-assortative model; cross-checked against part C.
  F. controls: atom-level Sinkhorn engine on genuinely union-closed
     families (gain must be <= 0) and on a mini slice+top CL family
     (gain must be > 0), plus engine cross-validation.

Usage: python3 uc_couplings.py [--part A|B|C|D|E|F|all] [--fast]
Deterministic.  Parts A,B,D,E are stdlib-only; parts C and F use numpy.
Results are checkpointed to ../attempts/union-closed/data/003_part*.json.
"""

import argparse
import json
import math
import os
import sys
import time
from math import lgamma, log, log2, sqrt, inf

LN2 = log(2.0)
PSI = (3.0 - sqrt(5.0)) / 2.0          # 0.381966...
RECORD = 0.38271                        # Liu, arXiv:2306.08824
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..", "attempts", "union-closed", "data")


def h(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * log2(p) - (1.0 - p) * log2(1.0 - p)


def log2C(n, w):
    if w < 0 or w > n:
        return -inf
    return (lgamma(n + 1) - lgamma(w + 1) - lgamma(n - w + 1)) / LN2


def logsumexp2(xs):
    m = max(xs)
    if m == -inf:
        return m
    return m + log2(sum(2.0 ** (x - m) for x in xs))


def save(name, obj):
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, name)
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=float)
    print(f"  [checkpoint: {os.path.relpath(path)}]")


# ---------------------------------------------------------------- part A
# Per-coordinate calculus.  A coupling of Bern(1-x) x Bern(1-y) (zero-probs
# x, y) is one number z = Pr[both zero] in [max(0,x+y-1), min(x,y)]; the
# odds ratio rho = z(1-x-y+z)/((x-z)(y-z)) parameterizes it (Plackett).
# Union bit: Pr[U_i = 1] = 1 - z, so per-coordinate union entropy = h(z).

def z_rho(x, y, rho):
    """Both-zero probability of the Plackett coupling with odds ratio rho."""
    if x <= 0 or y <= 0:
        return 0.0
    if x >= 1:
        return y
    if y >= 1:
        return x
    if abs(rho - 1.0) < 1e-12:
        return x * y
    a = rho - 1.0
    S = 1.0 + a * (x + y)
    disc = S * S - 4.0 * a * rho * x * y
    z = (S - sqrt(disc)) / (2.0 * a)
    return min(max(z, max(0.0, x + y - 1.0)), min(x, y))


def rho_star(p):
    """Odds ratio at which the product Bern(p)^n stops obstructing:
    z_rho(x,x) = 1-x at x = 1-p, i.e. the tilted union marginal hits 1-p
    (mirror of p) so h(union marg) = h(p).  Closed form p(3p-1)/(1-2p)^2."""
    return p * (3.0 * p - 1.0) / (1.0 - 2.0 * p) ** 2


def p_star(rho):
    """Ceiling of the FIXED odds-ratio rho tilt over Bernoulli products:
    the root of rho*(p) = rho in (psi, 1/2).  Closed form; = psi at rho=1."""
    if abs(rho - 1.0) < 1e-15:
        return PSI
    return (4.0 * rho - 1.0 - sqrt(4.0 * rho + 1.0)) / (2.0 * (4.0 * rho - 3.0))


def part_A():
    print("=" * 78)
    print("A. Per-coordinate overlap-tilt calculus on Bernoulli products")
    print("=" * 78)
    print("Coupling of Bern(p) with itself, odds ratio rho (rho>1 = overlap-")
    print("biased).  Union marginal m(rho) = 1 - z_rho(1-p,1-p); product")
    print("obstructs the tilt iff h(m) < h(p) for ALL rho, i.e. never once")
    print("rho can exceed rho*(p) = p(3p-1)/(1-2p)^2.\n")

    # verify z_rho: odds ratio round-trip + range, on a grid
    err = 0.0
    for x in [0.05, 0.3, 0.5, 0.6177, 0.8, 0.95]:
        for y in [0.1, 0.4, 0.6177, 0.9]:
            for rho in [0.05, 0.5, 1.0, 1.03, 2.0, 8.0, 100.0]:
                z = z_rho(x, y, rho)
                lo, hi = max(0.0, x + y - 1.0), min(x, y)
                assert lo - 1e-12 <= z <= hi + 1e-12
                if lo + 1e-9 < z < hi - 1e-9 and abs(rho - 1) > 1e-9:
                    r2 = z * (1 - x - y + z) / ((x - z) * (y - z))
                    err = max(err, abs(r2 - rho) / rho)
    print(f"z_rho validated: max odds-ratio round-trip rel. error = {err:.2e}")

    # verify the two closed forms
    e1 = e2 = 0.0
    for p in [0.383, 0.39, 0.40, 0.42, 0.45, 0.48, 0.499]:
        x = 1.0 - p
        z = z_rho(x, x, rho_star(p))
        e1 = max(e1, abs(z - (1.0 - x)))
        e2 = max(e2, abs(rho_star(p_star(rho_star(p))) - rho_star(p)))
    print(f"rho*(p) closed form:  max |z_rho*(x,x)-(1-x)| = {e1:.2e}")
    print(f"p*(rho) closed form:  max inversion error      = {e2:.2e}")
    print(f"continuity: p*(1) = {p_star(1.0):.6f} = psi;  "
          f"p*(1e6) = {p_star(1e6):.6f} -> 1/2\n")

    rows = []
    print(f"{'rho':>8}  {'p*(rho)':>9}   (products with p < p* never "
          f"obstruct tilt rho)")
    for rho in [1.0, 1.01, 1.0303, 1.1, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0, 100.0]:
        ps = p_star(rho)
        rows.append({"rho": rho, "p_star": ps})
        print(f"{rho:8.4f}  {ps:9.6f}")
    rho_rec = rho_star(RECORD)
    print(f"\nOdds ratio needed to clear Liu's record {RECORD}: "
          f"rho*({RECORD}) = {rho_rec:.6f}")
    print("(the tilt analogue of 002's c_needed = 0.007: a 3% odds-ratio")
    print(" tilt already moves the PRODUCT ceiling past the record --")
    print(" products are cheap; the mixture/slice adversaries decide.)")

    # adaptive-rho claim: for every p in (psi, 1/2) there is rho with
    # strictly positive per-coordinate gain; the maximizing m is 1/2.
    worst = inf
    for i in range(1, 400):
        p = PSI + (0.5 - PSI) * i / 400.0
        x = 1.0 - p
        # best achievable m = 1/2 (since x >= 1/2 >= x^2 here): gain 1-h(p)
        g = 1.0 - h(p)
        worst = min(worst, g)
    print(f"\nAdaptive tilt, per-coordinate: max_m h(m) - h(p) = 1 - h(p) > 0")
    print(f"for all p in (psi, 1/2); min over grid = {worst:.6f} (at p->1/2).")
    print("==> the Bernoulli-product obstruction of the adaptive overlap-")
    print("    tilt class sits at exactly 1/2 (vs psi for iid, c=0..1 KL).")
    save("003_partA.json",
         {"p_star_table": rows, "rho_needed_record": rho_rec,
          "z_rho_roundtrip_err": err, "closed_form_errs": [e1, e2]})
    return rho_rec


# ---------------------------------------------------------------- part B
# Pure slice, pure tilt: pi_rho(A,B) prop rho^{|A n B|} on slice x slice.
# Sum_B rho^{|A n B|} is the same for every slice set A, so both marginals
# are automatically Unif(slice): this IS a legal coupling for every rho.
# j = |A n B| has  P(j) prop C(w,j) C(n-w, w-j) rho^j,  and U = A u B is,
# given j, uniform on the slice of weight 2w - j (exchangeability).

def slice_tilt_eval(n, w0, lam2):
    """Exact H(U) - H(A) for the tilted slice coupling, lam2 = log2(rho)."""
    js = range(max(0, 2 * w0 - n), w0 + 1)
    lw = [log2C(w0, j) + log2C(n - w0, w0 - j) + lam2 * j for j in js]
    Z = logsumexp2(lw)
    HU = 0.0
    for j, l in zip(js, lw):
        pj = 2.0 ** (l - Z)
        if pj > 0.0:
            HU += pj * (log2C(n, 2 * w0 - j) - (l - Z))
    return HU - log2C(n, w0)


def part_B(fast):
    print("\n" + "=" * 78)
    print("B. Pure Chase-Lovett slice under the pure overlap tilt (exact)")
    print("=" * 78)
    print("mu = Unif(slice w0 = round(p n)).  This family kills the iid")
    print("functional for p > psi (h(2p-p^2) < h(p)) and has D = inf, so 002")
    print("could only reach it via smoothing.  The tilt coupling reaches it")
    print("directly.  Gain(lam) = H_pi(U) - log2 C(n,w0):\n")
    out = []
    cases = [(1000, 0.3823), (1000, 0.39), (10000, 0.3823), (10000, 0.42)]
    for n, p in cases:
        w0 = round(p * n)
        best, arg = -inf, None
        sweep = {}
        for lam2 in [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0,
                     1.25, 1.5, 2.0, 3.0]:
            g = slice_tilt_eval(n, w0, lam2)
            sweep[lam2] = g
            if g > best:
                best, arg = g, lam2
        pred = n * (1.0 - h(w0 / n))
        print(f"n={n:6d} p={w0/n:.4f}: gain(iid=lam 0) = {sweep[0.0]:+9.2f}"
              f"   max gain = {best:+9.2f} at lam2={arg:g}"
              f"   [n(1-h(p)) = {pred:.1f}]")
        out.append({"n": n, "p": w0 / n, "sweep": sweep, "max_gain": best,
                    "argmax_lam2": arg, "n_times_1mh": pred})
    print("\nReading: the iid row is the CL barrier (negative gain above psi);")
    print("a moderate overlap tilt makes the union sweep the deleted middle")
    print("band -- gain ~ n(1-h(p)) = Theta(n).  The tilt detects the slice's")
    print("non-closure by PROBABILITY, not log-likelihood.")
    save("003_partB.json", out)


def part_B2():
    print("\n" + "=" * 78)
    print("B2. Brute-force validation of the exchangeable bookkeeping (n=8)")
    print("=" * 78)
    n, w0, lam2 = 8, 3, 0.7
    rho = 2.0 ** lam2
    atoms = [m for m in range(1 << n) if bin(m).count("1") == w0]
    # tilted coupling, direct enumeration
    W = {}
    Zt = 0.0
    for a in atoms:
        for b in atoms:
            wgt = rho ** bin(a & b).count("1")
            Zt += wgt
            W[(a, b)] = wgt
    # marginal uniformity
    margs = []
    for a in atoms:
        margs.append(sum(W[(a, b)] for b in atoms) / Zt)
    mspread = (max(margs) - min(margs)) / max(margs)
    # U-law
    ulaw = {}
    for (a, b), wgt in W.items():
        u = a | b
        ulaw[u] = ulaw.get(u, 0.0) + wgt / Zt
    # within-slice uniformity of U
    byw = {}
    for u, pu in ulaw.items():
        byw.setdefault(bin(u).count("1"), []).append(pu)
    uspread = max((max(v) - min(v)) / max(v) for v in byw.values())
    HU_direct = -sum(p * log2(p) for p in ulaw.values())
    gain_direct = HU_direct - log2C(n, w0)
    gain_formula = slice_tilt_eval(n, w0, lam2)
    print(f"marginal-of-A spread over slice atoms : {mspread:.2e} (exact 0)")
    print(f"U within-slice nonuniformity          : {uspread:.2e} (exact 0)")
    print(f"gain direct enumeration               : {gain_direct:+.12f}")
    print(f"gain via part-B profile formula       : {gain_formula:+.12f}")
    ok = (mspread < 1e-12 and uspread < 1e-12
          and abs(gain_direct - gain_formula) < 1e-10)
    print("VALIDATED" if ok else "MISMATCH -- do not trust part B")
    save("003_partB2.json",
         {"marg_spread": mspread, "u_spread": uspread,
          "gain_direct": gain_direct, "gain_formula": gain_formula,
          "validated": ok})


# ---------------------------------------------------------------- part C
# Exchangeable Sinkhorn engine.  pi(A,B) = 2^(g_|A| + g_|B| + lam2 |AnB|).
# For an exchangeable mu (atom log-prob la_w at weight w) the marginal
# condition is g_w + M_w(g) = la_w with
#   M_w(g) = log2 sum_{j,m} C(w,j) C(n-w,m) 2^(g_{j+m} + lam2 j),
# (b = j+m: j common elements inside A, m elements of B outside A).
# The converged pi is a legal coupling with both marginals exactly mu.
# U = A u B is exchangeable with |U| = |A| + m, so the U-profile needs
# only (a, m) after a logsumexp over j.

def np_mod():
    import numpy as np
    return np


def lc_matrix(np, n):
    lg = np.array([lgamma(i + 1) for i in range(n + 1)])
    i = np.arange(n + 1)
    LC = (lg[:, None] - lg[None, :] - lg[np.maximum(i[:, None] - i[None, :],
                                                    0)]) / LN2
    LC[i[:, None] < i[None, :]] = -np.inf
    return LC


def sinkhorn_profile(np, n, la, lam2, g0=None, tol=1e-9, itmax=4000):
    """Fit g so that pi(A,B)=2^(g_|A|+g_|B|+lam2|AnB|) has profile la.
    Returns (g, residual, iterations)."""
    LC = sinkhorn_profile.LC
    sup = la > -1e300
    g = np.where(sup, la / 2.0, -np.inf) if g0 is None else g0.copy()

    def M_of(g):
        M = np.full(n + 1, -np.inf)
        for w in range(n + 1):
            if not sup[w]:
                continue
            J = np.arange(w + 1)
            Mm = np.arange(n - w + 1)
            mat = (LC[w, J][:, None] + LC[n - w, Mm][None, :]
                   + g[J[:, None] + Mm[None, :]] + lam2 * J[:, None])
            mx = mat.max()
            if mx == -np.inf:
                continue
            M[w] = mx + np.log2(np.exp2(mat - mx).sum())
        return M

    for it in range(itmax):
        M = M_of(g)
        gnew = np.where(sup, la - M, -np.inf)
        resid = np.abs(np.where(sup, g - gnew, 0.0)).max()
        g = np.where(sup, 0.5 * (g + gnew), -np.inf)
        if resid < tol:
            return g, resid, it + 1
    return g, resid, itmax


def profile_HA(np, n, la):
    LC = sinkhorn_profile.LC
    mass = np.exp2(np.where(la > -1e300, la + LC[n, :n + 1], -np.inf))
    return float(-(mass * np.where(mass > 0, la, 0.0)).sum()), mass


def tilt_gain_profile(np, n, la, lam2, g0=None):
    """Gain H_pi(U) - H(A) for the Sinkhorn tilt; also sanity numbers."""
    g, resid, its = sinkhorn_profile(np, n, la, lam2, g0)
    LC = sinkhorn_profile.LC
    sup = la > -1e300
    # accumulate q_u, u = a + m, after logsumexp over j
    lq = np.full(2 * n + 1, -np.inf)
    rows = []
    for a in range(n + 1):
        if not sup[a]:
            continue
        J = np.arange(a + 1)
        Mm = np.arange(n - a + 1)
        mat = (LC[a, J][:, None] + LC[n - a, Mm][None, :]
               + g[J[:, None] + Mm[None, :]] + lam2 * J[:, None])
        mx = mat.max()
        if mx == -np.inf:
            continue
        lrow = mx + np.log2(np.exp2(mat - mx).sum(axis=0))  # over j -> per m
        lrow += LC[n, a] + g[a]
        rows.append((a, lrow))
    for a, lrow in rows:
        u = a + np.arange(lrow.size)
        m1 = np.maximum(lq[u], lrow)
        both = np.exp2(np.where(np.isfinite(lq[u]), lq[u] - m1, -np.inf)) \
            + np.exp2(lrow - m1)
        lq[u] = m1 + np.log2(both)
    lq = lq[:n + 1]  # |U| <= n
    qu = np.exp2(lq)
    tot = qu.sum()
    HU = float((qu * (LC[n, :n + 1] - np.where(qu > 0, lq, 0.0))).sum())
    HA, _ = profile_HA(np, n, la)
    return {"gain": HU - HA, "HU": HU, "HA": HA, "mass": float(tot),
            "resid": float(resid), "iters": its, "g": g}


def sawin_profile(np, n, ubar, theta, tail=1e-28):
    K = 0
    while theta ** (K + 1) > tail:
        K += 1
    raw = [(1 - theta) * theta ** k for k in range(K + 1)]
    Z = sum(raw)
    Pk = [r / Z for r in raw]
    pk = [1.0 - (1.0 - ubar) ** (k + 1) for k in range(K + 1)]
    w = np.arange(n + 1)
    la = np.full(n + 1, -np.inf)
    for P, p in zip(Pk, pk):
        comp = log2(P) + w * log2(p) + (n - w) * log2(1.0 - p)
        la = np.logaddexp2(la, comp)
    marg = sum(P * p for P, p in zip(Pk, pk))
    return la, marg


def smoothed_slice_profile(np, n, p0, t_log2):
    """(1-t) Unif(slice w0) + t Bern(q)^n with q = 1-(1-m)^2, m = fixed
    point of the marginal equation (as in 002 part D)."""
    w0 = round(p0 * n)
    t = 2.0 ** t_log2
    m = w0 / n
    for _ in range(200):
        m = (1 - t) * (w0 / n) + t * (1.0 - (1.0 - m) ** 2)
    q = 1.0 - (1.0 - m) ** 2
    w = np.arange(n + 1)
    la = log2(t) + w * log2(q) + (n - w) * log2(1.0 - q)
    slice_atom = log2(1.0 - t) - log2C(n, w0)
    la[w0] = np.logaddexp2(la[w0], slice_atom)
    return la, m, w0, q


def part_C(fast):
    np = np_mod()
    print("\n" + "=" * 78)
    print("C. Sinkhorn overlap tilt on exchangeable adversaries (exact)")
    print("=" * 78)
    t0 = time.time()
    results = {}

    LAMS = [0.0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.45, 0.65, 1.0, 1.5, 2.5]
    NEG = [-0.05, -0.15, -0.4]

    def sweep(name, n, la, note, extra=None):
        sinkhorn_profile.LC = lc_matrix(np, n)
        HA, mass = profile_HA(np, n, la)
        marg = float((mass * np.arange(n + 1)).sum() / n)
        print(f"\n{name}  (n={n}, H(A)={HA:.2f}, mean marginal={marg:.4f}"
              f"{', ' + note if note else ''})")
        print(f"  {'lam2':>6} {'rho':>8} {'gain':>10} {'resid':>9} "
              f"{'mass':>8} {'iters':>6}")
        rows = []
        best = (-inf, None)
        g0 = None
        for lam2 in LAMS:
            r = tilt_gain_profile(np, n, la, lam2, g0)
            g0 = r["g"]
            rows.append({"lam2": lam2, "gain": r["gain"], "resid": r["resid"],
                         "mass": r["mass"], "iters": r["iters"]})
            if r["gain"] > best[0]:
                best = (r["gain"], lam2)
            print(f"  {lam2:6.2f} {2**lam2:8.3f} {r['gain']:+10.3f} "
                  f"{r['resid']:9.1e} {r['mass']:8.5f} {r['iters']:6d}")
        g0 = None
        for lam2 in NEG:
            r = tilt_gain_profile(np, n, la, lam2, g0)
            g0 = r["g"]
            rows.append({"lam2": lam2, "gain": r["gain"], "resid": r["resid"],
                         "mass": r["mass"], "iters": r["iters"]})
            if r["gain"] > best[0]:
                best = (r["gain"], lam2)
            print(f"  {lam2:6.2f} {2**lam2:8.3f} {r['gain']:+10.3f} "
                  f"{r['resid']:9.1e} {r['mass']:8.5f} {r['iters']:6d}")
        print(f"  ==> gain(iid) = {rows[0]['gain']:+.3f};  max gain = "
              f"{best[0]:+.3f} at lam2 = {best[1]}")
        results[name] = {"n": n, "HA": HA, "mean_marg": marg, "rows": rows,
                         "best_gain": best[0], "best_lam2": best[1]}
        if extra:
            results[name].update(extra)
        save("003_partC.json", results)
        return best

    # (i) engine validation: product Bern(p)^n must match part-A closed form
    n = 60 if fast else 80
    p = 0.40
    w = np.arange(n + 1)
    la = w * log2(p) + (n - w) * log2(1 - p)
    best = sweep(f"validation: Bern({p})^{n}", n, la,
                 "closed form: n(h(1-z_rho)-h(p))")
    x = 1 - p
    errs = []
    for lam2 in [0.1, 0.3, 0.65, 1.0]:
        pred = n * (h(z_rho(x, x, 2.0 ** lam2)) - h(p))
        got = [r["gain"] for r in results[f"validation: Bern({p})^{n}"]["rows"]
               if r["lam2"] == lam2][0]
        errs.append(abs(pred - got))
        print(f"  closed-form check lam2={lam2:g}: predicted {pred:+.6f} "
              f"got {got:+.6f}  (|err| = {abs(pred-got):.1e})")
    results[f"validation: Bern({p})^{n}"]["closedform_err"] = max(errs)

    # (ii) Sawin Prop-6 gadgets
    gadgets = [(200, 0.40, 0.05), (200, 0.42, 0.08), (300, 0.3823, 0.02)]
    if fast:
        gadgets = [(120, 0.40, 0.05), (200, 0.42, 0.08)]
    for n, ubar, theta in gadgets:
        la, marg = sawin_profile(np, n, ubar, theta)
        sweep(f"Sawin gadget n={n} ubar={ubar} theta={theta}", n, la,
              "the idea-B killer", {"ubar": ubar, "theta": theta})

    # (iii) smoothed slice (002 part D shape, moderate n)
    for (n, p0, tl2) in ([(200, 0.3927, -8.0)] if fast
                         else [(200, 0.3927, -8.0), (300, 0.40, -10.0)]):
        la, m, w0, q = smoothed_slice_profile(np, n, p0, tl2)
        sweep(f"smoothed slice n={n} w0={w0} t=2^{tl2:g}", n, la,
              f"max marg={m:.4f}", {"p0": p0, "t_log2": tl2, "q": q})

    # (iv) part-E model cross-check: bimodal mixture whose assortative
    # model is strongly negative at moderate tilt but whose iid endpoint
    # is negative too (a tail-killer of the fixed tilt).
    if not fast:
        n, p1, p2, P2 = 300, 0.0003, 0.862, 0.5684
        w = np.arange(n + 1)
        la = np.logaddexp2(
            log2(1 - P2) + w * log2(p1) + (n - w) * log2(1 - p1),
            log2(P2) + w * log2(p2) + (n - w) * log2(1 - p2))
        best = sweep(f"part-E model-killer n={n} bimodal", n, la,
                     f"marg={(1-P2)*p1+P2*p2:.3f}",
                     {"p_low": p1, "p_high": p2, "P_high": P2})
        blocks = [(1 - P2, p1), (P2, p2)]
        print("  model overlay (per-coord x n):")
        for lam2 in [0.0, 0.1, 0.3, 0.65, 1.0, 2.5]:
            mi = net_iid(blocks) * n
            ma = net_gain(blocks, 2.0 ** lam2) * n
            print(f"    lam2={lam2:4g}: iid-model {mi:+8.2f}  "
                  f"assort-model {ma:+8.2f}")

    print(f"\n[part C total time: {time.time()-t0:.1f}s]")
    return results


# ---------------------------------------------------------------- part D
# Block-adaptive couplings, exact large-n evaluation (stdlib).
# Sawin gadget: share the latent k (block-diagonal), then couple the two
# conditionally-iid Bern(p_k)^n copies per coordinate with Pr[both 1]=r_k:
#   m_k := union marginal = 2 p_k - r_k;  choose m_k = 1/2 if p_k <= 1/2
#   (legal and conditionally iid since then r_k = 2 p_k - 1/2 >= p_k^2
#   whenever q_k = 2p_k - p_k^2 >= 1/2, true for p_k >= 1 - 1/sqrt2 = .293),
#   else m_k = p_k (diagonal).  U | k ~ Bern(m_k)^n exactly.
# This coupling is a member of Liu's conditionally-iid class C3(mu).

def mixture_entropy(n, comps):
    """comps: list of (log2 weight, ('prod', p) or ('slice', w0)).
    Exact entropy of the exchangeable mixture, via weight profiles."""
    H = 0.0
    tot = 0.0
    for wgt in range(n + 1):
        lc = log2C(n, wgt)
        ls = []
        for lw, kind in comps:
            if kind[0] == "prod":
                p = kind[1]
                if 0.0 < p < 1.0:
                    ls.append(lw + wgt * log2(p) + (n - wgt) * log2(1 - p))
                elif (p == 0.0 and wgt == 0) or (p == 1.0 and wgt == n):
                    ls.append(lw)
            else:
                if wgt == kind[1]:
                    ls.append(lw - lc)
        if not ls:
            continue
        la = logsumexp2(ls)
        mass = 2.0 ** (lc + la)
        tot += mass
        H -= mass * la
    assert abs(tot - 1.0) < 1e-6, tot
    return H


def part_D():
    print("\n" + "=" * 78)
    print("D. Block-adaptive (conditionally-iid) couplings, exact large n")
    print("=" * 78)
    print("Exactly the instances whose c* certificates killed idea B in 002.")
    print("Coupling: share Sawin's latent k; per coordinate, common-")
    print("randomness mixture with union marginal m_k = 1/2 (p_k<1/2) or")
    print("diagonal (p_k>1/2).  Member of Liu's conditionally-iid class.\n")
    out = []
    print(f"{'n':>6} {'ubar':>7} {'theta':>6} {'max marg':>9} "
          f"{'gain_iid':>10} {'D':>7} {'gain_coupled':>13}")
    for n, ubar, theta in [(2000, 0.390, 0.05), (20000, 0.386, 0.02),
                           (60000, 0.3823, 0.001)]:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from uc_weighted_kl import sawin_eval
        r = sawin_eval(n, ubar, theta, [1.0])
        K = r["K"]
        raw = [(1 - theta) * theta ** k for k in range(K + 1)]
        Z = sum(raw)
        Pk = [x / Z for x in raw]
        pk = [1.0 - (1.0 - ubar) ** (k + 1) for k in range(K + 1)]
        HA = mixture_entropy(n, [(log2(P), ("prod", p))
                                 for P, p in zip(Pk, pk)])
        assert abs(HA - r["HA"]) < 1e-6 * max(1.0, abs(r["HA"])), \
            (HA, r["HA"])
        mk = [0.5 if p <= 0.5 else p for p in pk]
        rk = [2 * p - m for p, m in zip(pk, mk)]
        assert all(x >= p * p - 1e-12 for x, p in zip(rk, pk)), \
            "not conditionally iid"
        HU = mixture_entropy(n, [(log2(P), ("prod", m))
                                 for P, m in zip(Pk, mk)])
        gain = HU - HA
        print(f"{n:6d} {ubar:7.4f} {theta:6.3f} {r['marg']:9.6f} "
              f"{r['HU']-r['HA']:+10.2f} {r['D']:7.3f} {gain:+13.2f}")
        out.append({"n": n, "ubar": ubar, "theta": theta, "marg": r["marg"],
                    "gain_iid": r["HU"] - r["HA"], "D": r["D"],
                    "gain_coupled": gain, "m_k": mk[:4],
                    "cond_iid": True})
    print("\nSmoothed slice (002 part D instance), general dependent block")
    print("coupling: share the component; slice-slice pairs coupled at fixed")
    print("intersection s = 2 w0 - n/2 (U ~ Unif(slice n/2), exact);")
    print("Bernoulli-Bernoulli pairs diagonal.\n")
    print(f"{'n':>7} {'p0':>7} {'eps':>8} {'max marg':>9} {'F(1)_002':>9} "
          f"{'gain_coupled':>13}")
    for n, p0, eps, F1_002 in [(200000, 0.3927, 0.004, -2543.6),
                               (2000000, 0.3835, 0.0004, -3956.2)]:
        w0 = round(p0 * n)
        t = 2.0 ** (-eps * n)
        m = w0 / n
        for _ in range(200):
            m = (1 - t) * (w0 / n) + t * (1.0 - (1.0 - m) ** 2)
        q = 1.0 - (1.0 - m) ** 2
        HA = mixture_entropy(n, [(log2(1 - t), ("slice", w0)),
                                 (log2(t), ("prod", q))])
        HU = mixture_entropy(n, [(log2(1 - t), ("slice", n // 2)),
                                 (log2(t), ("prod", q))])
        gain = HU - HA
        print(f"{n:7d} {p0:7.4f} {eps:8.5f} {m:9.5f} {F1_002:9.1f} "
              f"{gain:+13.1f}")
        out.append({"n": n, "p0": p0, "eps": eps, "marg": m,
                    "gain_coupled": gain, "cond_iid": False,
                    "note": "slice part uses fixed-intersection coupling"})
    print("\nReading: every 002 certificate flips sign under a family-")
    print("adaptive coupling: the entropy deficit Theta(n) that KL could")
    print("only charge O(log 1/delta) for is recovered at probability scale.")
    save("003_partD.json", out)


# ---------------------------------------------------------------- part E
# Adversary hunt for the FIXED-global-tilt recipe, asymptotic model:
# for a mixture of homogeneous products (P_k, p_k), large n and fixed
# rho > 1, the Sinkhorn block coupling becomes assortative (same-k) and
# per-coordinate Plackett(rho) within blocks, so
#   net(rho) = sum_k P_k [ h(1 - z_rho(x_k, x_k)) - h(p_k) ],  x_k = 1-p_k.
# Killer = mixture with max marginal < 1/2 and net(rho) <= 0 for all rho.

def net_gain(blocks, rho):
    return sum(P * (h(z_rho(1 - p, 1 - p, rho)) - h(p)) for P, p in blocks)


def net_iid(blocks):
    """Exact lam = 0 member of the recipe, per coordinate, n -> inf:
    U is the mixture over INDEPENDENT block pairs (a,b) of products with
    union marginal 1 - x_a x_b."""
    tot = 0.0
    for Pa, pa in blocks:
        for Pb, pb in blocks:
            tot += Pa * Pb * h(1.0 - (1 - pa) * (1 - pb))
    return tot - sum(P * h(p) for P, p in blocks)


def tail_coeff(blocks):
    """Large-rho expansion: m_k - p_k ~ sqrt(p_k(1-p_k)/rho), so
    net(rho) ~ T / sqrt(rho) with T = sum_k P_k h'(p_k) sqrt(p_k(1-p_k)).
    T < 0 certifies that no LARGE tilt ever recovers a strict gain."""
    tot = 0.0
    for P, p in blocks:
        if 0.0 < p < 1.0:
            tot += P * log2((1 - p) / p) * sqrt(p * (1 - p))
    return tot


def net_anti0(blocks):
    """lam -> 0^- large-n limit for a 2-block mixture: the block coupling
    becomes extreme ANTI-assortative (cross mass c = min(P_l, P_h)) while
    per-coordinate pairs stay ~iid.  Cross unions are Bern(1 - x_a x_b)
    products.  (This is the channel that separates the bimodal instances,
    confirmed by the engine: part C, bimodal, lam2 = -0.05.)"""
    (Pl, pl), (Ph, ph) = blocks
    c = min(Pl, Ph)
    a, b = Pl - c, Ph - c
    xl, xh = 1 - pl, 1 - ph
    ent = (a * h(1 - xl * xl) + b * h(1 - xh * xh)
           + 2 * c * h(1 - xl * xh))
    return ent - (Pl * h(pl) + Ph * h(ph))


def J_obj(blocks, rhos):
    """Model gain of the single-lambda tilt recipe over three exactly-
    modeled regimes of the sweep: lam = 0 (iid, block-independent),
    lam > 0 fixed (block-assortative + Plackett(rho) within), and
    lam -> 0^- (block-anti-assortative + iid within).  A model-killer must
    defeat all three.  (Other regimes -- lam ~ gamma/n schedules, lam < 0
    fixed -- are not modeled; the true sup can only be larger, so killer
    status here is a necessary screen, not sufficient.)"""
    a = net_iid(blocks)
    b = max(net_gain(blocks, r) for r in rhos)
    c = net_anti0(blocks)
    return max(a, b, c), a, b, c


def hunt(pcap, rhos):
    worst = (inf, None)
    p1s = [0.001, 0.01] + [0.02 * i for i in range(1, 25)] \
        + [0.4 + 0.01 * i for i in range(10)] \
        + [0.49, 0.495, 0.498, 0.499, 0.4995]
    p2s = [0.51, 0.55] + [0.6 + 0.02 * i for i in range(20)]
    for p1 in p1s:
        if p1 >= min(pcap, 0.5):
            continue
        for p2 in p2s:
            if p2 >= 1.0:
                continue
            Pmax = (pcap - p1) / (p2 - p1)
            if Pmax <= 0:
                continue
            for frac in [0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.95, 1.0]:
                P2 = min(Pmax * frac, 0.999)
                blocks = [(1 - P2, p1), (P2, p2)]
                best, a, b, c = J_obj(blocks, rhos)
                if best < worst[0]:
                    worst = (best, (p1, p2, P2, a, b, c))
    # local refinement
    if worst[1]:
        p1, p2, P2 = worst[1][:3]
        step = 0.004
        for _ in range(8):
            improved = False
            for d1 in (-step, 0, step):
                for d2 in (-step, 0, step):
                    for dP in (-step / 2, 0, step / 2):
                        q1, q2, Q2 = p1 + d1, p2 + d2, P2 + dP
                        if not (0 < q1 < 0.5 < q2 < 1 and 0 < Q2 < 1):
                            continue
                        if (1 - Q2) * q1 + Q2 * q2 > pcap + 1e-12:
                            continue
                        blocks = [(1 - Q2, q1), (Q2, q2)]
                        best, a, b, c = J_obj(blocks, rhos)
                        if best < worst[0]:
                            worst = (best, (q1, q2, Q2, a, b, c))
                            p1, p2, P2 = q1, q2, Q2
                            improved = True
            if not improved:
                step /= 2
    return worst


def part_E():
    print("\n" + "=" * 78)
    print("E. Adversary hunt: 2-block product mixtures vs the single-lambda")
    print("   overlap-tilt recipe (asymptotic model)")
    print("=" * 78)
    print("Model of the recipe's sweep over lam (per coordinate, n -> inf):")
    print("  lam = 0:      net_iid  = sum_{a,b} P_a P_b h(1-x_a x_b) - E h(p)")
    print("  lam > 0 fixed: net(rho) = sum_k P_k [h(1-z_rho(x_k,x_k))-h(p_k)]")
    print("               (blocks become assortative since n*lam -> inf)")
    print("Killer candidate = both <= 0 for all rho.  NOTE: intermediate")
    print("lam ~ gamma/n schedules are not modeled; a model-killer is a")
    print("necessary screen only.\n")
    rhos = [2.0 ** (e / 8.0) for e in range(-48, 241)]  # 2^-6 .. 2^30
    out = []
    for pcap in [0.39, 0.42, 0.43, 0.44, 0.45, 0.47, 0.49, 0.499]:
        g, (p1, p2, P2, a, b, c) = hunt(pcap, rhos)
        marg = (1 - P2) * p1 + P2 * p2
        blocks = [(1 - P2, p1), (P2, p2)]
        T = tail_coeff(blocks)
        tag = "> 0 (no killer)" if g > 0 else "<= 0  MODEL-KILLER"
        print(f"cap {pcap:.3f}: hardest  p_low={p1:.4f} p_high={p2:.3f} "
              f"P_high={P2:.4f} (marg={marg:.4f}):\n"
              f"   J = {g:+.7f} b/coord (iid {a:+.5f}, tilt {b:+.7f}, "
              f"anti {c:+.5f}, tail T {T:+.4f})  {tag}")
        out.append({"pcap": pcap, "p_low": p1, "p_high": p2, "P_high": P2,
                    "marg": marg, "J": g, "net_iid": a, "net_tilt": b,
                    "net_anti0": c, "tail_coeff": T,
                    "model_killer": bool(g <= 0)})
    print("\nReading: bimodal mixtures (p_high ~ 0.8) defeat every FIXED")
    print("tilt rho > 1 (tail certificate T < 0), but the lam = 0 endpoint")
    print("(iid) separates them again through the cross-block unions.")
    print("Whether any 2-block mixture defeats the full sweep, and at which")
    print("marginal cap killers first appear, is read off the J column.")
    print("Model-killers are NOT counterexamples to (S-coup): the block-")
    print("adaptive coupling of part D separates every product mixture with")
    print("some marginal < 1/2 once n log-dominates the mixing entropy.")
    print("What they kill is the single family-oblivious tilt parameter:")
    print("adaptivity to the mixture structure is then forced.")
    save("003_partE.json", out)
    return out


# ---------------------------------------------------------------- part F
# Atom-level Sinkhorn engine (arbitrary mu on subsets of [n], small n),
# used for union-closed CONTROLS (gain must be <= 0 for every coupling,
# by the licensing lemma) and for cross-validation.

def atom_sinkhorn_gain(np, atoms, la, lam2, tol=1e-11, itmax=6000):
    """atoms: int bitmasks; la: log2 atom probs. Returns gain, resid."""
    A = np.array(atoms, dtype=np.uint64)
    pc = np.bitwise_count(A[:, None] & A[None, :]).astype(np.float64)
    K = lam2 * pc
    g = la / 2.0
    for it in range(itmax):
        Mrow = K + g[None, :]
        mx = Mrow.max(axis=1)
        M = mx + np.log2(np.exp2(Mrow - mx[:, None]).sum(axis=1))
        gnew = la - M
        resid = np.abs(g - gnew).max()
        g = 0.5 * (g + gnew)
        if resid < tol:
            break
    lp = g[:, None] + g[None, :] + K
    mx = lp.max()
    W = np.exp2(lp - mx)
    U = A[:, None] | A[None, :]
    uniq, inv = np.unique(U, return_inverse=True)
    pu = np.bincount(inv.ravel(), weights=W.ravel())
    pu = pu / pu.sum() * (np.exp2(mx) * W.sum() / W.sum())  # normalized below
    pu = pu / pu.sum()
    HU = float(-(pu * np.log2(np.where(pu > 0, pu, 1.0))).sum())
    pa = np.exp2(la)
    HA = float(-(pa * la).sum())
    escaped = float(pu[~np.isin(uniq, A)].sum()) if uniq.size else 0.0
    return {"gain": HU - HA, "resid": float(resid), "iters": it + 1,
            "escape_mass": escaped}


def close_family(gens):
    fam = set(gens)
    changed = True
    while changed:
        changed = False
        cur = list(fam)
        for a in cur:
            for b in cur:
                if a | b not in fam:
                    fam.add(a | b)
                    changed = True
    return sorted(fam)


def part_F(fast):
    np = np_mod()
    print("\n" + "=" * 78)
    print("F. Controls: union-closed families (gain must be <= 0) and a")
    print("   mini CL slice+top family (gain must be > 0); cross-checks")
    print("=" * 78)
    out = {}
    LAMS = [-0.5, 0.0, 0.3, 0.7, 1.2, 2.0, 3.0, 5.0]

    def uniform_la(np, fam):
        return np.full(len(fam), -log2(len(fam)))

    # union-closed controls
    fams = {
        "{0,{1}} n=1": [0b0, 0b1],
        "powerset n=3": list(range(8)),
        "upset w>=2 plus empty, n=6":
            [0] + [m for m in range(64) if bin(m).count("1") >= 2],
    }
    # a generated union-closed family at n=10
    gens = [0b0000101001, 0b0110000110, 0b1001010010, 0b0001100100,
            0b1010001000]
    fams["closure of 5 generators, n=10"] = close_family(gens)

    for name, fam in fams.items():
        la = uniform_la(np, fam)
        freqs = [sum(1 for m in fam if m >> i & 1) / len(fam)
                 for i in range(10)]
        best = (-inf, None)
        for lam2 in LAMS:
            r = atom_sinkhorn_gain(np, fam, la, lam2)
            if r["gain"] > best[0]:
                best = (r["gain"], lam2)
            assert r["escape_mass"] == 0.0, (name, "not union-closed?")
        ok = best[0] <= 1e-9
        print(f"union-closed control: {name:34s} |F|={len(fam):4d} "
              f"maxfreq={max(freqs):.3f}  max gain = {best[0]:+.2e}  "
              f"{'OK (<=0)' if ok else 'VIOLATION -- BUG'}")
        out[name] = {"size": len(fam), "maxfreq": max(freqs),
                     "max_gain": best[0], "ok": bool(ok)}

    # mini Chase-Lovett: slice w0=4 plus top w>=6, n=10
    n = 10
    F1 = [m for m in range(1 << n) if bin(m).count("1") == 4]
    F2 = [m for m in range(1 << n) if bin(m).count("1") >= 6]
    fam = F1 + F2
    la = uniform_la(np, fam)
    freqs = [sum(1 for m in fam if m >> i & 1) / len(fam) for i in range(n)]
    rows = {}
    best = (-inf, None)
    for lam2 in LAMS:
        r = atom_sinkhorn_gain(np, fam, la, lam2)
        rows[lam2] = {"gain": r["gain"], "escape_mass": r["escape_mass"]}
        if r["gain"] > best[0]:
            best = (r["gain"], lam2)
    print(f"CL mini slice+top n=10:                  |F|={len(fam):4d} "
          f"maxfreq={max(freqs):.3f}  max gain = {best[0]:+.3f} "
          f"at lam2={best[1]} (escape mass at opt: "
          f"{rows[best[1]]['escape_mass']:.3f})")
    out["CL mini slice+top n=10"] = {"size": len(fam),
                                     "maxfreq": max(freqs), "rows": rows,
                                     "max_gain": best[0],
                                     "best_lam2": best[1]}

    # cross-validation: pure slice n=10 w=4, atom engine vs part-B formula
    lam2 = 0.7
    r = atom_sinkhorn_gain(np, F1, uniform_la(np, F1), lam2)
    gform = slice_tilt_eval(10, 4, lam2)
    err = abs(r["gain"] - gform)
    print(f"cross-check atom engine vs slice formula (n=10,w=4,lam2=0.7): "
          f"|err| = {err:.2e}")
    out["crosscheck_slice"] = {"atom": r["gain"], "formula": gform,
                               "err": err}

    # cross-validation: small Sawin profile n=10, atom vs exchangeable engine
    ubar, theta = 0.40, 0.10
    n = 10
    la_prof, _ = sawin_profile(np, n, ubar, theta)
    sinkhorn_profile.LC = lc_matrix(np, n)
    rprof = tilt_gain_profile(np, n, la_prof, 0.5)
    atoms = list(range(1 << n))
    # sawin_profile's la is already the PER-ATOM log2-probability at weight w
    la_atoms = np.array([la_prof[bin(m).count("1")] for m in atoms])
    ratom = atom_sinkhorn_gain(np, atoms, la_atoms, 0.5)
    err2 = abs(rprof["gain"] - ratom["gain"])
    print(f"cross-check exchangeable vs atom engine (Sawin n=10, lam2=0.5): "
          f"gain {rprof['gain']:+.6f} vs {ratom['gain']:+.6f} "
          f"|err| = {err2:.2e}")
    out["crosscheck_engines"] = {"profile": rprof["gain"],
                                 "atom": ratom["gain"], "err": err2}
    save("003_partF.json", out)


# ----------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="all")
    ap.add_argument("--fast", action="store_true")
    args = ap.parse_args()
    todo = args.part.upper()

    if todo in ("A", "ALL"):
        part_A()
    if todo in ("B", "ALL"):
        part_B(args.fast)
        part_B2()
    if todo in ("D", "ALL"):
        part_D()
    if todo in ("E", "ALL"):
        part_E()
    if todo in ("F", "ALL"):
        part_F(args.fast)
    if todo in ("C", "ALL"):
        part_C(args.fast)


if __name__ == "__main__":
    main()
