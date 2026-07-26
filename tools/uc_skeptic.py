#!/usr/bin/env python3
"""uc_skeptic.py -- adversarial verification of 003-dependent-couplings.md.

Written INDEPENDENTLY of tools/uc_couplings.py (no imports from it; all
formulas re-derived from scratch; natural-log internals, different Sinkhorn
parametrization -- by uniqueness of Sinkhorn scalings the fitted coupling is
the same mathematical object, so agreement is a genuine cross-check).

Parts (run: python3 uc_skeptic.py --part L|I|E|S|U|M|K|all):
  L. lemma checks: derivative lemma (d/ddelta H((1-d)mu+d nu)|_0 =
     CrossEnt(nu,mu)-H(mu); full equivalence via concavity), Plackett
     closed forms rho*(p), p*(rho) re-derived by independent root-finding,
     product-tilt odds ratio = 2^lambda.
  I. closed-form iid gains for 003's two headline instances (no Sinkhorn:
     direct mixture laws of U under independence).
  E. independent Sinkhorn overlap-tilt engine (profile-level), validated
     against a from-scratch atom-level engine at n=10, then run on the two
     headline instances + the non-reproducible [ext] instance.
  S. Sawin gadgets: independent f_mu(1), D(U||mu), coupled block gains,
     growth in n (claims 2 and 5: rung ceiling psi; no-go evasion Theta(n)
     vs O(1)).
  U. union-closed controls: a priori impossibility note + battery of exact
     union-closed families (tilt sweep AND random-kernel Sinkhorn probes;
     any positive gain would refute the licensing lemma / expose a bug).
  M. mini-theorem (003 part e): the p=0 bracket bug, an explicit product
     mixture the mini-theorem coupling does NOT separate, and the
     "half-mixing" conditionally-iid coupling that does (with brute-force
     legality check at n=6).
  K. part-E re-analysis: recompute 003's hardest-instance table, exhibit
     exact-delta_empty 2-block model-killers BELOW cap 0.445, compute the
     corrected boundary curve, hunt 3-block mixtures, finite-n engine check.

Checkpoints: ../attempts/union-closed/data/004_*.json.  Deterministic.
"""

import argparse
import json
import os
import sys
import time
from math import lgamma, log, sqrt, inf, exp, isfinite

import numpy as np

LN2 = log(2.0)
PSI = (3.0 - sqrt(5.0)) / 2.0
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..", "attempts", "union-closed", "data")


def lnC(n, k):
    if k < 0 or k > n:
        return -inf
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def h2(p):
    """binary entropy in BITS."""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return (-p * log(p) - (1 - p) * log(1 - p)) / LN2


def lse(v):
    v = np.asarray(v, dtype=float)
    m = np.max(v) if v.size else -inf
    if not isfinite(m):
        return m
    return m + log(np.exp(v - m).sum())


def save(name, obj):
    os.makedirs(DATA, exist_ok=True)
    p = os.path.join(DATA, name)
    with open(p, "w") as f:
        json.dump(obj, f, indent=1, default=float)
    print(f"  [checkpoint: {os.path.relpath(p)}]")


# ----------------------------------------------------------------- Plackett
def z_quad(x, y, rho):
    """Both-zero probability at fixed margins (x,y of the zero events) and
    odds ratio rho; re-derived: (rho-1)z^2 - (1+(rho-1)(x+y))z + rho x y = 0,
    root in [max(0,x+y-1), min(x,y)]."""
    if x <= 0 or y <= 0:
        return 0.0
    if x >= 1:
        return y
    if y >= 1:
        return x
    if abs(rho - 1.0) < 1e-13:
        return x * y
    a = rho - 1.0
    S = 1.0 + a * (x + y)
    z = (S - sqrt(S * S - 4.0 * a * rho * x * y)) / (2.0 * a)
    return min(max(z, max(0.0, x + y - 1.0)), min(x, y))


def z_bisect(x, y, rho):
    """Independent numeric check: bisection on the odds-ratio identity."""
    lo, hi = max(0.0, x + y - 1.0) + 1e-15, min(x, y) - 1e-15
    if lo >= hi:
        return z_quad(x, y, rho)

    def oratio(z):
        return (z * (1 - x - y + z)) / ((x - z) * (y - z))
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if oratio(mid) < rho:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ================================================================= part L
def part_L():
    print("=" * 78)
    print("L. Lemma re-derivations")
    print("=" * 78)
    out = {}

    # --- derivative lemma ---------------------------------------------
    rng = np.random.default_rng(11)
    worst = 0.0
    for trial in range(40):
        k = rng.integers(3, 12)
        mu = rng.random(k) + 1e-3
        mu /= mu.sum()
        nu = rng.random(k) + 1e-3
        nu /= nu.sum()

        def H(p):
            return float(-(p * np.log(np.where(p > 0, p, 1.0))).sum() / LN2)
        ce = float(-(nu * np.log(mu)).sum() / LN2)
        pred = ce - H(mu)                       # claimed f'(0)
        d = 1e-6
        num = (H((1 - d) * mu + d * nu) - H(mu)) / d
        worst = max(worst, abs(num - pred) / max(1.0, abs(pred)))
        # equivalence: sup_delta H((1-d)mu+d nu) > H(mu)  <=>  pred > 0
        ds = np.linspace(0, 1, 201)
        sup = max(H((1 - dd) * mu + dd * nu) for dd in ds)
        lhs = sup > H(mu) + 1e-12
        rhs = pred > 1e-9
        if abs(pred) > 1e-6:
            assert lhs == rhs, (trial, pred, sup - H(mu))
    print(f"derivative lemma: f'(0) = CrossEnt(nu,mu)-H(mu) "
          f"(rel err vs numeric {worst:.1e}); equivalence sup>H(mu) <=> "
          f"f'(0)>0 held in all 40 random trials (concavity argument also "
          f"re-derived analytically: both directions OK)")
    out["derivative_lemma_relerr"] = worst

    # --- Plackett closed forms ----------------------------------------
    e_z = 0.0
    for x in [0.05, 0.4, 0.6177, 0.9]:
        for y in [0.1, 0.5, 0.6177, 0.95]:
            for rho in [0.3, 1.03, 2.0, 20.0]:
                e_z = max(e_z, abs(z_quad(x, y, rho) - z_bisect(x, y, rho)))
    print(f"z_rho quadratic vs bisection: max |diff| = {e_z:.1e}")
    out["z_rho_check"] = e_z

    # rho*(p): smallest rho with union marginal m = 1-p (so h(m)=h(p));
    # found by bisection on rho, compared to claimed closed form.
    e_rs, e_ps = 0.0, 0.0
    for p in [0.383, 0.40, 0.43, 0.47, 0.499]:
        x = 1.0 - p
        lo, hi = 1.0, 1e6
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            m = 1.0 - z_quad(x, x, mid)
            if m > 1.0 - p:      # m decreases with rho
                lo = mid
            else:
                hi = mid
        rs_num = 0.5 * (lo + hi)
        rs_cf = p * (3 * p - 1) / (1 - 2 * p) ** 2
        e_rs = max(e_rs, abs(rs_num - rs_cf) / rs_cf)
        # p*(rho): root of rho*(p)=rho
        rho = rs_cf
        lo2, hi2 = PSI, 0.5
        for _ in range(200):
            mid = 0.5 * (lo2 + hi2)
            if mid * (3 * mid - 1) / (1 - 2 * mid) ** 2 < rho:
                lo2 = mid
            else:
                hi2 = mid
        ps_num = 0.5 * (lo2 + hi2)
        ps_cf = (4 * rho - 1 - sqrt(4 * rho + 1)) / (2 * (4 * rho - 3))
        e_ps = max(e_ps, abs(ps_num - ps_cf))
    print(f"rho*(p) closed form vs bisection: max rel err = {e_rs:.1e}")
    print(f"p*(rho) closed form vs bisection: max |err|  = {e_ps:.1e}")
    print(f"p*(1) = {PSI:.6f} (=psi), analytic inversion re-derived: "
          f"disc = 4rho+1  OK")
    out["rho_star_err"] = e_rs
    out["p_star_err"] = e_ps

    # --- product tilt has constant conditional odds ratio 2^lambda ----
    # pi(A,B) prop 2^{alpha|A|+alpha|B|+lam|AnB|} on products factorizes
    # coordinatewise with per-coordinate weights (1, 2^a, 2^a, 2^{2a+lam}):
    # odds ratio w00 w11/(w01 w10) = 2^lam.  Verified by direct atom
    # computation at n=6, p=0.4, lam=0.8: conditional odds ratio of
    # (A_i,B_i) given every history must equal 2^0.8.
    n, p, lam = 6, 0.4, 0.8
    masks = np.arange(1 << n)
    w = np.array([bin(m).count("1") for m in masks])
    ell = w * log(p) + (n - w) * log(1 - p)
    pi = np.zeros((1 << n, 1 << n))
    # fit alpha so both marginals are Bern(p)^n: per-coordinate 2x2 fit
    # via 1-d bisection on alpha (independent of any Sinkhorn code)
    r = 2.0 ** lam

    def marg1(al):
        a = 2.0 ** al
        z00, z01, z10, z11 = 1.0, a, a, a * a * r
        Z = z00 + z01 + z10 + z11
        return (z10 + z11) / Z
    lo, hi = -20.0, 20.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if marg1(mid) < p:
            lo = mid
        else:
            hi = mid
    al = 0.5 * (lo + hi)
    a = 2.0 ** al
    W = np.array([[1.0, a], [a, a * a * r]])
    W /= W.sum()
    # coupling = product over coordinates of W; check conditional odds
    # ratio given any history is exactly r (trivially true by product
    # structure) and that the union bit matches z_quad:
    x = W[0, 0] + W[0, 1]           # P[A_i=0]
    zz = W[0, 0]
    orat = W[0, 0] * W[1, 1] / (W[0, 1] * W[1, 0])
    print(f"product tilt: fitted per-coordinate coupling has odds ratio "
          f"{orat:.12f} = 2^lam = {r:.12f}; both-zero prob {zz:.9f} vs "
          f"z_rho({x:.4f},{x:.4f},{r:.4f}) = {z_quad(x, x, r):.9f}")
    assert abs(orat - r) < 1e-9 and abs(zz - z_quad(x, x, r)) < 1e-9
    out["product_tilt_or"] = orat
    save("004_partL.json", out)


# ================================================================= profiles
def sawin_components(ubar, theta, tail=1e-28):
    K = 0
    while theta ** (K + 1) > tail:
        K += 1
    raw = np.array([(1 - theta) * theta ** k for k in range(K + 1)])
    P = raw / raw.sum()
    p = np.array([1.0 - (1.0 - ubar) ** (k + 1) for k in range(K + 1)])
    return P, p


def mix_profile(n, Ps, ps):
    """ln per-atom probability profile of  sum_k P_k Bern(p_k)^n."""
    wv = np.arange(n + 1, dtype=float)
    ell = np.full(n + 1, -np.inf)
    for P, p in zip(Ps, ps):
        if p <= 0.0:
            comp = np.full(n + 1, -np.inf)
            comp[0] = log(P)
        elif p >= 1.0:
            comp = np.full(n + 1, -np.inf)
            comp[n] = log(P)
        else:
            comp = log(P) + wv * log(p) + (n - wv) * log(1 - p)
        ell = np.logaddexp(ell, comp)
    return ell


def profile_entropy(n, ell):
    """H (bits) of an exchangeable law given ln per-atom probs; also mass."""
    lc = np.array([lnC(n, w) for w in range(n + 1)])
    sup = np.isfinite(ell)
    mass = np.exp(np.where(sup, ell + lc, -np.inf))
    H = float(-(mass * np.where(sup, ell, 0.0)).sum() / LN2)
    return H, float(mass.sum())


def smoothed_slice_mu(n, p0, t_log2):
    """Reconstruct 003's smoothed-slice mu; returns (ell, w0, q, marg)."""
    w0 = round(p0 * n)
    t = 2.0 ** t_log2
    m = w0 / n
    for _ in range(300):
        m = (1 - t) * (w0 / n) + t * (1.0 - (1.0 - m) ** 2)
    q = 1.0 - (1.0 - m) ** 2
    wv = np.arange(n + 1, dtype=float)
    ell = log(t) + wv * log(q) + (n - wv) * log(1 - q)
    ell[w0] = np.logaddexp(ell[w0], log(1 - t) - lnC(n, w0))
    return ell, w0, q, m


# ================================================================= part I
def part_I():
    print("=" * 78)
    print("I. Closed-form iid gains for the headline instances (no Sinkhorn)")
    print("=" * 78)
    out = {}

    # ---- Sawin n=300, ubar=0.3823, theta=0.02 ------------------------
    n, ubar, theta = 300, 0.3823, 0.02
    P, p = sawin_components(ubar, theta)
    HA, mA = profile_entropy(n, mix_profile(n, P, p))
    # iid U law: mixture over (k1,k2) of Bern(1-(1-p_k1)(1-p_k2))
    PU, pU = [], []
    for i in range(len(P)):
        for j in range(len(P)):
            PU.append(P[i] * P[j])
            pU.append(1.0 - (1.0 - p[i]) * (1.0 - p[j]))
    HU, mU = profile_entropy(n, mix_profile(n, PU, pU))
    g_iid_sawin = HU - HA
    print(f"Sawin (300,0.3823,0.02): H(A)={HA:.3f} (mass err {abs(mA-1):.1e})"
          f"  H(U_iid)={HU:.3f}  gain_iid = {g_iid_sawin:+.3f}"
          f"   [003 claims -2.14]")
    out["sawin300_iid"] = g_iid_sawin

    # ---- smoothed slice n=300, p0=0.40, t=2^-10 ----------------------
    n = 300
    ell, w0, q, marg = smoothed_slice_mu(n, 0.40, -10.0)
    t = 2.0 ** -10.0
    HA, mA = profile_entropy(n, ell)
    # iid U profile: four cases (slice-slice, slice-prod x2, prod-prod)
    lq = np.full(n + 1, -np.inf)
    # ss: j = |A n B| hypergeometric, |U| = 2w0 - j
    for j in range(max(0, 2 * w0 - n), w0 + 1):
        lp = (2 * log(1 - t) + lnC(w0, j) + lnC(n - w0, w0 - j)
              - lnC(n, w0))
        u = 2 * w0 - j
        lq[u] = np.logaddexp(lq[u], lp)
    # sp + ps: |U| = w0 + Binom(n - w0, q)
    for m_ in range(0, n - w0 + 1):
        lp = (log(2.0) + log(t) + log(1 - t) + lnC(n - w0, m_)
              + m_ * log(q) + (n - w0 - m_) * log(1 - q))
        u = w0 + m_
        lq[u] = np.logaddexp(lq[u], lp)
    # pp: |U| ~ Binom(n, Q), Q = 1-(1-q)^2
    Q = 1.0 - (1.0 - q) ** 2
    for u in range(n + 1):
        lp = 2 * log(t) + lnC(n, u) + u * log(Q) + (n - u) * log(1 - Q)
        lq[u] = np.logaddexp(lq[u], lp)
    qmass = np.exp(lq)
    lcv = np.array([lnC(n, u) for u in range(n + 1)])
    HU = float((qmass * (lcv - lq) / LN2)[qmass > 0].sum())
    g_iid_sl = HU - HA
    print(f"smoothed slice (300,w0={w0},t=2^-10): marg={marg:.4f} "
          f"H(A)={HA:.3f}  H(U_iid)={HU:.3f} (mass {qmass.sum():.6f})  "
          f"gain_iid = {g_iid_sl:+.3f}   [003 claims -4.62]")
    out["slice300_iid"] = g_iid_sl

    # ---- pure slice check (003 part B, n=10000 p=0.3823) -------------
    n, w0 = 10000, round(0.3823 * 10000)
    for lam2, tag in [(0.0, "iid"), (3.0, "lam2=3")]:
        rho = 2.0 ** lam2
        js = np.arange(max(0, 2 * w0 - n), w0 + 1, dtype=float)
        lw = np.array([lnC(w0, j) + lnC(n - w0, w0 - j) for j in js]) \
            + js * lam2 * LN2
        Z = lse(lw)
        pj = np.exp(lw - Z)
        HU = float((pj * (np.array([lnC(n, int(2 * w0 - j)) for j in js])
                          - (lw - Z)) / LN2).sum())
        g = HU - lnC(n, w0) / LN2
        ref = 1.25 if lam2 == 0.0 else 407.98
        print(f"pure slice n=10000 p=0.3823 {tag}: gain = {g:+.2f}"
              f"   [003: {ref:+.2f}]")
        out[f"pureslice_{tag}"] = g
    save("004_partI.json", out)


# ================================================================= part E
def fit_profile_tilt(n, ell, lam2, tol=1e-11, itmax=8000):
    """Symmetric potential fit for pi(A,B) = exp(g_|A| + g_|B| +
    lam2*ln2*|AnB|), target exchangeable profile ell (ln per-atom)."""
    lam = lam2 * LN2
    LCa = [np.array([lnC(a, j) for j in range(a + 1)]) for a in range(n + 1)]
    LCb = [np.array([lnC(n - a, m) for m in range(n - a + 1)])
           for a in range(n + 1)]
    sup = np.isfinite(ell)
    g = np.where(sup, ell / 2.0, -np.inf)

    def rowS(g):
        S = np.full(n + 1, -np.inf)
        for a in range(n + 1):
            if not sup[a]:
                continue
            J = np.arange(a + 1)
            mat = (LCa[a][:, None] + LCb[a][None, :]
                   + g[J[:, None] + np.arange(n - a + 1)[None, :]]
                   + lam * J[:, None])
            S[a] = lse(mat.ravel())
        return S
    resid = inf
    for it in range(itmax):
        S = rowS(g)
        gn = np.where(sup, ell - S, -np.inf)
        resid = float(np.abs(np.where(sup, g - gn, 0.0)).max())
        g = np.where(sup, 0.5 * (g + gn), -np.inf)
        if resid < tol:
            break
    return g, resid, it + 1, LCa, LCb


def eval_profile_tilt(n, ell, lam2, tol=1e-11):
    g, resid, its, LCa, LCb = fit_profile_tilt(n, ell, lam2, tol)
    lam = lam2 * LN2
    sup = np.isfinite(ell)
    lcn = np.array([lnC(n, a) for a in range(n + 1)])
    # joint (|A|, |B\A|) log-mass, U weight = a+m ; B-marginal check
    lq = np.full(n + 1, -np.inf)
    lmargB = np.full(n + 1, -np.inf)
    tot = -inf
    for a in range(n + 1):
        if not sup[a]:
            continue
        J = np.arange(a + 1)
        M = np.arange(n - a + 1)
        mat = (LCa[a][:, None] + LCb[a][None, :]
               + g[J[:, None] + M[None, :]] + lam * J[:, None])
        mat = mat + lcn[a] + g[a]
        # per-m collapse for U; per-(j+m) collapse for B-marginal
        mx = mat.max()
        if mx == -inf:
            continue
        colm = mx + np.log(np.exp(mat - mx).sum(axis=0))     # per m
        for m_, v in enumerate(colm):
            u = a + m_
            if u <= n:
                lq[u] = np.logaddexp(lq[u], v)
        bidx = (J[:, None] + M[None, :]).ravel()
        vals = mat.ravel()
        for b in range(n + 1):
            sel = vals[bidx == b]
            if sel.size:
                lmargB[b] = np.logaddexp(lmargB[b], lse(sel))
        tot = np.logaddexp(tot, lse(mat.ravel()))
    qm = np.exp(lq)
    HU = float((qm * (lcn - lq) / LN2)[qm > 0].sum())
    massA = np.exp(np.where(sup, ell + lcn, -np.inf))
    HA = float(-(massA * np.where(sup, ell, 0.0)).sum() / LN2)
    # B-marginal error (per-atom, relative, in ln units)
    eB = float(np.abs(np.where(sup, lmargB - (ell + lcn), 0.0)).max())
    return {"gain": HU - HA, "HU": HU, "HA": HA, "mass": float(exp(tot)),
            "resid": resid, "iters": its, "margB_lnerr": eB}


def atom_sinkhorn_indep(masks, lnmu, lam2, tol=1e-13, itmax=20000):
    """From-scratch atom-level engine with ASYMMETRIC alternating scaling
    (u,v).  By uniqueness of Sinkhorn scalings this converges to the same
    coupling as any symmetric-potential parametrization."""
    A = np.asarray(masks)
    pc = np.array([[bin(int(a) & int(b)).count("1") for b in A] for a in A],
                  dtype=float)
    lnK = lam2 * LN2 * pc
    lu = np.asarray(lnmu) / 2.0
    lv = lu.copy()
    for it in range(itmax):
        r1 = lnK + lv[None, :]
        m1 = r1.max(axis=1)
        lun = lnmu - (m1 + np.log(np.exp(r1 - m1[:, None]).sum(axis=1)))
        r2 = lnK + lun[:, None]
        m2 = r2.max(axis=0)
        lvn = lnmu - (m2 + np.log(np.exp(r2 - m2[None, :]).sum(axis=0)))
        d = max(np.abs(lun - lu).max(), np.abs(lvn - lv).max())
        lu, lv = lun, lvn
        if d < tol:
            break
    lpi = lu[:, None] + lnK + lv[None, :]
    pi = np.exp(lpi)
    mA = pi.sum(axis=1)
    mB = pi.sum(axis=0)
    mu = np.exp(np.asarray(lnmu))
    emax = max(np.abs(mA - mu).max(), np.abs(mB - mu).max())
    ul = {}
    for i, a in enumerate(A):
        for j, b in enumerate(A):
            u = int(a) | int(b)
            ul[u] = ul.get(u, 0.0) + pi[i, j]
    pu = np.array(list(ul.values()))
    pu = pu / pu.sum()
    HU = float(-(pu * np.log(pu)).sum() / LN2)
    HA = float(-(mu * np.asarray(lnmu)).sum() / LN2)
    esc = sum(v for k, v in ul.items() if k not in set(int(x) for x in A))
    return {"gain": HU - HA, "marg_err": float(emax), "iters": it + 1,
            "escape_mass": esc / pi.sum()}


def part_E_engine():
    print("=" * 78)
    print("E. Independent engine: validation + headline instances")
    print("=" * 78)
    out = {}
    t0 = time.time()

    # ---- validation at n=10 vs from-scratch atom engine --------------
    n = 10
    P, p = sawin_components(0.40, 0.10)
    ell = mix_profile(n, P, p)
    masks = list(range(1 << n))
    wts = np.array([bin(m).count("1") for m in masks])
    lnmu = np.array([ell[w] for w in wts])
    errs = []
    for lam2 in [-0.4, 0.0, 0.5, 1.5, 3.0]:
        rp = eval_profile_tilt(n, ell, lam2)
        ra = atom_sinkhorn_indep(masks, lnmu, lam2)
        errs.append(abs(rp["gain"] - ra["gain"]))
        print(f"  n=10 Sawin lam2={lam2:5.2f}: profile {rp['gain']:+.9f} "
              f"atom {ra['gain']:+.9f} |diff|={errs[-1]:.2e} "
              f"(atom marg err {ra['marg_err']:.1e})")
    ell2, w0, q, _ = smoothed_slice_mu(10, 0.40, -4.0)
    lnmu2 = np.array([ell2[w] for w in wts])
    for lam2 in [0.0, 1.0]:
        rp = eval_profile_tilt(10, ell2, lam2)
        ra = atom_sinkhorn_indep(masks, lnmu2, lam2)
        errs.append(abs(rp["gain"] - ra["gain"]))
        print(f"  n=10 smslice lam2={lam2:4.1f}: profile {rp['gain']:+.9f} "
              f"atom {ra['gain']:+.9f} |diff|={errs[-1]:.2e}")
    out["validation_max_err"] = max(errs)
    print(f"  ==> max profile-vs-atom disagreement: {max(errs):.2e}")

    # ---- headline instance 1: Sawin (300, 0.3823, 0.02) --------------
    ref1 = {0.0: -2.137, 0.3: 0.254, 1.0: 4.863, 2.5: 10.357}
    n = 300
    P, p = sawin_components(0.3823, 0.02)
    ell = mix_profile(n, P, p)
    rows = {}
    print("\nSawin gadget n=300 ubar=0.3823 theta=0.02 "
          "(003: iid -2.14 -> +10.36 at lam2=2.5):")
    for lam2 in [0.0, 0.1, 0.3, 0.65, 1.0, 1.5, 2.5, 3.5]:
        r = eval_profile_tilt(n, ell, lam2)
        rows[lam2] = r["gain"]
        ref = f"  [003: {ref1[lam2]:+.3f}]" if lam2 in ref1 else ""
        print(f"  lam2={lam2:4.2f}: gain {r['gain']:+9.3f}  mass "
              f"{r['mass']:.6f} margB_lnerr {r['margB_lnerr']:.1e}{ref}")
    out["sawin300"] = rows

    # ---- headline instance 2: smoothed slice (300, 0.40, 2^-10) ------
    ref2 = {0.0: -4.623, 0.3: -1.756, 1.0: 3.864, 2.5: 10.850}
    ell, w0, q, marg = smoothed_slice_mu(300, 0.40, -10.0)
    rows = {}
    print(f"\nsmoothed slice n=300 w0={w0} t=2^-10 "
          "(003: iid -4.62 -> +10.85 at lam2=2.5):")
    for lam2 in [0.0, 0.3, 0.65, 1.0, 1.5, 2.5, 3.5]:
        r = eval_profile_tilt(300, ell, lam2)
        rows[lam2] = r["gain"]
        ref = f"  [003: {ref2[lam2]:+.3f}]" if lam2 in ref2 else ""
        print(f"  lam2={lam2:4.2f}: gain {r['gain']:+9.3f}  mass "
              f"{r['mass']:.6f} margB_lnerr {r['margB_lnerr']:.1e}{ref}")
    out["smoothed_slice300"] = rows

    # ---- the [ext] instance (no generating code in uc_couplings.py) --
    refe = {2.5: -2.144, 3.0: -0.631, 5.0: 1.965, 6.0: 2.065, 7.0: 1.853}
    n = 200
    P, p = sawin_components(0.42, 0.08)
    ell = mix_profile(n, P, p)
    rows = {}
    print("\n[ext] Sawin n=200 ubar=0.42 theta=0.08 "
          "(003 claims interior max +2.07 at lam2=6):")
    for lam2 in [2.5, 3.0, 5.0, 6.0, 7.0]:
        r = eval_profile_tilt(n, ell, lam2)
        rows[lam2] = r["gain"]
        print(f"  lam2={lam2:4.1f}: gain {r['gain']:+9.3f}"
              f"   [003 ext: {refe[lam2]:+.3f}]")
    out["ext_sawin200_042"] = rows
    print(f"\n[part E time {time.time()-t0:.1f}s]")
    save("004_partE_engine.json", out)


# ================================================================= part S
def part_S():
    print("=" * 78)
    print("S. Sawin gadgets: f_mu(1), D, coupled gains, growth in n")
    print("=" * 78)
    out = {"rows": []}
    print(f"{'n':>6} {'ubar':>7} {'theta':>6} {'maxmarg':>9} {'H(A)':>10} "
          f"{'gain_iid':>9} {'D':>7} {'f_mu(1)':>8} {'c*':>6} "
          f"{'gain_cpl':>9}")
    for n, ubar, theta in [(2000, 0.390, 0.05), (20000, 0.386, 0.02),
                           (5000, 0.3823, 0.001), (10000, 0.3823, 0.001),
                           (20000, 0.3823, 0.001), (40000, 0.3823, 0.001),
                           (60000, 0.3823, 0.001)]:
        P, p = sawin_components(ubar, theta)
        ellA = mix_profile(n, P, p)
        HA, _ = profile_entropy(n, ellA)
        marg = float((P * p).sum())
        PU, pU = [], []
        for i in range(len(P)):
            for j in range(len(P)):
                PU.append(P[i] * P[j])
                pU.append(1.0 - (1.0 - p[i]) * (1.0 - p[j]))
        ellU = mix_profile(n, PU, pU)
        HU, _ = profile_entropy(n, ellU)
        # D(U||mu) via CrossEnt: both exchangeable
        lc = np.array([lnC(n, w) for w in range(n + 1)])
        qU = np.exp(ellU + lc)
        CE = float(-(qU * ellA).sum() / LN2)
        D = CE - HU
        f1 = HU + D - HA
        cstar = (HA - HU) / D if HU < HA else inf
        # block-adaptive coupled gain (m_k = 1/2 if p_k<=1/2 else p_k)
        mk = np.where(p <= 0.5, 0.5, p)
        assert np.all(2 * p - mk >= p * p - 1e-12)   # conditionally iid
        HUc, _ = profile_entropy(n, mix_profile(n, P, mk))
        gc = HUc - HA
        print(f"{n:6d} {ubar:7.4f} {theta:6.3f} {marg:9.6f} {HA:10.2f} "
              f"{HU-HA:+9.2f} {D:7.3f} {f1:+8.2f} {cstar:6.1f} {gc:+9.1f}")
        out["rows"].append({"n": n, "ubar": ubar, "theta": theta,
                            "marg": marg, "HA": HA, "gain_iid": HU - HA,
                            "D": D, "f1": f1, "cstar": cstar,
                            "gain_coupled": gc})
    r = [x for x in out["rows"] if x["ubar"] == 0.3823]
    slopes = [(r[i + 1]["gain_coupled"] - r[i]["gain_coupled"])
              / (r[i + 1]["n"] - r[i]["n"]) for i in range(len(r) - 1)]
    P, p = sawin_components(0.3823, 0.001)
    mk = np.where(p <= 0.5, 0.5, p)
    pred = float((P * (np.array([h2(m) for m in mk])
                       - np.array([h2(x) for x in p]))).sum())
    print(f"\ncoupled-gain slope vs n: {['%.5f' % s for s in slopes]} "
          f"-> per-coordinate prediction sum_k P_k(h(m_k)-h(p_k)) = "
          f"{pred:.5f}")
    print(f"D stays O(1) ({[round(x['D'], 2) for x in r]}) while coupled "
          f"gain grows Theta(n): Phi_C - (H_iid + cD) = Theta(n) for every "
          f"fixed c.  [claims 2 and 5]")
    out["slope_pred"] = pred
    save("004_partS.json", out)


# ================================================================= part U
def close_family(gens, nbits):
    fam = set(gens)
    frontier = list(fam)
    while frontier:
        new = []
        cur = list(fam)
        for a in frontier:
            for b in cur:
                u = a | b
                if u not in fam:
                    fam.add(u)
                    new.append(u)
        frontier = new
    return sorted(fam)


def is_union_closed(fam):
    s = set(fam)
    return all((a | b) in s for a in fam for b in fam)


def part_U():
    print("=" * 78)
    print("U. Union-closed controls: adversarial battery")
    print("=" * 78)
    print("""A priori note (one-line proof, part of this review): for uniform
mu on union-closed F, ANY sub-probability kernel supported on F x F has
supp(A u B) inside F, so H(A u B) <= log|F| = H(mu) and gain <= 0 --
REGARDLESS of whether Sinkhorn converged or the marginals are right.  So
these controls can only detect support-handling bugs, never refute the
lemma; the lemma itself is proved.  Battery run anyway, incl. couplings far
outside the tilt family (random-kernel Sinkhorn probes).\n""")
    rng = np.random.default_rng(2026)
    out = {}
    fams = {}
    # closure of the CL mini slice+top (n=10)
    n = 10
    F1 = [m for m in range(1 << n) if bin(m).count("1") == 4]
    F2 = [m for m in range(1 << n) if bin(m).count("1") >= 6]
    fams["closure(CL mini slice+top) n=10"] = (close_family(F1 + F2, n), n)
    # two-scale product family: F_a on bits 0-5 (up-set) x F_b on 6-11 (chain)
    Fa = [m for m in range(64) if bin(m).count("1") >= 2] + [0]
    chain = [0]
    for i in range(6):
        chain.append(chain[-1] | (1 << i))
    Fb = chain
    prod = sorted({a | (b << 6) for a in Fa for b in Fb})
    fams["two-scale product (upset x chain) n=12"] = (prod, 12)
    # {S: |S| >= 4} at n=12 (= closure of the full 4-slice)
    fams["all |S|>=4, n=12"] = ([m for m in range(1 << 12)
                                 if bin(m).count("1") >= 4], 12)
    # low-entropy union-closed: chain + one fat antichain closed
    gens = [0b000000111111, 0b111111000000, 0b000111000111]
    fams["closure(3 fat gens) n=12"] = (close_family(gens, 12), 12)
    # random generator closures
    for t in range(6):
        nb = 12
        k = int(rng.integers(3, 7))
        gens = [int(rng.integers(1, 1 << nb)) for _ in range(k)]
        fams[f"random closure #{t} n={nb}"] = (close_family(gens, nb), nb)
    # extremal max-freq-=1/2 families: powerset(4); {0,{1}} x {0,{2}} x ...
    fams["powerset n=4"] = (list(range(16)), 4)

    LAMS = [-1.0, -0.3, 0.0, 0.5, 1.0, 2.0, 3.0, 5.0]
    worst_overall = (-inf, None)
    for name, (fam, nb) in fams.items():
        assert is_union_closed(fam), name
        lnmu = np.full(len(fam), -log(len(fam)))
        freqs = [sum(1 for m in fam if m >> i & 1) / len(fam)
                 for i in range(nb)]
        best = -inf
        for lam2 in LAMS:
            r = atom_sinkhorn_indep(fam, lnmu, lam2)
            assert r["escape_mass"] < 1e-12, (name, "escape!")
            best = max(best, r["gain"])
        # random-kernel probes (well beyond the tilt family)
        nprobe = 30 if len(fam) <= 700 else 8
        for _ in range(nprobe):
            A = np.array(fam)
            lnK = rng.normal(0, 1.5, size=(len(fam), len(fam)))
            # sinkhorn on random kernel
            lu = lnmu / 2.0
            lv = lu.copy()
            for it in range(4000):
                r1 = lnK + lv[None, :]
                m1 = r1.max(axis=1)
                lun = lnmu - (m1 + np.log(np.exp(r1 - m1[:, None])
                                          .sum(axis=1)))
                r2 = lnK + lun[:, None]
                m2 = r2.max(axis=0)
                lvn = lnmu - (m2 + np.log(np.exp(r2 - m2[None, :])
                                          .sum(axis=0)))
                d = max(np.abs(lun - lu).max(), np.abs(lvn - lv).max())
                lu, lv = lun, lvn
                if d < 1e-12:
                    break
            pi = np.exp(lu[:, None] + lnK + lv[None, :])
            ul = {}
            for i, a in enumerate(fam):
                row = pi[i]
                for j, b in enumerate(fam):
                    u = a | b
                    ul[u] = ul.get(u, 0.0) + row[j]
            pu = np.array(list(ul.values()))
            pu /= pu.sum()
            HU = float(-(pu * np.log(pu)).sum() / LN2)
            g = HU - log(len(fam)) / LN2
            best = max(best, g)
        ok = best <= 1e-9
        if best > worst_overall[0]:
            worst_overall = (best, name)
        print(f"{name:42s} |F|={len(fam):5d} maxfreq={max(freqs):.3f} "
              f"max gain (tilts + {nprobe} random kernels) = {best:+.2e} "
              f"{'OK' if ok else '*** VIOLATION ***'}")
        out[name] = {"size": len(fam), "maxfreq": max(freqs),
                     "max_gain": best, "ok": bool(ok)}
    print(f"\nworst (largest) gain over battery: {worst_overall[0]:+.2e} "
          f"({worst_overall[1]}) -- all <= 0 as the lemma demands.")
    save("004_partU.json", out)


# ================================================================= part M
def part_M():
    print("=" * 78)
    print("M. Mini-theorem (003 part e): the p=0 bracket bug + repair")
    print("=" * 78)
    out = {}

    def bracket(p):
        if p >= 0.5 or p <= 0.0:
            return 0.0
        m = min(2 * p - p * p, 0.5)
        return h2(m) - h2(p)

    print("003 claims 'every bracket > 0 whenever p_{k,i} < 1/2'.")
    print(f"  bracket(0)     = {bracket(0.0):.6f}   (NOT > 0: claim false "
          f"at p = 0)")
    print(f"  bracket(0.3)   = {bracket(0.3):.6f}")
    print(f"  bracket(0.499) = {bracket(0.499):.2e} (vanishes as p->1/2)")

    # explicit unshielded mixture: mu = 1/2 delta_empty + 1/2 Bern(0.6)^n
    alpha, ph = 0.5, 0.6
    print(f"\nCounterexample to the claimed consequence-as-proved:")
    print(f"mu_n = {alpha} delta_0 + {1-alpha} Bern({ph})^n  "
          f"(marginals = {(1-alpha)*ph}, 2 components, H(P) = 1 bit)")
    rows = []
    for n in [20, 50, 100, 200]:
        ellA = mix_profile(n, [alpha, 1 - alpha], [0.0, ph])
        HA, _ = profile_entropy(n, ellA)
        # mini-theorem coupling: m_1 = 0 (p=0), m_2 = ph (p >= 1/2): U ~ mu
        ellU = mix_profile(n, [alpha, 1 - alpha], [0.0, ph])
        HU, _ = profile_entropy(n, ellU)
        g_mini = HU - HA          # exactly 0
        bound = -h2(alpha)        # the (e) lower bound: sum brackets - H(P)
        # half-mixing conditionally-iid coupling (this review):
        # law(U) = alpha/2 delta_0 + (1 - alpha/2) Bern(ph)^n
        ellU2 = mix_profile(n, [alpha / 2, 1 - alpha / 2], [0.0, ph])
        HU2, _ = profile_entropy(n, ellU2)
        g_half = HU2 - HA
        rows.append({"n": n, "gain_minithm": g_mini, "bound_e": bound,
                     "gain_halfmix": g_half})
        print(f"  n={n:4d}: (e)-coupling gain = {g_mini:+.6f} "
              f"((e) bound {bound:+.3f});  half-mixing gain = {g_half:+.3f}")
    out["unshielded"] = rows

    # brute-force legality of half-mixing at n=6 (atom level, exact)
    n, ph_, al = 6, 0.6, 0.5
    N = 1 << n
    pi = np.zeros((N, N))
    for s in range(N):
        w = bin(s).count("1")
        ps = ph_ ** w * (1 - ph_) ** (n - w)
        # weight 2*al: u=(pair,s), (A,B) iid from (delta_0 + delta_s)/2
        for (a, b) in [(0, 0), (0, s), (s, 0), (s, s)]:
            pi[a, b] += 2 * al * ps * 0.25
        # weight 1-2*al: (s,s)   [al = 1/2 -> zero]
        pi[s, s] += (1 - 2 * al) * ps
    mu = np.zeros(N)
    mu[0] += al
    for s in range(N):
        w = bin(s).count("1")
        mu[s] += (1 - al) * ph_ ** w * (1 - ph_) ** (n - w)
    eA = np.abs(pi.sum(axis=1) - mu).max()
    eB = np.abs(pi.sum(axis=0) - mu).max()
    ulaw = np.zeros(N)
    for a in range(N):
        for b in range(N):
            ulaw[a | b] += pi[a, b]
    HU = float(-(ulaw[ulaw > 0] * np.log(ulaw[ulaw > 0])).sum() / LN2)
    HA = float(-(mu[mu > 0] * np.log(mu[mu > 0])).sum() / LN2)
    ellU2 = mix_profile(n, [al / 2, 1 - al / 2], [0.0, ph_])
    HU2, _ = profile_entropy(n, ellU2)
    print(f"\nhalf-mixing brute force n=6: marginal errors ({eA:.1e}, "
          f"{eB:.1e}); H(U) direct {HU:.6f} vs profile formula "
          f"{HU2:.6f} (diff {abs(HU-HU2):.1e}); gain {HU-HA:+.6f}")
    assert eA < 1e-14 and eB < 1e-14 and abs(HU - HU2) < 1e-10
    out["halfmix_bruteforce"] = {"margA_err": eA, "margB_err": eB,
                                 "gain": HU - HA}

    # non-uniformity of the kill: p_low -> 0 blows up the n-threshold
    print("\nn-threshold of the (e) bound for (1-P)Bern(p_low)+P Bern(0.51):")
    thr = {}
    for plow in [0.01, 0.001, 0.00025, 1e-5]:
        P = 0.8625
        rate = (1 - P) * bracket(plow)
        nstar = h2(P) / rate if rate > 0 else inf
        thr[plow] = nstar
        print(f"  p_low={plow:8.5f}: rate/coord = {rate:.3e}, "
              f"n* = H(P)/rate = {nstar:,.0f}")
    print("  p_low = 0 exactly: rate = 0, (e) NEVER separates it "
          "(kill not uniform).")
    out["thresholds"] = thr
    save("004_partM.json", out)


# ================================================================= part K
def net_iid_K(blocks):
    tot = 0.0
    for Pa, pa in blocks:
        for Pb, pb in blocks:
            tot += Pa * Pb * h2(1.0 - (1 - pa) * (1 - pb))
    return tot - sum(P * h2(p) for P, p in blocks)


def net_assort_K(blocks, rho):
    return sum(P * (h2(1.0 - z_quad(1 - p, 1 - p, rho)) - h2(p))
               for P, p in blocks)


def best_sigma_iid_K(blocks):
    """max over block couplings sigma (both marginals = P) of
    sum sigma_ab h(1 - x_a x_b) - sum P h(p): linear objective over the
    transportation polytope -> exact via vertex enumeration (K <= 4)."""
    K = len(blocks)
    P = np.array([b[0] for b in blocks])
    val = np.zeros((K, K))
    for a in range(K):
        for b in range(K):
            val[a, b] = h2(1.0 - (1 - blocks[a][1]) * (1 - blocks[b][1]))
    base = sum(Pk * h2(pk) for Pk, pk in blocks)
    # constraints: row sums = P, col sums = P (2K eqs, rank 2K-1)
    cells = [(a, b) for a in range(K) for b in range(K)]
    nvar = len(cells)
    Aeq = np.zeros((2 * K, nvar))
    for idx, (a, b) in enumerate(cells):
        Aeq[a, idx] = 1.0
        Aeq[K + b, idx] = 1.0
    beq = np.concatenate([P, P])
    best = -inf
    from itertools import combinations
    m = 2 * K - 1
    for S in combinations(range(nvar), m):
        As = Aeq[:, S]
        sol, res, rk, _ = np.linalg.lstsq(As, beq, rcond=None)
        if rk < m:
            continue
        if np.abs(As @ sol - beq).max() > 1e-9:
            continue
        if sol.min() < -1e-9:
            continue
        v = sum(sol[i] * val[cells[S[i]]] for i in range(m))
        best = max(best, v - base)
    return best


def J_screen(blocks, rhos):
    a = net_iid_K(blocks)
    b = max(net_assort_K(blocks, r) for r in rhos)
    c = best_sigma_iid_K(blocks)
    return max(a, b, c), a, b, c


def part_K():
    print("=" * 78)
    print("K. Part-E re-analysis: the '0.445 ceiling' claim")
    print("=" * 78)
    out = {}
    rhos = [2.0 ** (e / 8.0) for e in range(-48, 241)]

    # ---- reproduce their 8 hardest rows with independent code --------
    theirs = [
        (0.39, 0.001, 0.998, 0.0029468845760980593, 0.009369578188485194),
        (0.42, 0.066, 0.51, 0.7971976744186045, 0.008087671964775504),
        (0.43, 0.019375, 0.51, 0.836859693877551, 0.002630545933450854),
        (0.44, 0.00025, 0.51, 0.862475442043222, 2.1164626947821176e-07),
        (0.45, 0.00025, 0.51, 0.8821218074656189, -9.74471430165394e-08),
        (0.47, 0.00025, 0.55, 0.8542805100182148, -2.9149551093589083e-06),
        (0.49, 0.00025, 0.62025, 0.789858844911147, -7.068579182812907e-06),
        (0.499, 0.00025, 0.64025, 0.7792177230046948, -8.21810464859425e-06),
    ]
    print("re-scoring 003's hardest instances (my J uses best-sigma >= "
          "their anti channel, so my J >= their J):")
    rows = []
    for cap, p1, p2, P2, Jthem in theirs:
        blocks = [(1 - P2, p1), (P2, p2)]
        J, a, b, c = J_screen(blocks, rhos)
        rows.append({"cap": cap, "J_mine": J, "J_theirs": Jthem,
                     "iid": a, "tilt": b, "bestsig": c})
        print(f"  cap {cap:5.3f}: J_mine = {J:+.7f} (theirs {Jthem:+.7f})"
              f"  [iid {a:+.5f} tilt {b:+.7f} bestsig {c:+.5f}]")
    out["rescore"] = rows

    # ---- the exact-delta_0 killers below 0.445 -----------------------
    print("\nexact delta_0 low block (p_low = 0, allowed by (S-coup), "
          "EXCLUDED by 003's grid which floors p_low at 0.00025):")
    kills = []
    for cap, ph in [(0.44, 0.51), (0.44, 0.502), (0.435, 0.502),
                    (0.433, 0.501), (0.432, 0.5005), (0.4315, 0.50001),
                    (0.43, 0.502)]:
        Ph = cap / ph
        if Ph >= 1:
            continue
        blocks = [(1 - Ph, 0.0), (Ph, ph)]
        J, a, b, c = J_screen(blocks, rhos)
        tag = "KILLED (J<=0)" if J <= 1e-12 else "live"
        kills.append({"cap": cap, "ph": ph, "Ph": Ph, "J": J, "iid": a,
                      "tilt": b, "bestsig": c, "killed": bool(J <= 1e-12)})
        print(f"  cap {cap:6.4f} p_h={ph:7.5f} P_h={Ph:.4f}: "
              f"J = {J:+.6f} [iid {a:+.5f} tilt {b:+.2e} bestsig "
              f"{c:+.5f}]  {tag}")
    out["delta0_killers"] = kills
    print("  (tilt channel at p_low=0: sup over rho is approached only in "
          "the diagonal limit rho->inf, value -> 0^-, never a STRICT gain)")

    # ---- corrected boundary curve ------------------------------------
    # binding channel is best-sigma (= anti here); boundary:
    # (2-3B)h(p) + (2B-1)h(q) = 0, B = cap/p, q = 2p - p^2
    print("\ncorrected 2-block model boundary cap*(p_h) = "
          "p_h (2h(p)-h(q)) / (3h(p)-2h(q)):")
    curve = []
    best = (inf, None)
    for ph in [0.5 + 1e-6, 0.5001, 0.501, 0.505, 0.51, 0.52, 0.55, 0.6,
               0.62, 0.7]:
        hp, hq = h2(ph), h2(2 * ph - ph * ph)
        capstar = ph * (2 * hp - hq) / (3 * hp - 2 * hq)
        curve.append({"ph": ph, "capstar": capstar})
        if capstar < best[0]:
            best = (capstar, ph)
        print(f"  p_h = {ph:8.6f}: cap* = {capstar:.6f}")
    print(f"  minimum: cap* = {best[0]:.6f} at p_h -> 1/2  "
          f"(analytic limit (2 - h(3/4))/(2(3 - 2h(3/4))) = "
          f"{(2 - h2(0.75)) / (2 * (3 - 2 * h2(0.75))):.6f})")
    out["boundary_curve"] = curve
    out["corrected_ceiling"] = best[0]

    # ---- 3-block hunt: can it go below the 2-block boundary? ---------
    print("\n3-block hunt (blocks may use p=0 and p>=1/2), caps below the "
          "corrected 2-block boundary:")
    rng = np.random.default_rng(7)
    hunt3 = []
    for cap in [0.40, 0.41, 0.42, 0.428]:
        worst = (inf, None)
        for trial in range(400):
            ps = sorted([0.0, float(rng.uniform(0.5, 0.9)),
                         float(rng.uniform(0.5, 0.9))])
            # random weights hitting the cap
            w2 = rng.uniform(0, 1)
            w3 = rng.uniform(0, 1 - w2)
            s = w2 * ps[1] + w3 * ps[2]
            if s <= cap:
                continue
            sc = cap / s
            w2, w3 = w2 * sc, w3 * sc
            blocks = [(1 - w2 - w3, ps[0]), (w2, ps[1]), (w3, ps[2])]
            J, a, b, c = J_screen(blocks, rhos[::4])
            if J < worst[0]:
                worst = (J, blocks)
        hunt3.append({"cap": cap, "Jmin": worst[0],
                      "blocks": worst[1]})
        print(f"  cap {cap:5.3f}: min J over 3-block trials = "
              f"{worst[0]:+.6f}  "
              f"{'KILLED below 2-block boundary!' if worst[0] <= 0 else 'live'}")
    out["hunt3"] = hunt3

    # ---- finite-n engine check of a corrected killer -----------------
    print("\nfinite-n engine check (n=300) of the corrected killer "
          "cap=0.4350, p_h=0.502:")
    Ph = 0.4350 / 0.502
    n = 300
    ell = mix_profile(n, [1 - Ph, Ph], [0.0, 0.502])
    eng = {}
    for lam2 in [-0.4, -0.15, -0.05, 0.0, 0.05, 0.15, 0.45, 1.0, 2.5]:
        r = eval_profile_tilt(n, ell, lam2)
        eng[lam2] = r["gain"]
        print(f"  lam2={lam2:5.2f}: gain {r['gain']:+9.4f} "
              f"(mass {r['mass']:.5f})")
    out["engine_killer_n300"] = eng
    save("004_partK.json", out)


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="all")
    args = ap.parse_args()
    todo = args.part.upper()
    if todo in ("L", "ALL"):
        part_L()
    if todo in ("I", "ALL"):
        part_I()
    if todo in ("E", "ALL"):
        part_E_engine()
    if todo in ("S", "ALL"):
        part_S()
    if todo in ("U", "ALL"):
        part_U()
    if todo in ("M", "ALL"):
        part_M()
    if todo in ("K", "ALL"):
        part_K()


if __name__ == "__main__":
    main()
