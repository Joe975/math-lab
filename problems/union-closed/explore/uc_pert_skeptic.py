#!/usr/bin/env python3
"""Attempt 012: skeptic review of 008 (perturbative assembly).  Default
stance: refute.

Independent implementations (nothing imported from uc_pert.py):
  * direct engine: alternating-scaling Sinkhorn written from scratch
    (stable z_rho root via 2*rho*x*y/(S+sqrt(disc)), both-margin residuals,
    symmetry of the fitted coupling checked, not assumed);
  * orbit engine: profile-space Sinkhorn exploiting the S_{n-2} symmetry of
    crash-direction mixtures (states (a1,a2,k), exact combinatorial census)
    -- structurally unrelated to 008's dense atom engine, and able to run
    the budget census to n = 32, far past 008's n <= 9.

Parts:
  SK0  validation: product identity, closed-form gain, crash identity,
       orbit-vs-direct cross-check at n = 6, 7.
  SK1  calibration claims: rho*(0.383), rho*(0.38271), g(1.031, 0.383) < 0,
       margin table, AHS coefficient h'(psi) phi^6 (closed form re-derived).
  SK2  Theorem P6' bound test: max |log2 OR - lam| <= 924 min(p,1-p)^{-3n}
       delta^2 on a stress grid incl. delta = 0.5; Lemma B first-order
       formula at an interior delta; Lemma C sup-bounds |vdot|, |vddot|.
  SK3  spot-verification of 008's 78-run assembly table + hunt for a
       marginal-capped negative instance 008's grid missed (new direction
       half-empty = 0.5 delta_empty + 0.5 Bern(2p): marginals exactly p,
       the 004 killer genre inside the class; rho = 1.06 cap-stress cells
       008 never ran; off-grid delta; rho = 1.045 near threshold).
  SK4  conditional-theorem constants: tau_half at box center vs corner
       (the box-uniformity claim), M_tau on [0,tau_half] vs sup|G_tau| on
       tau < 0 (where (B1) losses actually live), M2, L2, ||G||_box (needed
       by the stated delta_0 formula but never computed by 008), global
       inf G, corrected delta_0 variants.
  SK5  the kill shot: downward budget census avg_i E[(lam - log2 OR)_+] /
       delta^2, crash direction, both regimes, n = 4..16 (008: n <= 9),
       plus avg tax and E[(xt-x0)^2] trends (B2/B3 shapes).
  SK6  centered kernel: crash OR = lam(3-n)/n and mirror lam(n-1)/n by the
       potential-free identity AND by Sinkhorn fit; product band at n = 6, 8.

Standard library only.  Deterministic (no RNG anywhere).  Runtime ~90 s
for all parts (SK5B/SK5C/SK7 extend the census to n = 32).
Usage:
    python problems/union-closed/explore/uc_pert_skeptic.py [--part SK0|...|all]
Checkpoints: problems/union-closed/data/pertsk_SK*.json (refuse overwrite)
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

LN2 = math.log(2.0)
PSI = (3.0 - math.sqrt(5.0)) / 2.0
PHI = (math.sqrt(5.0) - 1.0) / 2.0  # 1 - PSI; PHI**2 == PSI


def h2(x: float) -> float:
    """Binary entropy in bits."""
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return (-x * math.log(x) - (1.0 - x) * math.log(1.0 - x)) / LN2


def zplack(x: float, y: float, rho: float) -> float:
    """Plackett both-zero probability: root of a z^2 - S z + rho x y = 0,
    a = rho - 1, S = 1 + a(x+y), via the stable conjugate form
    z = 2 rho x y / (S + sqrt(S^2 - 4 a rho x y)) (valid for both signs of a;
    no special-casing near rho = 1 needed)."""
    if x <= 0.0 or y <= 0.0:
        return 0.0
    if x >= 1.0:
        return y
    if y >= 1.0:
        return x
    a = rho - 1.0
    if a == 0.0:
        return x * y
    S = 1.0 + a * (x + y)
    disc = S * S - 4.0 * a * rho * x * y
    z = 2.0 * rho * x * y / (S + math.sqrt(max(disc, 0.0)))
    lo, hi = max(0.0, x + y - 1.0), min(x, y)
    return min(max(z, lo), hi)


def rho_star(p: float) -> float:
    return p * (3.0 * p - 1.0) / (1.0 - 2.0 * p) ** 2


def margin_g(rho: float, p: float) -> float:
    x = 1.0 - p
    return h2(zplack(x, x, rho)) - h2(p)


# ------------------------------------------------------------ direct engine

def prod_mu(n: int, p: float) -> dict[int, float]:
    out = {}
    for m in range(1 << n):
        k = bin(m).count("1")
        out[m] = p ** k * (1.0 - p) ** (n - k)
    return out


def mix(base: dict[int, float], nu: dict[int, float], d: float,
        n: int) -> dict[int, float]:
    out = {}
    for m in range(1 << n):
        out[m] = (1.0 - d) * base.get(m, 0.0) + d * nu.get(m, 0.0)
    return out


def crash_nu(n: int) -> dict[int, float]:
    full = (1 << n) - 1
    return {0b10: 0.33, full ^ 0b11: 0.33, full: 0.04, 0b01: 0.30}


def mirror_nu(n: int) -> dict[int, float]:
    full = (1 << n) - 1
    return {0b10 | (full ^ 0b11): 0.3, 0: 0.2, full: 0.3, 0b01: 0.2}


def emb0k2_nu(n: int) -> dict[int, float]:
    """008's zero-prefix embedded crash, k = 2 (re-built from its definition)."""
    zeros = 0b11
    tail_full = ((1 << n) - 1) ^ zeros
    s1, s4 = 1 << 3, 1 << 2
    s2 = tail_full ^ s1 ^ s4
    return {s1: 0.33, s2: 0.33, tail_full: 0.04, s4: 0.30}


def pm03_nu(n: int, p: float) -> dict[int, float]:
    lo, hi = prod_mu(n, p - 0.3), prod_mu(n, p + 0.3)
    return {m: 0.5 * lo[m] + 0.5 * hi[m] for m in range(1 << n)}


def halfempty_nu(n: int, p: float) -> dict[int, float]:
    """0.5 delta_empty + 0.5 Bern(2p)^n: element marginals exactly p.
    The 004 two-blocks-with-opposite-needs genre, inside the capped class."""
    hi = prod_mu(n, 2.0 * p)
    out = {m: 0.5 * hi[m] for m in range(1 << n)}
    out[0] = out.get(0, 0.0) + 0.5
    return out


def slice_nu(n: int, p: float) -> dict[int, float]:
    w0 = round(p * n)
    atoms = [m for m in range(1 << n) if bin(m).count("1") == w0]
    return {m: 1.0 / len(atoms) for m in atoms}


def elem_marg(mu: dict[int, float], n: int) -> list[float]:
    out = [0.0] * n
    for m, w in mu.items():
        for j in range(n):
            if m >> j & 1:
                out[j] += w
    return out


def Hbits(mu) -> float:
    vals = mu.values() if isinstance(mu, dict) else mu
    return -sum(w * math.log(w) for w in vals if w > 0.0) / LN2


def sinkhorn(mu: dict[int, float], lam: float, expo=None,
             tol: float = 1e-13, iters: int = 100000):
    """Alternating scaling, both-margin residual, from scratch.
    expo(a, b) defaults to lam * |a & b| (log2 kernel exponent).
    Returns (atoms, r, c, K, resid, sym_dev)."""
    atoms = sorted(m for m, w in mu.items() if w > 0.0)
    N = len(atoms)
    if expo is None:
        K = [[2.0 ** (lam * bin(a & b).count("1")) for b in atoms]
             for a in atoms]
    else:
        K = [[2.0 ** expo(a, b) for b in atoms] for a in atoms]
    tgt = [mu[a] for a in atoms]
    r = [1.0] * N
    c = [1.0] * N
    for it in range(iters):
        for j in range(N):
            s = 0.0
            Kj = K[j]  # K symmetric in |a&b|; rows double as columns
            for i in range(N):
                s += r[i] * Kj[i]
            c[j] = tgt[j] / s
        for i in range(N):
            s = 0.0
            Ki = K[i]
            for j in range(N):
                s += c[j] * Ki[j]
            r[i] = tgt[i] / s
        if it % 5 == 4 or it < 3:
            col_err = max(abs(c[j] * sum(r[i] * K[j][i] for i in range(N))
                              - tgt[j]) for j in range(N))
            if col_err < tol:
                break
    row_err = max(abs(r[i] * sum(c[j] * K[i][j] for j in range(N)) - tgt[i])
                  for i in range(N))
    col_err = max(abs(c[j] * sum(r[i] * K[j][i] for i in range(N)) - tgt[j])
                  for j in range(N))
    ratios = [r[i] / c[i] for i in range(N)]
    sym_dev = (max(ratios) - min(ratios)) / max(abs(ratios[0]), 1e-300)
    return atoms, r, c, K, max(row_err, col_err), sym_dev


def census_direct(atoms, r, c, K, n: int, lam: float, x0: float):
    """Per-coordinate stats of pi(A,B) = r_A c_B K(A,B).  Independent census."""
    N = len(atoms)
    out = []
    for i in range(1, n + 1):
        pref = (1 << (i - 1)) - 1
        bit = 1 << (i - 1)
        cells: dict[tuple[int, int], list[float]] = {}
        for ii in range(N):
            A = atoms[ii]
            rA = r[ii]
            Ki = K[ii]
            key_a = A & pref
            sa = 2 if A & bit else 0
            for jj in range(N):
                B = atoms[jj]
                cl = cells.setdefault((key_a, B & pref), [0.0] * 4)
                cl[sa + (1 if B & bit else 0)] += rA * c[jj] * Ki[jj]
        E_hz = E_hx = E_hy = E_xi2 = 0.0
        max_dev = down_sum = or_lt1_mass = 0.0
        min_l2 = float("inf")
        max_l2 = -float("inf")
        for cl in cells.values():
            mass = cl[0] + cl[1] + cl[2] + cl[3]
            if mass <= 0.0:
                continue
            z = cl[0] / mass
            xt = (cl[0] + cl[1]) / mass
            yt = (cl[0] + cl[2]) / mass
            E_hz += mass * h2(z)
            E_hx += mass * h2(xt)
            E_hy += mass * h2(yt)
            E_xi2 += mass * (xt - x0) ** 2
            if min(cl) > 0.0:
                l2 = math.log2(cl[3] * cl[0] / (cl[2] * cl[1]))
                dev = abs(l2 - lam)
                if dev > max_dev:
                    max_dev = dev
                if l2 < lam:
                    down_sum += (lam - l2) * mass
                if l2 < 0.0:
                    or_lt1_mass += mass
                min_l2 = min(min_l2, l2)
                max_l2 = max(max_l2, l2)
        out.append({"i": i, "E_hz": E_hz, "E_hx": E_hx, "E_hy": E_hy,
                    "E_xi2": E_xi2, "max_abs_dev": max_dev,
                    "down": down_sum, "or_lt1_mass": or_lt1_mass,
                    "min_l2": min_l2, "max_l2": max_l2})
    return out


def chain_H_direct(mu: dict[int, float], n: int) -> list[float]:
    out = []
    for i in range(1, n + 1):
        pref = (1 << (i - 1)) - 1
        bit = 1 << (i - 1)
        tot: dict[int, float] = {}
        zer: dict[int, float] = {}
        for m, w in mu.items():
            a = m & pref
            tot[a] = tot.get(a, 0.0) + w
            if not m & bit:
                zer[a] = zer.get(a, 0.0) + w
        out.append(sum(w * h2(zer.get(a, 0.0) / w)
                       for a, w in tot.items() if w > 0.0))
    return out


def union_gain_direct(atoms, r, c, K, mu):
    N = len(atoms)
    ul: dict[int, float] = {}
    for ii in range(N):
        A = atoms[ii]
        rA = r[ii]
        Ki = K[ii]
        for jj in range(N):
            u = A | atoms[jj]
            ul[u] = ul.get(u, 0.0) + rA * c[jj] * Ki[jj]
    return Hbits(ul) - Hbits(mu)


def run_instance(mu, n, lam, p_base):
    atoms, r, c, K, resid, sym = sinkhorn(mu, lam)
    st = census_direct(atoms, r, c, K, n, lam, 1.0 - p_base)
    Ha = chain_H_direct(mu, n)
    nets = [s["E_hz"] - Ha[s["i"] - 1] for s in st]
    taxes = [Ha[s["i"] - 1] - s["E_hx"] for s in st]
    return {
        "resid": resid, "sym_dev": sym,
        "min_net": min(nets), "gain_lb": sum(nets),
        "gain": union_gain_direct(atoms, r, c, K, mu),
        "max_tax": max(taxes), "avg_tax": sum(taxes) / n,
        "max_abs_dev": max(s["max_abs_dev"] for s in st),
        "or_lt1_mass": sum(s["or_lt1_mass"] for s in st),
        "max_E_xi2": max(s["E_xi2"] for s in st),
        "down_avg": sum(s["down"] for s in st) / n,
        "down_max": max(s["down"] for s in st),
    }


# ------------------------------------------------------------- orbit engine
# States s = (a1, a2, k): membership of coordinates 1, 2 and tail count
# k = |A cap {3..n}|.  Valid for measures invariant under permutations of
# {3..n} -- the crash-direction mixtures are.

def orbit_states(n: int):
    r = n - 2
    return [(a1, a2, k) for a1 in (0, 1) for a2 in (0, 1)
            for k in range(r + 1)]


def orbit_mu_crash(n: int, p: float, d: float):
    """Per-atom mass mu_atom(s) of (1-d) Bern(p)^n + d crash."""
    r = n - 2
    out = {}
    for (a1, a2, k) in orbit_states(n):
        w = a1 + a2 + k
        out[(a1, a2, k)] = (1.0 - d) * p ** w * (1.0 - p) ** (n - w)
    out[(0, 1, 0)] += d * 0.33
    out[(0, 0, r)] += d * 0.33
    out[(1, 1, r)] += d * 0.04
    out[(1, 0, 0)] += d * 0.30
    return out


def wtail(r: int, lam: float):
    """Wt[k][l] = sum over tail sets T_B of size l of 2^{lam |T_A cap T_B|},
    T_A fixed of size k: sum_j C(k,j) C(r-k, l-j) 2^{lam j}."""
    two = [2.0 ** (lam * j) for j in range(r + 1)]
    Wt = [[0.0] * (r + 1) for _ in range(r + 1)]
    for k in range(r + 1):
        for l in range(r + 1):
            s = 0.0
            for j in range(max(0, k + l - r), min(k, l) + 1):
                s += math.comb(k, j) * math.comb(r - k, l - j) * two[j]
            Wt[k][l] = s
    return Wt


def orbit_sinkhorn(n: int, p: float, d: float, lam: float,
                   tol: float = 1e-14, iters: int = 500000):
    r = n - 2
    states = orbit_states(n)
    idx = {s: q for q, s in enumerate(states)}
    tgt = orbit_mu_crash(n, p, d)
    Wt = wtail(r, lam)
    M = [[2.0 ** (lam * (s[0] * t[0] + s[1] * t[1])) * Wt[s[2]][t[2]]
          for t in states] for s in states]
    NS = len(states)
    tv = [tgt[s] for s in states]
    rv = [1.0] * NS
    cv = [1.0] * NS
    for it in range(iters):
        for j in range(NS):
            s = 0.0
            for i in range(NS):
                s += rv[i] * M[j][i]
            cv[j] = tv[j] / s
        for i in range(NS):
            s = 0.0
            Mi = M[i]
            for j in range(NS):
                s += cv[j] * Mi[j]
            rv[i] = tv[i] / s
        if it % 10 == 9:
            err = max(abs(cv[j] * sum(rv[i] * M[j][i] for i in range(NS))
                          - tv[j]) for j in range(NS))
            if err < tol:
                break
    resid = max(
        max(abs(rv[i] * sum(cv[j] * M[i][j] for j in range(NS)) - tv[i])
            for i in range(NS)),
        max(abs(cv[j] * sum(rv[i] * M[j][i] for i in range(NS)) - tv[j])
            for j in range(NS)))
    R = {s: rv[idx[s]] for s in states}
    C = {s: cv[idx[s]] for s in states}
    return R, C, resid


def tmat(f: int, lam: float):
    """T[s][t] = sum_w C(f,w) C(f-w,s-w) C(f-s,t-w) 2^{lam w}."""
    two = [2.0 ** (lam * w) for w in range(f + 1)]
    T = [[0.0] * (f + 1) for _ in range(f + 1)]
    for s in range(f + 1):
        for t in range(f + 1):
            acc = 0.0
            for w in range(max(0, s + t - f), min(s, t) + 1):
                acc += (math.comb(f, w) * math.comb(f - w, s - w)
                        * math.comb(f - s, t - w) * two[w])
            T[s][t] = acc
    return T


def orbit_census(n: int, p: float, d: float, lam: float, R, C, x0: float):
    """Exact per-coordinate stats for the crash mixture via orbit sums.
    Returns list of dicts with down, max_abs_dev, E_hx, E_xi2, tax, mass."""
    r = n - 2
    orb4 = [(0, 0), (0, 1), (1, 0), (1, 1)]
    out = []

    def cell_stats(cell_list):
        E_hx = E_xi2 = down = or_lt1 = 0.0
        max_dev = 0.0
        tot_mass = 0.0
        for c00, c01, c10, c11 in cell_list:
            mass = c00 + c01 + c10 + c11
            if mass <= 0.0:
                continue
            tot_mass += mass
            xt = (c00 + c01) / mass
            E_hx += mass * h2(xt)
            E_xi2 += mass * (xt - x0) ** 2
            if min(c00, c01, c10, c11) > 0.0:
                l2 = math.log2(c11 * c00 / (c10 * c01))
                dev = abs(l2 - lam)
                max_dev = max(max_dev, dev)
                if l2 < lam:
                    down += (lam - l2) * mass
                if l2 < 0.0:
                    or_lt1 += mass
        return E_hx, E_xi2, down, max_dev, or_lt1, tot_mass

    for i in range(1, n + 1):
        cells = []
        if i == 1:
            T = tmat(r, lam)
            cell = [0.0] * 4
            for al in (0, 1):
                for be in (0, 1):
                    acc = 0.0
                    for a2 in (0, 1):
                        for b2 in (0, 1):
                            fac = 2.0 ** (lam * (al * be + a2 * b2))
                            ss = 0.0
                            for s in range(r + 1):
                                Ra = R[(al, a2, s)]
                                Ts = T[s]
                                for t in range(r + 1):
                                    ss += Ra * C[(be, b2, t)] * Ts[t]
                            acc += fac * ss
                    cell[2 * al + be] = acc
            cells.append((cell[0], cell[1], cell[2], cell[3]))
        elif i == 2:
            T = tmat(r, lam)
            for a in (0, 1):
                for b in (0, 1):
                    cell = [0.0] * 4
                    for al in (0, 1):
                        for be in (0, 1):
                            fac = 2.0 ** (lam * (a * b + al * be))
                            ss = 0.0
                            for s in range(r + 1):
                                Ra = R[(a, al, s)]
                                Ts = T[s]
                                for t in range(r + 1):
                                    ss += Ra * C[(b, be, t)] * Ts[t]
                            cell[2 * al + be] = fac * ss
                    cells.append(tuple(cell))
        else:
            q = i - 3
            f = n - i
            T = tmat(f, lam)
            # S[(oa, ka0)][(ob, kb0)] partial sums
            S = {}
            for oa in orb4:
                for ka0 in range(q + 2):
                    for ob in orb4:
                        for kb0 in range(q + 2):
                            ss = 0.0
                            for s in range(f + 1):
                                Ra = R[(oa[0], oa[1], ka0 + s)]
                                Ts = T[s]
                                for t in range(f + 1):
                                    ss += Ra * C[(ob[0], ob[1], kb0 + t)] \
                                        * Ts[t]
                            S[(oa, ka0, ob, kb0)] = ss
            for oa in orb4:
                for ob in orb4:
                    base = 2.0 ** (lam * (oa[0] * ob[0] + oa[1] * ob[1]))
                    for ja in range(q + 1):
                        for jb in range(q + 1):
                            for j in range(max(0, ja + jb - q),
                                           min(ja, jb) + 1):
                                Ncl = (math.comb(q, j)
                                       * math.comb(q - j, ja - j)
                                       * math.comb(q - ja, jb - j))
                                if Ncl == 0:
                                    continue
                                pre = Ncl * base * 2.0 ** (lam * j)
                                cell = [0.0] * 4
                                for al in (0, 1):
                                    for be in (0, 1):
                                        v = pre * 2.0 ** (lam * al * be) \
                                            * S[(oa, ja + al, ob, jb + be)]
                                        cell[2 * al + be] = v
                                cells.append(tuple(cell))
        E_hx, E_xi2, down, max_dev, or_lt1, tot = cell_stats(cells)
        out.append({"i": i, "E_hx": E_hx, "E_xi2": E_xi2, "down": down,
                    "max_abs_dev": max_dev, "or_lt1_mass": or_lt1,
                    "mass": tot})
    return out


def orbit_chain_H(n: int, p: float, d: float) -> list[float]:
    """H(A_i | A_{<i}) for the crash mixture, per coordinate, exactly."""
    x = 1.0 - p
    out = []
    # crash atoms: (prefix pattern, has coordinate i for i in R, mass)
    # S1 = (0,1,tail empty) m=.33; S2 = (0,0,tail full) m=.33;
    # S3 = (1,1,tail full) m=.04; S4 = (1,0,tail empty) m=.30
    for i in range(1, n + 1):
        if i == 1:
            p0 = (1.0 - d) * x + d * (0.33 + 0.33)  # S1,S2 lack coord 1
            out.append(h2(p0))
            continue
        if i == 2:
            H = 0.0
            # prefix a1 = 0 carries S1, S2; a1 = 1 carries S3, S4
            for a1, atoms in ((0, ((0.33, 1), (0.33, 0))),
                              (1, ((0.04, 1), (0.30, 0)))):
                pp = p if a1 else x
                tot = (1.0 - d) * pp + d * sum(m for m, _ in atoms)
                zer = (1.0 - d) * pp * x + d * sum(
                    m for m, has in atoms if not has)
                H += tot * h2(zer / tot)
            out.append(H)
            continue
        q = i - 3
        H = 0.0
        crash_pref = {(0, 1, 0): (0.33, 0), (0, 0, q): (0.33, 1),
                      (1, 1, q): (0.04, 1), (1, 0, 0): (0.30, 0)}
        for a1 in (0, 1):
            for a2 in (0, 1):
                for ja in range(q + 1):
                    w = a1 + a2 + ja
                    t0 = (1.0 - d) * p ** w * x ** (i - 1 - w)
                    cnt = math.comb(q, ja)
                    key = (a1, a2, ja)
                    if key in crash_pref:
                        m, has = crash_pref[key]
                        # matching classes are singletons (C(q,0)=C(q,q)=1)
                        tot = t0 + d * m
                        zer = t0 * x + (0.0 if has else d * m)
                        H += tot * h2(zer / tot)
                        cnt -= 1
                    if cnt > 0:
                        H += cnt * t0 * h2(x)
        out.append(H)
    return out


# ---------------------------------------------------------------- utilities

def save(name: str, obj) -> None:
    DATA.mkdir(exist_ok=True)
    path = DATA / name
    if path.exists():
        raise SystemExit(f"refusing to overwrite existing {path}")
    path.write_text(json.dumps(obj, indent=1, sort_keys=True, default=float)
                    + "\n", encoding="utf-8")


# ------------------------------------------------------------------ parts

def part_sk0() -> None:
    print("=" * 78)
    print("SK0. Engine validation (both engines)")
    print("=" * 78)
    out = {}
    # direct: product identity + closed-form gain
    n, p, lam = 5, 0.30, 1.0
    mu = prod_mu(n, p)
    rep = run_instance(mu, n, lam, p)
    x = 1.0 - p
    gpred = n * (h2(zplack(x, x, 2.0 ** lam)) - h2(p))
    print(f"[SK0a] product Bern({p})^{n}: max|log2OR-lam| = "
          f"{rep['max_abs_dev']:.2e}, gain err vs closed form "
          f"{abs(rep['gain']-gpred):.1e}, tax {rep['max_tax']:.1e}, "
          f"sym_dev {rep['sym_dev']:.1e}")
    out["product"] = {"max_dev": rep["max_abs_dev"],
                      "gain_err": abs(rep["gain"] - gpred)}
    # direct: crash identity
    crash_ok = []
    for n2 in (5, 8):
        atoms, r, c, K, resid, sym = sinkhorn(crash_nu(n2), 1.0, tol=1e-14)
        st = census_direct(atoms, r, c, K, n2, 1.0, 0.7)
        # designated history i=2, a=0, b=1 -> recover from min_l2 at i=2
        l2 = st[1]["min_l2"]
        crash_ok.append({"n": n2, "log2_or": l2, "err": abs(l2 - (3 - n2))})
        print(f"[SK0b] crash n={n2}: min log2OR at i=2 = {l2:+.6f} "
              f"(predict {3-n2}), err {abs(l2-(3-n2)):.1e}")
    out["crash"] = crash_ok
    # orbit vs direct cross-check at n = 6, 7 (crash mixture)
    xrows = []
    for n2 in (6, 7):
        p2, lam2, d2 = 0.30, 1.0, 0.05
        mu2 = mix(prod_mu(n2, p2), crash_nu(n2), d2, n2)
        atoms, r, c, K, resid, sym = sinkhorn(mu2, lam2, tol=1e-14)
        stD = census_direct(atoms, r, c, K, n2, lam2, 1.0 - p2)
        HaD = chain_H_direct(mu2, n2)
        R, C, oresid = orbit_sinkhorn(n2, p2, d2, lam2)
        stO = orbit_census(n2, p2, d2, lam2, R, C, 1.0 - p2)
        HaO = orbit_chain_H(n2, p2, d2)
        worst = 0.0
        for k in range(n2):
            worst = max(worst,
                        abs(stD[k]["down"] - stO[k]["down"]),
                        abs(stD[k]["max_abs_dev"] - stO[k]["max_abs_dev"]),
                        abs(stD[k]["E_hx"] - stO[k]["E_hx"]),
                        abs(stD[k]["E_xi2"] - stO[k]["E_xi2"]),
                        abs(HaD[k] - HaO[k]),
                        abs(stO[k]["mass"] - 1.0))
        xrows.append({"n": n2, "max_abs_disagreement": worst,
                      "direct_resid": resid, "orbit_resid": oresid})
        print(f"[SK0c] orbit vs direct, n={n2} crash d=0.05: max "
              f"disagreement over all per-i stats = {worst:.2e}")
    out["orbit_vs_direct"] = xrows
    save("pertsk_SK0.json", out)


def part_sk1() -> None:
    print("=" * 78)
    print("SK1. Calibration claims (rho*, margins, AHS coefficient)")
    print("=" * 78)
    rs383 = rho_star(0.383)
    rs382 = rho_star(0.38271)
    g1031 = margin_g(1.031, 0.383)
    g106 = margin_g(1.06, 0.383)
    g12 = margin_g(1.20, 0.383)
    print(f"[SK1] rho*(0.383)   = {rs383:.6f}   (008 claims 1.042205)")
    print(f"[SK1] rho*(0.38271) = {rs382:.6f}   (008 claims 1.030222; "
          f"queue's 'rho ~ 1.03')")
    print(f"[SK1] g(1.031, 0.383) = {g1031:+.6f}  (008 claims -0.0004)")
    print(f"[SK1] g(1.06, 0.383) = {g106:+.6f}, g(1.20, 0.383) = {g12:+.6f}")
    # AHS coefficient.  Hand derivation: z(1-x-y+z) = rho(x-z)(y-z);
    # implicit diff at rho=1, z=xy gives zdot * 1 = x y (1-x)(1-y).  At
    # p = psi: x = phi, z = phi^2 = psi, dg/drho = h'(psi) phi^2 (1-phi)^2
    # = h'(psi) phi^6 since 1-phi = phi^2.
    hp = math.log2((1.0 - PSI) / PSI)
    closed = hp * PHI ** 6
    e = 1e-7
    num = (margin_g(1.0 + e, PSI) - margin_g(1.0 - e, PSI)) / (2 * e)
    print(f"[SK1] AHS dg/drho|_1: closed h'(psi) phi^6 = {closed:.7f}, "
          f"central numeric = {num:.7f}, |err| {abs(closed-num):.1e}")
    print(f"[SK1] identity check phi^2 - psi = {PHI**2 - PSI:.1e}")
    save("pertsk_SK1.json", {
        "rho_star_0383": rs383, "rho_star_038271": rs382,
        "g_1031_0383": g1031, "g_106": g106, "g_12": g12,
        "ahs_closed": closed, "ahs_numeric": num})


def part_sk2() -> None:
    print("=" * 78)
    print("SK2. Theorem P6' bound + Lemma B/C quantitative checks")
    print("=" * 78)
    out = {"bound_rows": []}
    grid = []
    for (n, p, lam) in ((4, 0.30, 1.0), (5, 0.30, 1.0), (5, 0.30, 2.0),
                        (5, 0.383, math.log2(1.2)), (4, 0.20, 1.5),
                        (6, 0.30, 1.0)):
        dirs = {"crash": crash_nu(n), "delta_empty": {0: 1.0},
                "delta_full": {(1 << n) - 1: 1.0},
                "bern005": prod_mu(n, 0.05),
                "halfempty": halfempty_nu(n, p)}
        for dname, nu in dirs.items():
            for d in (0.05, 0.1, 0.25, 0.5):
                grid.append((n, p, lam, dname, nu, d))
    worst_ratio = 0.0
    worst_at = None
    viol = 0
    for (n, p, lam, dname, nu, d) in grid:
        mu = mix(prod_mu(n, p), nu, d, n)
        atoms, r, c, K, resid, sym = sinkhorn(mu, lam)
        st = census_direct(atoms, r, c, K, n, lam, 1.0 - p)
        dev = max(s["max_abs_dev"] for s in st)
        bound = 924.0 * min(p, 1.0 - p) ** (-3 * n) * d * d
        ratio = dev / bound
        if ratio > worst_ratio:
            worst_ratio, worst_at = ratio, (n, p, lam, dname, d)
        if dev > bound:
            viol += 1
        out["bound_rows"].append({"n": n, "p": p, "lam": lam, "dir": dname,
                                  "delta": d, "dev": dev, "bound": bound,
                                  "ratio": ratio})
    print(f"[SK2a] P6' bound test: {len(grid)} instances (incl delta=0.5), "
          f"violations = {viol}; loosest gap: max dev/bound = "
          f"{worst_ratio:.2e} at {worst_at}")
    # Lemma B first-order formula at an interior delta (008 only checked
    # delta = 0), plus Lemma C sup-bounds.
    n, p, lam = 5, 0.30, 1.0
    m0 = min(p, 1 - p) ** n
    mu0 = prod_mu(n, p)
    nu = crash_nu(n)
    d0 = 0.25
    e = 0.005

    def vpot(d):
        mu = mix(mu0, nu, d, n)
        atoms, r, c, K, resid, sym = sinkhorn(mu, lam, tol=1e-14)
        return atoms, [0.5 * math.log(r[i] * c[i]) for i in range(len(atoms))]

    atoms, v0 = vpot(d0)
    _, vp = vpot(d0 + e)
    _, vm = vpot(d0 - e)
    vdot_fd = [(vp[i] - vm[i]) / (2 * e) for i in range(len(v0))]
    vddot_fd = [(vp[i] - 2 * v0[i] + vm[i]) / e ** 2 for i in range(len(v0))]
    # Lemma B formula at d0: (D_mu + Pi) vdot = nu - mu0
    mu = mix(mu0, nu, d0, n)
    atomsS, r, c, K, resid, sym = sinkhorn(mu, lam, tol=1e-14)
    N = len(atomsS)
    H = [[r[i] * c[j] * K[i][j] + (mu[atomsS[i]] if i == j else 0.0)
          for j in range(N)] for i in range(N)]
    rhs = [nu.get(atomsS[i], 0.0) - mu0[atomsS[i]] for i in range(N)]
    # Gaussian elimination
    Aug = [row[:] + [rhs[i]] for i, row in enumerate(H)]
    for col in range(N):
        piv = max(range(col, N), key=lambda rr: abs(Aug[rr][col]))
        Aug[col], Aug[piv] = Aug[piv], Aug[col]
        for rr in range(col + 1, N):
            fct = Aug[rr][col] / Aug[col][col]
            if fct:
                for cc in range(col, N + 1):
                    Aug[rr][cc] -= fct * Aug[col][cc]
    vdot_ift = [0.0] * N
    for rr in range(N - 1, -1, -1):
        s = Aug[rr][N] - sum(Aug[rr][cc] * vdot_ift[cc]
                             for cc in range(rr + 1, N))
        vdot_ift[rr] = s / Aug[rr][rr]
    err_ift = max(abs(vdot_fd[i] - vdot_ift[i]) for i in range(N))
    nvd = max(abs(x) for x in vdot_fd)
    nvdd = max(abs(x) for x in vddot_fd)
    print(f"[SK2b] Lemma B at interior delta={d0}: |FD vdot - IFT vdot|_inf "
          f"= {err_ift:.2e} (FD step {e})")
    print(f"[SK2b] Lemma C: |vdot|_inf = {nvd:.3f} vs bound 4/m0 = "
          f"{4/m0:.1f}; |vddot|_inf = {nvdd:.3f} vs bound 128/m0^3 = "
          f"{128/m0**3:.2e}")
    ok = nvd <= 4 / m0 and nvdd <= 128 / m0 ** 3
    print(f"[SK2b] Lemma C sup-bounds hold on this instance: {ok}")
    out["lemma_bc"] = {"n": n, "p": p, "lam": lam, "delta": d0,
                       "err_ift_vs_fd": err_ift, "vdot_inf": nvd,
                       "vddot_inf": nvdd, "bound_vdot": 4 / m0,
                       "bound_vddot": 128 / m0 ** 3, "ok": ok}
    save("pertsk_SK2.json", out)


def part_sk3() -> None:
    print("=" * 78)
    print("SK3. Spot-verify 008's assembly rows + hunt capped negatives")
    print("=" * 78)
    p = 0.383
    out = {"spot": [], "hunt": []}
    # spot checks against pert_full_run.log values
    spots = [
        ("crash", crash_nu, 6, math.log2(1.2), 0.10,
         {"min_net": 0.007407, "gain_lb": 0.05927, "gain": 0.09252}),
        ("pm03", lambda n: pm03_nu(n, p), 7, math.log2(1.2), 0.20,
         {"min_net": 0.005480, "gain_lb": 0.05060, "gain": 0.06340}),
        ("slice", lambda n: slice_nu(n, p), 7, math.log2(1.2), 0.20,
         {"min_net": -0.008118, "gain_lb": -0.04357, "gain": -0.00411}),
        ("crash", crash_nu, 6, math.log2(1.06), 0.10,
         {"min_net": 0.002865, "gain_lb": 0.03440, "gain": 0.06711}),
    ]
    for dname, mk, n, lam, d, ref in spots:
        mu = mix(prod_mu(n, p), mk(n), d, n)
        rep = run_instance(mu, n, lam, p)
        diffs = {k: abs(rep[k] - ref[k]) for k in ref}
        print(f"[SK3a] {dname} n={n} rho={2**lam:.3f} d={d}: "
              f"min_net {rep['min_net']:+.6f} gain_lb {rep['gain_lb']:+.5f} "
              f"gain {rep['gain']:+.5f}; max diff vs 008 log "
              f"{max(diffs.values()):.1e}")
        out["spot"].append({"dir": dname, "n": n, "rho": 2 ** lam,
                            "delta": d, "mine": {k: rep[k] for k in ref},
                            "ref_008": ref, "max_diff": max(diffs.values())})
    # the hunt
    t0 = time.time()
    hunts = []
    for rho in (1.06, 1.20, 1.045):
        lam = math.log2(rho)
        for n in (6, 7, 8, 9):
            dirs = {"pm03": pm03_nu(n, p), "halfempty": halfempty_nu(n, p),
                    "slice": slice_nu(n, p), "emb0_crash_k2": emb0k2_nu(n),
                    "crash": crash_nu(n)}
            for dname, nu in dirs.items():
                for d in (0.05, 0.10, 0.15, 0.20):
                    if n == 9 and (d not in (0.10, 0.20)
                                   or dname in ("slice", "crash")):
                        continue
                    if rho == 1.045 and (n > 7 or d < 0.10
                                         or dname == "emb0_crash_k2"):
                        continue
                    if rho == 1.20 and dname not in ("halfempty",):
                        continue  # 1.20 cells other than halfempty were run
                    hunts.append((rho, lam, n, dname, nu, d))
    neg = []
    for rho, lam, n, dname, nu, d in hunts:
        mu = mix(prod_mu(n, p), nu, d, n)
        capped = max(elem_marg(mu, n)) <= p + 1e-12
        rep = run_instance(mu, n, lam, p)
        row = {"rho": rho, "n": n, "dir": dname, "delta": d,
               "capped": capped, "min_net": rep["min_net"],
               "gain_lb": rep["gain_lb"], "gain": rep["gain"],
               "or_lt1_mass": rep["or_lt1_mass"], "resid": rep["resid"]}
        out["hunt"].append(row)
        flag = ""
        if capped and (rep["gain_lb"] <= 0 or rep["min_net"] <= 0):
            flag = "  <-- CAPPED NEGATIVE"
            neg.append(row)
        print(f"[SK3b] rho={rho:5.3f} n={n} nu={dname:13s} d={d:4.2f} "
              f"capped={str(capped):5s}: min_net {rep['min_net']:+.6f} "
              f"gain_lb {rep['gain_lb']:+.5f} gain {rep['gain']:+.5f} "
              f"OR<1 mass {rep['or_lt1_mass']:.1e}{flag}")
    print(f"[SK3b] {len(hunts)} hunt runs in {time.time()-t0:.0f}s; "
          f"capped instances with gain_lb <= 0 or min_net <= 0: {len(neg)}")
    out["n_capped_negative"] = len(neg)
    save("pertsk_SK3.json", out)


def part_sk4() -> None:
    print("=" * 78)
    print("SK4. Conditional-theorem constants (tau_half, M_tau, ||G||, d0)")
    print("=" * 78)
    p, rho, s0 = 0.383, 1.20, 0.10
    lam = math.log2(rho)
    x0 = 1.0 - p

    def G(x, y, tau):
        return h2(zplack(x, y, 2.0 ** (lam + tau))) - 0.5 * h2(x) \
            - 0.5 * h2(y)

    def bisect_tau_half(x, y):
        lo, hi = 0.0, 40.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if zplack(x, y, 2.0 ** (lam + mid)) < 0.5:
                lo = mid
            else:
                hi = mid
        return lo

    g = G(x0, x0, 0.0)
    th_center = bisect_tau_half(x0, x0)
    th_corner = bisect_tau_half(x0 + s0, x0 + s0)
    z_corner_at_center_tau = zplack(x0 + s0, x0 + s0, 2.0 ** (lam + th_center))
    print(f"[SK4a] g = {g:.6f}; tau_half(center) = {th_center:.4f} "
          f"(008: 3.017); tau_half(corner x0+s0) = {th_corner:.4f}")
    print(f"[SK4a] z at corner at tau_half(center): "
          f"{z_corner_at_center_tau:.4f}  ({'>' if z_corner_at_center_tau > 0.5 else '<='} 1/2 "
          f"-> box-uniform claim {'FALSE' if z_corner_at_center_tau > 0.5 else 'holds'})")
    # violation region: fraction of the 25x25 box grid with z > 1/2 at
    # tau just below tau_half(center); and worst negative dG/dtau there
    e = 1e-5
    grid = [x0 - s0 + 2 * s0 * k / 24 for k in range(25)]
    tprobe = th_center * 0.999
    nviol = 0
    worst_neg_slope = 0.0
    for x in grid:
        for y in grid:
            if zplack(x, y, 2.0 ** (lam + tprobe)) > 0.5:
                nviol += 1
                sl = (G(x, y, tprobe + e) - G(x, y, tprobe - e)) / (2 * e)
                worst_neg_slope = min(worst_neg_slope, sl)
    print(f"[SK4a] grid points (of 625) with z > 1/2 at tau = "
          f"0.999*tau_half(center): {nviol}; most negative dG/dtau there: "
          f"{worst_neg_slope:+.5f}")
    # constants over box x [0, tau_half(center)] (008's convention)
    taus = [th_center * k / 12 for k in range(13)]
    M_tau = M2 = L2 = 0.0
    dip = 0.0
    Gbox = 0.0
    for x in grid:
        for y in grid:
            for t in taus:
                gv = G(x, y, t)
                Gbox = max(Gbox, abs(gv))
                M_tau = max(M_tau, abs((G(x, y, t + e) - G(x, y, t - e))
                                       / (2 * e)))
                L2 = max(L2, abs((G(x + e, y, t) - G(x - e, y, t)) / (2 * e)))
            gxx = (G(x + e, y, 0) - 2 * G(x, y, 0) + G(x - e, y, 0)) / e ** 2
            gyy = (G(x, y + e, 0) - 2 * G(x, y, 0) + G(x, y - e, 0)) / e ** 2
            gxy = (G(x + e, y + e, 0) - G(x + e, y - e, 0)
                   - G(x - e, y + e, 0) + G(x - e, y - e, 0)) / (4 * e ** 2)
            tr, det = gxx + gyy, gxx * gyy - gxy * gxy
            lmin = 0.5 * (tr - math.sqrt(max(tr * tr - 4 * det, 0.0)))
            M2 = max(M2, max(0.0, -lmin))
            for t in (5.0, 10.0, 30.0, 60.0):
                dip = max(dip, -G(x, y, t))
    # sup |dG/dtau| for tau < 0 (where (B1) losses actually live)
    M_tau_neg = 0.0
    arg_neg = None
    for x in grid:
        for y in grid:
            for t in [-6.0 + 6.0 * k / 40 for k in range(41)]:
                sl = abs((G(x, y, t + e) - G(x, y, t - e)) / (2 * e))
                if sl > M_tau_neg:
                    M_tau_neg, arg_neg = sl, (x, y, t)
    print(f"[SK4b] M_tau on [0,tau_half] = {M_tau:.4f} (008: 0.0592); "
          f"sup|dG/dtau| on tau in [-6,0] = {M_tau_neg:.4f} at {arg_neg}")
    print(f"[SK4b] M2 = {M2:.4f} (008: 2.2518); L2 = {L2:.4f} (008: "
          f"1.0220); dip = {dip:.6f} (008: 0.007112); ||G||_box = {Gbox:.4f}"
          f" (008: not computed)")
    # global negative part of G over (0,1)^2 x [0, tau_half] (off-box charge)
    Gneg_glob = 0.0
    fine = [0.02 + 0.96 * k / 96 for k in range(97)]
    for x in fine:
        for y in fine:
            for t in (0.0, th_center / 2, th_center):
                Gneg_glob = max(Gneg_glob, -G(x, y, t))
    print(f"[SK4b] sup of (-G) over (0,1)^2 x [0,tau_half] = {Gneg_glob:.4f}"
          f"  (the constant the off-box Markov charge actually needs; "
          f"008 used ||G||_box)")
    c1, c2, c3 = 0.07, 0.01, 2.5
    d0_code = math.sqrt(0.5 * g / (M_tau * c1 + c2 + M2 * c3))
    d0_stmt = math.sqrt(g / (2 * (M_tau * c1 + c2
                                  + (M2 + Gbox / s0 ** 2) * c3)))
    d0_corr = math.sqrt(g / (2 * (M_tau_neg * c1 + c2
                                  + (M2 + Gneg_glob / s0 ** 2) * c3)))
    print(f"[SK4c] delta_0: as coded by 008 (drops ||G||_box term) = "
          f"{d0_code:.4f}; from 008's stated formula with ||G||_box = "
          f"{d0_stmt:.4f}; with corrected constants (M_tau over tau<0, "
          f"global G) = {d0_corr:.4f}")
    save("pertsk_SK4.json", {
        "g": g, "tau_half_center": th_center, "tau_half_corner": th_corner,
        "z_corner_at_center_tau": z_corner_at_center_tau,
        "n_box_grid_viol": nviol, "worst_neg_dGdtau_at_tau_half":
        worst_neg_slope, "M_tau": M_tau, "M_tau_neg_side": M_tau_neg,
        "M2": M2, "L2": L2, "dip": dip, "G_box": Gbox,
        "G_neg_global": Gneg_glob, "delta0_as_coded": d0_code,
        "delta0_stated_formula": d0_stmt, "delta0_corrected": d0_corr})


def part_sk5() -> None:
    print("=" * 78)
    print("SK5. Downward budget census extended to n = 16 (008: n <= 9)")
    print("=" * 78)
    out = {}
    d = 0.05
    for tag, p, lam in (("p030_lam1", 0.30, 1.0),
                        ("p0383_rho12", 0.383, math.log2(1.2))):
        rows = []
        t0 = time.time()
        for n in range(4, 17):
            R, C, resid = orbit_sinkhorn(n, p, d, lam)
            st = orbit_census(n, p, d, lam, R, C, 1.0 - p)
            Ha = orbit_chain_H(n, p, d)
            downs = [s["down"] for s in st]
            taxes = [Ha[k] - st[k]["E_hx"] for k in range(n)]
            xi2 = [s["E_xi2"] for s in st]
            mass_err = max(abs(s["mass"] - 1.0) for s in st)
            row = {"n": n, "delta": d, "resid": resid, "mass_err": mass_err,
                   "down_avg_over_d2": sum(downs) / n / d ** 2,
                   "down_max_over_d2": max(downs) / d ** 2,
                   "worst_i": downs.index(max(downs)) + 1,
                   "avg_tax_over_d2": sum(taxes) / n / d ** 2,
                   "max_tax": max(taxes),
                   "max_E_xi2_over_d2": max(xi2) / d ** 2,
                   "max_abs_dev": max(s["max_abs_dev"] for s in st),
                   "or_lt1_mass": sum(s["or_lt1_mass"] for s in st)}
            rows.append(row)
            print(f"[SK5] {tag} n={n:2d}: avg down/d^2 = "
                  f"{row['down_avg_over_d2']:.4f}  max_i/d^2 = "
                  f"{row['down_max_over_d2']:.4f} (i={row['worst_i']})  "
                  f"avg tax/d^2 = {row['avg_tax_over_d2']:.5f}  "
                  f"max E_xi2/d^2 = {row['max_E_xi2_over_d2']:.3f}  "
                  f"({time.time()-t0:.0f}s)")
        # increments of the average
        avgs = [r["down_avg_over_d2"] for r in rows]
        incs = [avgs[k + 1] - avgs[k] for k in range(len(avgs) - 1)]
        ratios = [incs[k + 1] / incs[k] for k in range(len(incs) - 1)]
        print(f"[SK5] {tag} avg-budget increments: "
              + " ".join(f"{v:.5f}" for v in incs))
        print(f"[SK5] {tag} increment ratios: "
              + " ".join(f"{v:.3f}" for v in ratios))
        out[tag] = {"rows": rows, "increments": incs,
                    "increment_ratios": ratios}
    # delta-stability spot check at n = 12
    for tag, p, lam in (("p0383_rho12", 0.383, math.log2(1.2)),):
        for d2 in (0.025,):
            R, C, resid = orbit_sinkhorn(12, p, d2, lam)
            st = orbit_census(12, p, d2, lam, R, C, 1.0 - p)
            v = sum(s["down"] for s in st) / 12 / d2 ** 2
            print(f"[SK5] {tag} n=12 delta={d2}: avg down/d^2 = {v:.4f} "
                  f"(delta^2-scaling check vs delta=0.05 row)")
            out.setdefault("delta_check", []).append(
                {"tag": tag, "n": 12, "delta": d2, "down_avg_over_d2": v})
    save("pertsk_SK5.json", out)


def part_sk5b() -> None:
    print("=" * 78)
    print("SK5b. Census pushed to n = 28 (shape of the average budgets)")
    print("=" * 78)
    out = {}
    d = 0.05
    for tag, p, lam in (("p030_lam1", 0.30, 1.0),
                        ("p0383_rho12", 0.383, math.log2(1.2))):
        rows = []
        t0 = time.time()
        for n in range(4, 29):
            R, C, resid = orbit_sinkhorn(n, p, d, lam)
            st = orbit_census(n, p, d, lam, R, C, 1.0 - p)
            Ha = orbit_chain_H(n, p, d)
            downs = [s["down"] for s in st]
            taxes = [Ha[k] - st[k]["E_hx"] for k in range(n)]
            rows.append({"n": n, "resid": resid,
                         "down_avg_over_d2": sum(downs) / n / d ** 2,
                         "down_max_over_d2": max(downs) / d ** 2,
                         "avg_tax_over_d2": sum(taxes) / n / d ** 2,
                         "max_E_xi2_over_d2":
                         max(s["E_xi2"] for s in st) / d ** 2,
                         "max_abs_dev": max(s["max_abs_dev"] for s in st),
                         "or_lt1_mass":
                         sum(s["or_lt1_mass"] for s in st)})
        avgs = [r["down_avg_over_d2"] for r in rows]
        incs = [avgs[k + 1] - avgs[k] for k in range(len(avgs) - 1)]
        ratios = [incs[k + 1] / incs[k] for k in range(len(incs) - 1)]
        # local power-law exponent of the increments: a_n from
        # d_{n+1}/d_n = ((n)/(n+1))^a  -> a > 1 needed for a bounded limit
        expo = []
        for k in range(len(incs) - 1):
            nn = rows[k + 1]["n"]
            expo.append(math.log(ratios[k]) / math.log(nn / (nn + 1.0))
                        if ratios[k] > 0 else float("nan"))
        print(f"[SK5b] {tag} ({time.time()-t0:.0f}s): avg down/d^2 at "
              f"n=9/16/22/28: {avgs[5]:.4f} {avgs[12]:.4f} {avgs[18]:.4f} "
              f"{avgs[24]:.4f}")
        print(f"[SK5b] {tag} increment ratios (n=17..28): "
              + " ".join(f"{v:.3f}" for v in ratios[12:]))
        print(f"[SK5b] {tag} local power-law exponent a_n (bounded limit "
              f"needs a > 1): " + " ".join(f"{v:.2f}" for v in expo[8:]))
        tx = [r["avg_tax_over_d2"] for r in rows]
        tinc = [tx[k + 1] - tx[k] for k in range(len(tx) - 1)]
        print(f"[SK5b] {tag} avg tax/d^2 at n=9/16/22/28: {tx[5]:.4f} "
              f"{tx[12]:.4f} {tx[18]:.4f} {tx[24]:.4f}; last 6 increments: "
              + " ".join(f"{v:.5f}" for v in tinc[-6:]))
        devs = [r["max_abs_dev"] for r in rows]
        print(f"[SK5b] {tag} max|log2OR-lam| at n=9/16/22/28: {devs[5]:.3f} "
              f"{devs[12]:.3f} {devs[18]:.3f} {devs[24]:.3f}; "
              f"OR<1 mass max: {max(r['or_lt1_mass'] for r in rows):.1e}")
        out[tag] = {"rows": rows, "increments": incs,
                    "increment_ratios": ratios, "local_exponent": expo}
    save("pertsk_SK5b.json", out)


def part_sk5c() -> None:
    """delta -> 0 scaling of the average budgets at large n: is the SK5b
    shape (turnover at p030, acceleration at p0383) a fixed-delta artifact
    of the shrinking asymptotic window (006 S7), or real in the delta^2
    coefficient?"""
    print("=" * 78)
    print("SK5c. delta-scaling of avg budgets at large n")
    print("=" * 78)
    out = {}
    for tag, p, lam in (("p030_lam1", 0.30, 1.0),
                        ("p0383_rho12", 0.383, math.log2(1.2))):
        rows = []
        for n in (20, 24, 28, 32):
            for d in (0.05, 0.025, 0.0125):
                R, C, resid = orbit_sinkhorn(n, p, d, lam)
                st = orbit_census(n, p, d, lam, R, C, 1.0 - p)
                Ha = orbit_chain_H(n, p, d)
                downs = [s["down"] for s in st]
                taxes = [Ha[k] - st[k]["E_hx"] for k in range(n)]
                rows.append({"n": n, "delta": d, "resid": resid,
                             "down_avg_over_d2": sum(downs) / n / d ** 2,
                             "avg_tax_over_d2": sum(taxes) / n / d ** 2})
            trip = [r for r in rows if r["n"] == n]
            print(f"[SK5c] {tag} n={n}: avg down/d^2 at d=0.05/0.025/0.0125"
                  f" = " + " ".join(f"{r['down_avg_over_d2']:.4f}"
                                    for r in trip)
                  + "   avg tax/d^2 = "
                  + " ".join(f"{r['avg_tax_over_d2']:.5f}" for r in trip))
        out[tag] = rows
    save("pertsk_SK5c.json", out)


def part_sk7() -> None:
    """Residual checks: (i) delta -> 0 stability of the n = 32 assembly-regime
    budget coefficients (are we in the quadratic window?); (ii) the analytic
    s0 ceiling for a genuinely box-uniform tau_half at (p, rho) =
    (0.383, 1.2): tau_half(corner) > 0 iff z_rho(x0+s0, x0+s0) < 1/2."""
    print("=" * 78)
    print("SK7. Residual checks (delta window at n=32; s0 ceiling)")
    print("=" * 78)
    out = {}
    p, lam = 0.383, math.log2(1.2)
    rows = []
    for d in (0.00625, 0.003125):
        R, C, resid = orbit_sinkhorn(32, p, d, lam)
        st = orbit_census(32, p, d, lam, R, C, 1.0 - p)
        Ha = orbit_chain_H(32, p, d)
        downs = [s["down"] for s in st]
        taxes = [Ha[k] - st[k]["E_hx"] for k in range(32)]
        rows.append({"n": 32, "delta": d, "resid": resid,
                     "down_avg_over_d2": sum(downs) / 32 / d ** 2,
                     "avg_tax_over_d2": sum(taxes) / 32 / d ** 2})
        print(f"[SK7a] p0383_rho12 n=32 delta={d}: avg down/d^2 = "
              f"{rows[-1]['down_avg_over_d2']:.4f}, avg tax/d^2 = "
              f"{rows[-1]['avg_tax_over_d2']:.5f}")
    out["n32_small_delta"] = rows
    # s0 ceiling
    x0 = 1.0 - p
    lo, hi = x0, 1.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if zplack(mid, mid, 1.2) < 0.5:
            lo = mid
        else:
            hi = mid
    xc = lo
    print(f"[SK7b] z_1.2(x, x) crosses 1/2 at x_c = {xc:.6f}; corner check "
          f"(x0+s0)^2 at 008's s0=0.10: {(x0+0.10)**2:.4f} (> 1/2, and "
          f"rho > 1 makes z > xy, so tau_half(corner) = 0 analytically); "
          f"max box-uniform s0 at (0.383, 1.2) = x_c - x0 = {xc-x0:.4f}")
    out["s0_ceiling"] = {"x_c": xc, "s0_max": xc - x0,
                         "corner_xy_at_s0_010": (x0 + 0.10) ** 2}
    save("pertsk_SK7.json", out)


def part_sk6() -> None:
    print("=" * 78)
    print("SK6. Centered kernel: crash/mirror closed forms + product band")
    print("=" * 78)
    lam = 1.0
    out = {}

    def cexp(n):
        return lambda a, b: lam * (bin(a & b).count("1")
                                   - bin(a).count("1") * bin(b).count("1") / n)

    rows = []
    for n in (5, 8):
        # potential-free identity (no Sinkhorn): OR over the 4 crash cells
        full = (1 << n) - 1
        S1, S2, S3, S4 = 0b10, full ^ 0b11, full, 0b01
        e = cexp(n)
        ident = e(S1, S3) + e(S2, S4) - e(S1, S4) - e(S2, S3)
        # Sinkhorn fit
        atoms, r, c, K, resid, sym = sinkhorn(crash_nu(n), lam,
                                              expo=cexp(n), tol=1e-14)
        st = census_direct(atoms, r, c, K, n, lam, 1.0 - 0.383)
        l2 = st[1]["min_l2"]
        pred = lam * (3 - n) / n
        rows.append({"n": n, "identity": ident, "sinkhorn": l2,
                     "closed_form": pred})
        print(f"[SK6a] centered crash n={n}: identity {ident:+.6f}, "
              f"Sinkhorn {l2:+.6f}, closed form lam(3-n)/n = {pred:+.6f}")
    out["crash"] = rows
    # mirror
    n = 5
    full = (1 << n) - 1
    M1, M2_, M3, M4 = 0b10 | (full ^ 0b11), 0, full, 0b01
    e = cexp(n)
    ident = e(M1, M3) + e(M2_, M4) - e(M1, M4) - e(M2_, M3)
    print(f"[SK6a] centered mirror n=5: identity {ident:+.6f} vs closed "
          f"form lam(n-1)/n = {lam*(n-1)/n:+.6f}")
    out["mirror"] = {"identity": ident, "closed": lam * (n - 1) / n}
    # product band
    prows = []
    for n in (6, 8):
        p = 0.383
        mu = prod_mu(n, p)
        atoms, r, c, K, resid, sym = sinkhorn(mu, lam, expo=cexp(n))
        st = census_direct(atoms, r, c, K, n, lam, 1.0 - p)
        mn = min(s["min_l2"] for s in st)
        mx = max(s["max_l2"] for s in st)
        gain = union_gain_direct(atoms, r, c, K, mu)
        prows.append({"n": n, "min_l2": mn, "max_l2": mx, "gain": gain})
        print(f"[SK6b] centered product n={n}: log2 OR in [{mn:+.4f}, "
              f"{mx:+.4f}] (008: n=6 [0.8265,0.8333], n=8 [0.8695,0.8750]); "
              f"gain {gain:+.5f}")
    out["product"] = prows
    save("pertsk_SK6.json", out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="all")
    args = ap.parse_args()
    todo = args.part.upper()
    t0 = time.time()
    for name, fn in (("SK0", part_sk0), ("SK1", part_sk1), ("SK2", part_sk2),
                     ("SK3", part_sk3), ("SK4", part_sk4), ("SK5", part_sk5),
                     ("SK5B", part_sk5b), ("SK5C", part_sk5c), ("SK6", part_sk6), ("SK7", part_sk7)):
        if todo in (name, "ALL"):
            fn()
    print(f"done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
