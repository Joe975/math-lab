#!/usr/bin/env python3
"""Attempt 008: perturbative assembly around product measures (Gap 3(c)).

Context (informed; 003/004/005/006): the dependent-couplings route needs, near
product measures mu_0 = Bern(p)^n and tilt rho = 2^lambda just above
rho*(p) = p(3p-1)/(1-2p)^2, a second-order assembly: per-coordinate margin
g(rho,p) = h(z_rho(x,x)) - h(p) > 0 (x = 1-p) vs losses that are O(delta^2)
by 005 Prop 6 (odds ratios) and O(delta^2) tax.  This tool measures every
term of that assembly exactly on small n and probes whether the budget
constants are uniform in n.

Parts:
  P0  engine validation: product identity OR == 2^lambda and closed-form gain;
      reproduction of 005 part G (checkpoint diff); crash-family identity.
  P1  the smoothness step of 005 Prop 6 (closed in the 008 record by a
      quantitative IFT): numerical check of the first-order formula
      v'(0) = (D_mu + Pi)^{-1} mu', quadratic remainder ratios.
  P2  the assembly, end to end at small n: margins g(rho,p) (closed form,
      incl. the AHS-point coefficient h'(psi) phi^6), then for
      mu_delta = (1-delta) Bern(p)^n + delta nu (marginal-capped nu) the exact
      per-coordinate net  N_i = E[h(z_i)] - H(A_i|A_{<i})  whose sum is the
      chain-rule gain lower bound; taxes; OR deviations.  This is 005 lead 3's
      falsifiable check run for real instead of by expansion.
  P3  budget growth in n: for adversarial directions nu (crash, embedded
      crash, point mass at empty set, Bern(0.7)^n) measure
      max_hist |log2 OR - lambda|, its mass-weighted mean, per-coordinate tax
      and E[(xtilde-x0)^2] at fixed deltas across n; extract c(n) = dev/delta^2.
  P4  005 lead 4's one-line kernel change: two-scale kernel
      2^{lam1 |A&B| + lam2 |A||B|/n} (lam2 = -lam1: "centered overlap") on the
      crash family; closed-form crash OR check; behavior on products.

Standard library only.  Deterministic.  Runtime ~6-10 min for --part all.
Usage:
    python problems/union-closed/explore/uc_pert.py [--part P0|P1|P2|P3|P4|all]
Checkpoints: problems/union-closed/data/pert_P[0-4]*.json
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

LOG2 = math.log(2.0)
PSI = (3.0 - math.sqrt(5.0)) / 2.0
PHI = 1.0 - PSI


def popcount(x: int) -> int:
    return bin(x).count("1")


def h(p: float) -> float:
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log2(p) - (1.0 - p) * math.log2(1.0 - p)


# ---------------------------------------------------------------- closed forms
# Plackett both-zero probability z_rho(x,y): root of
#   z(1-x-y+z) = rho (x-z)(y-z)   in [max(0,x+y-1), min(x,y)]
# (re-derived; same object as 003 part A).

def z_rho(x: float, y: float, rho: float) -> float:
    if x <= 0.0 or y <= 0.0:
        return 0.0
    if x >= 1.0:
        return y
    if y >= 1.0:
        return x
    if abs(rho - 1.0) < 1e-12:
        return x * y
    a = rho - 1.0
    S = 1.0 + a * (x + y)
    disc = max(S * S - 4.0 * a * rho * x * y, 0.0)  # exact disc >= 0; clamp fp noise
    z = (S - math.sqrt(disc)) / (2.0 * a)
    return min(max(z, max(0.0, x + y - 1.0)), min(x, y))


def rho_star(p: float) -> float:
    return p * (3.0 * p - 1.0) / (1.0 - 2.0 * p) ** 2


def margin(rho: float, p: float) -> float:
    """Per-coordinate gain of the Plackett(rho) coupling on Bern(p) vs its own
    entropy: g = h(z_rho(x,x)) - h(p), x = 1-p.  Positive iff rho > rho*(p)."""
    x = 1.0 - p
    return h(z_rho(x, x, rho)) - h(p)


def dz_drho_at1(x: float, y: float) -> float:
    """d z_rho / d rho at rho = 1 (closed form x y (1-x)(1-y))."""
    return x * y * (1.0 - x) * (1.0 - y)


# ---------------------------------------------------------------- measures

def product_measure(n: int, p: float) -> dict[int, float]:
    return {m: p ** popcount(m) * (1.0 - p) ** (n - popcount(m))
            for m in range(1 << n)}


def mix_measures(base: dict[int, float], nu: dict[int, float],
                 delta: float, n: int) -> dict[int, float]:
    return {m: (1.0 - delta) * base.get(m, 0.0) + delta * nu.get(m, 0.0)
            for m in range(1 << n)}


def crash_family(n: int) -> dict[int, float]:
    full = (1 << n) - 1
    return {0b10: 0.33, full ^ 0b11: 0.33, full: 0.04, 0b01: 0.30}


def embedded_crash(n: int, k: int) -> dict[int, float]:
    """Crash structure on tail coordinates k+1..n, all atoms AVOIDING the
    first k coordinates (zero-prefix embedding, so element marginals stay
    <= 0.37 and the direction is legal for the marginal-capped class).
    k = 0 reduces to crash_family(n).  The designated deep history is
    i = k+2, a = 0, b = 1<<k."""
    zeros = (1 << k) - 1  # coordinates 1..k, absent from every atom
    tail_full = ((1 << n) - 1) ^ zeros
    s1 = 1 << (k + 1)
    s4 = 1 << k
    s2 = tail_full ^ s1 ^ s4
    s3 = tail_full
    return {s1: 0.33, s2: 0.33, s3: 0.04, s4: 0.30}


def mirror_family(n: int) -> dict[int, float]:
    full = (1 << n) - 1
    return {0b10 | (full ^ 0b11): 0.3, 0b00: 0.2, full: 0.3, 0b01: 0.2}


def elem_marginals(mu: dict[int, float], n: int) -> list[float]:
    return [sum(w for m, w in mu.items() if m >> j & 1) for j in range(n)]


def entropy_bits(mu: dict[int, float]) -> float:
    return -sum(w * math.log2(w) for w in mu.values() if w > 0.0)


# ---------------------------------------------------------------- Sinkhorn

def fit_tilt(mu: dict[int, float], lam: float, kernel=None,
             tol: float = 1e-12, iters: int = 200000):
    """Alternating scaling for pi(A,B) = r(A) c(B) K(A,B), both marginals mu.
    K defaults to the overlap tilt 2^{lam |A&B|}.  Returns
    (atoms, pi_rows, r, c, resid) with pi_rows[i][j] the coupling matrix.
    Residual = max abs marginal error on both margins of the final coupling."""
    atoms = sorted(m for m, w in mu.items() if w > 0.0)
    N = len(atoms)
    if kernel is None:
        K = [[2.0 ** (lam * popcount(a & b)) for b in atoms] for a in atoms]
    else:
        K = [[kernel(a, b) for b in atoms] for a in atoms]
    tgt = [mu[a] for a in atoms]
    r = [1.0] * N
    c = [1.0] * N
    resid = float("inf")
    for it in range(iters):
        # column update (K symmetric: rows double as columns)
        c = [t / s for t, s in
             zip(tgt, (sum(ri * kij for ri, kij in zip(r, K[j]))
                       for j in range(N)))]
        # row update
        r = [t / s for t, s in
             zip(tgt, (sum(cj * kij for cj, kij in zip(c, K[i]))
                       for i in range(N)))]
        if it % 4 == 3 or it < 4:
            col_err = max(abs(c[j] * sum(ri * kij
                                         for ri, kij in zip(r, K[j])) - tgt[j])
                          for j in range(N))
            if col_err < tol:
                break
    row_err = max(abs(r[i] * sum(cj * kij for cj, kij in zip(c, K[i]))
                      - tgt[i]) for i in range(N))
    col_err = max(abs(c[j] * sum(ri * kij for ri, kij in zip(r, K[j]))
                      - tgt[j]) for j in range(N))
    resid = max(row_err, col_err)
    pi_rows = [[r[i] * c[j] * K[i][j] for j in range(N)] for i in range(N)]
    return atoms, pi_rows, r, c, resid


# ---------------------------------------------------------------- statistics

def chain_H(mu: dict[int, float], n: int) -> list[float]:
    """H(A_i | A_{<i}) per coordinate (bits), coordinate i = bit i-1."""
    out = []
    for i in range(1, n + 1):
        pref = (1 << (i - 1)) - 1
        bit = 1 << (i - 1)
        tot: dict[int, float] = {}
        zer: dict[int, float] = {}
        for m, w in mu.items():
            a = m & pref
            tot[a] = tot.get(a, 0.0) + w
            if not (m & bit):
                zer[a] = zer.get(a, 0.0) + w
        out.append(sum(w * h(zer.get(a, 0.0) / w) for a, w in tot.items()
                       if w > 0.0))
    return out


def coord_stats(atoms, pi_rows, n: int, lam: float, x0: float):
    """Exact per-coordinate assembly statistics of the coupling.

    Per i: E[h(z_i)], E[h(xt_i)], E[h(yt_i)] over histories (a,b);
    OR census stats over nondegenerate cells: max |log2 OR - lam|,
    mass-weighted mean of |log2 OR - lam|, mass with OR < 1, min/max log2 OR;
    E[(xt - x0)^2]; degenerate-cell mass."""
    N = len(atoms)
    stats = []
    for i in range(1, n + 1):
        pref = (1 << (i - 1)) - 1
        bit = 1 << (i - 1)
        cells: dict[tuple[int, int], list[float]] = {}
        for ii in range(N):
            A = atoms[ii]
            row = pi_rows[ii]
            ai = A & pref
            abit = 2 if A & bit else 0
            for jj in range(N):
                B = atoms[jj]
                key = (ai, B & pref)
                cell = cells.setdefault(key, [0.0, 0.0, 0.0, 0.0])
                cell[abit + (1 if B & bit else 0)] += row[jj]
        E_hz = E_hx = E_hy = 0.0
        max_dev = 0.0
        wsum_dev = 0.0
        def_mass = 0.0
        mass_or_lt1 = 0.0
        min_l2 = float("inf")
        max_l2 = -float("inf")
        E_xi2 = 0.0
        for (a, b), (p00, p01, p10, p11) in cells.items():
            mass = p00 + p01 + p10 + p11
            if mass <= 0.0:
                continue
            z = p00 / mass
            xt = (p00 + p01) / mass
            yt = (p00 + p10) / mass
            E_hz += mass * h(z)
            E_hx += mass * h(xt)
            E_hy += mass * h(yt)
            E_xi2 += mass * (xt - x0) ** 2
            if min(p00, p01, p10, p11) > 0.0:
                l2 = math.log2(p11 * p00 / (p10 * p01))
                dev = abs(l2 - lam)
                max_dev = max(max_dev, dev)
                wsum_dev += mass * dev
                def_mass += mass
                if l2 < 0.0:
                    mass_or_lt1 += mass
                min_l2 = min(min_l2, l2)
                max_l2 = max(max_l2, l2)
        stats.append({
            "i": i, "E_hz": E_hz, "E_hx": E_hx, "E_hy": E_hy,
            "max_abs_dev": max_dev,
            "wmean_abs_dev": (wsum_dev / def_mass) if def_mass > 0 else None,
            "defined_mass": def_mass, "mass_or_below_1": mass_or_lt1,
            "min_log2_or": (min_l2 if def_mass > 0 else None),
            "max_log2_or": (max_l2 if def_mass > 0 else None),
            "E_xi2": E_xi2,
        })
    return stats


def or_at_history(atoms, pi_rows, i: int, a: int, b: int):
    """log2 OR and mass at one history (i, a, b); None if degenerate."""
    pref = (1 << (i - 1)) - 1
    bit = 1 << (i - 1)
    cell = [0.0, 0.0, 0.0, 0.0]
    N = len(atoms)
    for ii in range(N):
        A = atoms[ii]
        if A & pref != a:
            continue
        row = pi_rows[ii]
        abit = 2 if A & bit else 0
        for jj in range(N):
            B = atoms[jj]
            if B & pref != b:
                continue
            cell[abit + (1 if B & bit else 0)] += row[jj]
    p00, p01, p10, p11 = cell
    mass = sum(cell)
    if min(cell) <= 0.0:
        return None, mass
    return math.log2(p11 * p00 / (p10 * p01)), mass


def union_gain(atoms, pi_rows, mu: dict[int, float]):
    N = len(atoms)
    ulaw: dict[int, float] = {}
    for ii in range(N):
        A = atoms[ii]
        row = pi_rows[ii]
        for jj in range(N):
            u = A | atoms[jj]
            ulaw[u] = ulaw.get(u, 0.0) + row[jj]
    HU = -sum(w * math.log2(w) for w in ulaw.values() if w > 0.0)
    return HU - entropy_bits(mu)


def assembly_report(mu, n, lam, p_base):
    """Fit + all assembly numbers for one instance."""
    atoms, pi_rows, _r, _c, resid = fit_tilt(mu, lam)
    Ha = chain_H(mu, n)
    st = coord_stats(atoms, pi_rows, n, lam, 1.0 - p_base)
    nets = [s["E_hz"] - Ha[s["i"] - 1] for s in st]
    taxA = [Ha[s["i"] - 1] - s["E_hx"] for s in st]
    taxB = [Ha[s["i"] - 1] - s["E_hy"] for s in st]
    gain_lb = sum(nets)
    gain = union_gain(atoms, pi_rows, mu)
    return {
        "resid": resid, "gain_lb": gain_lb, "gain": gain,
        "min_net": min(nets), "nets": nets,
        "max_taxA": max(taxA), "sum_taxA": sum(taxA),
        "max_taxB": max(taxB),
        "max_abs_dev": max(s["max_abs_dev"] for s in st),
        "max_wmean_dev": max((s["wmean_abs_dev"] or 0.0) for s in st),
        "max_E_xi2": max(s["E_xi2"] for s in st),
        "mass_or_below_1": sum(s["mass_or_below_1"] for s in st),
        "per_coord": st,
    }, atoms, pi_rows


def save(name: str, obj) -> None:
    DATA.mkdir(exist_ok=True)
    (DATA / name).write_text(json.dumps(obj, indent=1, sort_keys=True,
                                        default=float) + "\n",
                             encoding="utf-8")


# ---------------------------------------------------------------- part P0

def part_p0() -> None:
    print("=" * 78)
    print("P0. Engine validation")
    print("=" * 78)
    out = {}
    # (a) product: OR == 2^lambda at every history; gain closed form
    n, p, lam = 5, 0.30, 1.0
    mu = product_measure(n, p)
    rep, atoms, pi_rows = assembly_report(mu, n, lam, p)
    x = 1.0 - p
    g_pred = n * (h(z_rho(x, x, 2.0 ** lam)) - h(p))
    print(f"[P0a] product Bern({p})^{n} lam={lam}: max|log2OR-lam| = "
          f"{rep['max_abs_dev']:.2e}, gain = {rep['gain']:+.9f} vs closed form "
          f"{g_pred:+.9f} (|err| {abs(rep['gain']-g_pred):.1e}), "
          f"gain_lb = {rep['gain_lb']:+.9f}, resid {rep['resid']:.1e}")
    ok_a = rep["max_abs_dev"] < 1e-7 and abs(rep["gain"] - g_pred) < 1e-6
    out["product"] = {"max_dev": rep["max_abs_dev"], "gain": rep["gain"],
                      "gain_pred": g_pred, "gain_lb": rep["gain_lb"],
                      "pass": ok_a}
    # taxes must vanish on product
    print(f"[P0a] max tax_i on product = {rep['max_taxA']:.2e} (must be ~0)")
    # (b) reproduce 005 part G
    ref = json.loads((DATA / "005_partG.json").read_text())
    rows = []
    worst = 0.0
    for row in ref["rows"]:
        d = row["delta"]
        mu = mix_measures(product_measure(5, 0.30), product_measure(5, 0.05),
                          d, 5)
        atoms, pi_rows, _r, _c, resid = fit_tilt(mu, 1.0, tol=1e-14)
        st = coord_stats(atoms, pi_rows, 5, 1.0, 0.70)
        dev = max(s["max_abs_dev"] for s in st)
        diff = abs(dev - row["max_overall_dev"])
        worst = max(worst, diff)
        rows.append({"delta": d, "max_dev": dev,
                     "ref_005": row["max_overall_dev"], "abs_diff": diff})
    print(f"[P0b] 005 part-G reproduction (6 deltas): max abs diff vs "
          f"checkpoint = {worst:.2e} ({'PASS' if worst < 1e-6 else 'FAIL'})")
    out["partG_repro"] = {"rows": rows, "max_diff": worst}
    # (c) crash identity
    crash_rows = []
    for n2 in (5, 8):
        mu = crash_family(n2)
        atoms, pi_rows, _r, _c, resid = fit_tilt(mu, 1.0, tol=1e-14)
        l2, mass = or_at_history(atoms, pi_rows, 2, 0, 1)
        err = abs(l2 - (3 - n2))
        crash_rows.append({"n": n2, "log2_or": l2, "predicted": 3 - n2,
                           "err": err, "mass": mass})
        print(f"[P0c] crash n={n2}: log2 OR = {l2:+.6f} (predicted {3-n2}), "
              f"err {err:.1e}, history mass {mass:.4f}")
    out["crash"] = crash_rows
    save("pert_P0_validation.json", out)


# ---------------------------------------------------------------- part P1

def solve_linear(M: list[list[float]], b: list[float]) -> list[float]:
    """Gaussian elimination with partial pivoting (small N only)."""
    N = len(b)
    A = [row[:] + [b[i]] for i, row in enumerate(M)]
    for col in range(N):
        piv = max(range(col, N), key=lambda r: abs(A[r][col]))
        A[col], A[piv] = A[piv], A[col]
        d = A[col][col]
        for r in range(col + 1, N):
            f = A[r][col] / d
            if f != 0.0:
                for cc in range(col, N + 1):
                    A[r][cc] -= f * A[col][cc]
    x = [0.0] * N
    for r in range(N - 1, -1, -1):
        s = A[r][N] - sum(A[r][cc] * x[cc] for cc in range(r + 1, N))
        x[r] = s / A[r][r]
    return x


def part_p1() -> None:
    print("=" * 78)
    print("P1. IFT first-order formula and quadratic remainder (T1 numerics)")
    print("=" * 78)
    n, p, lam = 5, 0.30, 1.0
    mu0 = product_measure(n, p)
    nu = product_measure(n, 0.05)
    atoms, pi_rows, r0, c0, resid0 = fit_tilt(mu0, lam, tol=1e-14)
    N = len(atoms)
    v0 = [0.5 * math.log(r0[i] * c0[i]) for i in range(N)]
    # Hessian at the solution: diag(mu) + Pi
    H = [[pi_rows[i][j] + (mu0[atoms[i]] if i == j else 0.0)
          for j in range(N)] for i in range(N)]
    mudot = [nu.get(atoms[i], 0.0) - mu0[atoms[i]] for i in range(N)]
    vdot = solve_linear(H, mudot)
    # residual of the solve
    sres = max(abs(sum(H[i][j] * vdot[j] for j in range(N)) - mudot[i])
               for i in range(N))
    rows = []
    for delta in (0.01, 0.02, 0.04, 0.08):
        mu = mix_measures(mu0, nu, delta, n)
        a2, p2, r2, c2, res2 = fit_tilt(mu, lam, tol=1e-14)
        v = [0.5 * math.log(r2[i] * c2[i]) for i in range(N)]
        # gauge: scaling r -> t r, c -> c/t leaves pi fixed but shifts
        # sqrt(rc) by nothing (product rc is gauge-invariant).  v is gauge-free.
        err = max(abs(v[i] - v0[i] - delta * vdot[i]) for i in range(N))
        lin = max(abs(v[i] - v0[i]) for i in range(N))
        rows.append({"delta": delta, "remainder_inf": err,
                     "remainder_over_delta2": err / delta ** 2,
                     "step_inf": lin, "resid": res2})
        print(f"[P1] delta={delta:5.2f}: |v - v0 - delta*vdot|_inf = "
              f"{err:.3e}  /delta^2 = {err/delta**2:.4f}   "
              f"|v-v0|_inf = {lin:.3e}")
    ratios = [rows[k + 1]["remainder_over_delta2"] /
              rows[k]["remainder_over_delta2"] for k in range(len(rows) - 1)]
    print(f"[P1] remainder/delta^2 ratios between successive deltas: "
          + " ".join(f"{q:.3f}" for q in ratios)
          + "  (→ 1 means clean quadratic remainder; solve resid "
          f"{sres:.1e})")
    save("pert_P1_ift.json", {"n": n, "p": p, "lambda": lam,
                              "rows": rows, "solve_resid": sres})


# ---------------------------------------------------------------- part P2

def part_p2() -> None:
    print("=" * 78)
    print("P2. Assembly margins and end-to-end per-coordinate nets (lead 3)")
    print("=" * 78)
    out = {}
    # (a) margin closed forms
    p = 0.383
    rs = rho_star(p)
    print(f"[P2a] rho*(0.383) = {rs:.6f}; rho*(0.38271) = "
          f"{rho_star(0.38271):.6f}")
    mtab = []
    for rho in (1.031, 1.05, 1.06, 1.10, 1.20, 1.40, 2.00, 4.00):
        g = margin(rho, p)
        mtab.append({"rho": rho, "margin_bits": g})
        print(f"[P2a]   g(rho={rho:5.3f}, p=0.383) = {g:+.6f} bits/coord")
    # AHS-point coefficient: dg/drho at rho=1, p=psi equals h'(psi) phi^6
    x = PHI
    hp = math.log2((1.0 - PSI) / PSI)
    coef_closed = hp * x ** 6  # dz/drho|_1 = x^2 (1-x)^2 = phi^2 phi^4
    eps = 1e-6
    coef_num = (margin(1.0 + eps, PSI) - margin(1.0, PSI)) / eps
    print(f"[P2a] AHS point dg/drho|_1 = {coef_num:.6f} numeric vs "
          f"{coef_closed:.6f} closed form h'(psi)*phi^6 "
          f"(|err| {abs(coef_num-coef_closed):.1e})")
    out["margins"] = {"rho_star_0383": rs, "table": mtab,
                      "ahs_coef_closed": coef_closed,
                      "ahs_coef_numeric": coef_num}
    # (b) end-to-end nets on the delta-ball, marginal-capped directions
    p = 0.383
    dirs = {
        "crash": crash_family,
        "delta_empty": lambda n: {0: 1.0},
        "bern005": lambda n: product_measure(n, 0.05),
    }
    runs = []
    grid = [(math.log2(1.20), nn, dname, d)
            for nn in (5, 6, 7, 8)
            for dname in dirs
            for d in (0.05, 0.10, 0.20)]
    grid += [(math.log2(1.06), 6, dname, d)
             for dname in dirs for d in (0.02, 0.05, 0.10)]
    t0 = time.time()
    for lam, nn, dname, d in grid:
        nu = dirs[dname](nn)
        mu = mix_measures(product_measure(nn, p), nu, d, nn)
        capped = max(elem_marginals(mu, nn)) <= p + 1e-12
        rep, _a, _p = assembly_report(mu, nn, lam, p)
        del rep["per_coord"], rep["nets"]
        rep.update({"rho": 2.0 ** lam, "n": nn, "dir": dname, "delta": d,
                    "marginal_capped": capped,
                    "margin_at_base": margin(2.0 ** lam, p)})
        runs.append(rep)
        print(f"[P2b] rho={2.0**lam:5.3f} n={nn} nu={dname:11s} "
              f"delta={d:5.3f}: min net_i {rep['min_net']:+.6f} "
              f"(base margin {rep['margin_at_base']:+.6f}), gain_lb "
              f"{rep['gain_lb']:+.5f}, gain {rep['gain']:+.5f}, "
              f"max tax {rep['max_taxA']:.2e}, max|dev| "
              f"{rep['max_abs_dev']:.2e}, OR<1 mass {rep['mass_or_below_1']:.1e}")
    print(f"[P2b] total time {time.time()-t0:.1f}s")
    out["runs"] = runs
    # (c) cap-saturating and deep-history directions: the marginal-slack
    # first-order term is ~0 (pm03, slice) or the direction plants a
    # nu-dominated deep history (zero-prefix embedded crash), so the
    # second-order losses are exposed.
    def slice_dir(nn):
        w0 = round(p * nn)
        atoms = [m for m in range(1 << nn) if popcount(m) == w0]
        return {m: 1.0 / len(atoms) for m in atoms}

    dirs2 = {
        "pm03": lambda nn: mix_measures(product_measure(nn, p - 0.3),
                                        product_measure(nn, p + 0.3),
                                        0.5, nn),
        "slice": slice_dir,
        "emb0_crash_k2": lambda nn: embedded_crash(nn, 2),
    }
    runs2 = []
    lam = math.log2(1.20)
    for nn in (5, 6, 7, 8):
        for dname, mk in dirs2.items():
            if dname == "emb0_crash_k2" and nn < 6:
                continue
            nu = mk(nn)
            for d in (0.05, 0.10, 0.20):
                mu = mix_measures(product_measure(nn, p), nu, d, nn)
                capped = max(elem_marginals(mu, nn)) <= p + 1e-9
                rep, _a, _p2 = assembly_report(mu, nn, lam, p)
                del rep["per_coord"], rep["nets"]
                rep.update({"rho": 1.20, "n": nn, "dir": dname, "delta": d,
                            "marginal_capped": capped,
                            "margin_at_base": margin(1.20, p)})
                runs2.append(rep)
                print(f"[P2c] rho=1.200 n={nn} nu={dname:13s} "
                      f"delta={d:5.3f}: min net_i {rep['min_net']:+.6f} "
                      f"(base {rep['margin_at_base']:+.6f}), gain_lb "
                      f"{rep['gain_lb']:+.5f}, gain {rep['gain']:+.5f}, "
                      f"max tax {rep['max_taxA']:.2e}, max|dev| "
                      f"{rep['max_abs_dev']:.2e}, OR<1 mass "
                      f"{rep['mass_or_below_1']:.1e}, capped={capped}")
    out["runs_cap_saturating"] = runs2
    allruns = runs + runs2
    neg = [r for r in allruns if r["gain_lb"] <= 0]
    print(f"[P2] runs with chain-rule bound gain_lb <= 0: {len(neg)}"
          + ("" if not neg else "  <-- assembly loses at: "
             + ", ".join(f"(rho={r['rho']:.2f},n={r['n']},{r['dir']},"
                         f"d={r['delta']})" for r in neg)))
    negnet = [r for r in allruns if r["min_net"] <= 0]
    print(f"[P2] runs with some per-coordinate net_i <= 0: {len(negnet)}"
          + ("" if not negnet else "  at: "
             + ", ".join(f"(n={r['n']},{r['dir']},d={r['delta']})"
                         for r in negnet)))
    save("pert_P2_assembly.json", out)


def part_p2d() -> None:
    """Numeric constants for the conditional assembly theorem at
    (p, rho, s0) = (0.383, 1.2, 0.10): G(x,y,tau) = h(z_{2^{lam+tau}}(x,y))
    - h(x)/2 - h(y)/2; sups of derivatives over the good box
    [x0-s0, x0+s0]^2 by dense finite differences."""
    print("=" * 78)
    print("P2d. Explicit constants for the conditional theorem")
    print("=" * 78)
    p, rho, s0 = 0.383, 1.20, 0.10
    lam = math.log2(rho)
    x0 = 1.0 - p

    def G(x, y, tau):
        return h(z_rho(x, y, 2.0 ** (lam + tau))) - 0.5 * h(x) - 0.5 * h(y)

    g = G(x0, x0, 0.0)
    e = 1e-5
    Gx = (G(x0 + e, x0, 0) - G(x0 - e, x0, 0)) / (2 * e)
    # tau_half: tau with z = 1/2 at the box center (beyond it G can decrease)
    lo_t, hi_t = 0.0, 20.0
    for _ in range(60):
        mid = 0.5 * (lo_t + hi_t)
        if z_rho(x0, x0, 2.0 ** (lam + mid)) < 0.5:
            lo_t = mid
        else:
            hi_t = mid
    tau_half = lo_t
    grid = [x0 - s0 + 2 * s0 * k / 24 for k in range(25)]
    taus = [tau_half * k / 12 for k in range(13)]
    M_tau = 0.0      # sup |dG/dtau| on box x [0, tau_half]
    M2 = 0.0         # sup negative curvature of (x,y) -> G(x,y,tau), tau=0
    L2 = 0.0         # sup |dG/dx| on box at tau in [0, tau_half]
    dip = 0.0        # sup of (-G) over box x [0, 60] (worst huge-tau dip)
    for x in grid:
        for y in grid:
            for t in taus:
                M_tau = max(M_tau, abs((G(x, y, t + e) - G(x, y, t - e))
                                       / (2 * e)))
                L2 = max(L2, abs((G(x + e, y, t) - G(x - e, y, t)) / (2 * e)))
            for t in (0.0,):
                gxx = (G(x + e, y, t) - 2 * G(x, y, t) + G(x - e, y, t)) / e**2
                gyy = (G(x, y + e, t) - 2 * G(x, y, t) + G(x, y - e, t)) / e**2
                gxy = (G(x + e, y + e, t) - G(x + e, y - e, t)
                       - G(x - e, y + e, t) + G(x - e, y - e, t)) / (4 * e**2)
                tr, det = gxx + gyy, gxx * gyy - gxy * gxy
                lam_min = 0.5 * (tr - math.sqrt(max(tr * tr - 4 * det, 0.0)))
                M2 = max(M2, max(0.0, -lam_min))
            for t in (5.0, 10.0, 30.0, 60.0):
                dip = max(dip, -G(x, y, t))
    print(f"[P2d] p={p} rho={rho} s0={s0}: g = {g:.6f}, G_x(center) = "
          f"{Gx:+.4f}, tau_half = {tau_half:.3f}")
    print(f"[P2d] M_tau = sup|dG/dtau| = {M_tau:.4f};  M2 = sup neg "
          f"curvature (tau=0) = {M2:.4f};  L2 = sup|dG/dx| = {L2:.4f}")
    print(f"[P2d] worst huge-tau dip on box: sup(-G) = {dip:.6f}")
    # instantiated delta_0 with the budget constants measured in P3/P3c
    c1, c2, c3 = 0.07, 0.01, 2.5   # avg downward OR, avg tax, xi-variance
    # per-coordinate loss <= M_tau*c1*d^2 + c2*d^2 + M2*c3*d^2 + dip-term
    denom = M_tau * c1 + c2 + M2 * c3
    d0 = math.sqrt(0.5 * g / denom)
    print(f"[P2d] with measured budgets (c1,c2,c3) = ({c1},{c2},{c3}): "
          f"loss <= {denom:.3f}*delta^2; delta_0 = sqrt(g/2 / loss-rate) = "
          f"{d0:.4f}")
    save("pert_P2d_constants.json", {
        "p": p, "rho": rho, "s0": s0, "g": g, "G_x_center": Gx,
        "tau_half": tau_half, "M_tau": M_tau, "M2": M2, "L2": L2,
        "dip": dip, "c_budgets": [c1, c2, c3], "delta0": d0})


# ---------------------------------------------------------------- part P3

def part_p3() -> None:
    print("=" * 78)
    print("P3. Budget growth in n (is any budget constant n-uniform?)")
    print("=" * 78)
    p, lam = 0.30, 1.0
    dirs = {
        "crash": lambda n: crash_family(n),
        "emb0_crash_k2": lambda n: embedded_crash(n, 2),
        "delta_empty": lambda n: {0: 1.0},
        "bern07": lambda n: product_measure(n, 0.70),
    }
    deltas = (0.0125, 0.025, 0.05, 0.1)
    rows = []
    t0 = time.time()
    for nn in (4, 5, 6, 7, 8, 9):
        for dname, mk in dirs.items():
            if dname == "emb0_crash_k2" and nn < 6:
                continue
            if nn == 9 and dname not in ("crash", "emb0_crash_k2"):
                continue
            nu = mk(nn)
            for d in deltas:
                mu = mix_measures(product_measure(nn, p), nu, d, nn)
                atoms, pi_rows, _r, _c, resid = fit_tilt(mu, lam)
                st = coord_stats(atoms, pi_rows, nn, lam, 1.0 - p)
                Ha = chain_H(mu, nn)
                taxes = [Ha[s["i"] - 1] - s["E_hx"] for s in st]
                row = {
                    "n": nn, "dir": dname, "delta": d, "resid": resid,
                    "max_abs_dev": max(s["max_abs_dev"] for s in st),
                    "max_wmean_dev": max((s["wmean_abs_dev"] or 0.0)
                                         for s in st),
                    "max_tax": max(taxes), "sum_tax": sum(taxes),
                    "max_E_xi2": max(s["E_xi2"] for s in st),
                    "mass_or_below_1": sum(s["mass_or_below_1"] for s in st),
                }
                if dname == "emb0_crash_k2":
                    l2, mass = or_at_history(atoms, pi_rows, 4, 0b000, 0b100)
                    row["emb_hist_log2_or"] = l2
                    row["emb_hist_mass"] = mass
                rows.append(row)
        done = [r for r in rows if r["n"] == nn]
        print(f"[P3] n={nn} done ({len(done)} rows, "
              f"{time.time()-t0:.0f}s elapsed)")
    # digest: c(n) = dev/delta^2 at the two smallest deltas; growth at fixed d
    print(f"\n[P3] c_hat(n) = max_abs_dev/delta^2 at delta=0.025, and "
          f"max_abs_dev at fixed delta=0.1:")
    print(f"[P3] {'dir':>13s} " + " ".join(f"n={nn}: c_hat / dev@0.1"
                                           for nn in (4, 5, 6, 7, 8, 9)))
    for dname in dirs:
        cells = []
        for nn in (4, 5, 6, 7, 8, 9):
            r_small = [r for r in rows if r["n"] == nn and r["dir"] == dname
                       and r["delta"] == 0.025]
            r_big = [r for r in rows if r["n"] == nn and r["dir"] == dname
                     and r["delta"] == 0.1]
            if not r_small:
                cells.append("      --      ")
                continue
            cells.append(f"{r_small[0]['max_abs_dev']/0.025**2:8.2f} / "
                         f"{r_big[0]['max_abs_dev']:7.4f}")
        print(f"[P3] {dname:>13s} " + " ".join(cells))
    print(f"\n[P3] mass-weighted analogue c_hat_w(n) = "
          f"max_wmean_dev/delta^2 at delta=0.025:")
    for dname in dirs:
        cells = []
        for nn in (4, 5, 6, 7, 8, 9):
            r_small = [r for r in rows if r["n"] == nn and r["dir"] == dname
                       and r["delta"] == 0.025]
            cells.append("    --  " if not r_small else
                         f"{r_small[0]['max_wmean_dev']/0.025**2:8.3f}")
        print(f"[P3] {dname:>13s} " + " ".join(cells))
    print(f"\n[P3] per-coordinate tax and E_xi2 growth at delta=0.1:")
    for dname in dirs:
        cells = []
        for nn in (4, 5, 6, 7, 8, 9):
            rr = [r for r in rows if r["n"] == nn and r["dir"] == dname
                  and r["delta"] == 0.1]
            cells.append("      --      " if not rr else
                         f"{rr[0]['max_tax']:.2e}/{rr[0]['max_E_xi2']:.2e}")
        print(f"[P3] {dname:>13s} " + " ".join(cells))
    emb = [r for r in rows if r["dir"] == "emb0_crash_k2"]
    if emb:
        print(f"\n[P3] embedded-crash designated history (i=4, a=0011, "
              f"b=0111): log2 OR by (n, delta):")
        for r in emb:
            print(f"[P3]   n={r['n']} delta={r['delta']:6.4f}: "
                  f"log2 OR = {r['emb_hist_log2_or']:+9.4f} "
                  f"(pure-family value {1.0*(3-(r['n']-2)):+.1f}), "
                  f"mass {r['emb_hist_mass']:.2e}")
    save("pert_P3_budgets.json", {"p": p, "lambda": lam, "rows": rows})


def part_p3b() -> None:
    """Anatomy of the worst deviations: which histories, what mass, what sign."""
    print("=" * 78)
    print("P3b. Anatomy of the max deviations (crash direction)")
    print("=" * 78)
    p, lam = 0.30, 1.0
    out = {}
    for nn, d in ((8, 0.05), (9, 0.05)):
        nu = crash_family(nn)
        mu = mix_measures(product_measure(nn, p), nu, d, nn)
        atoms, pi_rows, _r, _c, resid = fit_tilt(mu, lam)
        N = len(atoms)
        rows = []
        below_lam_mass = 0.0
        below_1_mass = 0.0
        tot_def = 0.0
        for i in range(1, nn + 1):
            pref = (1 << (i - 1)) - 1
            bit = 1 << (i - 1)
            cells: dict[tuple[int, int], list[float]] = {}
            for ii in range(N):
                A = atoms[ii]
                row = pi_rows[ii]
                ai = A & pref
                abit = 2 if A & bit else 0
                for jj in range(N):
                    B = atoms[jj]
                    key = (ai, B & pref)
                    cell = cells.setdefault(key, [0.0, 0.0, 0.0, 0.0])
                    cell[abit + (1 if B & bit else 0)] += row[jj]
            for (a, b), cl in cells.items():
                if min(cl) <= 0.0:
                    continue
                mass = sum(cl)
                l2 = math.log2(cl[3] * cl[0] / (cl[2] * cl[1]))
                tot_def += mass
                if l2 < lam - 1e-12:
                    below_lam_mass += mass
                if l2 < 0.0:
                    below_1_mass += mass
                rows.append((abs(l2 - lam), l2, i, a, b, mass))
        rows.sort(reverse=True)
        print(f"[P3b] crash dir, n={nn}, delta={d}: top-5 deviation histories "
              f"(dev, log2OR, i, a, b, mass):")
        for dev, l2, i, a, b, mass in rows[:5]:
            print(f"[P3b]   dev={dev:7.4f} log2OR={l2:+8.4f} i={i} "
                  f"a={a:0{nn}b} b={b:0{nn}b} mass={mass:.3e}")
        # the DOWNWARD side (what a one-sided averaged budget must control)
        down = sorted((r for r in rows if r[1] < lam), key=lambda t: t[1])
        worst_down = down[0] if down else None
        down_wsum = sum((lam - r[1]) * r[5] for r in down)
        if worst_down:
            print(f"[P3b]   worst BELOW-lam history: log2OR="
                  f"{worst_down[1]:+8.4f} i={worst_down[2]} "
                  f"mass={worst_down[5]:.3e}")
        print(f"[P3b]   mass-weighted sum of (lam - log2OR)_+ over all "
              f"coords: {down_wsum:.3e}  (per coordinate "
              f"{down_wsum/nn:.3e}, /delta^2 = {down_wsum/nn/d**2:.4f})")
        print(f"[P3b]   mass with log2OR < lam: {below_lam_mass:.3e}; "
              f"with OR < 1: {below_1_mass:.3e} (defined mass "
              f"{tot_def:.4f}, per-coordinate avg {tot_def/nn:.4f})")
        out[f"n{nn}"] = {
            "delta": d, "top": [
                {"dev": t[0], "log2_or": t[1], "i": t[2], "a": t[3],
                 "b": t[4], "mass": t[5]} for t in rows[:10]],
            "worst_below_lam": (None if not worst_down else
                                {"log2_or": worst_down[1],
                                 "i": worst_down[2],
                                 "mass": worst_down[5]}),
            "downward_wsum": down_wsum,
            "mass_below_lam": below_lam_mass, "mass_below_1": below_1_mass}
    save("pert_P3b_anatomy.json", out)


def part_p3c() -> None:
    """n-trend of the ONE-SIDED (downward) mass-weighted OR budget --
    the perturbative form of gap a-prime -- in both parameter regimes."""
    print("=" * 78)
    print("P3c. Downward budget c_minus(n) = E[(lam - log2OR)_+]/delta^2")
    print("=" * 78)
    out = {}
    d = 0.05
    for tag, p, lam in (("p030_lam1", 0.30, 1.0),
                        ("p0383_rho12", 0.383, math.log2(1.2))):
        rows = []
        for nn in (4, 5, 6, 7, 8, 9):
            nu = crash_family(nn)
            mu = mix_measures(product_measure(nn, p), nu, d, nn)
            atoms, pi_rows, _r, _c, resid = fit_tilt(mu, lam)
            N = len(atoms)
            down_wsum = 0.0
            worst_i = None
            per_i = []
            for i in range(1, nn + 1):
                pref = (1 << (i - 1)) - 1
                bit = 1 << (i - 1)
                cells: dict[tuple[int, int], list[float]] = {}
                for ii in range(N):
                    A = atoms[ii]
                    row = pi_rows[ii]
                    ai = A & pref
                    abit = 2 if A & bit else 0
                    for jj in range(N):
                        B = atoms[jj]
                        key = (ai, B & pref)
                        cell = cells.setdefault(key, [0.0, 0.0, 0.0, 0.0])
                        cell[abit + (1 if B & bit else 0)] += row[jj]
                di = 0.0
                for cl in cells.values():
                    if min(cl) <= 0.0:
                        continue
                    l2 = math.log2(cl[3] * cl[0] / (cl[2] * cl[1]))
                    if l2 < lam:
                        di += (lam - l2) * sum(cl)
                per_i.append(di)
                down_wsum += di
            worst_i = max(range(nn), key=lambda k: per_i[k]) + 1
            rows.append({"n": nn, "delta": d,
                         "down_per_coord_max": max(per_i),
                         "down_per_coord_max_over_d2": max(per_i) / d ** 2,
                         "down_total": down_wsum,
                         "down_avg_over_d2": down_wsum / nn / d ** 2,
                         "worst_i": worst_i, "per_i": per_i, "resid": resid})
            print(f"[P3c] {tag}: n={nn}: max_i E[(lam-log2OR)_+] = "
                  f"{max(per_i):.3e}  /delta^2 = {max(per_i)/d**2:.4f} "
                  f"(worst at i={worst_i}); per-coord AVG /delta^2 = "
                  f"{down_wsum/nn/d**2:.4f}")
        out[tag] = rows
    save("pert_P3c_downward.json", out)


# ---------------------------------------------------------------- part P4

def part_p4() -> None:
    print("=" * 78)
    print("P4. Lead 4's one-line kernel change: two-scale / centered kernel")
    print("=" * 78)
    out = {}
    lam1 = 1.0
    rows = []
    for nn in (5, 8):
        # centered-overlap kernel: 2^{lam1 (|A&B| - |A||B|/n)}
        kern = (lambda n_: (lambda a, b: 2.0 ** (
            lam1 * (popcount(a & b) - popcount(a) * popcount(b) / n_))))(nn)
        mu = crash_family(nn)
        atoms, pi_rows, _r, _c, resid = fit_tilt(mu, lam1, kernel=kern,
                                                 tol=1e-14)
        l2, mass = or_at_history(atoms, pi_rows, 2, 0, 1)
        pred = lam1 * (3 - nn) / nn
        st = coord_stats(atoms, pi_rows, nn, lam1, 1.0 - 0.383)
        min_l2 = min(s["min_log2_or"] for s in st
                     if s["min_log2_or"] is not None)
        max_l2 = max(s["max_log2_or"] for s in st
                     if s["max_log2_or"] is not None)
        rows.append({"n": nn, "crash_log2_or": l2, "predicted": pred,
                     "err": abs(l2 - pred), "min_log2_or": min_l2,
                     "max_log2_or": max_l2, "resid": resid})
        print(f"[P4] centered kernel, crash n={nn}: crash-history log2 OR = "
              f"{l2:+.6f} vs closed form lam(3-n)/n = {pred:+.6f} "
              f"(err {abs(l2-pred):.1e}); census min/max log2 OR = "
              f"{min_l2:+.3f}/{max_l2:+.3f}  [plain tilt: {3-nn:+d}]")
    out["crash_centered"] = rows
    # mirror family under centered kernel
    nn = 5
    kern = (lambda a, b: 2.0 ** (
        lam1 * (popcount(a & b) - popcount(a) * popcount(b) / nn)))
    mu = mirror_family(nn)
    atoms, pi_rows, _r, _c, resid = fit_tilt(mu, lam1, kernel=kern, tol=1e-14)
    l2, mass = or_at_history(atoms, pi_rows, 2, 0, 1)
    print(f"[P4] centered kernel, mirror n={nn}: history log2 OR = {l2:+.4f} "
          f"[plain tilt: {nn-1:+d}]")
    out["mirror_centered"] = {"n": nn, "log2_or": l2}
    # products under centered kernel: does OR stay in a band around lam1?
    prows = []
    for nn in (6, 8):
        kern = (lambda n_: (lambda a, b: 2.0 ** (
            lam1 * (popcount(a & b) - popcount(a) * popcount(b) / n_))))(nn)
        p = 0.383
        mu = product_measure(nn, p)
        atoms, pi_rows, _r, _c, resid = fit_tilt(mu, lam1, kernel=kern)
        st = coord_stats(atoms, pi_rows, nn, lam1, 1.0 - p)
        min_l2 = min(s["min_log2_or"] for s in st
                     if s["min_log2_or"] is not None)
        max_l2 = max(s["max_log2_or"] for s in st
                     if s["max_log2_or"] is not None)
        Ha = chain_H(mu, nn)
        nets = [s["E_hz"] - Ha[s["i"] - 1] for s in st]
        gain = union_gain(atoms, pi_rows, mu)
        prows.append({"n": nn, "min_log2_or": min_l2, "max_log2_or": max_l2,
                      "gain": gain, "gain_lb": sum(nets), "resid": resid})
        print(f"[P4] centered kernel, product Bern(0.383)^{nn}: log2 OR in "
              f"[{min_l2:+.4f}, {max_l2:+.4f}] (plain tilt: == {lam1}); "
              f"gain {gain:+.5f}, gain_lb {sum(nets):+.5f}")
    out["product_centered"] = prows
    save("pert_P4_twoscale.json", out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="all")
    args = ap.parse_args()
    todo = args.part.upper()
    t0 = time.time()
    if todo in ("P0", "ALL"):
        part_p0()
    if todo in ("P1", "ALL"):
        part_p1()
    if todo in ("P2", "ALL"):
        part_p2()
    if todo in ("P2D", "ALL"):
        part_p2d()
    if todo in ("P3", "ALL"):
        part_p3()
    if todo in ("P3B", "ALL"):
        part_p3b()
    if todo in ("P3C", "ALL"):
        part_p3c()
    if todo in ("P4", "ALL"):
        part_p4()
    print(f"done in {time.time()-t0:.1f}s; checkpoints in {DATA}")


if __name__ == "__main__":
    main()
