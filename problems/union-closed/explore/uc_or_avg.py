#!/usr/bin/env python3
"""Attempt 007: the restated Gap 1 -- history-averaged odds-ratio control.

Object under test (005 lead 1, restated after 006's bookkeeping note):

    M_i(mu, lambda) = E_{(a,b) ~ pi_lambda}[ log2 OR_i(a,b) | nondegenerate ]
                    >= lambda   for all mu, all i, all lambda >= 0 ?

Well-posed definition used throughout (part W provides the computational side
of the forcing argument; the record has the proofs):

  * P_i = reachable (i-1)-bit prefixes; N_i = prefixes in P_i whose next bit
    is NOT deterministic (both extensions reachable in supp(mu)).
  * For the tilt coupling, a positive-mass history (a,b) has all four
    conditional cells positive  <=>  a in N_i and b in N_i  (degeneracy
    dichotomy: cells inherit positivity from the strictly positive kernel, so
    a zero cell is always a zero conditional MARGIN, never a mixed pattern;
    the odds ratio there is unidentifiable, not 0 or infinity).
  * M_i = sum over N_i x N_i of mass * log2 OR / (mass of N_i x N_i);
    undefined (vacuous) iff N_i is empty.

Parts:
  W   well-posedness computations: (W1) degeneracy dichotomy assertion over
      random sparse supports; (W2) no-continuous-extension: at a degenerate
      history of the 4-atom crash family, the eps->0 limit of log2 OR depends
      on the mixing direction, so no value assignment is continuous.
  C   cross-checks of this (third) engine against 006_partS8 (crash/mirror
      M_i) and 005_partE (mass-weighted means) checkpoints.
  K   kill census: crash + mirror (several n, off-grid lambda), eps-mixed
      crash at i >= 3 (histories 006 flagged as zero-mass become positive-
      mass), delta_0 + Bern(1/2+eps) genre, Chase-Lovett slice and smoothed
      slice, Sawin-genre geometric mixtures, random sparse supports.
  L   structure probe: per instance, the matrix L_ab = log2 OR_i(a,b) - lambda
      over N_i x N_i -- minimum eigenvalue (Jacobi), plus the identity
      M_i - lambda = <m, L>_F / <m, J>_F with m the history-mass matrix
      (m is PSD: Gram in the 2^{lambda|A^B|} inner product -- proved in the
      record; L PSD would imply M_i >= lambda, but L is NOT always PSD).
  S   adversarial search: multiplicative hill-climb on potentials over sparse
      supports (5..12 atoms, sizes where the point-mass theorem does not
      apply), minimizing min_{i<n} (M_i - lambda), with a penalty steering
      max elementwise marginal below the record 0.38271; plus a second
      search minimizing the min eigenvalue of L and re-testing M_i there.
      (Finds nothing: generic supports do not violate.)
  T   THE HEADLINE: explicit 10-atom in-regime counterexample (n = 7,
      lambda = 3.5, M_5 - lambda = -0.122, all marginals <= 0.318,
      H(mu) = 2.72 bits) with its verification battery: Sinkhorn
      round-trip, dilution invariance, 3% perturbations, 2-digit rounding,
      lambda profile (violation for all lambda >~ 0.03; positive below --
      Theorem C), light-slice weight sweep.  See the part-T banner for the
      float-underflow artifact this replaced and how it was caught.
  U   underflow-clamped hill-climbs re-discovering in-regime violations
      independently (single-i mode) and hunting the i-aggregated form
      (mass-weighted average over i < n).

Standard library only.  Deterministic (fixed seeds).  Runtime ~4-6 min.

Usage:
    python problems/union-closed/explore/uc_or_avg.py
Checkpoints: problems/union-closed/data/or_avg_part[WCKLSTU].json
             problems/union-closed/data/or_avg_run.log (tee by hand)
"""

from __future__ import annotations

import json
import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

RECORD = 0.38271


def popcount(x: int) -> int:
    return bin(x).count("1")


def save(name: str, obj) -> None:
    DATA.mkdir(exist_ok=True)
    (DATA / name).write_text(json.dumps(obj, indent=1, sort_keys=True) + "\n",
                             encoding="utf-8")


# ----------------------------------------------------------------------
# core (written fresh for 007; cross-checked against 005/006 checkpoints
# in part C)
# ----------------------------------------------------------------------

def tilt_of_potential(u: dict[int, float], lam: float) -> dict[tuple[int, int], float]:
    """pi(A,B) prop u(A)u(B)2^{lam|A&B|}.  By uniqueness of Sinkhorn scaling
    this IS the tilt coupling of its own marginal (005's trick)."""
    z = 0.0
    pi = {}
    for a, ua in u.items():
        for b, ub in u.items():
            v = ua * ub * 2.0 ** (lam * popcount(a & b))
            pi[(a, b)] = v
            z += v
    return {k: v / z for k, v in pi.items()}


def sinkhorn(mu: dict[int, float], lam: float, tol: float = 1e-13,
             iters: int = 60000) -> tuple[dict[tuple[int, int], float], float]:
    """Alternating scaling for pi(A,B)=r(A)c(B)2^{lam|A&B|} with both
    marginals mu.  Returns (pi, max marginal residual over both margins)."""
    supp = sorted(mu)
    N = len(supp)
    K = [[2.0 ** (lam * popcount(a & b)) for b in supp] for a in supp]
    t = [mu[a] for a in supp]
    r = [1.0] * N
    c = [1.0] * N
    for _ in range(iters):
        for i in range(N):
            r[i] = t[i] / sum(c[j] * K[i][j] for j in range(N))
        for j in range(N):
            c[j] = t[j] / sum(r[i] * K[i][j] for i in range(N))
        row = max(abs(sum(r[i] * c[j] * K[i][j] for j in range(N)) - t[i])
                  for i in range(N))
        if row < tol:
            break
    col = max(abs(sum(r[i] * c[j] * K[i][j] for i in range(N)) - t[j])
              for j in range(N))
    pi = {(supp[i], supp[j]): r[i] * c[j] * K[i][j]
          for i in range(N) for j in range(N)}
    return pi, max(row, col)


def marginal_elements(mu: dict[int, float], n: int) -> list[float]:
    return [sum(w for m, w in mu.items() if m >> j & 1) for j in range(n)]


def mu_of(pi: dict[tuple[int, int], float]) -> dict[int, float]:
    mu: dict[int, float] = {}
    for (a, _b), w in pi.items():
        mu[a] = mu.get(a, 0.0) + w
    return mu


def census(pi: dict[tuple[int, int], float], n: int) -> list[dict]:
    """Per coordinate i: the well-posed M_i and its bookkeeping.

    Returns one dict per i with:
      n_prefixes (|P_i| reachable), n_active (|N_i|), defined_mass (mass of
      N_i x N_i), M_i (None if N_i empty), min/max log2 OR over N_i x N_i,
      rows (list of (a, b, mass, log2or) for the L-matrix), and the
      degeneracy-dichotomy check flag.
    """
    out = []
    supp = {a for (a, _b) in pi}
    for i in range(1, n + 1):
        pmask = (1 << (i - 1)) - 1
        bit = 1 << (i - 1)
        # N_i from the support alone (this is the point: membership is a
        # property of supp(mu), not of the table masses)
        prefixes = {A & pmask for A in supp}
        active = {c for c in prefixes
                  if any(A & pmask == c and A & bit for A in supp)
                  and any(A & pmask == c and not A & bit for A in supp)}
        tables: dict[tuple[int, int], list[float]] = {}
        for (A, B), w in pi.items():
            key = (A & pmask, B & pmask)
            tab = tables.setdefault(key, [0.0] * 4)
            tab[(2 if A & bit else 0) + (1 if B & bit else 0)] += w
        dichotomy_ok = True
        rows = []
        def_mass = 0.0
        for (a, b), (p00, p01, p10, p11) in tables.items():
            nondeg = min(p00, p01, p10, p11) > 0.0
            in_nn = a in active and b in active
            if nondeg != in_nn:
                dichotomy_ok = False
            if nondeg:
                mass = p00 + p01 + p10 + p11
                # sum of logs: robust to underflow of cell products
                l2or = (math.log2(p11) + math.log2(p00)
                        - math.log2(p10) - math.log2(p01))
                rows.append((a, b, mass, l2or))
                def_mass += mass
        mi = (sum(m * v for (_a, _b, m, v) in rows) / def_mass
              if def_mass > 0 else None)
        out.append({
            "i": i, "n_prefixes": len(prefixes), "n_active": len(active),
            "defined_mass": def_mass, "M_i": mi,
            "min_log2_or": min((v for (_a, _b, _m, v) in rows), default=None),
            "max_log2_or": max((v for (_a, _b, _m, v) in rows), default=None),
            "rows": rows, "dichotomy_ok": dichotomy_ok,
        })
    return out


def min_excess(cens: list[dict], lam: float, lt_n: bool = False):
    """min over defined i of M_i - lambda (None if nothing defined).

    lt_n=True restricts to i < n: at i = n, M_n = lambda identically
    (Prop 1), so the all-i minimum is floored at ~0 by float noise there and
    carries no information; the i < n margin is the real quantity."""
    n = len(cens)
    vals = [c["M_i"] - lam for c in cens
            if c["M_i"] is not None and (not lt_n or c["i"] < n)]
    return min(vals) if vals else None


# ----------------------------------------------------------------------
# instance builders
# ----------------------------------------------------------------------

def bern(n: int, p: float) -> dict[int, float]:
    return {m: p ** popcount(m) * (1 - p) ** (n - popcount(m))
            for m in range(1 << n)}


def mixture(n: int, comps) -> dict[int, float]:
    mu = {m: 0.0 for m in range(1 << n)}
    for w, p in comps:
        for m, v in bern(n, p).items():
            mu[m] += w * v
    return {m: v for m, v in mu.items() if v > 0}


def crash_mu(n: int) -> dict[int, float]:
    full = (1 << n) - 1
    return {0b10: 0.33, full ^ 0b11: 0.33, full: 0.04, 0b01: 0.30}


def mirror_mu(n: int) -> dict[int, float]:
    full = (1 << n) - 1
    return {0b10 | (full ^ 0b11): 0.3, 0b00: 0.2, full: 0.3, 0b01: 0.2}


def slice_mu(n: int, w0: int) -> dict[int, float]:
    atoms = [m for m in range(1 << n) if popcount(m) == w0]
    return {m: 1.0 / len(atoms) for m in atoms}


def smoothed_slice_mu(n: int, w0: int, t: float) -> dict[int, float]:
    mu = {m: (1 - t) * v for m, v in slice_mu(n, w0).items()}
    for m, v in bern(n, w0 / n).items():
        mu[m] = mu.get(m, 0.0) + t * v
    return mu


# ----------------------------------------------------------------------
# Jacobi eigenvalues (for the L-matrix probe; pure stdlib)
# ----------------------------------------------------------------------

def jacobi_min_eig(mat: list[list[float]]) -> float:
    n = len(mat)
    if n == 1:
        return mat[0][0]
    a = [row[:] for row in mat]
    for _sweep in range(200):
        off = math.sqrt(sum(a[i][j] ** 2 for i in range(n)
                            for j in range(n) if i != j))
        if off < 1e-13:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-18:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = (1 if theta >= 0 else -1) / (abs(theta) + math.sqrt(theta * theta + 1))
                cth = 1 / math.sqrt(t * t + 1)
                sth = t * cth
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = cth * akp - sth * akq
                    a[k][q] = sth * akp + cth * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = cth * apk - sth * aqk
                    a[q][k] = sth * apk + cth * aqk
    return min(a[i][i] for i in range(n))


def l_matrix_probe(cens_row: dict, pi, lam: float):
    """Build L over N_i x N_i from census rows; return (min_eig, frob_check)
    where frob_check = |<m,L>/<m,J> - (M_i - lam)| (identity sanity)."""
    rows = cens_row["rows"]
    if not rows:
        return None, None
    idx = sorted({a for (a, _b, _m, _v) in rows} | {b for (_a, b, _m, _v) in rows})
    pos = {c: k for k, c in enumerate(idx)}
    A = len(idx)
    L = [[0.0] * A for _ in range(A)]
    M = [[0.0] * A for _ in range(A)]
    for (a, b, m, v) in rows:
        L[pos[a]][pos[b]] = v - lam
        M[pos[a]][pos[b]] = m
    frob = sum(M[x][y] * L[x][y] for x in range(A) for y in range(A))
    tot = sum(M[x][y] for x in range(A) for y in range(A))
    frob_check = abs(frob / tot - (cens_row["M_i"] - lam))
    return jacobi_min_eig(L), frob_check


# ----------------------------------------------------------------------
# part W -- well-posedness computations
# ----------------------------------------------------------------------

def part_w() -> dict:
    rng = random.Random(70001)
    # W1: degeneracy dichotomy over random sparse supports
    checks = 0
    all_ok = True
    for _ in range(200):
        n = rng.choice([4, 5, 6])
        k = rng.randrange(3, 10)
        atoms = rng.sample(range(1 << n), min(k, 1 << n))
        u = {a: math.exp(rng.uniform(-2, 2)) for a in atoms}
        lam = rng.choice([0.3, 1.0, 2.7])
        pi = tilt_of_potential(u, lam)
        for c in census(pi, n):
            checks += 1
            if not c["dichotomy_ok"]:
                all_ok = False
    print(f"[W1] degeneracy dichotomy (cell zero <=> margin zero <=> prefix "
          f"outside N_i): {checks} coordinate-censuses on 200 random sparse "
          f"supports -- {'PASS' if all_ok else 'FAIL'}")

    # W2: no continuous extension at a degenerate history.
    # Crash family n=5: at i=3 the history a=b=(0,1) pins atom S1={2} on both
    # sides (positive mass, fully degenerate table).  Mix mu with eps * nu for
    # two different nu; if lim_{eps->0} log2 OR differ, no assignment of a
    # value to the degenerate history is continuous in mu.
    n, lam = 5, 1.0
    base = crash_mu(n)
    hist = (0b10, 0b10)   # prefixes (bit1,bit2) = (0,1) both sides
    def or_at(mu):
        pi, resid = sinkhorn(mu, lam)
        assert resid < 1e-9, resid
        pmask, bit = 0b11, 0b100
        tab = [0.0] * 4
        for (A, B), w in pi.items():
            if (A & pmask, B & pmask) == hist:
                tab[(2 if A & bit else 0) + (1 if B & bit else 0)] += w
        return math.log2(tab[3] * tab[0] / (tab[2] * tab[1]))
    dirs = {
        "unif": {m: 1.0 / (1 << n) for m in range(1 << n)},
        "bern015": bern(n, 0.15),
        # two atoms with prefix (0,1) and both values of bit 3, so the
        # degenerate table is repopulated along a very sparse direction
        "two_atoms": {0b00110: 0.5, 0b01010: 0.5},
    }
    lims = {}
    for name, nu in dirs.items():
        vals = []
        for eps in (1e-2, 1e-3, 1e-4):
            mu = {m: (1 - eps) * base.get(m, 0.0) + eps * nu.get(m, 0.0)
                  for m in range(1 << n)}
            mu = {m: v for m, v in mu.items() if v > 0}
            vals.append(or_at(mu))
        lims[name] = vals
        print(f"[W2] direction {name:9s}: log2 OR at degenerate history "
              f"(i=3, a=b=01) for eps=1e-2,1e-3,1e-4 -> "
              + ", ".join(f"{v:+.4f}" for v in vals))
    spread = max(v[-1] for v in lims.values()) - min(v[-1] for v in lims.values())
    print(f"[W2] spread of eps=1e-4 values across directions = {spread:.4f} "
          f"({'direction-DEPENDENT limit -> no continuous extension' if spread > 0.05 else 'inconclusive'})")
    out = {"dichotomy_checks": checks, "dichotomy_all_ok": all_ok,
           "w2_limits": lims, "w2_spread": spread}
    save("or_avg_partW.json", out)
    return out


# ----------------------------------------------------------------------
# part C -- cross-checks against 005/006 checkpoints
# ----------------------------------------------------------------------

def part_c() -> dict:
    # against 006_partS8 (crash + mirror M_i at lambda = 1)
    theirs = json.loads((DATA / "006_partS8.json").read_text())
    max_diff_s8 = 0.0
    for rec in theirs:
        mk = crash_mu if rec["family"] == "crash" else mirror_mu
        n, lam = rec["n"], rec["lambda"]
        pi, resid = sinkhorn(mk(n), lam)
        cens = census(pi, n)
        for p in rec["per_i"]:
            mine = cens[p["i"] - 1]["M_i"]
            if p["M_i"] is None:
                assert mine is None, (rec["family"], n, p["i"])
            else:
                max_diff_s8 = max(max_diff_s8, abs(mine - p["M_i"]))
    print(f"[C] vs 006_partS8 (crash+mirror M_i): max |diff| = {max_diff_s8:.2e} "
          f"({'PASS' if max_diff_s8 < 1e-7 else 'FAIL'})")

    # against 005_partE (mass-weighted means on the realistic mixture)
    theirs_e = json.loads((DATA / "005_partE.json").read_text())
    n = theirs_e["n"]
    mu = mixture(n, [(0.7, 0.30), (0.3, 0.05)])
    max_diff_e = 0.0
    for run in theirs_e["runs"]:
        lam = run["lambda"]
        pi, _res = sinkhorn(mu, lam)
        cens = census(pi, n)
        for rec in run["per_coordinate"]:
            mine = cens[rec["i"] - 1]["M_i"]
            max_diff_e = max(max_diff_e,
                             abs(mine - rec["mass_weighted_mean_log2_or"]))
    print(f"[C] vs 005_partE (mass-weighted mean log2 OR): max |diff| = "
          f"{max_diff_e:.2e} ({'PASS' if max_diff_e < 1e-7 else 'FAIL'})")
    out = {"max_diff_vs_006_S8": max_diff_s8, "max_diff_vs_005_E": max_diff_e}
    save("or_avg_partC.json", out)
    return out


# ----------------------------------------------------------------------
# part K -- the kill census
# ----------------------------------------------------------------------

LAMS_MAIN = (0.05, 0.31, 0.73, 1.0, 1.37, 2.6)


def census_instance(name: str, mu: dict[int, float], n: int, lam: float,
                    results: list, use_sinkhorn: bool = True,
                    pot: dict[int, float] | None = None) -> list[dict]:
    if pot is not None:
        pi = tilt_of_potential(pot, lam)
        mu = mu_of(pi)
        resid = 0.0
    else:
        pi, resid = sinkhorn(mu, lam)
    marg = marginal_elements(mu, n)
    cens = census(pi, n)
    exc = min_excess(cens, lam)
    exc_lt = min_excess(cens, lam, lt_n=True)
    rec = {
        "name": name, "n": n, "lambda": lam, "max_marginal": max(marg),
        "in_regime": max(marg) < RECORD, "sinkhorn_resid": resid,
        "per_i": [{k: v for k, v in c.items() if k != "rows"} for c in cens],
        "min_excess": exc, "min_excess_lt_n": exc_lt,
    }
    results.append(rec)
    return cens


def part_k() -> dict:
    results: list[dict] = []

    # K1: crash + mirror, several n, off-grid lambda
    for n in (4, 5, 6, 8):
        for lam in LAMS_MAIN:
            census_instance(f"crash_n{n}", crash_mu(n), n, lam, results)
            census_instance(f"mirror_n{n}", mirror_mu(n), n, lam, results)

    # K2: eps-mixed crash -- the i >= 3 histories 006 flagged as zero-mass
    # now carry positive mass; this is the first genuinely new probe.
    for n in (6, 8):
        for lam in (0.31, 1.0, 2.6):
            for eps in (1e-2, 1e-3, 1e-4, 1e-5):
                base = crash_mu(n)
                mu = {m: (1 - eps) * base.get(m, 0.0) + eps / (1 << n)
                      for m in range(1 << n)}
                census_instance(f"crash_epsmix_n{n}_eps{eps:g}", mu, n, lam,
                                results)

    # K3: delta_0 + Bern(1/2+eps) genre (004's recipe killers / half-mixing)
    for n in (6, 7):
        for (pl, ph) in ((0.30, 0.51), (0.25, 0.51), (0.30, 0.55),
                         (0.35, 0.62), (0.50, 0.75)):
            mu = mixture(n, [(pl, 0.0), (1 - pl, ph)])
            for lam in (0.31, 1.0, 2.6):
                census_instance(f"d0_bern_n{n}_p{ph}", mu, n, lam, results)

    # K4: Chase-Lovett slice + smoothed slice (marginals w0/n < 0.38271)
    for (n, w0) in ((8, 3), (7, 2), (6, 2)):
        for lam in (0.31, 1.0, 2.6):
            census_instance(f"slice_n{n}_w{w0}", slice_mu(n, w0), n, lam,
                            results)
            census_instance(f"smslice_n{n}_w{w0}",
                            smoothed_slice_mu(n, w0, 2.0 ** -8), n, lam,
                            results)

    # K5: Sawin-genre geometric mixtures (high components, geometric weights)
    saw = [
        ("sawin_geo_a", [(0.80, 0.28), (0.12, 0.55), (0.06, 0.75), (0.02, 0.92)]),
        ("sawin_geo_b", [(0.70, 0.32), (0.20, 0.50), (0.08, 0.70), (0.02, 0.90)]),
        ("sawin_geo_c", [(0.60, 0.30), (0.25, 0.45), (0.10, 0.60), (0.05, 0.80)]),
        ("partE_mix", [(0.7, 0.30), (0.3, 0.05)]),
    ]
    for name, comps in saw:
        for n in (6,):
            mu = mixture(n, comps)
            for lam in LAMS_MAIN:
                census_instance(f"{name}_n{n}", mu, n, lam, results)

    # K6: random sparse supports (sizes 3..10 -- above the 4-atom theorem),
    # potential trick, marginals recorded; off-grid lambdas
    rng = random.Random(70002)
    for _ in range(300):
        n = rng.choice([4, 5, 6])
        k = rng.randrange(3, 11)
        atoms = rng.sample(range(1 << n), min(k, 1 << n))
        pot = {a: math.exp(rng.uniform(-2.5, 2.5)) for a in atoms}
        lam = rng.choice(LAMS_MAIN)
        census_instance(f"rand_sparse_n{n}_k{k}", {}, n, lam, results,
                        pot=pot)

    # summary.  Violations counted over ALL i (a kill anywhere is a kill);
    # margins reported over i < n, where the i = n identity M_n = lambda
    # does not floor the minimum at 0.
    viol = [r for r in results if r["min_excess"] is not None
            and r["min_excess"] < -1e-9]
    with_lt = [r for r in results if r["min_excess_lt_n"] is not None]
    worst_rec = min(with_lt, key=lambda r: r["min_excess_lt_n"])
    n_inst = len(results)
    print(f"[K] census: {n_inst} instances "
          f"({sum(1 for r in results if r['in_regime'])} in-regime), "
          f"violations of M_i >= lambda (any i, tol 1e-9): {len(viol)}")
    print(f"[K] global min over i<n of (M_i - lambda) = "
          f"{worst_rec['min_excess_lt_n']:+.6f} at {worst_rec['name']} "
          f"lam={worst_rec['lambda']}")
    # closest approaches among structured families (i < n margin)
    by_family: dict[str, float] = {}
    for r in with_lt:
        fam = r["name"].split("_n")[0]
        by_family[fam] = min(by_family.get(fam, float("inf")),
                             r["min_excess_lt_n"])
    for fam in sorted(by_family, key=lambda f: by_family[f])[:10]:
        print(f"[K]   closest approach {fam:22s}: min i<n excess "
              f"{by_family[fam]:+.6f}")
    save("or_avg_partK.json", {"n_instances": n_inst,
                               "n_violations": len(viol),
                               "global_min_excess_lt_n": worst_rec["min_excess_lt_n"],
                               "worst_instance": {k: v for k, v in worst_rec.items()},
                               "by_family_min_excess_lt_n": by_family,
                               "results": results})
    return {"results": results, "violations": viol}


# ----------------------------------------------------------------------
# part L -- L-matrix probe on census + random instances
# ----------------------------------------------------------------------

def part_l() -> dict:
    rng = random.Random(70003)
    min_eig_global = float("inf")
    min_eig_where = None
    max_frob_err = 0.0
    neg_eig_instances = []
    n_probed = 0

    def probe(name, pi, n, lam):
        nonlocal min_eig_global, min_eig_where, max_frob_err, n_probed
        for c in census(pi, n):
            if c["M_i"] is None or c["n_active"] < 2:
                continue
            ev, ferr = l_matrix_probe(c, pi, lam)
            n_probed += 1
            max_frob_err = max(max_frob_err, ferr)
            if ev < min_eig_global:
                min_eig_global = ev
                min_eig_where = (name, lam, c["i"])
            if ev < -1e-9:
                neg_eig_instances.append(
                    {"name": name, "lambda": lam, "i": c["i"],
                     "min_eig": ev, "M_i_minus_lam": c["M_i"] - lam})

    # structured instances
    for n in (5, 6):
        for lam in (0.31, 1.0, 2.6):
            pi, _ = sinkhorn(crash_mu(n), lam)
            probe(f"crash_n{n}", pi, n, lam)
            pi, _ = sinkhorn(mirror_mu(n), lam)
            probe(f"mirror_n{n}", pi, n, lam)
            base = crash_mu(n)
            eps = 1e-3
            mu = {m: (1 - eps) * base.get(m, 0.0) + eps / (1 << n)
                  for m in range(1 << n)}
            pi, _ = sinkhorn(mu, lam)
            probe(f"crash_epsmix_n{n}", pi, n, lam)
            pi, _ = sinkhorn(mixture(n, [(0.3, 0.0), (0.7, 0.51)]), lam)
            probe(f"d0_bern_n{n}", pi, n, lam)
    # random sparse
    for _ in range(300):
        n = rng.choice([5, 6])
        k = rng.randrange(4, 12)
        atoms = rng.sample(range(1 << n), k)
        pot = {a: math.exp(rng.uniform(-2.5, 2.5)) for a in atoms}
        lam = rng.choice((0.31, 1.0, 2.6, 6.0))
        pi = tilt_of_potential(pot, lam)
        probe(f"rand_n{n}_k{k}", pi, n, lam)

    print(f"[L] L-matrix probe: {n_probed} (instance, i) pairs; Frobenius "
          f"identity max err = {max_frob_err:.2e} "
          f"({'PASS' if max_frob_err < 1e-8 else 'FAIL'})")
    print(f"[L] min eigenvalue of L over probe = {min_eig_global:+.6f} at "
          f"{min_eig_where}; instances with L not PSD: {len(neg_eig_instances)}")
    if neg_eig_instances:
        worst_m = min(r["M_i_minus_lam"] for r in neg_eig_instances)
        print(f"[L]   ...yet min (M_i - lambda) among non-PSD instances = "
              f"{worst_m:+.6f} (mass weighting rescues the average)")
    out = {"n_probed": n_probed, "max_frob_err": max_frob_err,
           "global_min_eig": min_eig_global,
           "min_eig_where": min_eig_where,
           "n_non_psd": len(neg_eig_instances),
           "non_psd_instances": neg_eig_instances[:50]}
    save("or_avg_partL.json", out)
    return out


# ----------------------------------------------------------------------
# part S -- adversarial searches
# ----------------------------------------------------------------------

def eval_potential(pot: dict[int, float], n: int, lam: float):
    """Returns (min i<n excess, max_marginal, census).  i = n is excluded
    from the objective: M_n = lambda identically (Prop 1) floors the all-i
    minimum at ~0 and would blind the search."""
    pi = tilt_of_potential(pot, lam)
    mu = mu_of(pi)
    cens = census(pi, n)
    return min_excess(cens, lam, lt_n=True), max(marginal_elements(mu, n)), cens


def part_s() -> dict:
    rng = random.Random(70004)

    # S-a: minimize min_{i<n} (M_i - lambda) with marginal penalty.
    # Small lambdas (0.05, 0.15) included deliberately: the first-order
    # coefficient is a perfect square ||E[D_a]||^2, so a small-lambda kill
    # would need E[D] ~ 0 with a negative second order -- let the search try.
    best = []
    for restart in range(36):
        n = rng.choice([5, 6])
        k = rng.randrange(5, 13)          # sizes above the 4-atom theorem
        atoms = rng.sample(range(1 << n), k)
        pot = {a: math.exp(rng.uniform(-1.5, 1.5)) for a in atoms}
        lam = rng.choice((0.05, 0.15, 0.31, 0.73, 1.0, 2.6))

        def obj(p):
            exc, mm, _ = eval_potential(p, n, lam)
            if exc is None:
                return float("inf"), exc, mm
            return exc + 20.0 * max(0.0, mm - RECORD), exc, mm

        cur, cur_exc, cur_mm = obj(pot)
        for step in range(400):
            sigma = 1.2 * (1 - step / 400) + 0.05
            cand = dict(pot)
            a = rng.choice(list(cand))
            cand[a] *= math.exp(rng.gauss(0, sigma))
            v, exc, mm = obj(cand)
            if v < cur:
                pot, cur, cur_exc, cur_mm = cand, v, exc, mm
        best.append({"restart": restart, "n": n, "k": k, "lambda": lam,
                     "min_excess_lt_n": cur_exc, "max_marginal": cur_mm,
                     "objective": cur})
    worst = min(best, key=lambda b: b["min_excess_lt_n"])
    print(f"[S-a] hill-climb on min_(i<n) (M_i - lambda), 36 restarts x 400 "
          f"steps (5..12 atoms, marginal-penalized):")
    print(f"[S-a] best (most negative) i<n excess found = "
          f"{worst['min_excess_lt_n']:+.6f} at n={worst['n']} k={worst['k']} "
          f"lam={worst['lambda']} (max marginal {worst['max_marginal']:.4f})")
    approach = sorted(b["min_excess_lt_n"] for b in best)[:5]
    print(f"[S-a] five most negative endpoints: "
          + " ".join(f"{v:+.6f}" for v in approach))

    # S-b: minimize the min eigenvalue of any L_i, then check M_i there
    lbest = []
    for restart in range(12):
        n = 6
        k = rng.randrange(6, 12)
        atoms = rng.sample(range(1 << n), k)
        pot = {a: math.exp(rng.uniform(-1.5, 1.5)) for a in atoms}
        lam = rng.choice((1.0, 2.6, 6.0))

        def lobj(p):
            pi = tilt_of_potential(p, lam)
            cens = census(pi, n)
            ev_min = float("inf")
            exc = min_excess(cens, lam, lt_n=True)
            for c in cens:
                if c["M_i"] is None or c["n_active"] < 2:
                    continue
                ev, _ = l_matrix_probe(c, pi, lam)
                ev_min = min(ev_min, ev)
            return (ev_min if ev_min < float("inf") else 0.0), exc

        cur_ev, cur_exc = lobj(pot)
        for step in range(250):
            sigma = 1.0 * (1 - step / 250) + 0.05
            cand = dict(pot)
            a = rng.choice(list(cand))
            cand[a] *= math.exp(rng.gauss(0, sigma))
            ev, exc = lobj(cand)
            if ev < cur_ev:
                pot, cur_ev, cur_exc = cand, ev, exc
        lbest.append({"restart": restart, "n": n, "k": k, "lambda": lam,
                      "min_L_eig": cur_ev, "min_excess_lt_n_at_opt": cur_exc})
    most_neg = min(lbest, key=lambda b: b["min_L_eig"])
    worst_exc_at_negL = min(lbest, key=lambda b: b["min_excess_lt_n_at_opt"])
    print(f"[S-b] hill-climb on min eig of L, 12 restarts x 250 steps:")
    print(f"[S-b] most negative L eigenvalue found = {most_neg['min_L_eig']:+.4f} "
          f"(lam={most_neg['lambda']}); i<n excess at that optimum = "
          f"{most_neg['min_excess_lt_n_at_opt']:+.6f}")
    print(f"[S-b] min i<n excess across all L-optimized instances = "
          f"{worst_exc_at_negL['min_excess_lt_n_at_opt']:+.6f}")
    out = {"search_a": best, "search_b": lbest,
           "global_min_excess_lt_n_search": worst["min_excess_lt_n"]}
    save("or_avg_partS.json", out)
    return out

# ----------------------------------------------------------------------
# part T -- the HEADLINE: an explicit in-regime counterexample to the
# restated Gap 1 (M_i >= lambda), with its full verification battery.
#
# Construction (found by reducing the |N_i| = 2 case of M_i to free
# slices, then stripping marginal-1 coordinates and adding a fresh-prefix
# dilution block, which leaves M_i at the response coordinate EXACTLY
# invariant while driving all elementwise marginals below the record):
#
#   coords: 1 = block marker, 2..4 = dilution markers, 5 = response i,
#   6..7 = futures.  Potential u (its own marginal's Sinkhorn potential):
#     {1,6}: 3.8027..., {1,7}: 8.0947..., {1,5,7}: 0.004793...,
#     {}: 0.25267..., {6}: 27.286..., {7}: 12.530..., {5,6}: 0.36092...,
#     {2}: 20, {3}: 20, {4}: 20.
#   At lambda = 3.5: M_5 - lambda = -0.122033, all marginals <= 0.3178,
#   H(mu) = 2.72 bits.  A 2-significant-digit rounding still violates.
#
# HISTORY / methodological lesson: the first "witnesses" this cycle's
# hill-climbs produced were FLOAT ARTIFACTS -- potentials whose lightest
# atom (~1e-299) made u(A)u(B) underflow to exact 0 in one diagonal cell,
# silently reclassifying a hugely positive diagonal history as degenerate
# and dragging the average down.  BOTH engines shared the underflow, so
# cross-engine agreement did not catch it; an epsilon-sweep of the light
# atom mass did (the honest family was positive).  The witness below has
# atom weights spanning 4.8e-3..27 -- no underflow anywhere -- and its
# violation is stable under 3% perturbations, 2-digit rounding, and the
# light-slice weight swept from 5e-2 down to 1e-6.
# ----------------------------------------------------------------------

def entropy_of(mu: dict[int, float]) -> float:
    return -sum(w * math.log2(w) for w in mu.values() if w > 0)


MARG_TARGET = 0.375     # search penalty target, comfortably below 0.38271


def rand_slice(rng, dim: int):
    if rng.random() < 0.5:
        v = [0.0] * dim
        for _ in range(rng.randrange(1, 3)):
            v[rng.randrange(dim)] = math.exp(rng.uniform(-2, 2))
        if sum(v) == 0:
            v[rng.randrange(dim)] = 1.0
        return v
    return [math.exp(rng.uniform(-3, 3)) if rng.random() < 0.7 else 0.0
            for _ in range(dim)]


WITNESS_LAM = 3.5
WITNESS_N = 7
WITNESS_I = 5
WITNESS_U = {
    0b0100001: 3.802705744653898,      # {1,6}
    0b1000001: 8.094727225702503,      # {1,7}
    0b1010001: 0.004792968160237178,   # {1,5,7}  (the light slice)
    0b0000000: 0.2526660836353666,     # {}
    0b0100000: 27.28641199464732,      # {6}
    0b1000000: 12.529780963546393,     # {7}
    0b0110000: 0.3609234435898425,     # {5,6}
    0b0000010: 20.0,                   # {2}  dilution
    0b0000100: 20.0,                   # {3}  dilution
    0b0001000: 20.0,                   # {4}  dilution
}
WITNESS_U_TIDY = {k: float(f"{v:.2g}") for k, v in WITNESS_U.items()}


def witness_eval(u: dict[int, float], lam: float,
                 n: int = WITNESS_N, i: int = WITNESS_I):
    """(M_i - lam, max marginal, aggregate over i<n, census, mu)."""
    pi = tilt_of_potential(u, lam)
    mu = mu_of(pi)
    cens = census(pi, n)
    exc = (None if cens[i - 1]["M_i"] is None
           else cens[i - 1]["M_i"] - lam)
    rows = [(c["defined_mass"], c["M_i"]) for c in cens
            if c["M_i"] is not None and c["i"] < n]
    agg = (sum(w * (v - lam) for w, v in rows) / sum(w for w, _v in rows)
           if rows else None)
    return exc, max(marginal_elements(mu, n)), agg, cens, mu


def cell_floor(pi, n: int, i: int):
    """min over nondegenerate histories at coordinate i of
    (smallest cell) / (history mass) -- the underflow-artifact detector."""
    pmask = (1 << (i - 1)) - 1
    bit = 1 << (i - 1)
    tabs: dict[tuple[int, int], list[float]] = {}
    for (A, B), w in pi.items():
        t = tabs.setdefault((A & pmask, B & pmask), [0.0] * 4)
        t[(2 if A & bit else 0) + (1 if B & bit else 0)] += w
    vals = [min(t) / sum(t) for t in tabs.values() if min(t) > 0]
    return min(vals) if vals else None


def part_t() -> dict:
    lam, n, i = WITNESS_LAM, WITNESS_N, WITNESS_I
    exc, mm, agg, cens, mu = witness_eval(WITNESS_U, lam)
    pi = tilt_of_potential(WITNESS_U, lam)
    gamma = cell_floor(pi, n, i)
    print(f"[T] EXPLICIT WITNESS: n={n}, 10 atoms, lambda={lam}, "
          f"response coordinate i={i}")
    print(f"[T]   M_i - lambda = {exc:+.6f}  (VIOLATION)  |  max marginal = "
          f"{mm:.4f} < 0.38271  |  H(mu) = {entropy_of(mu):.4f} bits")
    print(f"[T]   all M_i - lambda: "
          + " ".join("--" if c["M_i"] is None else f"{c['M_i'] - lam:+.4f}"
                     for c in cens))
    print(f"[T]   defined mass at i = {cens[i - 1]['defined_mass']:.4f}; "
          f"cell floor at i = {gamma:.2e} (no underflow: atoms span "
          f"{min(WITNESS_U.values()):.2g}..{max(WITNESS_U.values()):.2g})")
    print(f"[T]   aggregate over i<n at witness = {agg:+.4f} "
          f"(the i-averaged form is NOT killed by this instance)")

    # (a) Sinkhorn round-trip: refit the tilt from mu alone
    pi2, resid = sinkhorn(mu, lam)
    cens2 = census(pi2, n)
    rt = max(abs((c1["M_i"] or 0.0) - (c2["M_i"] or 0.0))
             for c1, c2 in zip(cens, cens2))
    print(f"[T]   (a) Sinkhorn round-trip from mu: resid {resid:.1e}, "
          f"max |M_i diff| = {rt:.2e}")

    # (b) dilution invariance: the same slices without the 3 dilution atoms
    u0 = {k: v for k, v in WITNESS_U.items()
          if k not in (0b0000010, 0b0000100, 0b0001000)}
    exc0, mm0, _a0, _c0, _m0 = witness_eval(u0, lam)
    print(f"[T]   (b) dilution invariance: without dilution atoms "
          f"M_i - lambda = {exc0:+.6f} (diff {abs(exc0 - exc):.1e}), "
          f"max marginal {mm0:.4f} (out of regime -> dilution is what buys "
          f"in-regime, changing M_i by nothing)")

    # (c) perturbation robustness: 20 x 3% multiplicative noise
    rng = random.Random(70007)
    nviol = 0
    least = float("inf")
    for _ in range(20):
        up = {k: v * math.exp(rng.gauss(0, 0.03))
              for k, v in WITNESS_U.items()}
        e2, m2, _a, _c, _m = witness_eval(up, lam)
        if e2 is not None and e2 < -1e-9 and m2 < RECORD:
            nviol += 1
            least = min(least, -e2)
    print(f"[T]   (c) perturbation: 3% noise, 20 trials -> {nviol}/20 still "
          f"in-regime violations (weakest {-least:+.6f})")

    # (d) rounding: 2 significant digits
    exct, mmt, _at, _ct, _mt = witness_eval(WITNESS_U_TIDY, lam)
    print(f"[T]   (d) 2-significant-digit witness: M_i - lambda = "
          f"{exct:+.6f}, max marginal {mmt:.4f}")

    # (e) lambda profile of the SAME mu (Sinkhorn refits), incl. small
    # lambda where Theorem C forces the excess positive
    profile = []
    for lam2 in (0.005, 0.01, 0.02, 0.05, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5,
                 3.0, 3.5, 4.0, 5.0):
        pi3, res3 = sinkhorn(mu, lam2, tol=1e-15)
        c3 = census(pi3, n)
        profile.append({"lambda": lam2,
                        "excess": c3[i - 1]["M_i"] - lam2,
                        "resid": res3})
    print(f"[T]   (e) lambda profile (same mu): "
          + " ".join(f"{r['lambda']:g}:{r['excess']:+.2e}" for r in profile))
    small = [r for r in profile if r["lambda"] <= 0.02]
    print(f"[T]       small-lambda excesses all positive: "
          f"{all(r['excess'] > 0 for r in small)} (Theorem C consistency; "
          f"the violation is second-order in lambda, onset ~0.03)")

    # (f) light-slice weight sweep: violation is a stable limit family
    sweep = []
    for eps in (5e-2, 1e-2, 1e-3, 1e-4, 1e-6):
        ue = dict(WITNESS_U)
        ue[0b1010001] = eps
        e4, m4, _a, _c, _m = witness_eval(ue, lam)
        sweep.append({"light_weight": eps, "excess": e4, "max_marginal": m4})
    print(f"[T]   (f) light-slice weight sweep: "
          + " ".join(f"{r['light_weight']:g}:{r['excess']:+.5f}"
                     for r in sweep))

    out = {"lambda": lam, "n": n, "i": i,
           "potential_u": {format(k, "b").zfill(n): v
                           for k, v in WITNESS_U.items()},
           "mu_atoms": {format(k, "b").zfill(n): v for k, v in mu.items()},
           "marginals": marginal_elements(mu, n),
           "max_marginal": mm, "H_mu_bits": entropy_of(mu),
           "witness_excess": exc, "cell_floor": gamma,
           "defined_mass_at_i": cens[i - 1]["defined_mass"],
           "per_i": [{"i": c["i"], "n_active": c["n_active"],
                      "defined_mass": c["defined_mass"], "M_i": c["M_i"]}
                     for c in cens],
           "aggregate_lt_n": agg,
           "roundtrip": {"resid": resid, "max_Mi_diff": rt},
           "dilution_invariance_diff": abs(exc0 - exc),
           "perturbation": {"n_trials": 20, "n_violating": nviol},
           "tidy_2sig": {"excess": exct, "max_marginal": mmt,
                         "potential_u": {format(k, "b").zfill(n): v
                                         for k, v in WITNESS_U_TIDY.items()}},
           "lambda_profile": profile,
           "light_slice_sweep": sweep}
    save("or_avg_partT.json", out)
    return out


# ----------------------------------------------------------------------
# part U -- underflow-clamped searches: (Ua) rediscover in-regime
# single-i violations independently of the hard-coded witness; (Ub) hunt
# the i-AGGREGATED form (mass-weighted average over defined i < n).
# The clamp zeroes slice entries below 1e-9 of the largest entry, so no
# candidate can reach the underflow channel that produced the earlier
# artifacts; zeroed entries are honest support changes handled by the
# N_i bookkeeping.
# ----------------------------------------------------------------------

def build_two_prefix_dil(F, G, m: int, da: int, db: int, wdil: float,
                         K: int = 3):
    """Two-prefix instance plus K fresh-prefix dilution atoms (singleton
    marker coordinates, response bit 0, no future)."""
    p = da + db + K
    a_mask = (1 << da) - 1
    b_mask = ((1 << db) - 1) << da
    u: dict[int, float] = {}
    for alpha in (0, 1):
        for x in range(1 << m):
            if F[alpha][x] > 0:
                u[a_mask | (alpha << p) | (x << (p + 1))] = F[alpha][x]
            if G[alpha][x] > 0:
                key = b_mask | (alpha << p) | (x << (p + 1))
                u[key] = u.get(key, 0.0) + G[alpha][x]
    if wdil > 0:
        for j in range(K):
            u[1 << (da + db + j)] = wdil
    return u, p + 1 + m, p


def clamp_slices(vs, rel: float = 1e-9):
    mx = max(max(v) for v in vs)
    if mx <= 0:
        return vs
    return [[x if x >= mx * rel else 0.0 for x in v] for v in vs]


def part_u() -> dict:
    rng = random.Random(70008)
    out: dict = {}
    for mode in ("single_i", "aggregate"):
        endpoints = []
        best = None
        for restart in range(20):
            m = rng.choice([2, 3])
            dim = 1 << m
            lam = rng.choice([1.0, 2.0, 3.5, 5.0])
            da, db = rng.choice([(1, 0), (0, 1), (1, 1)])
            vs = [rand_slice(rng, dim) for _ in range(4)]
            wdil = [5.0]

            def obj(vv, wd):
                if sum(vv[0]) + sum(vv[1]) == 0 or sum(vv[2]) + sum(vv[3]) == 0:
                    return float("inf"), None, None
                u, n, p = build_two_prefix_dil((vv[0], vv[1]),
                                               (vv[2], vv[3]),
                                               m, da, db, wd[0])
                pi = tilt_of_potential(u, lam)
                mu = mu_of(pi)
                cens = census(pi, n)
                i = p + 1
                if mode == "single_i":
                    if cens[i - 1]["M_i"] is None:
                        return float("inf"), None, None
                    val = cens[i - 1]["M_i"] - lam
                    dmass = cens[i - 1]["defined_mass"]
                else:
                    rows = [(c["defined_mass"], c["M_i"]) for c in cens
                            if c["M_i"] is not None and c["i"] < n]
                    if not rows:
                        return float("inf"), None, None
                    dmass = sum(w for w, _v in rows)
                    val = sum(w * (v - lam) for w, v in rows) / dmass
                mm = max(marginal_elements(mu, n))
                pen = (25.0 * max(0.0, mm - MARG_TARGET)
                       + 10.0 * max(0.0, 0.05 - dmass)
                       + 2.0 * max(0.0, 0.5 - entropy_of(mu)))
                return val + pen, val, mm

            cur, cval, cmm = obj(vs, wdil)
            if not math.isfinite(cur):
                continue
            for step in range(700):
                sigma = 1.0 * (1 - step / 700) + 0.05
                cand = [v[:] for v in vs]
                wd2 = wdil[:]
                if rng.random() < 0.1:
                    wd2[0] = max(wd2[0], 1e-3) * math.exp(rng.gauss(0, sigma))
                else:
                    vi = rng.randrange(4)
                    xi = rng.randrange(dim)
                    if rng.random() < 0.15:
                        cand[vi][xi] = (0.0 if cand[vi][xi] > 0
                                        else math.exp(rng.uniform(-2, 2)))
                    else:
                        cand[vi][xi] = (max(cand[vi][xi], 1e-300)
                                        * math.exp(rng.gauss(0, sigma)))
                cand = clamp_slices(cand)
                v2, val2, mm2 = obj(cand, wd2)
                if v2 < cur:
                    vs, wdil, cur, cval, cmm = cand, wd2, v2, val2, mm2
            endpoints.append({"restart": restart, "m": m, "lambda": lam,
                              "da": da, "db": db, "value": cval,
                              "max_marginal": cmm, "wdil": wdil[0]})
            if (cval is not None and cval < -1e-9 and cmm < RECORD
                    and (best is None or cval < best["value"])):
                best = {"m": m, "lambda": lam, "da": da, "db": db,
                        "value": cval, "max_marginal": cmm,
                        "slices": [v[:] for v in vs], "wdil": wdil[0]}
        nv = sum(1 for e in endpoints if e["value"] is not None
                 and e["value"] < -1e-9 and e["max_marginal"] < RECORD)
        vals = sorted(e["value"] for e in endpoints if e["value"] is not None)
        print(f"[U] clamped search, mode={mode}: 20 restarts x 700 steps; "
              f"{nv} in-regime violations; three lowest endpoints: "
              + " ".join(f"{v:+.5f}" for v in vals[:3]))
        out[mode] = {"endpoints": endpoints,
                     "n_in_regime_violations": nv,
                     "best": ({k: v for k, v in best.items() if k != "slices"}
                              | {"slices": best["slices"]}) if best else None}
    save("or_avg_partU.json", out)
    return out


def main() -> None:
    part_w()
    part_c()
    part_k()
    part_l()
    part_s()
    part_t()
    part_u()
    print("done; checkpoints in", DATA)


if __name__ == "__main__":
    main()
