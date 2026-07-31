#!/usr/bin/env python3
"""Attempt 006: skeptic engine for 005 (odds-ratio control refuted).

Independent re-implementation of everything numerical in
`explore/uc_odds_ratio.py` / attempt 005, sharing NO code with it:

  - Sinkhorn: column-first alternating scaling (005 is row-first), residual
    measured on BOTH marginals of the final coupling, not the row marginal
    mid-iteration; symmetry of the fitted potentials (r = t*c) is CHECKED,
    not assumed, because 005's reduction takes equal potentials.
  - OR census: separate implementation; agreement is meaningful because the
    coupling defined by the scaling of a fixed positive kernel is unique.
  - The 005 checks are re-run at DIFFERENT n, lambda, seeds, and (new) on
    sparse random supports, which 005's random sweeps never touched.

Parts (S1-S8), one per load-bearing claim of 005:

  S1  the reduction OR_i(a,b) = 2^lambda * R, R the cross-ratio of
      <F,WG> future-slice inner products: verified numerically against the
      census at every history, random dense AND sparse supports.
  S2  Prop 2 (diagonal Cauchy-Schwarz OR >= 2^lambda): random potentials at
      n=6 incl. sparse supports; equality at a constructed proportional-slice
      instance; the lambda = 0 edge case of the "iff".
  S3  Prop 3 (4-atom crash family): own Sinkhorn fit, n in {4,5,8,11},
      lambda in {0.7, 1.0, 1.9} (005 only ran lambda = 1); marginals,
      entropy, crash-history mass; full-support eps-mixes at two lambdas.
  S4  Prop 4 (sharp range [2^{lam(1-m)}, 2^{lam(1+m)}]): random sweep with
      sparse supports at n=6; attainment of both ends re-checked.
  S5  part E (realistic mixture census): independent recomputation, diffed
      against data/005_partE.json to ~1e-9.
  S6  Prop 5 (MTP2): (a) direct random test of the Karlin-Rinott
      marginalization step 005 cites from memory; (b) OR >= 2^lambda for
      random log-supermodular potentials at n=5 (005: n=4); (c) refit of
      005's worst_case_instance with this Sinkhorn -- is mu really MTP2,
      does the fitted potential really break log-supermodularity at 0.9922?
  S7  Prop 6 (log2 OR = lambda + O(delta^2)): own delta sweep on 005's
      family plus a SECOND family with a >1/2 component; log-log slopes.
  S8  (bonus, 005's lead-1 falsifiable step) the mass-weighted mean
      M_i = E[log2 OR_i] on the crash family and its mirror.

Standard library only.  Deterministic (fixed seeds).  Runtime ~1 min.

Usage:
    python problems/union-closed/explore/uc_or_skeptic.py
Checkpoints: problems/union-closed/data/006_partS[1-8].json
"""

from __future__ import annotations

import json
import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

POW2 = math.pow


def popcount(x: int) -> int:
    return bin(x).count("1")


# ----------------------------------------------------------------------
# core, written fresh
# ----------------------------------------------------------------------

def fit_tilt(mu: dict[int, float], lam: float, tol: float = 1e-13,
             iters: int = 100000):
    """Alternating scaling, COLUMN-first.  Returns (pi, u, resid, sym_dev):
    pi the coupling, u the symmetric potential sqrt(r*c), resid the max
    abs marginal error of pi on BOTH margins, sym_dev the max relative
    deviation of r/c from constant (should be ~0: pi must be symmetric)."""
    supp = sorted(mu)
    K = [[POW2(2.0, lam * popcount(a & b)) for b in supp] for a in supp]
    tgt = [mu[a] for a in supp]
    N = len(supp)
    r = [1.0] * N
    c = [1.0] * N
    for _ in range(iters):
        for j in range(N):
            c[j] = tgt[j] / sum(r[i] * K[i][j] for i in range(N))
        for i in range(N):
            r[i] = tgt[i] / sum(c[j] * K[i][j] for j in range(N))
        col_err = max(abs(sum(r[i] * c[j] * K[i][j] for i in range(N)) - tgt[j])
                      for j in range(N))
        if col_err < tol:
            break
    # residuals of the final coupling, both margins
    row_err = max(abs(sum(r[i] * c[j] * K[i][j] for j in range(N)) - tgt[i])
                  for i in range(N))
    col_err = max(abs(sum(r[i] * c[j] * K[i][j] for i in range(N)) - tgt[j])
                  for j in range(N))
    ratios = [r[i] / c[i] for i in range(N)]
    sym_dev = max(ratios) / min(ratios) - 1.0
    pi = {(supp[i], supp[j]): r[i] * c[j] * K[i][j]
          for i in range(N) for j in range(N)}
    u = {supp[i]: math.sqrt(r[i] * c[i]) for i in range(N)}
    return pi, u, max(row_err, col_err), sym_dev


def tilt_from_potential(u: dict[int, float], lam: float):
    """pi(A,B) prop u(A)u(B)2^{lam|A&B|}; by Sinkhorn uniqueness this IS the
    tilt coupling of its own marginal (005's counterexample-builder trick)."""
    z = 0.0
    pi = {}
    for a, ua in u.items():
        for b, ub in u.items():
            v = ua * ub * POW2(2.0, lam * popcount(a & b))
            pi[(a, b)] = v
            z += v
    return {k: v / z for k, v in pi.items()}


def or_census(pi: dict[tuple[int, int], float], n: int):
    """Rows (i, a, b, mass, or, diag) for every nondegenerate history."""
    out = []
    for i in range(1, n + 1):
        mask = (1 << (i - 1)) - 1
        bit = 1 << (i - 1)
        tables: dict[tuple[int, int], list[float]] = {}
        for (A, B), w in pi.items():
            t = tables.setdefault((A & mask, B & mask), [0.0] * 4)
            t[(1 if A & bit else 0) * 2 + (1 if B & bit else 0)] += w
        for (a, b), t in sorted(tables.items()):
            p00, p01, p10, p11 = t
            if p00 > 0 and p01 > 0 and p10 > 0 and p11 > 0:
                out.append((i, a, b, sum(t), p11 * p00 / (p01 * p10), a == b))
    return out


def slice_inner(u: dict[int, float], n: int, i: int, a: int, b: int,
                lam: float):
    """The four <F_alpha, W G_beta> of 005's reduction, computed directly."""
    m = n - i
    bit = 1 << (i - 1)

    def sl(prefix, alpha):
        return [u.get(prefix | (alpha * bit) | (x << i), 0.0)
                for x in range(1 << m)]

    F = [sl(a, 0), sl(a, 1)]
    G = [sl(b, 0), sl(b, 1)]
    inn = [[0.0] * 2 for _ in range(2)]
    for al in (0, 1):
        for be in (0, 1):
            s = 0.0
            for x in range(1 << m):
                fx = F[al][x]
                if fx == 0.0:
                    continue
                for y in range(1 << m):
                    s += fx * G[be][y] * POW2(2.0, lam * popcount(x & y))
            inn[al][be] = s
    return inn


def bin_entropy_bits(mu: dict[int, float]) -> float:
    return -sum(w * math.log2(w) for w in mu.values() if w > 0)


def bern_product(n: int, p: float) -> dict[int, float]:
    return {m: p ** popcount(m) * (1 - p) ** (n - popcount(m))
            for m in range(1 << n)}


def mix(n: int, comps) -> dict[int, float]:
    mu = {m: 0.0 for m in range(1 << n)}
    for w, p in comps:
        for m, v in bern_product(n, p).items():
            mu[m] += w * v
    return mu


def rand_potential(rng, n: int, spread: float, sparse: bool):
    atoms = list(range(1 << n))
    if sparse:
        k = rng.randrange(6, (1 << n) // 2 + 1)
        atoms = rng.sample(atoms, k)
    return {a: math.exp(rng.uniform(-spread, spread)) for a in atoms}


def log_supermod_worst(u: dict[int, float], n: int) -> float:
    """min over pairs of u(S|jk)u(S)/(u(S|j)u(S|k)); >= 1 iff MTP2.
    Only meaningful for full-support u."""
    worst = float("inf")
    for j in range(n):
        for k in range(j + 1, n):
            bj, bk = 1 << j, 1 << k
            for s in range(1 << n):
                if s & (bj | bk):
                    continue
                num = u[s | bj | bk] * u[s]
                den = u[s | bj] * u[s | bk]
                worst = min(worst, num / den)
    return worst


def rand_ferro(rng, n: int, jmax: float) -> dict[int, float]:
    """exp(random field + NONNEG pair couplings): log-supermodular by design.
    Written fresh; jmax differs from 005's to sample a different region."""
    J = [[rng.uniform(0.0, jmax) for _ in range(n)] for _ in range(n)]
    h = [rng.uniform(-1.5, 1.5) for _ in range(n)]
    u = {}
    for m in range(1 << n):
        e = sum(h[j] for j in range(n) if m >> j & 1)
        e += sum(J[j][k] for j in range(n) for k in range(j + 1, n)
                 if m >> j & 1 and m >> k & 1)
        u[m] = math.exp(e)
    return u


def save(name: str, obj) -> None:
    DATA.mkdir(exist_ok=True)
    (DATA / name).write_text(json.dumps(obj, indent=1, sort_keys=True) + "\n",
                             encoding="utf-8")


# ----------------------------------------------------------------------
# S1 -- the reduction itself
# ----------------------------------------------------------------------

def part_s1():
    rng = random.Random(60001)
    worst = 0.0
    checked = 0
    for n in (4, 5, 6):
        for lam in (0.6, 1.7):
            for sparse in (False, True):
                for _ in range(4):
                    u = rand_potential(rng, n, 2.5, sparse)
                    pi = tilt_from_potential(u, lam)
                    for (i, a, b, _mass, orr, _d) in or_census(pi, n):
                        inn = slice_inner(u, n, i, a, b, lam)
                        if min(inn[0][0], inn[0][1], inn[1][0], inn[1][1]) <= 0:
                            continue
                        pred = POW2(2.0, lam) * (inn[1][1] * inn[0][0]) / (inn[1][0] * inn[0][1])
                        worst = max(worst, abs(pred / orr - 1.0))
                        checked += 1
    out = {"histories_checked": checked, "max_rel_err": worst}
    print(f"[S1] reduction OR = 2^lam * cross-ratio: {checked} histories "
          f"(dense+sparse, n=4..6), max rel err {worst:.2e} "
          f"({'PASS' if worst < 1e-9 else 'FAIL'})")
    save("006_partS1.json", out)


# ----------------------------------------------------------------------
# S2 -- Prop 2, diagonal Cauchy-Schwarz
# ----------------------------------------------------------------------

def part_s2():
    rng = random.Random(60002)
    n = 6
    min_lt, min_all = float("inf"), float("inf")
    trials = 0
    for lam in (0.25, 1.0, 3.0):
        for sparse in (False, True):
            for _ in range(10):
                u = rand_potential(rng, n, 2.0, sparse)
                pi = tilt_from_potential(u, lam)
                s = POW2(2.0, lam)
                for (i, a, b, _m, orr, diag) in or_census(pi, n):
                    if not diag:
                        continue
                    min_all = min(min_all, orr / s)
                    if i < n:
                        min_lt = min(min_lt, orr / s)
                trials += 1
    # equality construction: proportional future slices at one prefix
    n2, lam2, i2 = 4, 1.3, 2
    u = {m: math.exp(rng.uniform(-1, 1)) for m in range(1 << n2)}
    bit = 1 << (i2 - 1)
    for x in range(1 << (n2 - i2)):          # force F1 = 3 * F0 at prefix a=0
        u[(x << i2) | bit] = 3.0 * u[x << i2]
    pi = tilt_from_potential(u, lam2)
    eq_rows = [r for r in or_census(pi, n2) if r[0] == i2 and r[1] == 0 and r[2] == 0]
    eq_dev = abs(eq_rows[0][4] / POW2(2.0, lam2) - 1.0)
    # lambda = 0: OR == 1 with NON-proportional slices (the iff fails at 0)
    u0 = rand_potential(rng, 4, 2.0, False)
    pi0 = tilt_from_potential(u0, 0.0)
    lam0_dev = max(abs(r[4] - 1.0) for r in or_census(pi0, 4))
    out = {"n": n, "trials": trials,
           "min_diag_or_over_2lam_i_lt_n": min_lt,
           "min_diag_or_over_2lam_all": min_all,
           "proportional_slice_equality_dev": eq_dev,
           "lambda0_max_abs_or_minus_1": lam0_dev}
    print(f"[S2] diagonal C-S, {trials} random potentials (n=6, dense+sparse): "
          f"min OR/2^lam = {min_all:.12f} all, {min_lt:.12f} on i<n "
          f"({'PASS' if min_all > 1 - 1e-9 else 'FAIL'})")
    print(f"[S2] proportional-slice equality dev {eq_dev:.2e}; "
          f"lambda=0 max|OR-1| = {lam0_dev:.2e} (equality WITHOUT proportional "
          f"slices -> Prop 2's 'iff' needs lambda > 0)")
    save("006_partS2.json", out)


# ----------------------------------------------------------------------
# S3 -- Prop 3, the crash family
# ----------------------------------------------------------------------

def crash_mu(n: int) -> dict[int, float]:
    full = (1 << n) - 1
    return {0b10: 0.33, full ^ 0b11: 0.33, full: 0.04, 0b01: 0.30}


def crash_row(pi, n):
    rows = [r for r in or_census(pi, n) if r[0] == 2 and r[1] == 0 and r[2] == 1]
    assert len(rows) == 1
    return rows[0]


def part_s3():
    rows_out = []
    worst = 0.0
    for n in (4, 5, 8, 11):
        mu = crash_mu(n)
        marg = [sum(w for m, w in mu.items() if m >> j & 1) for j in range(n)]
        for lam in (0.7, 1.0, 1.9):
            pi, _u, resid, sym = fit_tilt(mu, lam)
            _i, _a, _b, mass, orr, _d = crash_row(pi, n)
            err = abs(math.log2(orr) - lam * (3 - n))
            worst = max(worst, err)
            rows_out.append({"n": n, "lambda": lam, "max_marginal": max(marg),
                             "H_bits": bin_entropy_bits(mu), "mass": mass,
                             "log2_or": math.log2(orr),
                             "predicted": lam * (3 - n), "abs_err": err,
                             "resid": resid, "sym_dev": sym})
        print(f"[S3] n={n:2d}: max marginal {max(marg):.4f}, "
              f"H={bin_entropy_bits(mu):.4f} bits, crash mass "
              f"{rows_out[-1]['mass']:.4f}, log2 OR errs at lam .7/1/1.9 "
              f"<= {max(r['abs_err'] for r in rows_out[-3:]):.1e}")
    # full-support eps-mix, two lambdas (005 ran lambda=1 only)
    fs = []
    n = 8
    for lam in (1.0, 1.5):
        for eps in (1e-2, 1e-3, 1e-4):
            base = crash_mu(n)
            mu = {m: (1 - eps) * base.get(m, 0.0) + eps / (1 << n)
                  for m in range(1 << n)}
            pi, _u, resid, _s = fit_tilt(mu, lam)
            _i, _a, _b, mass, orr, _d = crash_row(pi, n)
            fs.append({"n": n, "lambda": lam, "eps": eps,
                       "log2_or": math.log2(orr), "limit": lam * (3 - n),
                       "resid": resid})
            print(f"[S3] full support lam={lam} eps={eps:g}: "
                  f"log2 OR = {math.log2(orr):+.4f} (limit {lam * (3 - n):+.1f})")
    out = {"exact_rows": rows_out, "full_support": fs,
           "max_abs_err_exact": worst}
    print(f"[S3] crash identity log2 OR = lam(3-n): max |err| = {worst:.1e} "
          f"({'PASS' if worst < 1e-9 else 'FAIL'})")
    save("006_partS3.json", out)


# ----------------------------------------------------------------------
# S4 -- Prop 4, sharp range
# ----------------------------------------------------------------------

def part_s4():
    rng = random.Random(60004)
    n, lam = 6, 0.8
    lo_slack = hi_slack = float("inf")
    for sparse in (False, True):
        for _ in range(15):
            u = rand_potential(rng, n, 3.0, sparse)
            pi = tilt_from_potential(u, lam)
            for (i, _a, _b, _m, orr, _d) in or_census(pi, n):
                m = n - i
                l2 = math.log2(orr)
                lo_slack = min(lo_slack, l2 - lam * (1 - m))
                hi_slack = min(hi_slack, lam * (1 + m) - l2)
    # attainment: crash family (lower), mirror family (upper), several n
    att = []
    for n2 in (5, 7):
        full = (1 << n2) - 1
        mu = crash_mu(n2)
        pi, _u, _r, _s = fit_tilt(mu, 1.0)
        lo = math.log2(crash_row(pi, n2)[4])
        mirror = {0b10 | (full ^ 0b11): 0.3, 0b00: 0.2, full: 0.3, 0b01: 0.2}
        pi2 = tilt_from_potential(mirror, 1.0)
        hi = math.log2(crash_row(pi2, n2)[4])
        att.append({"n": n2, "lower_attained": lo, "lower_bound": 1.0 * (3 - n2),
                    "upper_attained": hi, "upper_bound": 1.0 * (n2 - 1)})
    ok = all(abs(a["lower_attained"] - a["lower_bound"]) < 1e-9
             and abs(a["upper_attained"] - a["upper_bound"]) < 1e-9 for a in att)
    out = {"n": n, "lambda": lam, "min_slack_lower": lo_slack,
           "min_slack_upper": hi_slack, "attainment": att}
    print(f"[S4] range sweep (30 potentials, dense+sparse, n=6): min slack "
          f"lo/hi = {lo_slack:.3f}/{hi_slack:.3f} "
          f"({'PASS' if min(lo_slack, hi_slack) > -1e-9 else 'FAIL'}); "
          f"both ends attained at n=5,7: {'PASS' if ok else 'FAIL'}")
    save("006_partS4.json", out)


# ----------------------------------------------------------------------
# S5 -- part E mixture census, diffed against 005's checkpoint
# ----------------------------------------------------------------------

def part_s5():
    n = 6
    mu = mix(n, [(0.7, 0.30), (0.3, 0.05)])
    theirs = json.loads((DATA / "005_partE.json").read_text())
    max_diff = 0.0
    summary = []
    for run in theirs["runs"]:
        lam = run["lambda"]
        pi, _u, resid, _s = fit_tilt(mu, lam)
        rows = or_census(pi, n)
        for rec in run["per_coordinate"]:
            i = rec["i"]
            ri = [r for r in rows if r[0] == i]
            mass = sum(r[3] for r in ri)
            mine = {
                "min_log2_or": min(math.log2(r[4]) for r in ri),
                "max_log2_or": max(math.log2(r[4]) for r in ri),
                "mass_weighted_mean_log2_or":
                    sum(r[3] * math.log2(r[4]) for r in ri) / mass,
                "mass_or_below_1": sum(r[3] for r in ri if r[4] < 1.0),
                "mass_or_above_2lam":
                    sum(r[3] for r in ri if r[4] > POW2(2.0, lam) * (1 + 1e-9)),
            }
            for k, v in mine.items():
                max_diff = max(max_diff, abs(v - rec[k]))
            summary.append({"lambda": lam, "i": i, **mine})
    out = {"max_abs_diff_vs_005": max_diff, "mine": summary}
    print(f"[S5] mixture census vs 005_partE.json: max abs diff over every "
          f"reported field = {max_diff:.2e} ({'PASS' if max_diff < 1e-8 else 'FAIL'})")
    save("006_partS5.json", out)


# ----------------------------------------------------------------------
# S6 -- Prop 5 MTP2: marginalization step, OR bound, Sinkhorn-breaks-MTP2
# ----------------------------------------------------------------------

def marginalize(f: dict[int, float], n: int, keep: list[int]) -> dict[int, float]:
    g: dict[int, float] = {}
    for m, v in f.items():
        key = 0
        for pos, j in enumerate(keep):
            if m >> j & 1:
                key |= 1 << pos
        g[key] = g.get(key, 0.0) + v
    return g


def part_s6():
    rng = random.Random(60006)
    # (a) Karlin-Rinott marginalization on the binary cube, tested directly
    worst_marg = float("inf")
    tests = 0
    for _ in range(120):
        n = rng.choice([5, 6, 7])
        f = rand_ferro(rng, n, 1.5)
        k = rng.randrange(2, n)
        keep = sorted(rng.sample(range(n), k))
        g = marginalize(f, n, keep)
        worst_marg = min(worst_marg, log_supermod_worst(g, k))
        tests += 1
    print(f"[S6a] MTP2 marginalization (Karlin-Rinott Prop 3.2), {tests} random "
          f"ferro f, n=5..7: worst log-supermod ratio of marginal = "
          f"{worst_marg:.12f} ({'PASS' if worst_marg >= 1 - 1e-9 else 'FAIL'})")
    # (b) log-supermodular potential => OR >= 2^lambda at every history
    n = 5
    min_ratio = float("inf")
    for lam in (0.5, 2.0):
        for _ in range(10):
            u = rand_ferro(rng, n, 2.0)
            pi = tilt_from_potential(u, lam)
            s = POW2(2.0, lam)
            min_ratio = min(min_ratio, min(r[4] / s for r in or_census(pi, n)))
    print(f"[S6b] 20 random log-supermodular potentials (n=5, lam 0.5/2.0): "
          f"min OR/2^lam over ALL histories = {min_ratio:.12f} "
          f"({'PASS' if min_ratio >= 1 - 1e-9 else 'FAIL'})")
    # (c) 005's worst_case_instance: refit and re-check
    inst = json.loads((DATA / "005_partF.json").read_text())["worst_case_instance"]
    mu = {int(k): v for k, v in inst["mu"].items()}
    lam = inst["lambda"]
    n4 = 4
    mu_worst = log_supermod_worst(mu, n4)
    pi, u, resid, sym = fit_tilt(mu, lam)
    u_worst = log_supermod_worst(u, n4)
    min_or = min(r[4] / POW2(2.0, lam) for r in or_census(pi, n4))
    out = {"marginalization_tests": tests, "worst_marginal_ratio": worst_marg,
           "min_or_ratio_logsupermod": min_ratio,
           "instance": {"lambda": lam, "mu_is_mtp2": mu_worst >= 1 - 1e-9,
                        "mu_worst_ratio": mu_worst,
                        "u_worst_ratio": u_worst, "claimed": 0.9921520205873654,
                        "resid": resid, "sym_dev": sym,
                        "min_or_over_2lam": min_or}}
    print(f"[S6c] worst_case_instance refit: mu MTP2 ratio {mu_worst:.6f} "
          f"({'is MTP2' if mu_worst >= 1 - 1e-9 else 'NOT MTP2!'}); fitted-u "
          f"worst ratio {u_worst:.9f} (claimed 0.992152021); "
          f"min OR/2^lam {min_or:.9f}")
    save("006_partS6.json", out)


# ----------------------------------------------------------------------
# S7 -- Prop 6 delta^2 scaling, 005's family + a second one
# ----------------------------------------------------------------------

def part_s7():
    out = {}
    for name, mk, deltas in (
        ("005_family_Bern30_Bern05_n5",
         lambda d: mix(5, [(1 - d, 0.30), (d, 0.05)]),
         [0.005, 0.01, 0.02, 0.04, 0.08]),
        # a >1/2 component: larger delta^3 coefficient, so the asymptotic
        # regime needs smaller delta -- the point of the extended range.
        ("second_family_Bern35_Bern70_n5",
         lambda d: mix(5, [(1 - d, 0.35), (d, 0.70)]),
         [0.00015625, 0.0003125, 0.000625, 0.00125, 0.0025, 0.005, 0.02, 0.08]),
    ):
        lam = 1.0
        devs = []
        for d in deltas:
            pi, _u, resid, _s = fit_tilt(mk(d), lam, tol=1e-15)
            dev = max(abs(math.log2(r[4]) - lam) for r in or_census(pi, 5))
            devs.append(dev)
        slopes = [math.log(devs[k + 1] / devs[k]) / math.log(deltas[k + 1] / deltas[k])
                  for k in range(len(devs) - 1)]
        out[name] = {"deltas": deltas, "max_dev": devs, "loglog_slopes": slopes}
        print(f"[S7] {name}: max|log2 OR - lam| = "
              + " ".join(f"{v:.2e}" for v in devs)
              + f"; log-log slopes {' '.join(f'{s:.3f}' for s in slopes)}")
    save("006_partS7.json", out)


# ----------------------------------------------------------------------
# S8 -- 005 lead 1's falsifiable step: M_i on crash family and mirror
# ----------------------------------------------------------------------

def part_s8():
    out = []
    for name, mk in (("crash", crash_mu),
                     ("mirror", lambda n: {0b10 | (((1 << n) - 1) ^ 0b11): 0.3,
                                           0b00: 0.2, (1 << n) - 1: 0.3,
                                           0b01: 0.2})):
        for n in (5, 8):
            for lam in (1.0,):
                mu = mk(n)
                pi, _u, _r, _s = fit_tilt(mu, lam)
                rows = or_census(pi, n)
                per_i = []
                for i in range(1, n + 1):
                    ri = [r for r in rows if r[0] == i]
                    mass = sum(r[3] for r in ri)
                    mi = (sum(r[3] * math.log2(r[4]) for r in ri) / mass
                          if mass > 0 else None)
                    per_i.append({"i": i, "defined_mass": mass, "M_i": mi})
                out.append({"family": name, "n": n, "lambda": lam,
                            "per_i": per_i})
                bad = [p for p in per_i if p["M_i"] is not None
                       and p["M_i"] < lam - 1e-9]
                print(f"[S8] {name} n={n} lam={lam}: M_i = "
                      + " ".join("--" if p["M_i"] is None else f"{p['M_i']:+.3f}"
                                 for p in per_i)
                      + ("  << M_i < lambda at i=" + ",".join(str(p["i"]) for p in bad)
                         if bad else "  (all M_i >= lambda)"))
    save("006_partS8.json", out)


def main() -> None:
    part_s1()
    part_s2()
    part_s3()
    part_s4()
    part_s5()
    part_s6()
    part_s7()
    part_s8()
    print("done; checkpoints in", DATA)


if __name__ == "__main__":
    main()
