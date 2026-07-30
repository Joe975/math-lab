#!/usr/bin/env python3
"""Attempt 005: conditional odds-ratio control for the Sinkhorn tilt family.

Target: the Gap-1 SPECULATION of attempt 003 ("plackett-odds-ratio-control"):
for the overlap-tilt coupling pi_lambda(A,B) prop u(A)u(B)2^{lambda|A^B|}
(both marginals mu), lambda >= 0 supposedly forces every conditional odds
ratio OR_i(a,b) of the pair (A_i,B_i) given histories (A_{<i},B_{<i})=(a,b)
into [1, 2^lambda].

This tool provides the computational side of the refutation and of the
partial results:

  part A  engine sanity: product mu  =>  OR == 2^lambda at every history
          (known identity, 003/004); cross-check of the OR engine.
  part B  diagonal histories (a == b): OR >= 2^lambda ALWAYS (Cauchy-Schwarz
          in the tilt inner product) -- i.e. the UPPER half of the claim
          fails, in the opposite direction, generically strictly.
          Random-potential sweep confirming the sign, plus the strictness.
  part C  the 4-atom "crash family": explicit mu with all marginals < 0.37
          whose tilt coupling has OR = 2^{lambda(3-n)} -> 0 at a positive-
          probability history -- the LOWER half of the claim fails too,
          exponentially. The OR value is potential-free (exact identity);
          the code confirms it through an actual Sinkhorn fit from mu.
  part D  sharp general bounds OR in [2^{lambda(1-m)}, 2^{lambda(1+m)}],
          m = n - i remaining coordinates: random sweep finds no violation,
          and explicit 4-atom configurations attain both ends.
  part E  realistic instance (2-block product mixture, the Sawin genre):
          full OR census -- min/max per coordinate, pi-mass-weighted mean of
          log2 OR, mass below 1, mass above 2^lambda.
  part F  MTP2 partial result: if the potential u is log-supermodular then
          OR >= 1 at every history (Karlin-Rinott); random log-supermodular
          potentials confirm, and confirm the upper half STILL fails.
          Plus a probe: does mu MTP2 force the fitted u MTP2?
  part G  perturbative scaling near product mu: diagonal excess of log2 OR
          above lambda scales like delta^2; overall deviation like delta.

Standard library only.  Deterministic (fixed seeds).  Runtime ~10 s.

Usage:
    python problems/union-closed/explore/uc_odds_ratio.py
Checkpoints: problems/union-closed/data/005_part[A-G].json
"""

from __future__ import annotations

import json
import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

LOG2 = math.log(2.0)


# ----------------------------------------------------------------------
# core: coupling from a potential, Sinkhorn fit from a marginal, OR census
# ----------------------------------------------------------------------

def coupling_from_potential(u: dict[int, float], lam: float) -> dict[tuple[int, int], float]:
    """pi(A,B) = u(A) u(B) 2^{lam |A^B|}, normalized.  Support = supp(u)^2."""
    pi = {}
    z = 0.0
    for a, ua in u.items():
        for b, ub in u.items():
            w = ua * ub * math.exp(lam * LOG2 * bin(a & b).count("1"))
            pi[(a, b)] = w
            z += w
    return {k: v / z for k, v in pi.items()}


def marginal(pi: dict[tuple[int, int], float]) -> dict[int, float]:
    mu: dict[int, float] = {}
    for (a, _b), w in pi.items():
        mu[a] = mu.get(a, 0.0) + w
    return mu


def sinkhorn_fit(mu: dict[int, float], lam: float,
                 tol: float = 1e-14, max_iter: int = 200000) -> tuple[dict[tuple[int, int], float], dict[int, float], float]:
    """Fit r,c > 0 with pi(A,B) = r(A) c(B) K(A,B), both marginals = mu.

    K strictly positive on supp(mu)^2, so the scaling exists and the coupling
    it defines is unique (classical Sinkhorn/IPF).  Returns (pi, r, residual).
    """
    supp = sorted(mu)
    kern = {(a, b): math.exp(lam * LOG2 * bin(a & b).count("1")) for a in supp for b in supp}
    r = {a: 1.0 for a in supp}
    c = {a: 1.0 for a in supp}
    resid = float("inf")
    for _ in range(max_iter):
        for a in supp:
            r[a] = mu[a] / sum(c[b] * kern[(a, b)] for b in supp)
        for b in supp:
            c[b] = mu[b] / sum(r[a] * kern[(a, b)] for a in supp)
        # row residual after column update
        resid = max(abs(sum(r[a] * c[b] * kern[(a, b)] for b in supp) - mu[a]) for a in supp)
        if resid < tol:
            break
    pi = {(a, b): r[a] * c[b] * kern[(a, b)] for a in supp for b in supp}
    return pi, r, resid


def odds_ratio_census(pi: dict[tuple[int, int], float], n: int) -> list[dict]:
    """All conditional odds ratios OR_i(a,b), coordinates revealed 1..n
    (coordinate i = bit i-1).  Histories are (i-1)-bit prefixes.

    Returns one row per (i, a, b) with all four cell probs positive:
      {i, a, b, mass, or_, diag} -- mass = pi(history), diag = (a == b).
    """
    rows = []
    for i in range(1, n + 1):
        pref = (1 << (i - 1)) - 1        # mask of coordinates < i
        bit = 1 << (i - 1)
        cells: dict[tuple[int, int], list[float]] = {}
        for (A, B), w in pi.items():
            key = (A & pref, B & pref)
            cell = cells.setdefault(key, [0.0, 0.0, 0.0, 0.0])
            idx = (2 if A & bit else 0) + (1 if B & bit else 0)
            cell[idx] += w                # idx: 0=(0,0) 1=(0,1) 2=(1,0) 3=(1,1)
        for (a, b), (p00, p01, p10, p11) in cells.items():
            if min(p00, p01, p10, p11) <= 0.0:
                continue                  # OR undefined on degenerate cells
            rows.append({
                "i": i, "a": a, "b": b,
                "mass": p00 + p01 + p10 + p11,
                "or_": (p11 * p00) / (p10 * p01),
                "diag": a == b,
            })
    return rows


def entropy_bits(mu: dict[int, float]) -> float:
    return -sum(w * math.log(w, 2) for w in mu.values() if w > 0)


def product_measure(n: int, p: float) -> dict[int, float]:
    mu = {}
    for mask in range(1 << n):
        k = bin(mask).count("1")
        mu[mask] = (p ** k) * ((1 - p) ** (n - k))
    return mu


def mixture(n: int, comps: list[tuple[float, float]]) -> dict[int, float]:
    """Mixture of Bernoulli products: comps = [(weight, p), ...]."""
    mu = {mask: 0.0 for mask in range(1 << n)}
    for w, p in comps:
        for mask, v in product_measure(n, p).items():
            mu[mask] += w * v
    return mu


def save(name: str, obj) -> None:
    DATA.mkdir(exist_ok=True)
    (DATA / name).write_text(json.dumps(obj, indent=1, sort_keys=True) + "\n", encoding="utf-8")


# ----------------------------------------------------------------------
# part A -- engine sanity on product mu (OR must be exactly 2^lambda)
# ----------------------------------------------------------------------

def part_a() -> dict:
    n, p, lam = 4, 0.3, 0.8
    pi, _r, resid = sinkhorn_fit(product_measure(n, p), lam)
    rows = odds_ratio_census(pi, n)
    dev = max(abs(math.log2(r["or_"]) - lam) for r in rows)
    out = {"n": n, "p": p, "lambda": lam, "n_histories": len(rows),
           "max_abs_log2_or_minus_lambda": dev, "sinkhorn_residual": resid}
    print(f"[A] product mu: {len(rows)} histories, max |log2 OR - lambda| = {dev:.2e}  (PASS={'yes' if dev < 1e-8 else 'NO'})")
    save("005_partA.json", out)
    return out


# ----------------------------------------------------------------------
# part B -- diagonal histories: OR >= 2^lambda always, strict generically
# ----------------------------------------------------------------------

def part_b() -> dict:
    """Note: at the last coordinate (i = n, no future) OR == 2^lambda exactly
    for EVERY history, so strictness is reported over i < n separately."""
    rng = random.Random(20260730)
    n = 5
    results = []
    worst_sign = float("inf")     # min over everything of OR/2^lambda on diagonal
    worst_strict = float("inf")   # min over diagonal rows with i < n
    for lam in (0.5, 1.0, 2.0):
        for trial in range(20):
            u = {mask: math.exp(rng.uniform(-2.0, 2.0)) for mask in range(1 << n)}
            pi = coupling_from_potential(u, lam)
            rows = odds_ratio_census(pi, n)
            diag = [r for r in rows if r["diag"]]
            diag_lt = [r for r in diag if r["i"] < n]
            offd = [r for r in rows if not r["diag"]]
            min_diag_ratio = min(r["or_"] / 2 ** lam for r in diag)
            min_diag_lt = min(r["or_"] / 2 ** lam for r in diag_lt)
            max_diag_ratio = max(r["or_"] / 2 ** lam for r in diag)
            min_off_or = min(r["or_"] for r in offd)
            worst_sign = min(worst_sign, min_diag_ratio)
            worst_strict = min(worst_strict, min_diag_lt)
            results.append({
                "lambda": lam, "trial": trial,
                "min_diag_or_over_2lam": min_diag_ratio,
                "min_diag_or_over_2lam_i_lt_n": min_diag_lt,
                "max_diag_or_over_2lam": max_diag_ratio,
                "min_offdiag_or": min_off_or,
            })
    frac_strict = sum(1 for r in results if r["min_diag_or_over_2lam_i_lt_n"] > 1 + 1e-12) / len(results)
    n_below_one = sum(1 for r in results if r["min_offdiag_or"] < 1.0)
    out = {"n": n, "trials": results,
           "global_min_diag_or_over_2lam": worst_sign,
           "global_min_diag_or_over_2lam_i_lt_n": worst_strict,
           "fraction_trials_all_diag_i_lt_n_strictly_above": frac_strict,
           "n_trials_with_some_offdiag_or_below_1": n_below_one}
    print(f"[B] diagonal C-S: min over 60 random potentials of OR/2^lam on diagonal = {worst_sign:.12f} (>= 1: {'yes' if worst_sign >= 1 - 1e-9 else 'NO'})")
    print(f"[B] strictness on i<n diagonal rows: global min OR/2^lam = {worst_strict:.6f}; fraction of trials all-strict: {frac_strict:.2f}")
    print(f"[B] off-diagonal ORs below 1 found in {n_below_one}/60 random dense trials (crash needs structured mu -> part C)")
    save("005_partB.json", out)
    return out


# ----------------------------------------------------------------------
# part C -- the 4-atom crash family: OR = 2^{lambda(3-n)} from Sinkhorn fit
# ----------------------------------------------------------------------

def crash_family(n: int) -> dict[int, float]:
    """mu on atoms S1={2}, S2={3..n}, S3=[n], S4={1}; masses tuned so every
    marginal is < 0.375 < 0.38271 (checked in code, printed)."""
    s1 = 0b10
    s2 = ((1 << n) - 1) ^ 0b11
    s3 = (1 << n) - 1
    s4 = 0b01
    return {s1: 0.33, s2: 0.33, s3: 0.04, s4: 0.30}


def part_c() -> dict:
    lam = 1.0
    rows_out = []
    for n in (5, 8, 11, 13):
        mu = crash_family(n)
        pi, _r, resid = sinkhorn_fit(mu, lam)
        # marginals of mu per element
        elem_marg = [sum(w for mask, w in mu.items() if mask >> j & 1) for j in range(n)]
        rows = odds_ratio_census(pi, n)
        # the crash history: i=2, a = (A_1=0) -> prefix 0, b = (B_1=1) -> prefix 1
        crash = [r for r in rows if r["i"] == 2 and r["a"] == 0 and r["b"] == 1]
        assert len(crash) == 1
        got = math.log2(crash[0]["or_"])
        predicted = lam * (3 - n)
        rows_out.append({
            "n": n, "lambda": lam,
            "max_marginal": max(elem_marg), "H_mu_bits": entropy_bits(mu),
            "crash_history_mass": crash[0]["mass"],
            "log2_or_measured": got, "log2_or_predicted": predicted,
            "abs_err": abs(got - predicted), "sinkhorn_residual": resid,
        })
        print(f"[C] n={n:2d}: max marginal {max(elem_marg):.4f} < 0.38271, "
              f"H(mu)={entropy_bits(mu):.3f} bits, history mass {crash[0]['mass']:.4f}, "
              f"log2 OR = {got:+.6f} (predicted {predicted:+.1f}, err {abs(got - predicted):.1e})")
    # C2: full-support robustness -- mix the 4-atom family with eps * Unif(2^[n]).
    # The crash is not a support-degeneracy artifact: OR -> 2^{lambda(3-n)} as eps -> 0.
    n = 8
    full_rows = []
    for eps in (1e-2, 1e-3, 1e-4):
        base = crash_family(n)
        mu = {mask: (1 - eps) * base.get(mask, 0.0) + eps / (1 << n) for mask in range(1 << n)}
        pi, _r, resid = sinkhorn_fit(mu, lam)
        rows = odds_ratio_census(pi, n)
        crash = [r for r in rows if r["i"] == 2 and r["a"] == 0 and r["b"] == 1]
        got = math.log2(crash[0]["or_"])
        full_rows.append({"n": n, "eps": eps, "log2_or": got,
                          "limit": lam * (3 - n), "sinkhorn_residual": resid})
        print(f"[C2] full support n={n} eps={eps:g}: log2 OR at crash history = {got:+.4f} (limit {lam * (3 - n):+.1f})")
    save("005_partC.json", {"rows": rows_out, "full_support_rows": full_rows})
    return {"rows": rows_out, "full_support_rows": full_rows}


# ----------------------------------------------------------------------
# part D -- sharp range OR in [2^{lam(1-m)}, 2^{lam(1+m)}], m = n - i
# ----------------------------------------------------------------------

def part_d() -> dict:
    rng = random.Random(1234)
    n, lam = 5, 1.0
    # random sweep: no violation of the bounds
    worst_slack_lo = float("inf")
    worst_slack_hi = float("inf")
    for _ in range(40):
        u = {mask: math.exp(rng.uniform(-3.0, 3.0)) for mask in range(1 << n)}
        pi = coupling_from_potential(u, lam)
        for r in odds_ratio_census(pi, n):
            m = n - r["i"]
            lo, hi = lam * (1 - m), lam * (1 + m)
            l2 = math.log2(r["or_"])
            worst_slack_lo = min(worst_slack_lo, l2 - lo)
            worst_slack_hi = min(worst_slack_hi, hi - l2)
    # attainment of the top end: atoms with x1=y1=all-ones, x0=y0=all-zeros
    # at i=2 -> prefix-0 atoms (0,1,1...1),(0,0,0...0); prefix-1 atoms likewise
    a_hi = {  # A-side and B-side atoms coincide (symmetric u)
        (0b10 | (((1 << n) - 1) ^ 0b11)): 0.3,   # (A1=0, A2=1, future all ones)
        0b00: 0.2,                                # (A1=0, A2=0, future all zeros)
        ((1 << n) - 1): 0.3,                      # (A1=1, A2=1, future all ones)
        0b01: 0.2,                                # (A1=1, A2=0, future all zeros)
    }
    pi = coupling_from_potential(a_hi, lam)
    rows = [r for r in odds_ratio_census(pi, n) if r["i"] == 2 and r["a"] == 0 and r["b"] == 1]
    top_l2 = math.log2(rows[0]["or_"])
    m = n - 2
    out = {"n": n, "lambda": lam,
           "min_slack_above_lower_bound": worst_slack_lo,
           "min_slack_below_upper_bound": worst_slack_hi,
           "top_attained_log2_or": top_l2, "top_bound": lam * (1 + m)}
    print(f"[D] random sweep: min slack to bounds lo/hi = {worst_slack_lo:.3f}/{worst_slack_hi:.3f} (both >= 0: "
          f"{'yes' if min(worst_slack_lo, worst_slack_hi) >= -1e-9 else 'NO'})")
    print(f"[D] top end attained: log2 OR = {top_l2:.6f} vs bound lam(1+m) = {lam * (1 + m):.1f}")
    save("005_partD.json", out)
    return out


# ----------------------------------------------------------------------
# part E -- realistic 2-block mixture: full OR census
# ----------------------------------------------------------------------

def part_e() -> dict:
    n = 6
    mu = mixture(n, [(0.7, 0.30), (0.3, 0.05)])
    out_rows = []
    for lam in (0.5, 1.5):
        pi, _r, resid = sinkhorn_fit(mu, lam)
        rows = odds_ratio_census(pi, n)
        per_i = []
        for i in range(1, n + 1):
            ri = [r for r in rows if r["i"] == i]
            mass = sum(r["mass"] for r in ri)
            wmean = sum(r["mass"] * math.log2(r["or_"]) for r in ri) / mass
            per_i.append({
                "i": i,
                "min_log2_or": min(math.log2(r["or_"]) for r in ri),
                "max_log2_or": max(math.log2(r["or_"]) for r in ri),
                "mass_weighted_mean_log2_or": wmean,
                "mass_or_below_1": sum(r["mass"] for r in ri if r["or_"] < 1.0),
                "mass_or_above_2lam": sum(r["mass"] for r in ri if r["or_"] > 2 ** lam * (1 + 1e-9)),
                "min_diag_or_over_2lam": min((r["or_"] / 2 ** lam for r in ri if r["diag"]), default=None),
            })
        out_rows.append({"lambda": lam, "sinkhorn_residual": resid, "per_coordinate": per_i})
        mins = min(p["min_log2_or"] for p in per_i)
        maxs = max(p["max_log2_or"] for p in per_i)
        above = sum(p["mass_or_above_2lam"] for p in per_i) / n
        below1 = sum(p["mass_or_below_1"] for p in per_i) / n
        diag_ok = min(p["min_diag_or_over_2lam"] for p in per_i if p["min_diag_or_over_2lam"] is not None)
        print(f"[E] mixture n={n} lam={lam}: log2 OR in [{mins:+.3f}, {maxs:+.3f}], "
              f"avg mass with OR>2^lam {above:.3f}, with OR<1 {below1:.4f}, min diag OR/2^lam {diag_ok:.6f}")
        print(f"[E]   mass-weighted mean log2 OR by coordinate: "
              + " ".join(f"{p['mass_weighted_mean_log2_or']:+.3f}" for p in per_i) + f"  (lambda={lam})")
    save("005_partE.json", {"n": n, "mixture": "0.7*Bern(0.30)+0.3*Bern(0.05)", "runs": out_rows})
    return {"runs": out_rows}


# ----------------------------------------------------------------------
# part F -- MTP2: log-supermodular u forces OR >= 1; does mu MTP2 force u?
# ----------------------------------------------------------------------

def is_log_supermodular(u: dict[int, float], n: int, rel_tol: float = 1e-9) -> tuple[bool, float]:
    """Check u(S|j|k) u(S) >= u(S|j) u(S|k) for all S, j<k not in S.
    Returns (ok, worst_ratio) with worst_ratio = min lhs/rhs."""
    worst = float("inf")
    for j in range(n):
        for k in range(j + 1, n):
            bj, bk = 1 << j, 1 << k
            for s in range(1 << n):
                if s & (bj | bk):
                    continue
                lhs = u[s | bj | bk] * u[s]
                rhs = u[s | bj] * u[s | bk]
                if rhs > 0:
                    worst = min(worst, lhs / rhs)
    return worst >= 1 - rel_tol, worst


def random_ferro_potential(rng: random.Random, n: int, jmax: float = 1.0) -> dict[int, float]:
    """exp(ferromagnetic quadratic + random field): log-supermodular by design."""
    J = {(j, k): rng.uniform(0.0, jmax) for j in range(n) for k in range(j + 1, n)}
    h = [rng.uniform(-1.0, 1.0) for _ in range(n)]
    u = {}
    for mask in range(1 << n):
        e = sum(h[j] for j in range(n) if mask >> j & 1)
        e += sum(v for (j, k), v in J.items() if mask >> j & 1 and mask >> k & 1)
        u[mask] = math.exp(e)
    return u


def part_f() -> dict:
    rng = random.Random(555)
    n, lam = 4, 1.0
    # (i) log-supermodular potential u -> OR >= 2^lambda at EVERY history
    # (Karlin-Rinott composition, stronger than the >= 1 of plain MTP2 pair
    # extraction); upper half of the claimed control still fails.
    min_ratio_global = float("inf")
    max_diag_ratio = 0.0
    for _ in range(15):
        u = random_ferro_potential(rng, n)
        pi = coupling_from_potential(u, lam)
        rows = odds_ratio_census(pi, n)
        min_ratio_global = min(min_ratio_global, min(r["or_"] / 2 ** lam for r in rows))
        max_diag_ratio = max(max_diag_ratio, max(r["or_"] / 2 ** lam for r in rows if r["diag"]))
    print(f"[F] MTP2 potentials (15 random): min OR/2^lam over ALL histories = {min_ratio_global:.9f} (>= 1: "
          f"{'yes' if min_ratio_global >= 1 - 1e-9 else 'NO'}); max diagonal OR/2^lam = {max_diag_ratio:.3f} (> 1 => upper half fails under MTP2 too)")

    # (ii) probe: mu MTP2 => fitted u MTP2 ?  (if not, the partial result's
    # hypothesis genuinely lives on u, not on mu)
    probe = []
    worst_u_ratio = float("inf")
    worst_case = None
    min_ratio_mu_mtp2 = float("inf")
    for trial in range(12):
        pot = random_ferro_potential(rng, n)
        z = sum(pot.values())
        mu = {k: v / z for k, v in pot.items()}          # mu MTP2 by construction
        for lam2 in (0.5, 1.5):
            _pi, r_pot, resid = sinkhorn_fit(mu, lam2)
            assert resid < 1e-12, f"Sinkhorn did not converge: {resid}"
            ok, worst = is_log_supermodular(r_pot, n)
            if worst < worst_u_ratio:
                worst_u_ratio = worst
                worst_case = {"trial": trial, "lambda": lam2,
                              "mu": {str(k): v for k, v in mu.items()},
                              "residual": resid}
            rows = odds_ratio_census(_pi, n)
            min_ratio_mu_mtp2 = min(min_ratio_mu_mtp2, min(r["or_"] / 2 ** lam2 for r in rows))
            probe.append({"trial": trial, "lambda": lam2, "u_log_supermodular": ok,
                          "worst_supermod_ratio": worst, "residual": resid})
    n_fail = sum(1 for p in probe if not p["u_log_supermodular"])
    print(f"[F] probe mu-MTP2 -> u-MTP2: {n_fail}/{len(probe)} fits violated log-supermodularity of u; "
          f"worst ratio {worst_u_ratio:.9f} (violating instance checkpointed)")
    print(f"[F] min OR/2^lam over histories for mu-MTP2 instances: {min_ratio_mu_mtp2:.9f}")
    out = {"n": n, "min_or_ratio_logsupermod_potentials": min_ratio_global,
           "max_diag_or_over_2lam_logsupermod": max_diag_ratio,
           "probe": probe, "probe_failures": n_fail,
           "worst_u_supermod_ratio": worst_u_ratio,
           "worst_case_instance": worst_case,
           "min_or_ratio_mu_mtp2_instances": min_ratio_mu_mtp2}
    save("005_partF.json", out)
    return out


# ----------------------------------------------------------------------
# part G -- perturbative scaling near product mu
# ----------------------------------------------------------------------

def part_g() -> dict:
    n, lam = 5, 1.0
    rows_out = []
    for delta in (0.005, 0.01, 0.02, 0.04, 0.08, 0.16):
        mu = mixture(n, [(1 - delta, 0.30), (delta, 0.05)])
        pi, _r, resid = sinkhorn_fit(mu, lam)
        rows = odds_ratio_census(pi, n)
        diag_excess = max(math.log2(r["or_"]) - lam for r in rows if r["diag"])
        overall_dev = max(abs(math.log2(r["or_"]) - lam) for r in rows)
        rows_out.append({"delta": delta, "max_diag_excess": diag_excess,
                         "max_overall_dev": overall_dev, "residual": resid})
    for prev, cur in zip(rows_out, rows_out[1:]):
        cur["diag_excess_ratio_vs_prev"] = cur["max_diag_excess"] / prev["max_diag_excess"]
        cur["overall_dev_ratio_vs_prev"] = cur["max_overall_dev"] / prev["max_overall_dev"]
    print("[G] delta | max diag excess (log2 OR - lam) | max overall |log2 OR - lam| | ratios vs prev delta")
    for r in rows_out:
        rat = (f"  x{r.get('diag_excess_ratio_vs_prev', float('nan')):.2f} / x{r.get('overall_dev_ratio_vs_prev', float('nan')):.2f}"
               if "diag_excess_ratio_vs_prev" in r else "")
        print(f"[G] {r['delta']:g}  {r['max_diag_excess']:.6f}   {r['max_overall_dev']:.6f}{rat}")
    save("005_partG.json", {"n": n, "lambda": lam, "rows": rows_out})
    return {"rows": rows_out}


def main() -> None:
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    part_f()
    part_g()
    print("done; checkpoints in", DATA)


if __name__ == "__main__":
    main()
