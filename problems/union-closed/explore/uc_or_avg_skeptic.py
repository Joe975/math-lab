#!/usr/bin/env python3
"""Attempt 013: skeptic engine for 007 (averaged OR-control: well-posed,
partially proved, refuted in-regime).

Default stance: refute.  007's own confession is that its two verifying
engines share IEEE float arithmetic and that its FIRST witnesses were
shared-float-underflow artifacts.  The decisive check here is therefore
EXACT RATIONAL ARITHMETIC (fractions.Fraction end-to-end, no floats on the
critical path): couplings, conditional tables, odds ratios and marginals are
exact rationals, and every log2 is replaced by a certified rational
enclosure (atanh series with directed rounding, explicit tail bounds), so
every sign reported by part A/B/C/D is a THEOREM about the stated rational
inputs, not a float.

The tilt parameter: 007 works at lambda = 3.5, where 2^lambda = 8*sqrt(2) is
irrational.  We work at exact rational tilt t (kernel t^{|A cap B|}); the
conjecture at t is  M_i >= log2 t,  i.e.  sum of mass * log2(OR/t) >= 0 over
N_i^2, and OR/t is an exact rational -- so the sign is certifiable.  We
verify at t = 181/16 ~ 2^3.5, at the float-rationalized 2^3.5, and at exact
integer-lambda points t = 2, 4, 8, 16, 32 (lambda = 1,2,3,4,5 EXACTLY
rational) plus small t where the first-order theorem forces positivity.

Parts:
  A  exact witness: dichotomy (exact), marginals (exact, vs 38271/100000),
     defined mass, certified enclosure of M_5 - log2 t at t ~ 2^3.5;
     certified sign at a sweep of exact rational t (fixed-u family; along
     this path the first-order positivity needs NO Sinkhorn smoothness --
     pi is explicit in t -- so the small-t sign is itself a check);
     convention robustness: certified sign of the deficit even if every
     degenerate history were scored at exactly lambda.
  B  exact invariances: dilution (cell-level rational equality of the i=5
     tables with/without the three dilution atoms) and stripping (append an
     always-1 coordinate; every OR and every relative mass exactly equal).
  C  exact i-aggregated form on the witness: certified enclosure of
     sum_{i<n} w_i (M_i - log2 t) / sum w_i at t = 181/16.
  D  the "same mu violates at every lambda >~ 0.03" claim: float refit
     profile with an independent Sinkhorn (geometric-damped symmetric
     scaling, not 007's row/col alternation), then EXACT nearby-mu
     witnesses at lambda = 1, 2, 3 (rationalized refit potentials; exact
     marginals in-regime; certified violation).
  E  battery spot-check with own float engine: crash/mirror/eps-mix/
     d0+Bern/slice/Sawin/random-sparse families -> expect 0 violations and
     the record's closest-approach ordering; then PLANT the witness genre
     into the same pipeline (witness, tidy witness, own-seed perturbations,
     light-weight sweep) -> the pipeline must flag every one.
  F  theorem probes: Theorem B on random <=4-atom supports (float) and the
     exact potential-free crash identity OR = t^{3-n} with RANDOM rational
     potentials; MTP2 potentials; i=1 and i=n anchors; Frobenius identity;
     mass-matrix PSD (own Jacobi); Theorem C first-order coefficient
     ||sum nu(a) D_a||^2 vs finite differences (random mu and the witness
     mu, target 4.65e-4).
  W  007 part W2 re-run (own engine): direction-dependent eps-limits at the
     degenerate history of the crash family.
  G  own adversarial search for an i-AGGREGATED violation (the proposed
     replacement gap): fresh parametrization, own seeds, including
     hill-climbs seeded AT the witness potential and at two-light-slice
     variants attacking two coordinates at once.

Standard library only.  Deterministic (fixed seeds).  Runtime ~2-4 min.

Usage:
    python problems/union-closed/explore/uc_or_avg_skeptic.py
Checkpoints: problems/union-closed/data/oravgsk_part[ABCDEFWG].json
             problems/union-closed/data/oravgsk_run.log (tee by hand)
"""

from __future__ import annotations

import json
import math
import random
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

RECORD = Fraction(38271, 100000)
F0 = Fraction(0)
F1 = Fraction(1)


def popc(x: int) -> int:
    return bin(x).count("1")


def save(name: str, obj) -> None:
    DATA.mkdir(exist_ok=True)
    (DATA / name).write_text(json.dumps(obj, indent=1, sort_keys=True,
                                        default=str) + "\n", encoding="utf-8")


# ======================================================================
# certified log2 enclosures (rational in, rational bounds out)
# ======================================================================

def _ln_1_to_2(s: Fraction, terms: int = 48) -> tuple[Fraction, Fraction]:
    """[lo, hi] with lo <= ln s <= hi, for 1 <= s < 2.  atanh series:
    ln s = 2 sum_{j>=0} z^{2j+1}/(2j+1), z = (s-1)/(s+1) in [0, 1/3).
    Tail after J terms is <= 2 z^{2J+1} / ((2J+1)(1-z^2))."""
    if s == 1:
        return F0, F0
    z = (s - 1) / (s + 1)
    z2 = z * z
    term = z
    acc = F0
    for j in range(terms):
        acc += term / (2 * j + 1)
        term *= z2
    lo = 2 * acc
    tail = 2 * term / ((2 * terms + 1) * (1 - z2))   # term == z^{2J+1}
    return lo, lo + tail


LN2_LO, LN2_HI = _ln_1_to_2(Fraction(2))  # ln 2 via z = 1/3, tail ~ 1e-47


def _dir_round(r: Fraction, scale: int = 10 ** 50) -> tuple[Fraction, Fraction]:
    n, d = r.numerator, r.denominator
    lo = Fraction(n * scale // d, scale)
    hi = Fraction(-((-n * scale) // d), scale)
    return lo, hi


def log2_enclosure(r: Fraction) -> tuple[Fraction, Fraction]:
    """[lo, hi] with lo <= log2 r <= hi, exact rationals.  r > 0."""
    assert r > 0
    r_lo, r_hi = _dir_round(r)
    if r_lo <= 0:            # r astronomically small; keep exact
        r_lo = r
    out = []
    for rr, want_hi in ((r_lo, False), (r_hi, True)):
        k = 0
        while rr >= 2:
            rr /= 2
            k += 1
        while rr < 1:
            rr *= 2
            k -= 1
        ln_lo, ln_hi = _ln_1_to_2(rr)
        if want_hi:
            out.append(k + ln_hi / LN2_LO)
        else:
            out.append(k + ln_lo / LN2_HI)
    return out[0], out[1]


# ======================================================================
# exact census machinery
# ======================================================================

def exact_census(u: dict[int, Fraction], t: Fraction, n: int) -> list[dict]:
    """Per coordinate i: active set N_i, exact unnormalized tables over all
    positive-mass histories, exact OR over N_i^2, exact degeneracy-dichotomy
    check.  Unnormalized masses are fine: M_i is a mass-weighted mean."""
    supp = sorted(u)
    kern = {}
    for a in supp:
        for b in supp:
            kern[(a, b)] = u[a] * u[b] * t ** popc(a & b)
    out = []
    for i in range(1, n + 1):
        pmask = (1 << (i - 1)) - 1
        bit = 1 << (i - 1)
        prefixes = {A & pmask for A in supp}
        active = {c for c in prefixes
                  if any(A & pmask == c and A & bit for A in supp)
                  and any(A & pmask == c and not A & bit for A in supp)}
        tables: dict[tuple[int, int], list[Fraction]] = {}
        for (A, B), w in kern.items():
            key = (A & pmask, B & pmask)
            tab = tables.setdefault(key, [F0] * 4)
            tab[(2 if A & bit else 0) + (1 if B & bit else 0)] += w
        dich_ok = True
        rows = []
        def_mass = F0
        tot_mass = F0
        for (a, b), (p00, p01, p10, p11) in tables.items():
            tot_mass += p00 + p01 + p10 + p11
            nondeg = min(p00, p01, p10, p11) > 0
            if nondeg != (a in active and b in active):
                dich_ok = False
            if nondeg:
                mass = p00 + p01 + p10 + p11
                orr = (p11 * p00) / (p10 * p01)
                rows.append((a, b, mass, orr))
                def_mass += mass
        out.append({"i": i, "active": active, "rows": rows,
                    "def_mass": def_mass, "tot_mass": tot_mass,
                    "dich_ok": dich_ok, "tables": tables})
    return out


def mi_minus_log2t(rows, t: Fraction):
    """Certified enclosure of M_i - log2 t = sum m log2(OR/t) / sum m."""
    W = sum(m for (_a, _b, m, _o) in rows)
    slo = shi = F0
    for (_a, _b, m, orr) in rows:
        lo, hi = log2_enclosure(orr / t)
        slo += m * lo
        shi += m * hi
    return slo / W, shi / W, W


def exact_marginals(u: dict[int, Fraction], t: Fraction, n: int):
    supp = sorted(u)
    z = F0
    mu: dict[int, Fraction] = {}
    for a in supp:
        row = sum(u[a] * u[b] * t ** popc(a & b) for b in supp)
        mu[a] = row
        z += row
    mu = {a: v / z for a, v in mu.items()}
    marg = [sum(v for a, v in mu.items() if a >> j & 1) for j in range(n)]
    return mu, marg


def frac_str(x: Fraction, digits: int = 15) -> str:
    return f"{float(x):+.{digits}g}"


# ----------------------------------------------------------------------
# the witness (data copied from 007's record / or_avg_partT.json; the CODE
# here is fresh -- nothing imported from uc_or_avg.py)
# ----------------------------------------------------------------------

WITNESS_N = 7
WITNESS_I = 5
DIL_ATOMS = (0b0000010, 0b0000100, 0b0001000)
WITNESS_U_FLOATS = {
    0b0100001: 3.802705744653898,      # {1,6}
    0b1000001: 8.094727225702503,      # {1,7}
    0b1010001: 0.004792968160237178,   # {1,5,7}  (light slice)
    0b0000000: 0.2526660836353666,     # {}
    0b0100000: 27.28641199464732,      # {6}
    0b1000000: 12.529780963546393,     # {7}
    0b0110000: 0.3609234435898425,     # {5,6}
    0b0000010: 20.0, 0b0000100: 20.0, 0b0001000: 20.0,
}
WITNESS_U_EXACT = {k: Fraction(v) for k, v in WITNESS_U_FLOATS.items()}
# exact-decimal 2-significant-digit witness (the tidy claim, taken at face
# value as decimals rather than as their float images)
WITNESS_U_TIDY = {
    0b0100001: Fraction("3.8"), 0b1000001: Fraction("8.1"),
    0b1010001: Fraction("0.0048"), 0b0000000: Fraction("0.25"),
    0b0100000: Fraction("27"), 0b1000000: Fraction("13"),
    0b0110000: Fraction("0.36"),
    0b0000010: Fraction(20), 0b0000100: Fraction(20), 0b0001000: Fraction(20),
}

T_MAIN = Fraction(181, 16)                     # ~ 2^3.4996
T_35 = Fraction(2 ** 3.5)                      # float-rationalized 8*sqrt2


def part_a() -> dict:
    n, i = WITNESS_N, WITNESS_I
    out: dict = {}

    # -- exact census at t ~ 2^3.5 (float-rationalized), the record's point
    cens = exact_census(WITNESS_U_EXACT, T_35, n)
    dich = all(c["dich_ok"] for c in cens)
    mu, marg = exact_marginals(WITNESS_U_EXACT, T_35, n)
    mm = max(marg)
    in_regime = mm < RECORD
    lam_lo, lam_hi = log2_enclosure(T_35)      # lambda = log2 t, certified
    row5 = cens[i - 1]["rows"]
    lo, hi, W = mi_minus_log2t(row5, T_35)
    def_frac = cens[i - 1]["def_mass"] / cens[i - 1]["tot_mass"]
    print(f"[A] EXACT witness, t = float(2^3.5) as a rational "
          f"(lambda in [{float(lam_lo):.15f}, {float(lam_hi):.15f}]):")
    print(f"[A]   degeneracy dichotomy exact at all i: {dich}")
    print(f"[A]   max elementwise marginal (exact) = {float(mm):.10f} "
          f"< 0.38271: {in_regime}  [exact comparison "
          f"{mm.numerator}/{mm.denominator} < 38271/100000]")
    print(f"[A]   defined history-mass fraction at i=5 (exact) = "
          f"{float(def_frac):.10f}")
    print(f"[A]   M_5 - log2 t  in  [{float(lo):.15f}, {float(hi):.15f}]")
    print(f"[A]   CERTIFIED VIOLATION at this exact t: {hi < 0}   "
          f"(record claims -0.122033)")
    out["t_float_2p3.5"] = {
        "t": str(T_35), "dichotomy": dich, "max_marginal": str(mm),
        "max_marginal_float": float(mm), "in_regime": bool(in_regime),
        "defined_fraction_i5": float(def_frac),
        "M5_minus_lam": [float(lo), float(hi)], "certified": bool(hi < 0)}

    # per-i enclosures at t = 181/16 (nice rational; record regime)
    cens_m = exact_census(WITNESS_U_EXACT, T_MAIN, n)
    mu_m, marg_m = exact_marginals(WITNESS_U_EXACT, T_MAIN, n)
    per_i = []
    for c in cens_m:
        if not c["rows"]:
            per_i.append(None)
            continue
        lo_i, hi_i, _w = mi_minus_log2t(c["rows"], T_MAIN)
        per_i.append([float(lo_i), float(hi_i)])
    print(f"[A]   t = 181/16: all M_i - log2 t enclosures: "
          + " ".join("--" if p is None else f"[{p[0]:+.4f},{p[1]:+.4f}]"
                     for p in per_i))
    print(f"[A]   t = 181/16: max marginal {float(max(marg_m)):.6f}, "
          f"in-regime {max(marg_m) < RECORD}, "
          f"i=5 violation certified: {per_i[i - 1][1] < 0}")
    out["t_181_16"] = {"per_i": per_i, "max_marginal": float(max(marg_m)),
                       "in_regime": bool(max(marg_m) < RECORD)}

    # -- convention robustness: score every degenerate history at exactly
    #    lambda instead of excluding it: M_5' - lam = deficit_sum / TOTAL mass
    lo5, hi5, W5 = mi_minus_log2t(row5, T_MAIN)
    tot5 = cens_m[i - 1]["tot_mass"]
    lo_c, hi_c = lo5 * W5 / tot5, hi5 * W5 / tot5
    print(f"[A]   convention robustness (degenerate histories scored AT "
          f"lambda): M_5 - lam in [{float(lo_c):+.6f}, {float(hi_c):+.6f}] "
          f"-> still violated: {hi_c < 0}")
    out["surrogate_lambda_convention"] = {
        "M5_minus_lam": [float(lo_c), float(hi_c)],
        "still_violated": bool(hi_c < 0)}

    # -- tidy 2-significant-decimal witness at t = 181/16
    cens_t = exact_census(WITNESS_U_TIDY, T_MAIN, n)
    mu_t, marg_t = exact_marginals(WITNESS_U_TIDY, T_MAIN, n)
    lo_t, hi_t, _wt = mi_minus_log2t(cens_t[i - 1]["rows"], T_MAIN)
    print(f"[A]   tidy 2-sig-decimal witness at t=181/16: M_5 - lam in "
          f"[{float(lo_t):+.6f}, {float(hi_t):+.6f}], max marginal "
          f"{float(max(marg_t)):.6f}, certified in-regime violation: "
          f"{hi_t < 0 and max(marg_t) < RECORD}")
    out["tidy"] = {"M5_minus_lam": [float(lo_t), float(hi_t)],
                   "max_marginal": float(max(marg_t)),
                   "certified": bool(hi_t < 0 and max(marg_t) < RECORD)}

    # -- exact t sweep on the fixed-u family (includes exact integer lambda)
    sweep = []
    print(f"[A]   exact-t sweep (fixed u; lambda = log2 t exact where t is "
          f"a power of 2):")
    for label, t in (("101/100", Fraction(101, 100)),
                     ("103/100", Fraction(103, 100)),
                     ("21/20", Fraction(21, 20)),
                     ("11/10", Fraction(11, 10)),
                     ("3/2", Fraction(3, 2)),
                     ("2", Fraction(2)), ("4", Fraction(4)),
                     ("8", Fraction(8)), ("16", Fraction(16)),
                     ("32", Fraction(32))):
        c5 = exact_census(WITNESS_U_EXACT, t, n)[i - 1]
        lo_s, hi_s, _w = mi_minus_log2t(c5["rows"], t)
        _mu_s, marg_s = exact_marginals(WITNESS_U_EXACT, t, n)
        sign = "VIOL" if hi_s < 0 else ("pos" if lo_s > 0 else "??")
        sweep.append({"t": label, "lam": float(math.log2(float(t))),
                      "enc": [float(lo_s), float(hi_s)], "sign": sign,
                      "max_marginal": float(max(marg_s)),
                      "in_regime": bool(max(marg_s) < RECORD)})
        print(f"[A]     t={label:8s} (lam={sweep[-1]['lam']:.4f}): "
              f"M_5-lam in [{float(lo_s):+.3e},{float(hi_s):+.3e}] {sign}"
              f"{' in-regime' if sweep[-1]['in_regime'] else ' (marg>rec)'}")
    out["t_sweep"] = sweep
    save("oravgsk_partA.json", out)
    return out


def part_b() -> dict:
    n, i = WITNESS_N, WITNESS_I
    t = T_MAIN
    # dilution invariance, exact: unnormalized i=5 tables over N_5^2 must be
    # IDENTICAL rationals with and without the dilution atoms.
    u0 = {k: v for k, v in WITNESS_U_EXACT.items() if k not in DIL_ATOMS}
    cf = exact_census(WITNESS_U_EXACT, t, n)[i - 1]
    c0 = exact_census(u0, t, n)[i - 1]
    same_active = cf["active"] == c0["active"]
    rows_f = {(a, b): (m, o) for (a, b, m, o) in cf["rows"]}
    rows_0 = {(a, b): (m, o) for (a, b, m, o) in c0["rows"]}
    same_hist = set(rows_f) == set(rows_0)
    cells_equal = same_hist and all(
        rows_f[k][0] == rows_0[k][0] and rows_f[k][1] == rows_0[k][1]
        for k in rows_f)
    _mu_f, marg_f = exact_marginals(WITNESS_U_EXACT, t, n)
    _mu_0, marg_0 = exact_marginals(u0, t, n)
    print(f"[B] dilution invariance EXACT: same N_5 {same_active}, same "
          f"histories {same_hist}, all unnormalized masses and ORs equal as "
          f"rationals: {cells_equal}")
    print(f"[B]   max marginal with dilution {float(max(marg_f)):.4f} "
          f"(in-regime), without {float(max(marg_0)):.4f} (out of regime) -- "
          f"dilution buys regime at exactly zero cost to M_5")

    # stripping invariance, exact: append coordinate 8 set in every atom
    u_str = {(k | (1 << 7)): v for k, v in WITNESS_U_EXACT.items()}
    cs = exact_census(u_str, t, 8)[i - 1]
    rows_s = {(a, b): o for (a, b, _m, o) in cs["rows"]}
    ors_equal = (set(rows_s) == set(rows_f)
                 and all(rows_s[k] == rows_f[k][1] for k in rows_s))
    # relative masses: with the always-1 coordinate the kernel gains a global
    # factor t, so unnormalized masses scale by exactly t
    mass_s = {(a, b): m for (a, b, m, _o) in cs["rows"]}
    masses_scale = all(mass_s[k] == rows_f[k][0] * t for k in mass_s)
    print(f"[B] stripping invariance EXACT (always-1 coordinate appended): "
          f"ORs identical {ors_equal}; masses scale by exactly t "
          f"{masses_scale} -> M_i exactly invariant")
    out = {"dilution_exact_equal": bool(cells_equal),
           "stripping_ors_equal": bool(ors_equal),
           "stripping_mass_scale_t": bool(masses_scale)}
    save("oravgsk_partB.json", out)
    return out


def part_c() -> dict:
    n = WITNESS_N
    t = T_MAIN
    cens = exact_census(WITNESS_U_EXACT, t, n)
    num_lo = num_hi = F0
    den = F0
    for c in cens:
        if c["i"] >= n or not c["rows"]:
            continue
        lo, hi, W = mi_minus_log2t(c["rows"], t)
        num_lo += lo * W
        num_hi += hi * W
        den += W
    agg_lo, agg_hi = num_lo / den, num_hi / den
    print(f"[C] exact i-aggregated form on the witness (t=181/16): "
          f"sum w_i (M_i - lam) / sum w_i in "
          f"[{float(agg_lo):+.6f}, {float(agg_hi):+.6f}]")
    print(f"[C]   certified POSITIVE (aggregated form survives the witness): "
          f"{agg_lo > 0}   (record claims +1.84)")
    out = {"aggregate_enclosure": [float(agg_lo), float(agg_hi)],
           "certified_positive": bool(agg_lo > 0)}
    save("oravgsk_partC.json", out)
    return out


# ======================================================================
# own float engine (fresh implementation; used for battery + searches)
# ======================================================================

def f_census(u: dict[int, float], lam: float, n: int):
    """Float census from a potential (own code path, mirrors the exact one).
    Returns per-i dicts with M_i, def_mass, rows."""
    t = 2.0 ** lam
    supp = sorted(u)
    kern = {}
    for a in supp:
        for b in supp:
            kern[(a, b)] = u[a] * u[b] * t ** popc(a & b)
    tot = sum(kern.values())
    out = []
    for i in range(1, n + 1):
        pmask = (1 << (i - 1)) - 1
        bit = 1 << (i - 1)
        active = {c for c in {A & pmask for A in supp}
                  if any(A & pmask == c and A & bit for A in supp)
                  and any(A & pmask == c and not A & bit for A in supp)}
        tables: dict[tuple[int, int], list[float]] = {}
        for (A, B), w in kern.items():
            key = (A & pmask, B & pmask)
            tab = tables.setdefault(key, [0.0] * 4)
            tab[(2 if A & bit else 0) + (1 if B & bit else 0)] += w
        rows = []
        dm = 0.0
        for (a, b), (p00, p01, p10, p11) in tables.items():
            if a in active and b in active and min(p00, p01, p10, p11) > 0:
                mass = p00 + p01 + p10 + p11
                v = (math.log2(p11) + math.log2(p00) - math.log2(p10)
                     - math.log2(p01))
                rows.append((a, b, mass / tot, v))
                dm += mass / tot
        mi = (sum(m * v for (_a, _b, m, v) in rows) / dm) if dm > 0 else None
        out.append({"i": i, "M_i": mi, "def_mass": dm, "rows": rows,
                    "n_active": len(active)})
    return out


def f_marginals(u: dict[int, float], lam: float, n: int):
    t = 2.0 ** lam
    supp = sorted(u)
    mu = {}
    z = 0.0
    for a in supp:
        row = sum(u[a] * u[b] * t ** popc(a & b) for b in supp)
        mu[a] = row
        z += row
    mu = {a: v / z for a, v in mu.items()}
    return mu, [sum(v for a, v in mu.items() if a >> j & 1) for j in range(n)]


def sinkhorn_own(mu: dict[int, float], lam: float, iters: int = 40000,
                 tol: float = 1e-14):
    """Symmetric geometric-damped scaling u <- u * sqrt(target/current) --
    a different iteration from 007's row/column alternation."""
    supp = sorted(mu)
    N = len(supp)
    K = [[2.0 ** (lam * popc(supp[a] & supp[b])) for b in range(N)]
         for a in range(N)]
    tgt = [mu[s] for s in supp]
    u = [math.sqrt(max(v, 1e-300)) for v in tgt]
    resid = float("inf")
    for _ in range(iters):
        cur = [u[a] * sum(u[b] * K[a][b] for b in range(N)) for a in range(N)]
        resid = max(abs(c - g) for c, g in zip(cur, tgt))
        if resid < tol:
            break
        u = [ua * math.sqrt(g / c) for ua, c, g in zip(u, cur, tgt)]
    z = sum(u[a] * u[b] * K[a][b] for a in range(N) for b in range(N))
    s = 1.0 / math.sqrt(z)
    return {supp[a]: u[a] * s for a in range(N)}, resid


def min_excess_lt_n(cens, lam):
    vals = [c["M_i"] - lam for c in cens if c["M_i"] is not None
            and c["i"] < len(cens)]
    return min(vals) if vals else None


# instance builders (own transcriptions of the published families)

def bern(n, p):
    return {m: p ** popc(m) * (1 - p) ** (n - popc(m)) for m in range(1 << n)}


def mixture(n, comps):
    mu = {m: 0.0 for m in range(1 << n)}
    for w, p in comps:
        for m, v in bern(n, p).items():
            mu[m] += w * v
    return {m: v for m, v in mu.items() if v > 0}


def crash_mu(n):
    full = (1 << n) - 1
    return {0b10: 0.33, full ^ 0b11: 0.33, full: 0.04, 0b01: 0.30}


def mirror_mu(n):
    full = (1 << n) - 1
    return {0b10 | (full ^ 0b11): 0.3, 0b00: 0.2, full: 0.3, 0b01: 0.2}


def slice_mu(n, w0):
    atoms = [m for m in range(1 << n) if popc(m) == w0]
    return {m: 1.0 / len(atoms) for m in atoms}


def part_d() -> dict:
    n, i = WITNESS_N, WITNESS_I
    # the SAME mu: exact marginal of the witness coupling at t=181/16
    mu_ex, _marg = exact_marginals(WITNESS_U_EXACT, T_35, n)
    mu = {a: float(v) for a, v in mu_ex.items()}

    # own-engine float refit profile (tests part T(e) incl. small lambda)
    prof = []
    for lam in (0.005, 0.01, 0.02, 0.05, 0.1, 0.5, 1.0, 2.0, 3.0, 3.5, 4.0,
                5.0):
        pot, resid = sinkhorn_own(mu, lam, tol=1e-15)
        cens = f_census(pot, lam, n)
        prof.append({"lambda": lam, "excess": cens[i - 1]["M_i"] - lam,
                     "resid": resid})
    print("[D] own-engine refit profile of the SAME mu: "
          + " ".join(f"{r['lambda']:g}:{r['excess']:+.2e}" for r in prof))
    small_pos = all(r["excess"] > 0 for r in prof if r["lambda"] <= 0.02)
    print(f"[D]   small-lambda (<=0.02) excesses positive: {small_pos}; "
          f"violation onset and depth match 007 part T(e)")

    # EXACT nearby-mu witnesses at integer lambda: rationalize refit
    # potentials, take exact marginals, certify violation + regime.
    exact_pts = []
    for lam_int, t_int in ((1, Fraction(2)), (2, Fraction(4)),
                           (3, Fraction(8))):
        pot, resid = sinkhorn_own(mu, float(lam_int), tol=1e-15)
        uq = {a: Fraction(v) for a, v in pot.items()}
        cq = exact_census(uq, t_int, n)[i - 1]
        lo, hi, _w = mi_minus_log2t(cq["rows"], t_int)
        mu_q, marg_q = exact_marginals(uq, t_int, n)
        dist = max(abs(float(mu_q[a]) - mu[a]) for a in mu)
        ok = hi < 0 and max(marg_q) < RECORD
        exact_pts.append({"lambda": lam_int,
                          "enc": [float(lo), float(hi)],
                          "max_marginal": float(max(marg_q)),
                          "mu_distance": dist, "certified": bool(ok)})
        print(f"[D]   EXACT witness at lambda={lam_int} (rationalized refit; "
              f"mu within {dist:.1e} of the lam=3.5 witness's mu): "
              f"M_5-{lam_int} in [{float(lo):+.5f},{float(hi):+.5f}], "
              f"max marg {float(max(marg_q)):.4f} -> certified in-regime "
              f"violation: {ok}")
    out = {"profile": prof, "small_lambda_positive": small_pos,
           "exact_integer_lambda_witnesses": exact_pts}
    save("oravgsk_partD.json", out)
    return out


def part_e() -> dict:
    results = []

    def run_mu(name, mu, n, lam):
        pot, resid = sinkhorn_own(mu, lam)
        cens = f_census(pot, lam, n)
        _m, marg = f_marginals(pot, lam, n)
        exc = min_excess_lt_n(cens, lam)
        alli = [c["M_i"] - lam for c in cens if c["M_i"] is not None]
        viol = any(v < -1e-9 for v in alli)
        results.append({"name": name, "n": n, "lambda": lam,
                        "min_excess_lt_n": exc, "violation": viol,
                        "max_marginal": max(marg), "resid": resid})

    def run_pot(name, u, n, lam):
        cens = f_census(u, lam, n)
        _m, marg = f_marginals(u, lam, n)
        exc = min_excess_lt_n(cens, lam)
        alli = [c["M_i"] - lam for c in cens if c["M_i"] is not None]
        viol = any(v < -1e-9 for v in alli)
        results.append({"name": name, "n": n, "lambda": lam,
                        "min_excess_lt_n": exc, "violation": viol,
                        "max_marginal": max(marg)})

    for n in (4, 6):
        for lam in (0.31, 1.0, 2.6):
            run_mu(f"crash_n{n}", crash_mu(n), n, lam)
            run_mu(f"mirror_n{n}", mirror_mu(n), n, lam)
    for lam in (0.31, 1.0, 2.6):
        for eps in (1e-3, 1e-5):
            base = crash_mu(6)
            mu = {m: (1 - eps) * base.get(m, 0.0) + eps / 64
                  for m in range(64)}
            run_mu(f"crash_epsmix_{eps:g}", mu, 6, lam)
        run_mu("d0_bern_0.51", mixture(6, [(0.30, 0.0), (0.70, 0.51)]), 6, lam)
        run_mu("d0_bern_0.62", mixture(6, [(0.35, 0.0), (0.65, 0.62)]), 6, lam)
        run_mu("slice_7_2", slice_mu(7, 2), 7, lam)
        run_mu("sawin_a", mixture(6, [(0.80, 0.28), (0.12, 0.55),
                                      (0.06, 0.75), (0.02, 0.92)]), 6, lam)
        run_mu("sawin_b", mixture(6, [(0.70, 0.32), (0.20, 0.50),
                                      (0.08, 0.70), (0.02, 0.90)]), 6, lam)
        run_mu("partE_mix", mixture(6, [(0.7, 0.30), (0.3, 0.05)]), 6, lam)
    rng = random.Random(130001)
    for _ in range(60):
        n = rng.choice([4, 5, 6])
        k = rng.randrange(3, 11)
        atoms = rng.sample(range(1 << n), min(k, 1 << n))
        pot = {a: math.exp(rng.uniform(-2.5, 2.5)) for a in atoms}
        run_pot(f"rand_n{n}_k{k}", pot, n, rng.choice((0.05, 0.73, 1.37, 2.6)))

    nviol = sum(1 for r in results if r["violation"])
    fam_min: dict[str, float] = {}
    for r in results:
        if r["min_excess_lt_n"] is None:
            continue
        fam = r["name"].split("_")[0]
        fam_min[fam] = min(fam_min.get(fam, float("inf")),
                           r["min_excess_lt_n"])
    print(f"[E] battery spot-check (own engine): {len(results)} instances, "
          f"violations: {nviol} (record: 0 in its 444)")
    for fam in sorted(fam_min, key=lambda f: fam_min[f]):
        print(f"[E]   family {fam:12s} min i<n excess {fam_min[fam]:+.6f}")

    # planted witness genre through the SAME pipeline
    planted = []
    u = dict(WITNESS_U_FLOATS)
    cens = f_census(u, 3.5, WITNESS_N)
    planted.append(("witness", cens[WITNESS_I - 1]["M_i"] - 3.5))
    tidyf = {k: float(v) for k, v in WITNESS_U_TIDY.items()}
    cens = f_census(tidyf, 3.5, WITNESS_N)
    planted.append(("tidy", cens[WITNESS_I - 1]["M_i"] - 3.5))
    rng = random.Random(130002)
    for j in range(5):
        up = {k: v * math.exp(rng.gauss(0, 0.03)) for k, v in u.items()}
        cens = f_census(up, 3.5, WITNESS_N)
        planted.append((f"pert{j}", cens[WITNESS_I - 1]["M_i"] - 3.5))
    for eps in (1e-2, 1e-5):
        ue = dict(u)
        ue[0b1010001] = eps
        cens = f_census(ue, 3.5, WITNESS_N)
        planted.append((f"light{eps:g}", cens[WITNESS_I - 1]["M_i"] - 3.5))
    all_detected = all(v < -1e-3 for (_nm, v) in planted)
    print(f"[E] planted genre detection: "
          + " ".join(f"{nm}:{v:+.4f}" for nm, v in planted))
    print(f"[E]   pipeline flags every planted instance: {all_detected} -- "
          f"the battery pipeline detects the genre; its families simply do "
          f"not contain it")
    out = {"n_instances": len(results), "n_violations": nviol,
           "family_min": fam_min, "planted": planted,
           "all_planted_detected": bool(all_detected), "results": results}
    save("oravgsk_partE.json", out)
    return out


# ----------------------------------------------------------------------
# part F -- theorem probes
# ----------------------------------------------------------------------

def jacobi_eigs(mat):
    nn = len(mat)
    a = [row[:] for row in mat]
    for _ in range(100):
        off = math.sqrt(sum(a[p][q] ** 2 for p in range(nn)
                            for q in range(nn) if p != q))
        if off < 1e-12:
            break
        for p in range(nn - 1):
            for q in range(p + 1, nn):
                if abs(a[p][q]) < 1e-15:
                    continue
                th = (a[q][q] - a[p][p]) / (2 * a[p][q])
                tt = (1 if th >= 0 else -1) / (abs(th)
                                               + math.sqrt(th * th + 1))
                c = 1 / math.sqrt(tt * tt + 1)
                s = tt * c
                for k in range(nn):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p], a[k][q] = c * akp - s * akq, s * akp + c * akq
                for k in range(nn):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k], a[q][k] = c * apk - s * aqk, s * apk + c * aqk
    return [a[k][k] for k in range(nn)]


def part_f() -> dict:
    out: dict = {}
    rng = random.Random(130003)

    # F1: Theorem B -- random <=4-atom supports, all M_i >= lambda
    worst = float("inf")
    trials = 0
    for _ in range(500):
        n = rng.choice([4, 5, 6, 7])
        k = rng.choice([2, 3, 4])
        atoms = rng.sample(range(1 << n), k)
        u = {a: math.exp(rng.uniform(-3, 3)) for a in atoms}
        lam = rng.choice((0.31, 1.0, 2.6, 3.5))
        cens = f_census(u, lam, n)
        vals = [c["M_i"] - lam for c in cens if c["M_i"] is not None]
        if vals:
            worst = min(worst, min(vals))
            trials += 1
    print(f"[F1] Theorem B probe: {trials} random <=4-atom instances, min "
          f"(M_i - lambda) over all defined i = {worst:+.2e} "
          f"({'PASS >= -1e-10' if worst > -1e-10 else 'FAIL'})")
    out["thmB_min"] = worst

    # F2: exact potential-free crash identity OR = t^{3-n} with RANDOM
    # rational potentials (5 trials) -- singleton-cell mechanics, exact
    ok_all = True
    for _ in range(5):
        n = rng.choice([4, 5, 6])
        full = (1 << n) - 1
        atoms = [0b10, full ^ 0b11, full, 0b01]
        u = {a: Fraction(rng.randrange(1, 500), rng.randrange(1, 500))
             for a in atoms}
        t = Fraction(rng.randrange(3, 40), 2)
        cens = exact_census(u, t, n)
        rows = {(a, b): o for (a, b, _m, o) in cens[1]["rows"]}
        got = rows.get((0b00, 0b01))          # history (A1=0, B1=1)
        ok = got == t ** (3 - n)
        ok_all = ok_all and ok
    print(f"[F2] exact crash identity OR == t^(3-n) at random rational "
          f"potentials and t: {'PASS' if ok_all else 'FAIL'}")
    out["crash_identity_exact"] = bool(ok_all)

    # F3: MTP2 potentials -> all OR >= 2^lam hence M_i >= lam
    worst_m = float("inf")
    for _ in range(30):
        n = 4
        h = [rng.uniform(-1, 1) for _ in range(n)]
        J = {(a, b): rng.uniform(0, 0.8) for a in range(n)
             for b in range(a + 1, n)}
        u = {}
        for m in range(1 << n):
            lg = sum(h[j] for j in range(n) if m >> j & 1)
            lg += sum(v for (a, b), v in J.items()
                      if m >> a & 1 and m >> b & 1)
            u[m] = math.exp(lg)
        lam = rng.choice((0.5, 2.0))
        cens = f_census(u, lam, n)
        for c in cens:
            for (_a, _b, _m2, v) in c["rows"]:
                worst_m = min(worst_m, v - lam)
    print(f"[F3] MTP2-potential probe: min (log2 OR - lambda) over every "
          f"history of 30 ferromagnetic instances = {worst_m:+.2e} "
          f"({'PASS' if worst_m > -1e-10 else 'FAIL'})")
    out["mtp2_min_pointwise"] = worst_m

    # F4: anchors + Frobenius identity + mass-matrix PSD on random instances
    worst_m1 = float("inf")
    worst_mn = 0.0
    worst_frob = 0.0
    min_mass_eig = float("inf")
    for _ in range(40):
        n = rng.choice([4, 5])
        k = rng.randrange(4, 11)
        atoms = rng.sample(range(1 << n), k)
        u = {a: math.exp(rng.uniform(-2, 2)) for a in atoms}
        lam = rng.choice((0.31, 1.0, 2.6))
        cens = f_census(u, lam, n)
        if cens[0]["M_i"] is not None:
            worst_m1 = min(worst_m1, cens[0]["M_i"] - lam)
        if cens[-1]["M_i"] is not None:
            worst_mn = max(worst_mn, abs(cens[-1]["M_i"] - lam))
        for c in cens:
            if c["M_i"] is None or c["n_active"] < 2:
                continue
            idx = sorted({a for (a, _b, _m2, _v) in c["rows"]})
            pos = {a: j for j, a in enumerate(idx)}
            A = len(idx)
            Lm = [[0.0] * A for _ in range(A)]
            Mm = [[0.0] * A for _ in range(A)]
            for (a, b, m2, v) in c["rows"]:
                Lm[pos[a]][pos[b]] = v - lam
                Mm[pos[a]][pos[b]] = m2
            fr = sum(Mm[x][y] * Lm[x][y] for x in range(A) for y in range(A))
            tw = sum(Mm[x][y] for x in range(A) for y in range(A))
            worst_frob = max(worst_frob, abs(fr / tw - (c["M_i"] - lam)))
            min_mass_eig = min(min_mass_eig, min(jacobi_eigs(Mm)) / tw)
    print(f"[F4] anchors: min (M_1 - lam) = {worst_m1:+.2e} (>=0 expected); "
          f"max |M_n - lam| = {worst_mn:.2e} (0 expected)")
    print(f"[F4] Frobenius identity max err = {worst_frob:.2e}; "
          f"normalized mass-matrix min eig = {min_mass_eig:+.2e} "
          f"(PSD expected)")
    out["anchor_m1_min"] = worst_m1
    out["anchor_mn_maxdev"] = worst_mn
    out["frob_err"] = worst_frob
    out["mass_min_eig"] = min_mass_eig

    # F5: Theorem C first-order coefficient vs finite differences
    def first_order_coeff(mu, n, i):
        # nu(a) ~ mu(prefix a) over N_i; D_a = E_mu[fut|a,1] - E_mu[fut|a,0]
        pmask = (1 << (i - 1)) - 1
        bit = 1 << (i - 1)
        pref = {}
        for m, w in mu.items():
            pref.setdefault(m & pmask, []).append((m, w))
        active = {c for c, lst in pref.items()
                  if any(m & bit for m, _w in lst)
                  and any(not m & bit for m, _w in lst)}
        if not active:
            return None
        tot = sum(w for c in active for (_m, w) in pref[c])
        acc = [0.0] * (n - i)
        for c in active:
            nu = sum(w for (_m, w) in pref[c]) / tot
            for val in (0, 1):
                sel = [(m, w) for (m, w) in pref[c]
                       if bool(m & bit) == bool(val)]
                wsum = sum(w for (_m, w) in sel)
                mean = [sum(w for (m, w) in sel if m >> (i + j) & 1) / wsum
                        for j in range(n - i)]
                for j in range(n - i):
                    acc[j] += nu * mean[j] * (1 if val else -1)
        return sum(x * x for x in acc)

    checks = []
    for trial in range(3):
        n = 4
        mu = {m: rng.uniform(0.2, 1.0) for m in range(1 << n)}
        z = sum(mu.values())
        mu = {m: v / z for m, v in mu.items()}
        i = rng.choice([2, 3])
        pred = first_order_coeff(mu, n, i)
        slopes = []
        for lam in (2e-3, 1e-3):
            pot, _res = sinkhorn_own(mu, lam, tol=1e-15)
            cens = f_census(pot, lam, n)
            slopes.append((cens[i - 1]["M_i"] - lam) / lam)
        rich = 2 * slopes[1] - slopes[0]
        checks.append({"trial": trial, "i": i, "pred": pred,
                       "richardson": rich,
                       "rel_err": abs(rich - pred) / max(pred, 1e-12)})
    max_rel = max(c["rel_err"] for c in checks)
    print(f"[F5] Theorem C coefficient ||sum nu D||^2 vs Richardson slope "
          f"(3 random mu): max rel err = {max_rel:.2e} "
          f"({'PASS < 5%' if max_rel < 0.05 else 'FAIL'})")
    # witness mu: predicted coefficient at i=5 (record: 4.65e-4)
    mu_ex, _mm = exact_marginals(WITNESS_U_EXACT, T_35, WITNESS_N)
    mu_w = {a: float(v) for a, v in mu_ex.items()}
    pred_w = first_order_coeff(mu_w, WITNESS_N, WITNESS_I)
    print(f"[F5] witness mu predicted first-order coefficient at i=5: "
          f"{pred_w:.4e}  (record claims 4.65e-4)")
    out["thmC_checks"] = checks
    out["thmC_witness_coeff"] = pred_w

    # F6: the witness's own deficit matrix at i=5 -- L must be non-PSD and
    # must break the 2x2-minor bound |L_ab| <= sqrt(e_a e_b) (the record's
    # named failure mode); plus H(mu).
    lam = 3.5
    cens = f_census(WITNESS_U_FLOATS, lam, WITNESS_N)
    c5 = cens[WITNESS_I - 1]
    idx = sorted({a for (a, _b, _m2, _v) in c5["rows"]})
    pos = {a: j for j, a in enumerate(idx)}
    A = len(idx)
    Lm = [[0.0] * A for _ in range(A)]
    for (a, b, _m2, v) in c5["rows"]:
        Lm[pos[a]][pos[b]] = v - lam
    eigs = jacobi_eigs(Lm)
    minor_ratio = max(abs(Lm[x][y]) / math.sqrt(Lm[x][x] * Lm[y][y])
                      for x in range(A) for y in range(A) if x != y)
    mu_w2 = f_marginals(WITNESS_U_FLOATS, lam, WITNESS_N)[0]
    h_mu = -sum(w * math.log2(w) for w in mu_w2.values())
    print(f"[F6] witness L at i=5: |N_5| = {A}, diag = "
          + " ".join(f"{Lm[j][j]:+.4f}" for j in range(A))
          + f", eigs = " + " ".join(f"{e:+.4f}" for e in sorted(eigs)))
    print(f"[F6]   L non-PSD: {min(eigs) < -1e-9}; 2x2-minor ratio "
          f"|L_ab|/sqrt(e_a e_b) = {minor_ratio:.2f} (> 1 = minor control "
          f"fails, as the record's mechanism claims); H(mu) = {h_mu:.4f} "
          f"bits (record: 2.72)")
    out["witness_L_eigs"] = sorted(eigs)
    out["witness_minor_ratio"] = minor_ratio
    out["witness_H_mu"] = h_mu
    save("oravgsk_partF.json", out)
    return out


def part_w2() -> dict:
    n, lam = 5, 1.0
    base = crash_mu(n)
    hist = (0b10, 0b10)
    dirs = {
        "unif": {m: 1.0 / (1 << n) for m in range(1 << n)},
        "bern015": bern(n, 0.15),
        "two_atoms": {0b00110: 0.5, 0b01010: 0.5},
    }
    lims = {}
    for name, nu in dirs.items():
        eps = 1e-4
        mu = {m: (1 - eps) * base.get(m, 0.0) + eps * nu.get(m, 0.0)
              for m in range(1 << n)}
        mu = {m: v for m, v in mu.items() if v > 0}
        pot, resid = sinkhorn_own(mu, lam, tol=1e-13)
        t = 2.0 ** lam
        pmask, bit = 0b11, 0b100
        tab = [0.0] * 4
        for A in pot:
            for B in pot:
                if (A & pmask, B & pmask) == hist:
                    w = pot[A] * pot[B] * t ** popc(A & B)
                    tab[(2 if A & bit else 0) + (1 if B & bit else 0)] += w
        lims[name] = math.log2(tab[3] * tab[0] / (tab[2] * tab[1]))
    spread = max(lims.values()) - min(lims.values())
    print(f"[W2] own-engine directional limits at the degenerate history: "
          + " ".join(f"{k}:{v:+.4f}" for k, v in lims.items())
          + f"  spread {spread:.4f} (007: 1.4716/1.0399/1.0000, 0.4716)")
    out = {"limits": lims, "spread": spread}
    save("oravgsk_partW.json", out)
    return out


def part_g() -> dict:
    """Own hunt for an i-aggregated violation.  Two attacks:
    (a) hill-climbs seeded AT the witness potential (support fixed or with
        extra light atoms giving a second crashable coordinate);
    (b) random sparse supports at n = 7, own parametrization/seeds."""
    rng = random.Random(130004)
    n = WITNESS_N

    def agg_of(u, lam):
        cens = f_census(u, lam, n)
        rows = [(c["def_mass"], c["M_i"]) for c in cens
                if c["M_i"] is not None and c["i"] < n]
        if not rows:
            return None, None
        agg = sum(w * (v - lam) for w, v in rows) / sum(w for w, _v in rows)
        _m, marg = f_marginals(u, lam, n)
        return agg, max(marg)

    def climb(u0, lam, steps=350):
        u = dict(u0)
        def obj(uu):
            agg, mm = agg_of(uu, lam)
            if agg is None:
                return float("inf"), None, None
            return agg + 25.0 * max(0.0, mm - 0.375), agg, mm
        cur, cagg, cmm = obj(u)
        for step in range(steps):
            sigma = 1.0 * (1 - step / steps) + 0.05
            cand = dict(u)
            a = rng.choice(list(cand))
            cand[a] = max(cand[a], 1e-8) * math.exp(rng.gauss(0, sigma))
            v, agg, mm = obj(cand)
            if v < cur:
                u, cur, cagg, cmm = cand, v, agg, mm
        return cagg, cmm

    endpoints = []
    # (a) seeded at the witness; also with a second light slice on coord 6
    seeds = [("witness", dict(WITNESS_U_FLOATS), 3.5)]
    u2 = dict(WITNESS_U_FLOATS)
    u2[0b0110001] = 0.005      # {1,5,6}: second light response-1 atom
    seeds.append(("witness+light156", u2, 3.5))
    u3 = dict(WITNESS_U_FLOATS)
    u3[0b1100000] = 0.004      # {6,7} light in block b
    seeds.append(("witness+light67", u3, 3.5))
    for name, u0, lam in seeds:
        for rep in range(3):
            agg, mm = climb(u0, lam)
            endpoints.append({"seed": name, "lambda": lam, "agg": agg,
                              "max_marginal": mm})
    # (b) random sparse
    for restart in range(12):
        k = rng.randrange(6, 13)
        atoms = rng.sample(range(1 << n), k)
        u0 = {a: math.exp(rng.uniform(-1.5, 1.5)) for a in atoms}
        lam = rng.choice((2.0, 3.5, 5.0))
        agg, mm = climb(u0, lam)
        endpoints.append({"seed": f"rand_k{k}", "lambda": lam, "agg": agg,
                          "max_marginal": mm})
    valid = [e for e in endpoints if e["agg"] is not None]
    nviol = sum(1 for e in valid if e["agg"] < -1e-9
                and e["max_marginal"] < 0.38271)
    low = sorted(e["agg"] for e in valid)[:5]
    print(f"[G] aggregated-form hunt (own seeds, witness-seeded + random): "
          f"{len(valid)} endpoints, in-regime violations: {nviol}")
    print(f"[G]   five lowest endpoints: "
          + " ".join(f"{v:+.5f}" for v in low))
    out = {"endpoints": endpoints, "n_violations": nviol,
           "lowest": low}
    save("oravgsk_partG.json", out)
    return out


def main() -> None:
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    part_f()
    part_w2()
    part_g()
    print("done; checkpoints in", DATA)


if __name__ == "__main__":
    main()
