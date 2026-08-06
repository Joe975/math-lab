#!/usr/bin/env python3
"""Attempt 018: run the MU(n, r) ladder adversary against BOTH surviving
Gap-1 candidates (016 leads 1 and 4; STATUS queue items 1 and 20).

After 016/017 killed the plain i-aggregated odds-ratio control at large n,
two restatements survived.  Each gets a falsifiable probe here, against the
one adversary genre known to kill at scale (unit replication with
adversarial shared weights), with the adversary allowed to re-optimize its
weights against the NEW statement (a positive without that step would be
worthless).

Candidate I -- margin-modulated control (007 lead 2).  The chain-rule cost
of an OR-deviation at a history (a, b) is the change it makes to the
history's Plackett gain h(z_rho(x, y)), where x, y are the conditional
zero-margins of the response bit at that history.  Primary form tested
(assembly-exact, secant form):

    MM(mu, lam) = E_hist[ h2(z~) - h2(z_{2^lam}(x, y)) ]  >=? 0

where z~ = realized conditional both-zero probability (so h2(z~) =
h2(z_OR(x, y)) exactly), and the expectation is mass-weighted over
nondegenerate histories of all i < n.  By the mean-value theorem this IS
E[sigma * (log2 OR - lam)] with sigma the Plackett sensitivity of 007
lead 2, evaluated on the secant; two pointwise-weighted variants
(signed derivative at the target, |derivative| at the target) are computed
alongside.  Degenerate-margin histories are excluded by the same N_i
bookkeeping as the plain control (their sensitivity weight -> 0 anyway --
that is the entire point of the candidate).

Candidate II -- lambda-restricted control (017 C2; 016 lead 4).  The plain
aggregate A(mu, lam) >= 0 restricted to the workable window
lam <= LAMWIN(n) = 4.847/(n-3) (the 009/011 lambda-window law).  The kill
in 016 lives at lam in [2, 2.5], ~40x above the window at n = 96; this
probe measures the ladder INSIDE the window, with theta re-optimized at
window lambda, and fits the small-lambda expansion
A ~ a1(n)*lam + a2(n)*lam^2 to see whether first order (Theorem C's
perfect square, a1 >= 0) protects the window asymptotically.

Parts:
  A  anchors: 007 witness aggregate; 016's theta* kill at n = 96 reproduced;
     Plackett round-trip identities on real census rows; padding invariance.
  B  MM on the witness + per-history table; MM along the raw-weight and
     theta*(lam=2) ladders, n up to 160; lambda profiles at n = 96.
  C  adversarial theta re-optimization AGAINST MM at n = 48 (lam 2 and 3.5),
     transfer to n = 64..160; free-support spot climbs at n = 8 against MM.
  D  window probe: per-n ladder at lam = LAMWIN(n), LAMWIN/2, LAMWIN/4 with
     theta re-optimized at window lambda (n = 48 and 96 opts, transfers to
     320); small-lambda coefficient fit a1(n), a2(n); free (theta, lam)
     climb with lam <= LAMWIN as a constraint.
  E  (--ext) raw-witness-weight ladder extension past n = 128 (016 lead 3),
     honest dilution-swept minimum, lam = 2.
  F  (--cert) exact rational certification at rational tilt t of the
     decisive points found by B-D: plain aggregate via 016's exact path
     (013's engine), MM via a new certified interval path (dyadic bisection
     for z_t, directed log2 enclosures for h2), with (x, y)-level caching.

Standard library only.  Deterministic (fixed seeds, fixed step counts).
Default run (A-D) ~20-40 min; --ext and --cert are separate longer passes.

Usage:
    python problems/union-closed/explore/uc_gap1_candidates.py [--fast]
    python problems/union-closed/explore/uc_gap1_candidates.py --ext
    python problems/union-closed/explore/uc_gap1_candidates.py --cert
Checkpoints: problems/union-closed/data/gap1c_part[A-F]*.json
"""

from __future__ import annotations

import json
import math
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
sys.path.insert(0, str(HERE))

import uc_or_avg as ORA                      # 007 census (float)
import uc_or_avg_skeptic as ORS              # 013 exact machinery
from uc_agg_ctrl_probe import (agg_eval, mu_ladder, exact_aggregate,
                               dec_bound, MARG_TARGET)
from uc_agg_ctrl_probe2 import (mu_ladder_theta, tune_dil, THETA_WITNESS,
                                marginals_of_u)

FAST = "--fast" in sys.argv
T0 = time.monotonic()

LAMWIN_C = 4.847        # 009/011 lambda-window law: lam_max ~ LAMWIN_C/(n-3)


def lamwin(n: int) -> float:
    return LAMWIN_C / (n - 3)


# theta* from 016 P2 (data/aggprobe2_partP2.json), the weights that force
# the certified kill at n >= 96, lam = 2
THETA_STAR_L2 = {
    "A1": 14.026917090492047, "A2": 22.622718089566355,
    "B0": 13.573329717299568, "B1": 43.00957928858568,
    "B2": 24.62436309125563, "dil": 21.044172077650057,
    "light": 0.007943082225456255, "resp": 0.0017259962695652103,
}


def save(name: str, obj) -> None:
    DATA.mkdir(exist_ok=True)
    (DATA / name).write_text(json.dumps(obj, indent=1, sort_keys=True,
                                        default=str) + "\n", encoding="utf-8")


def log(msg: str) -> None:
    print(f"[{time.monotonic() - T0:7.0f}s] {msg}", flush=True)


# ======================================================================
# Plackett helpers (zero-margin convention as in uc_pert/uc_couplings:
# x = P(A_i = 0 | hist), y = P(B_i = 0 | hist), z = P(both 0 | hist),
# rho = z(1-x-y+z) / ((x-z)(y-z)))
# ======================================================================

def h2(z: float) -> float:
    """binary entropy in bits; 0 at the endpoints"""
    if z <= 0.0 or z >= 1.0:
        return 0.0
    return -z * math.log2(z) - (1.0 - z) * math.log2(1.0 - z)


def z_plackett(x: float, y: float, rho: float) -> float:
    """Both-zero probability of the Plackett coupling with zero-margins
    x, y and odds ratio rho (stable root of a z^2 - S z + rho x y = 0)."""
    if rho == 1.0:
        return x * y
    a = rho - 1.0
    S = 1.0 + a * (x + y)
    disc = S * S - 4.0 * a * rho * x * y
    return 2.0 * rho * x * y / (S + math.sqrt(disc))


def dh_dlam(x: float, y: float, rho: float) -> float:
    """d/d lam of h2(z_{2^lam}(x, y)) at 2^lam = rho (signed).  Vanishes as
    a margin degenerates: dz/drho carries the factor (x-z)(y-z)."""
    z = z_plackett(x, y, rho)
    if z <= 0.0 or z >= 1.0 or z >= min(x, y) or z <= max(0.0, x + y - 1.0):
        return 0.0
    denom = (1.0 - x - y + 2.0 * z) + rho * (x + y - 2.0 * z)
    dzdrho = (x - z) * (y - z) / denom
    return math.log2((1.0 - z) / z) * dzdrho * rho * math.log(2.0)


# ======================================================================
# the margin-modulated census (float path)
# ======================================================================

def census_mm(u: dict[int, float], n: int, lam: float) -> dict:
    """Per-history margin-modulated bookkeeping over i < n.

    Returns the three MM aggregates plus the plain aggregate recomputed from
    the same tables (anchor against agg_eval), max marginal, entropy,
    dichotomy flag, and the per-i extreme MM contributions.
    """
    pi = ORA.tilt_of_potential(u, lam)
    supp = {a for (a, _b) in pi}
    rho_t = 2.0 ** lam

    num_sec = 0.0          # sum m * (h2(z~) - h2(z_target))
    num_der = 0.0          # sum m * dh_dlam * (l2or - lam)
    num_abs = 0.0          # sum m * |dh_dlam| * (l2or - lam)
    den_mass = 0.0
    den_absw = 0.0
    num_plain = 0.0
    per_i_sec: dict[int, float] = {}
    dich_all = True
    n_rows = 0

    for i in range(1, n):                       # i < n only
        pmask = (1 << (i - 1)) - 1
        bit = 1 << (i - 1)
        prefixes = {A & pmask for A in supp}
        active = {c for c in prefixes
                  if any(A & pmask == c and A & bit for A in supp)
                  and any(A & pmask == c and not A & bit for A in supp)}
        tables: dict[tuple[int, int], list[float]] = {}
        for (A, B), w in pi.items():
            key = (A & pmask, B & pmask)
            tab = tables.setdefault(key, [0.0] * 4)
            tab[(2 if A & bit else 0) + (1 if B & bit else 0)] += w
        i_sec = 0.0
        for (a, b), (p00, p01, p10, p11) in tables.items():
            nondeg = min(p00, p01, p10, p11) > 0.0
            if nondeg != (a in active and b in active):
                dich_all = False
            if not nondeg:
                continue
            mass = p00 + p01 + p10 + p11
            l2or = (math.log2(p11) + math.log2(p00)
                    - math.log2(p10) - math.log2(p01))
            x = (p00 + p01) / mass              # P(A_i = 0 | hist)
            y = (p00 + p10) / mass              # P(B_i = 0 | hist)
            z_rlz = p00 / mass
            z_t = z_plackett(x, y, rho_t)
            w_der = dh_dlam(x, y, rho_t)
            d = l2or - lam
            sec = mass * (h2(z_rlz) - h2(z_t))
            num_sec += sec
            i_sec += sec
            num_der += mass * w_der * d
            num_abs += mass * abs(w_der) * d
            den_mass += mass
            den_absw += mass * abs(w_der)
            num_plain += mass * d
            n_rows += 1
        per_i_sec[i] = i_sec

    if den_mass == 0.0:
        return None
    mu = ORA.mu_of(pi)
    worst_i = min(per_i_sec, key=lambda k: per_i_sec[k])
    return {
        "mm_sec": num_sec / den_mass,
        "mm_der": num_der / (den_absw if den_absw > 0 else 1.0),
        "mm_abs": num_abs / (den_absw if den_absw > 0 else 1.0),
        "agg_plain": num_plain / den_mass,
        "num_sec": num_sec, "num_der": num_der, "num_abs": num_abs,
        "den_mass": den_mass,
        "max_marginal": max(ORA.marginal_elements(mu, n)),
        "H": ORA.entropy_of(mu),
        "dichotomy_ok": dich_all, "n_rows": n_rows,
        "per_i_sec_min": per_i_sec[worst_i], "per_i_sec_argmin": worst_i,
    }


# ======================================================================
# part A -- anchors
# ======================================================================

def part_a() -> dict:
    out = {}
    wit = ORA.WITNESS_U

    # 1. plain aggregate reproduces 013/016
    ev = agg_eval(wit, 7, 3.5)
    log(f"[A] witness plain aggregate (lam 3.5): {ev['agg']:+.6f} "
        f"(013 certified +1.844669)")
    out["witness_agg"] = ev["agg"]
    ok_wit = abs(ev["agg"] - 1.844669) < 1e-3

    # 2. census_mm's internal plain aggregate agrees with agg_eval
    mm = census_mm(wit, 7, 3.5)
    log(f"[A] census_mm plain-agg cross-check: {mm['agg_plain']:+.6f} "
        f"(diff {abs(mm['agg_plain'] - ev['agg']):.2e})")
    ok_cross = abs(mm["agg_plain"] - ev["agg"]) < 1e-9

    # 3. 016's theta* kill at n = 96, lam = 2 reproduces
    th = dict(THETA_STAR_L2)
    th = tune_dil(96, 90, 2.0, th)
    ev96 = agg_eval(mu_ladder_theta(96, 90, th), 96, 2.0)
    log(f"[A] theta* ladder n=96 lam=2 plain aggregate: {ev96['agg']:+.9f} "
        f"(016 certified -0.000759865)")
    ok_kill = ev96["agg"] < 0 and abs(ev96["agg"] + 0.000759865) < 5e-5

    # 4. Plackett round-trip on real census rows: z_plackett(x, y, OR)
    #    must equal the realized z~ (the 2x2 table is determined by its
    #    margins and OR)
    pi = ORA.tilt_of_potential(wit, 3.5)
    supp = {a for (a, _b) in pi}
    worst = 0.0
    n_checked = 0
    for i in range(1, 7):
        pmask, bit = (1 << (i - 1)) - 1, 1 << (i - 1)
        tables: dict[tuple[int, int], list[float]] = {}
        for (A, B), w in pi.items():
            key = (A & pmask, B & pmask)
            tab = tables.setdefault(key, [0.0] * 4)
            tab[(2 if A & bit else 0) + (1 if B & bit else 0)] += w
        for (p00, p01, p10, p11) in tables.values():
            if min(p00, p01, p10, p11) <= 0:
                continue
            mass = p00 + p01 + p10 + p11
            x, y, z = (p00 + p01) / mass, (p00 + p10) / mass, p00 / mass
            orr = (p00 * p11) / (p01 * p10)
            worst = max(worst, abs(z_plackett(x, y, orr) - z))
            n_checked += 1
    log(f"[A] Plackett round-trip on {n_checked} witness rows: "
        f"max |z_rho(x,y,OR) - z~| = {worst:.2e}")
    ok_rt = worst < 1e-12

    # 5. secant-vs-derivative consistency: on one row-like configuration,
    #    h2(z_rho2) - h2(z_rho1) ~ integral of dh_dlam
    x, y = 0.62, 0.55
    l1, l2 = 1.0, 1.001
    sec = h2(z_plackett(x, y, 2 ** l2)) - h2(z_plackett(x, y, 2 ** l1))
    der = dh_dlam(x, y, 2 ** l1) * (l2 - l1)
    log(f"[A] secant vs derivative (dl=1e-3): {sec:.9f} vs {der:.9f} "
        f"(rel diff {abs(sec - der) / abs(der):.1e})")
    ok_sd = abs(sec - der) / abs(der) < 1e-2

    # 6. MM padding invariance: MU(20,1) == witness on the MM aggregates
    mm20 = census_mm(mu_ladder(20, 1), 20, 3.5)
    d_pad = abs(mm20["mm_sec"] - mm["mm_sec"])
    log(f"[A] MM padding invariance MU(20,1) vs witness: diff {d_pad:.2e}")
    ok_pad = d_pad < 1e-10

    out.update({"ok_witness": ok_wit, "ok_cross": ok_cross,
                "ok_kill_repro": ok_kill, "kill_repro_agg": ev96["agg"],
                "ok_roundtrip": ok_rt, "ok_secant_deriv": ok_sd,
                "ok_padding": ok_pad, "witness_mm_sec": mm["mm_sec"],
                "witness_mm_der": mm["mm_der"], "witness_mm_abs": mm["mm_abs"]})
    ok = all(v for k, v in out.items() if k.startswith("ok_"))
    log(f"[A] anchors {'ALL PASS' if ok else '*** FAILURE ***'}")
    out["all_ok"] = ok
    save("gap1c_partA.json", out)
    if not ok:
        raise SystemExit("anchor failure -- do not trust anything downstream")
    return out


# ======================================================================
# part B -- MM on witness and ladders
# ======================================================================

def theta_ladder_at(n: int, lam: float, th_base: dict) -> dict[int, float]:
    r = n - 6
    th = tune_dil(n, r, lam, th_base)
    return mu_ladder_theta(n, r, th)


def part_b() -> dict:
    out_rows = []
    wit = ORA.WITNESS_U

    mmw = census_mm(wit, 7, 3.5)
    log(f"[B] witness (n=7, lam=3.5): plain {mmw['agg_plain']:+.4f}  "
        f"MM_sec {mmw['mm_sec']:+.6f}  MM_der {mmw['mm_der']:+.6f}  "
        f"MM_abs {mmw['mm_abs']:+.6f}")
    out_rows.append({"name": "witness", "n": 7, "lambda": 3.5, **{
        k: mmw[k] for k in ("mm_sec", "mm_der", "mm_abs", "agg_plain",
                            "max_marginal", "dichotomy_ok")}})

    # ladders: raw witness weights and theta*(lam=2), lam = 2
    ns = (16, 32, 48, 64, 96, 128, 160) if not FAST else (16, 32)
    for name, th_base in (("raw", THETA_WITNESS), ("theta*", THETA_STAR_L2)):
        for n in ns:
            u = theta_ladder_at(n, 2.0, th_base)
            mm = census_mm(u, n, 2.0)
            log(f"[B] {name:7s} ladder n={n:3d} lam=2: "
                f"plain {mm['agg_plain']:+.6f}  MM_sec {mm['mm_sec']:+.6e}  "
                f"MM_der {mm['mm_der']:+.6e}  MM_abs {mm['mm_abs']:+.6e}  "
                f"(mm {mm['max_marginal']:.3f})")
            out_rows.append({"name": f"{name}_ladder", "n": n, "lambda": 2.0,
                             **{k: mm[k] for k in
                                ("mm_sec", "mm_der", "mm_abs", "agg_plain",
                                 "max_marginal", "dichotomy_ok",
                                 "per_i_sec_min", "per_i_sec_argmin")}})

    # lambda profile of the theta* ladder at n = 96
    if not FAST:
        for lam in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5):
            u = theta_ladder_at(96, lam, THETA_STAR_L2)
            mm = census_mm(u, 96, lam)
            log(f"[B] theta* ladder n=96 lam={lam:3.1f}: "
                f"plain {mm['agg_plain']:+.6f}  MM_sec {mm['mm_sec']:+.6e}  "
                f"MM_abs {mm['mm_abs']:+.6e}")
            out_rows.append({"name": "theta*_lamprof", "n": 96, "lambda": lam,
                             **{k: mm[k] for k in
                                ("mm_sec", "mm_der", "mm_abs", "agg_plain",
                                 "max_marginal", "dichotomy_ok")}})

    save("gap1c_partB.json", {"rows": out_rows})
    return {"rows": out_rows}


# ======================================================================
# part C -- adversarial theta re-optimization against MM
# ======================================================================

def mm_objective(th: dict, n: int, r: int, lam: float, key: str = "mm_abs"):
    """Objective = the |sigma|-weighted variant: parts A/B showed the two
    SIGNED variants are already negative on the witness itself, so mm_abs is
    the only form left standing as a candidate control -- it is what the
    adversary must now attack."""
    u = mu_ladder_theta(n, r, th)
    mm = census_mm(u, n, lam)
    if mm is None:
        return float("inf"), None
    pen = 300.0 * max(0.0, mm["max_marginal"] - MARG_TARGET)
    return mm[key] + pen, mm


def part_c() -> dict:
    out_rows = []
    n_opt = 48 if not FAST else 16
    r_opt = n_opt - 6
    steps = 260 if not FAST else 30
    best_by_lam = {}
    for lam in (2.0, 3.5):
        rng = random.Random(180000 + int(lam * 10))
        th = dict(tune_dil(n_opt, r_opt, lam, THETA_STAR_L2))
        cur, mm = mm_objective(th, n_opt, r_opt, lam)
        for step in range(steps):
            sigma = 0.8 * (1 - step / steps) + 0.05
            cand = dict(th)
            k = rng.choice(list(cand))
            cand[k] = max(cand[k], 1e-8) * math.exp(rng.gauss(0, sigma))
            v2, mm2 = mm_objective(cand, n_opt, r_opt, lam)
            if v2 < cur:
                th, cur, mm = cand, v2, mm2
        log(f"[C] MM theta-opt n={n_opt} lam={lam}: "
            f"MM_abs {mm['mm_abs']:+.6e} MM_sec {mm['mm_sec']:+.6e} "
            f"plain {mm['agg_plain']:+.4f} (mm {mm['max_marginal']:.3f})")
        best_by_lam[lam] = dict(th)
        out_rows.append({"stage": "opt", "n": n_opt, "lambda": lam,
                         "theta": th, **{k: mm[k] for k in
                                         ("mm_sec", "mm_der", "mm_abs",
                                          "agg_plain", "max_marginal")}})
        for n in ((64, 96, 128, 160) if not FAST else (24,)):
            r = n - 6
            th_n = tune_dil(n, r, lam, th)
            mm_n = census_mm(mu_ladder_theta(n, r, th_n), n, lam)
            log(f"[C]   MM theta({lam}) at n={n}: "
                f"MM_abs {mm_n['mm_abs']:+.6e} MM_sec {mm_n['mm_sec']:+.6e} "
                f"plain {mm_n['agg_plain']:+.4f} "
                f"(mm {mm_n['max_marginal']:.3f})")
            out_rows.append({"stage": "transfer", "n": n, "lambda": lam,
                             "theta": th_n, **{k: mm_n[k] for k in
                                               ("mm_sec", "mm_der", "mm_abs",
                                                "agg_plain", "max_marginal")}})

    # free-support climbs against MM_abs at n = 8 (the ladder hides its
    # deficit at near-degenerate margins, which |sigma| suppresses -- a
    # violation of mm_abs needs deficits at HEALTHY margins, a geometry the
    # theta-ladder cannot express; free atoms can)
    rng = random.Random(181111)
    best_free = None
    wit8 = {a: w for a, w in ORA.WITNESS_U.items()}     # n=7 embeds in n=8

    def free_obj(m):
        return (m["mm_abs"]
                + 300.0 * max(0.0, m["max_marginal"] - MARG_TARGET))

    for trial in range(8 if not FAST else 2):
        lam = (1.0, 2.0, 3.5, 2.0, 3.5, 1.0, 2.0, 3.5)[trial % 8]
        if trial < 3:
            u = dict(wit8)
        else:
            u = {a: math.exp(rng.uniform(-1.5, 1.5))
                 for a in rng.sample(range(256), 10)}
        mm = census_mm(u, 8, lam)
        if mm is None:
            continue
        cur = free_obj(mm)
        for step in range(600 if not FAST else 30):
            cand = dict(u)
            a = rng.choice(list(cand))
            if rng.random() < 0.15 and len(cand) < 14:
                cand[rng.getrandbits(8)] = math.exp(rng.uniform(-3, 1))
            else:
                cand[a] = max(cand[a] * math.exp(rng.gauss(0, 0.5)), 1e-9)
            mm2 = census_mm(cand, 8, lam)
            if mm2 is None:
                continue
            v2 = free_obj(mm2)
            if v2 < cur:
                u, cur, mm = cand, v2, mm2
        row = {"trial": trial, "lambda": lam, "mm_abs": mm["mm_abs"],
               "mm_sec": mm["mm_sec"], "max_marginal": mm["max_marginal"],
               "u": {format(a, "08b"): w for a, w in u.items()}}
        if (mm["max_marginal"] <= MARG_TARGET
                and (best_free is None or mm["mm_abs"] < best_free["mm_abs"])):
            best_free = row
        log(f"[C] free climb n=8 lam={lam}: endpoint MM_abs "
            f"{mm['mm_abs']:+.6e} MM_sec {mm['mm_sec']:+.6e} "
            f"(mm {mm['max_marginal']:.3f})")
    out = {"rows": out_rows, "best_free": best_free,
           "best_theta_by_lam": {str(k): v for k, v in best_by_lam.items()}}
    save("gap1c_partC.json", out)
    return out


# ======================================================================
# part D -- the lambda-window-restricted candidate
# ======================================================================

def part_d() -> dict:
    out_rows = []

    # window scan with theta* and raw weights, lam = LAMWIN(n) etc.
    ns = (48, 96, 160, 224, 320) if not FAST else (24, 48)
    for name, th_base in (("theta*", THETA_STAR_L2), ("raw", THETA_WITNESS)):
        for n in ns:
            lw = lamwin(n)
            row = {"name": f"{name}_window", "n": n, "lamwin": lw, "points": []}
            lams = [lw, lw / 2] + ([lw / 4] if n < 224 else [])
            for lam in lams:
                u = theta_ladder_at(n, lam, th_base)
                ev = agg_eval(u, n, lam)
                row["points"].append({"lambda": lam, "agg": ev["agg"],
                                      "max_marginal": ev["max_marginal"],
                                      "per_i_min": ev["per_i_min"]})
            # small-lambda fit A ~ a1 lam + a2 lam^2 from the two smallest
            # lambda points; with 3 points the lw row is an out-of-fit check
            p0 = row["points"][0]
            p1, p2 = row["points"][-1], row["points"][-2]
            l1, l2 = p1["lambda"], p2["lambda"]
            a1 = (p1["agg"] * l2 * l2 - p2["agg"] * l1 * l1) / \
                 (l1 * l2 * (l2 - l1))
            a2 = (p2["agg"] * l1 - p1["agg"] * l2) / (l1 * l2 * (l2 - l1))
            pred = a1 * p0["lambda"] + a2 * p0["lambda"] ** 2
            row.update({"a1": a1, "a2": a2, "fit_check_at_lw":
                        {"pred": pred, "actual": p0["agg"],
                         "in_fit": len(lams) == 2}})
            log(f"[D] {name:7s} n={n:3d} window lam={lw:.4f}: "
                f"A(lw) {p0['agg']:+.6e}  A(lw/2) {row['points'][1]['agg']:+.6e}  "
                f"a1 {a1:+.5e}  a2 {a2:+.5e}  a1*n {a1 * n:+.3f}")
            out_rows.append(row)

    # theta re-optimization AT window lambda (the adversary adapts to the
    # window), then per-n transfer where each n uses its own LAMWIN(n)
    best_theta = None
    for n_opt in ((48, 96) if not FAST else (24,)):
        r_opt = n_opt - 6
        lam = lamwin(n_opt)
        rng = random.Random(190000 + n_opt)
        th = dict(tune_dil(n_opt, r_opt, lam, THETA_STAR_L2))
        ev = agg_eval(mu_ladder_theta(n_opt, r_opt, th), n_opt, lam)
        cur = ev["agg"] + 100.0 * max(0.0, ev["max_marginal"] - MARG_TARGET)
        steps = 300 if not FAST else 30
        for step in range(steps):
            sigma = 0.8 * (1 - step / steps) + 0.05
            cand = dict(th)
            k = rng.choice(list(cand))
            cand[k] = max(cand[k], 1e-8) * math.exp(rng.gauss(0, sigma))
            ev2 = agg_eval(mu_ladder_theta(n_opt, r_opt, cand), n_opt, lam)
            if ev2 is None:
                continue
            v2 = (ev2["agg"]
                  + 100.0 * max(0.0, ev2["max_marginal"] - MARG_TARGET))
            if v2 < cur:
                th, cur, ev = cand, v2, ev2
        log(f"[D] window theta-opt n={n_opt} lam={lam:.4f}: "
            f"agg {ev['agg']:+.6e} (mm {ev['max_marginal']:.3f})")
        out_rows.append({"stage": "window_opt", "n": n_opt, "lambda": lam,
                         "theta": th, "agg": ev["agg"],
                         "max_marginal": ev["max_marginal"]})
        if best_theta is None or ev["agg"] < best_theta[1]:
            best_theta = (dict(th), ev["agg"], n_opt)
        for n in ((96, 160, 224, 320) if not FAST else (48,)):
            if n == n_opt:
                continue
            r = n - 6
            lw = lamwin(n)
            th_n = tune_dil(n, r, lw, th)
            ev_n = agg_eval(mu_ladder_theta(n, r, th_n), n, lw)
            log(f"[D]   window theta(n={n_opt}) at n={n} lam={lw:.4f}: "
                f"agg {ev_n['agg']:+.6e} (mm {ev_n['max_marginal']:.3f})")
            out_rows.append({"stage": "window_transfer", "n": n, "lambda": lw,
                             "from_n": n_opt, "theta": th_n,
                             "agg": ev_n["agg"],
                             "max_marginal": ev_n["max_marginal"]})

    # free (theta, lam <= LAMWIN) climb: adversary also picks lambda inside
    # the window
    n_free = 96 if not FAST else 24
    r_free = n_free - 6
    lw = lamwin(n_free)
    rng = random.Random(191111)
    th = dict(tune_dil(n_free, r_free, lw, THETA_STAR_L2))
    lam = lw
    ev = agg_eval(mu_ladder_theta(n_free, r_free, th), n_free, lam)
    cur = ev["agg"] + 100.0 * max(0.0, ev["max_marginal"] - MARG_TARGET)
    for step in range(200 if not FAST else 20):
        sigma = 0.8 * (1 - step / 200) + 0.05
        cand, cl = dict(th), lam
        if rng.random() < 0.25:
            cl = min(lw, max(lw / 64, lam * math.exp(rng.gauss(0, sigma))))
        else:
            k = rng.choice(list(cand))
            cand[k] = max(cand[k], 1e-8) * math.exp(rng.gauss(0, sigma))
        ev2 = agg_eval(mu_ladder_theta(n_free, r_free, cand), n_free, cl)
        if ev2 is None:
            continue
        v2 = ev2["agg"] + 100.0 * max(0.0, ev2["max_marginal"] - MARG_TARGET)
        if v2 < cur:
            th, lam, cur, ev = cand, cl, v2, ev2
    log(f"[D] free (theta, lam<=window) climb n={n_free}: endpoint "
        f"agg {ev['agg']:+.6e} at lam={lam:.4f} (window {lw:.4f})")
    out_rows.append({"stage": "free_window_climb", "n": n_free,
                     "lambda": lam, "lamwin": lw, "theta": th,
                     "agg": ev["agg"], "max_marginal": ev["max_marginal"]})

    out = {"rows": out_rows}
    save("gap1c_partD.json", out)
    return out


# ======================================================================
# part E (--ext) -- raw-weight ladder extension (016 lead 3)
# ======================================================================

def part_e() -> dict:
    out_rows = []
    for n in (160, 192, 224, 256, 288, 320):
        r = n - 6
        best = None
        for k in range(6):                      # 6-point dilution sweep
            th = dict(THETA_WITNESS)
            th["dil"] = 10.0 * (2.0 ** k)
            marg = marginals_of_u(mu_ladder_theta(n, r, th), n, 2.0)
            if max(marg) > MARG_TARGET:
                continue
            ev = agg_eval(mu_ladder_theta(n, r, th), n, 2.0)
            if best is None or ev["agg"] < best["agg"]:
                best = {"agg": ev["agg"], "dil": th["dil"],
                        "max_marginal": ev["max_marginal"]}
        th = tune_dil(n, r, 2.0)                # plus the bisection tuner
        ev = agg_eval(mu_ladder_theta(n, r, th), n, 2.0)
        if best is None or ev["agg"] < best["agg"]:
            best = {"agg": ev["agg"], "dil": th["dil"],
                    "max_marginal": ev["max_marginal"]}
        log(f"[E] raw ladder n={n}: honest min agg {best['agg']:+.6f} "
            f"(dil {best['dil']:.1f}, mm {best['max_marginal']:.3f})")
        out_rows.append({"n": n, **best})
    save("gap1c_partE.json", {"rows": out_rows})
    return {"rows": out_rows}


# ======================================================================
# part F (--cert) -- exact certification
# ======================================================================

F0, F1 = Fraction(0), Fraction(1)


def h2_enclosure(z: Fraction) -> tuple[Fraction, Fraction]:
    """[lo, hi] for h2(z), z rational in (0, 1)."""
    if z <= 0 or z >= 1:
        return F0, F0
    l1lo, l1hi = ORS.log2_enclosure(z)
    l2lo, l2hi = ORS.log2_enclosure(1 - z)
    return (-z * l1hi - (1 - z) * l2hi), (-z * l1lo - (1 - z) * l2lo)


def z_target_enclosure(x: Fraction, y: Fraction, t: Fraction,
                       steps: int = 80) -> tuple[Fraction, Fraction]:
    """Dyadic bisection enclosure of the Plackett root z_t(x, y):
    F(z) = z(1-x-y+z) - t(x-z)(y-z), unique root in the coupling range.
    For t > 1 the root lies in (xy, min(x,y)); t < 1 in (max(0,x+y-1), xy)."""
    if t == 1:
        return x * y, x * y

    def bigF(z: Fraction) -> Fraction:
        return z * (1 - x - y + z) - t * (x - z) * (y - z)

    if t > 1:
        lo, hi = x * y, min(x, y)
    else:
        lo, hi = max(F0, x + y - 1), x * y
    # F(lo) < 0 < F(hi) for t > 1 (and reversed signs handled uniformly)
    flo = bigF(lo)
    for _ in range(steps):
        mid = (lo + hi) / 2
        # keep denominators dyadic-bounded
        mid = Fraction(mid.numerator, mid.denominator)
        fm = bigF(mid)
        if fm == 0:
            return mid, mid
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return lo, hi


def h2_over_interval(zlo: Fraction, zhi: Fraction) -> tuple[Fraction, Fraction]:
    """[lo, hi] of h2 over [zlo, zhi] (h2 unimodal, peak 1 at z = 1/2)."""
    alo, ahi = h2_enclosure(zlo)
    blo, bhi = h2_enclosure(zhi)
    lo = min(alo, blo)
    hi = max(ahi, bhi)
    if zlo < Fraction(1, 2) < zhi:
        hi = F1
    return lo, hi


def exact_mm(u_float: dict[int, float], n: int, t: Fraction,
             rationalize=None):
    """Certified enclosure of the MM_sec numerator (and the aggregate) at
    rational tilt t for the exact measure whose potential is the
    rationalization of u_float.  Caches by (x, y) -- ladder units are
    exchangeable so the distinct-margin count is tiny compared to rows."""
    if rationalize is None:
        u = {a: Fraction(w) for a, w in u_float.items()}
    else:
        u = {a: rationalize(w) for a, w in u_float.items()}
        u = {a: w for a, w in u.items() if w > 0}
    cens = ORS.exact_census(u, t, n)
    zcache: dict[tuple[Fraction, Fraction], tuple] = {}
    num_lo = num_hi = F0
    plain_lo = plain_hi = F0
    wtot = F0
    dich = True
    n_rows = 0
    for c in cens:
        if c["i"] >= n:
            continue
        dich = dich and c["dich_ok"]
        active = c["active"]
        for (a, b), tab in c["tables"].items():
            p00, p01, p10, p11 = tab
            if min(p00, p01, p10, p11) <= 0:
                continue
            assert a in active and b in active
            mass = p00 + p01 + p10 + p11
            x, y = (p00 + p01) / mass, (p00 + p10) / mass
            z_rlz = p00 / mass
            key = (x, y)
            if key not in zcache:
                zlo, zhi = z_target_enclosure(x, y, t)
                zcache[key] = h2_over_interval(zlo, zhi)
            htlo, hthi = zcache[key]
            hrlo, hrhi = h2_enclosure(z_rlz)
            num_lo += mass * (hrlo - hthi)
            num_hi += mass * (hrhi - htlo)
            orr = (p11 * p00) / (p10 * p01)
            llo, lhi = ORS.log2_enclosure(orr / t)
            plain_lo += mass * llo
            plain_hi += mass * lhi
            wtot += mass
            n_rows += 1
    _mu, marg = ORS.exact_marginals(u, t, n)
    return {"mm_num_lo": num_lo, "mm_num_hi": num_hi,
            "mm_lo": num_lo / wtot, "mm_hi": num_hi / wtot,
            "agg_lo": plain_lo / wtot, "agg_hi": plain_hi / wtot,
            "max_marg": max(marg), "in_regime": max(marg) < ORS.RECORD,
            "dich_ok": dich, "n_rows": n_rows,
            "n_distinct_margins": len(zcache)}


def tidy(w: float) -> Fraction:
    return Fraction(w).limit_denominator(10 ** 4)


def part_f() -> dict:
    """Job list is fixed from the B-D probe outcomes (see the attempt
    record): certify the decisive MM and window points."""
    out_rows = []
    jobs = []

    # (1) MM on the witness (007 lead 2's first ask, now certified)
    jobs.append(("mm_witness_t181/16", ORA.WITNESS_U, 7,
                 Fraction(181, 16), "mm", None))

    # (2) MM on the theta* kill instance: does the margin-modulated
    #     aggregate survive on the exact family that kills the plain one?
    th = tune_dil(96, 90, 2.0, dict(THETA_STAR_L2))
    u96 = mu_ladder_theta(96, 90, th)
    jobs.append(("mm_theta*_n96_t4", u96, 96, Fraction(4), "mm", None))

    # (3) same at n = 160
    th160 = tune_dil(160, 154, 2.0, dict(THETA_STAR_L2))
    u160 = mu_ladder_theta(160, 154, th160)
    jobs.append(("mm_theta*_n160_t4", u160, 160, Fraction(4), "mm", tidy))

    # (4) window points: plain aggregate at a rational tilt inside the
    #     window at n = 96 and 160 (LAMWIN(96) ~ 0.0521, log2(26/25) ~
    #     0.0566 > it, log2(51/50) ~ 0.0286 < it -> use 51/50 at 96;
    #     LAMWIN(160) ~ 0.0309 -> 51/50 works there too)
    for n, u in ((96, None), (160, None)):
        lw = lamwin(n)
        t_win = Fraction(51, 50)
        lam_t = math.log2(51 / 50)
        assert lam_t < lw, (n, lam_t, lw)
        r = n - 6
        thw = tune_dil(n, r, lam_t, dict(THETA_STAR_L2))
        uw = mu_ladder_theta(n, r, thw)
        jobs.append((f"window_theta*_n{n}_t51/50", uw, n, t_win, "plain",
                     tidy if n >= 160 else None))

    for name, u, n, t, kind, rat in jobs:
        t_start = time.monotonic()
        if kind == "mm":
            r = exact_mm(u, n, t, rationalize=rat)
            lo, hi = float(r["mm_lo"]), float(r["mm_hi"])
            verdict = ("CERTIFIED NEGATIVE (KILL)" if r["mm_hi"] < 0 else
                       "CERTIFIED POSITIVE" if r["mm_lo"] > 0 else
                       "STRADDLES 0")
            log(f"[F] exact MM {name}: MM_sec in [{lo:+.9e}, {hi:+.9e}] "
                f"{verdict}; plain in [{float(r['agg_lo']):+.6f}, "
                f"{float(r['agg_hi']):+.6f}]; in-regime {r['in_regime']}; "
                f"dich {r['dich_ok']}; rows {r['n_rows']} "
                f"({r['n_distinct_margins']} distinct margins); "
                f"{time.monotonic() - t_start:.0f}s")
            out_rows.append({
                "name": name, "n": n, "t": str(t), "kind": kind,
                "mm_lo": dec_bound(r["mm_lo"]),
                "mm_hi": dec_bound(r["mm_hi"], up=True),
                "mm_lo_float": lo, "mm_hi_float": hi,
                "agg_lo_float": float(r["agg_lo"]),
                "agg_hi_float": float(r["agg_hi"]),
                "in_regime": r["in_regime"],
                "max_marg_float": float(r["max_marg"]),
                "dich_ok": r["dich_ok"], "verdict": verdict})
        else:
            u_in = ({a: float(tidy(w)) for a, w in u.items()}
                    if rat is not None else u)
            r = exact_aggregate(u_in, n, t)
            lo, hi = float(r["agg_lo"]), float(r["agg_hi"])
            verdict = ("CERTIFIED NEGATIVE (KILL)" if r["agg_hi"] < 0 else
                       "CERTIFIED POSITIVE" if r["agg_lo"] > 0 else
                       "STRADDLES 0")
            log(f"[F] exact plain {name}: agg in [{lo:+.9e}, {hi:+.9e}] "
                f"{verdict}; in-regime {r['in_regime']}; dich {r['dich_ok']}; "
                f"{time.monotonic() - t_start:.0f}s")
            out_rows.append({
                "name": name, "n": n, "t": str(t), "kind": kind,
                "agg_lo": dec_bound(r["agg_lo"]),
                "agg_hi": dec_bound(r["agg_hi"], up=True),
                "agg_lo_float": lo, "agg_hi_float": hi,
                "in_regime": r["in_regime"],
                "max_marg_float": float(r["max_marg"]),
                "dich_ok": r["dich_ok"], "verdict": verdict})
        save("gap1c_partF.json", {"rows": out_rows})
    return {"rows": out_rows}


# ======================================================================
# part R (--robust) -- robustness battery for the witness MM_sec kill,
# plus MM on friendly genres (is the negativity generic or structural?)
# ======================================================================

def part_r() -> dict:
    out = {}
    wit = ORA.WITNESS_U

    # product anchors: MM_sec must vanish identically (OR = 2^lam at every
    # history, so z~ = z_target row by row)
    prows = []
    for n, p, lam in ((6, 0.3, 2.0), (8, 0.38, 1.0)):
        u = {a: (p / (1 - p)) ** ORA.popcount(a) for a in range(1 << n)}
        mm = census_mm(u, n, lam)
        prows.append({"n": n, "p": p, "lambda": lam, "mm_sec": mm["mm_sec"]})
        log(f"[R] product Bern({p})^{n} lam={lam}: MM_sec {mm['mm_sec']:+.1e}")
    out["product_anchors"] = prows

    # friendly non-product genres: MM_sec sign is NOT generically negative
    frows = []
    n = 8
    genres = {}
    eps = 0.04
    genres["d0_bern"] = {0: 1.0} | {
        a: 0.55 * (0.5 + eps) ** ORA.popcount(a)
        * (0.5 - eps) ** (n - ORA.popcount(a)) for a in range(1, 1 << n)}
    genres["bern_mix"] = {
        a: 0.6 * 0.2 ** ORA.popcount(a) * 0.8 ** (n - ORA.popcount(a))
        + 0.4 * 0.45 ** ORA.popcount(a) * 0.55 ** (n - ORA.popcount(a))
        for a in range(1 << n)}
    genres["smoothed_slice"] = {
        a: (1.0 if ORA.popcount(a) == 3 else 1e-3) for a in range(1 << n)}
    rng = random.Random(7)
    for tr in range(3):
        genres[f"random_dense_{tr}"] = {
            a: math.exp(rng.gauss(0, 0.35)) for a in range(1 << 6)}
    for name, u in genres.items():
        nn = 6 if name.startswith("random") else n
        for lam in (1.0, 2.0):
            mm = census_mm(u, nn, lam)
            frows.append({"name": name, "n": nn, "lambda": lam,
                          **{k: mm[k] for k in ("mm_sec", "mm_der", "mm_abs",
                                                "agg_plain", "max_marginal")}})
            log(f"[R] {name:16s} lam={lam}: MM_sec {mm['mm_sec']:+.5f} "
                f"MM_der {mm['mm_der']:+.5f} (mm {mm['max_marginal']:.3f})")
    out["friendly_genres"] = frows

    # witness kill robustness: 20 x 3% perturbations, 2-digit tidy, lambda
    # profile
    rng = random.Random(2024)
    pert = []
    for _ in range(20):
        u = {a: w * math.exp(rng.uniform(-0.03, 0.03))
             for a, w in wit.items()}
        mm = census_mm(u, 7, 3.5)
        pert.append({"mm_sec": mm["mm_sec"],
                     "max_marginal": mm["max_marginal"]})
    n_fail = sum(1 for r in pert
                 if r["mm_sec"] >= 0 or r["max_marginal"] > 0.38271)
    log(f"[R] witness 20 x 3% perturbations: {n_fail} failures; MM_sec in "
        f"[{min(r['mm_sec'] for r in pert):+.5f}, "
        f"{max(r['mm_sec'] for r in pert):+.5f}]")
    out["perturbations"] = {"rows": pert, "n_fail": n_fail}

    u2 = {a: float(f"{w:.2g}") for a, w in wit.items()}
    mm2 = census_mm(u2, 7, 3.5)
    out["tidy_2digit"] = {"mm_sec": mm2["mm_sec"],
                          "max_marginal": mm2["max_marginal"]}
    log(f"[R] 2-digit tidy witness: MM_sec {mm2['mm_sec']:+.5f} "
        f"(mm {mm2['max_marginal']:.4f})")

    prof = []
    for l in (0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0):
        mm = census_mm(wit, 7, l)
        prof.append({"lambda": l, **{k: mm[k] for k in
                                     ("mm_sec", "mm_der", "mm_abs")}})
        log(f"[R] witness lam={l:4.2f}: MM_sec {mm['mm_sec']:+.6f} "
            f"MM_der {mm['mm_der']:+.6f} MM_abs {mm['mm_abs']:+.6f}")
    out["witness_lambda_profile"] = prof
    save("gap1c_partR.json", out)
    return out


def main() -> None:
    if "--ext" in sys.argv:
        part_e()
        return
    if "--cert" in sys.argv:
        part_a()
        part_f()
        return
    if "--robust" in sys.argv:
        part_a()
        part_r()
        return
    part_a()
    part_b()
    part_c()
    part_d()
    log("done (parts A-D)")


if __name__ == "__main__":
    main()
