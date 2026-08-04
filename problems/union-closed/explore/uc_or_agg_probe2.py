#!/usr/bin/env python3
"""Attempt 014, deep phase: push the aggregated-control probe past n = 32.

Round 1 (uc_or_agg_probe.py part B) found the multi-unit witness ladder
MU(n, n-6) has aggregate decaying like a - b log n over n = 7..32 (best
2-parameter fit; max residual ~0.05), which EXTRAPOLATES to a sign change at
n ~ 78 (lambda 5), ~140 (lambda 3.5), ~320 (lambda 2).  Extrapolation is not
evidence (012's lesson cuts both ways), so this phase measures instead:

  P1  ladder extension: MU(n, n-6), dilution tuned on marginals only,
      n up to 128, lambda in {2, 3.5} (+5 at spot points).  Sign check.
  P2  unit-symmetric ladder optimization: the ladder is determined by 8
      shared weights (A1, A2, light, B0, B1, B2, resp, dil); climb THETA at
      moderate n (in-regime penalty), then evaluate the optimized theta at
      large n.  This is the orbit-symmetrization idea in search form: the
      symmetric parametrization keeps the eval cheap and the search
      low-dimensional at n far beyond round 1's free-support climbs.
  P3  deep annealed free-support climbs at n = 24, 32, seeded from round
      1's best endpoints (checkpoints on disk) + fresh ladders; 2500 steps,
      stiffer marginal penalty, occasional uphill acceptance.
  P4  lambda profiles of the extreme instances (violations must not hide
      between grid points).
  P5  exact rational certification (013's engine) of the decisive points:
      any float violation found, else the tightest in-regime margins at the
      largest n reached, in tidy-rational form (limit_denominator) so the
      certification stays tractable at large n.

Standard library only.  Deterministic (fixed seeds).  Runtime ~30-60 min.

Usage:
    python problems/union-closed/explore/uc_or_agg_probe2.py
Checkpoints: problems/union-closed/data/aggprobe2_partP[12345]*.json
             problems/union-closed/data/aggprobe2_run.log (tee by hand)
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

import uc_or_avg as ORA
import uc_or_avg_skeptic as ORS
from uc_or_agg_probe import (agg_eval, clamp_u, mu_ladder, objective,
                             exact_aggregate, dec_bound, RECORD, MARG_TARGET)

T0 = time.monotonic()


def save(name: str, obj) -> None:
    DATA.mkdir(exist_ok=True)
    (DATA / name).write_text(json.dumps(obj, indent=1, sort_keys=True,
                                        default=str) + "\n", encoding="utf-8")


def log(msg: str) -> None:
    print(msg, flush=True)


# ----------------------------------------------------------------------
# fast marginal-only evaluation (dilution tuning at large n must not pay
# for the census)
# ----------------------------------------------------------------------

def marginals_of_u(u: dict[int, float], n: int, lam: float):
    supp = list(u)
    z = 0.0
    mu = {}
    for a in supp:
        row = sum(u[a] * u[b] * 2.0 ** (lam * ORA.popcount(a & b))
                  for b in supp)
        mu[a] = row
        z += row
    marg = [0.0] * n
    for a, w in mu.items():
        for j in range(n):
            if a >> j & 1:
                marg[j] += w / z
    return marg


THETA_WITNESS = {
    "A1": 3.802705744653898, "A2": 8.094727225702503,
    "light": 0.004792968160237178, "B0": 0.2526660836353666,
    "B1": 27.28641199464732, "B2": 12.529780963546393,
    "resp": 0.3609234435898425, "dil": 20.0,
}


def mu_ladder_theta(n: int, r: int, th: dict[str, float]) -> dict[int, float]:
    assert 1 <= r <= n - 6
    b_hi, b_lo = 1 << (n - 1), 1 << (n - 2)
    u = {0b1 | b_lo: th["A1"], 0b1 | b_hi: th["A2"], 0b0: th["B0"],
         b_lo: th["B1"], b_hi: th["B2"]}
    for j in range(r):
        c = 1 << (4 + j)
        u[0b1 | c | b_hi] = th["light"]
        u[c | b_lo] = th["resp"]
    for d in range(3):
        u[1 << (1 + d)] = th["dil"]
    return u


def tune_dil(n: int, r: int, lam: float, th=None):
    """Bisect the dilution weight on marginals only."""
    th = dict(th or THETA_WITNESS)
    lo, hi = 0.5, 2000.0
    for _ in range(18):
        w = math.sqrt(lo * hi)
        th["dil"] = w
        marg = marginals_of_u(mu_ladder_theta(n, r, th), n, lam)
        dil_m = max(marg[1:4])
        oth_m = max(m for j, m in enumerate(marg) if j not in (1, 2, 3))
        if oth_m > MARG_TARGET and dil_m < MARG_TARGET:
            lo = w
        elif dil_m > MARG_TARGET:
            hi = w
        else:
            break
    return th


# ----------------------------------------------------------------------
# P1 -- ladder extension
# ----------------------------------------------------------------------

def part_p1() -> dict:
    rows = []
    for lam, ns in ((2.0, (40, 48, 64, 80, 96, 128)),
                    (3.5, (40, 48, 64, 80, 96, 128)),
                    (5.0, (48, 96))):
        for n in ns:
            r = n - 6
            th = tune_dil(n, r, lam)
            t0 = time.monotonic()
            ev = agg_eval(mu_ladder_theta(n, r, th), n, lam)
            rows.append({"n": n, "r": r, "lambda": lam, "wdil": th["dil"],
                         "agg": ev["agg"], "max_marginal": ev["max_marginal"],
                         "per_i_min": ev["per_i_min"], "H": ev["H"],
                         "argmin_i": ev["argmin_i"],
                         "secs": round(time.monotonic() - t0, 1)})
            log(f"[P1] MU({n},{r}) lam={lam}: agg {ev['agg']:+.4f}, "
                f"max marg {ev['max_marginal']:.3f}, per-i min "
                f"{ev['per_i_min']:+.4f} ({rows[-1]['secs']}s)")
    save("aggprobe2_partP1.json", {"rows": rows})
    return {"rows": rows}


# ----------------------------------------------------------------------
# P2 -- unit-symmetric theta optimization
# ----------------------------------------------------------------------

def theta_objective(th, n, r, lam):
    ev = agg_eval(clamp_u(mu_ladder_theta(n, r, th)), n, lam)
    if ev is None:
        return float("inf"), None
    pen = (100.0 * max(0.0, ev["max_marginal"] - MARG_TARGET)
           + 10.0 * max(0.0, 0.05 - ev["wsum"])
           + 2.0 * max(0.0, 0.5 - ev["H"]))
    return ev["agg"] + pen, ev


def part_p2() -> dict:
    out_rows = []
    best_by_lam = {}
    n_opt, r_opt = 48, 42
    for lam in (2.0, 3.5):
        rng = random.Random(150000 + int(lam * 10))
        th = dict(tune_dil(n_opt, r_opt, lam))
        cur, ev = theta_objective(th, n_opt, r_opt, lam)
        for step in range(260):
            sigma = 0.8 * (1 - step / 260) + 0.05
            cand = dict(th)
            k = rng.choice(list(cand))
            cand[k] = max(cand[k], 1e-8) * math.exp(rng.gauss(0, sigma))
            v2, ev2 = theta_objective(cand, n_opt, r_opt, lam)
            if v2 < cur:
                th, cur, ev = cand, v2, ev2
        log(f"[P2] theta-opt at n={n_opt} lam={lam}: agg {ev['agg']:+.4f} "
            f"(mm {ev['max_marginal']:.3f}), theta="
            + " ".join(f"{k}:{v:.4g}" for k, v in sorted(th.items())))
        best_by_lam[lam] = th
        out_rows.append({"stage": "opt", "n": n_opt, "lambda": lam,
                         "theta": th, "agg": ev["agg"],
                         "max_marginal": ev["max_marginal"],
                         "per_i_min": ev["per_i_min"]})
        for n in (64, 96, 128):
            r = n - 6
            th_n = tune_dil(n, r, lam, th)
            ev_n = agg_eval(mu_ladder_theta(n, r, th_n), n, lam)
            log(f"[P2]   theta({lam}) at n={n}: agg {ev_n['agg']:+.4f}, "
                f"mm {ev_n['max_marginal']:.3f}, per-i min "
                f"{ev_n['per_i_min']:+.4f}")
            out_rows.append({"stage": "transfer", "n": n, "lambda": lam,
                             "theta": th_n, "agg": ev_n["agg"],
                             "max_marginal": ev_n["max_marginal"],
                             "per_i_min": ev_n["per_i_min"]})
    save("aggprobe2_partP2.json", {"rows": out_rows})
    return {"rows": out_rows, "best_by_lam": best_by_lam}


# ----------------------------------------------------------------------
# P3 -- deep annealed free-support climbs at n = 24, 32
# ----------------------------------------------------------------------

def anneal(u0, n, lam, rng, steps=2500, k_max=72):
    u = clamp_u(dict(u0))
    ev = agg_eval(u, n, lam)

    def obj(e):
        if e is None:
            return float("inf")
        return (e["agg"] + 100.0 * max(0.0, e["max_marginal"] - MARG_TARGET)
                + 10.0 * max(0.0, 0.05 - e["wsum"])
                + 2.0 * max(0.0, 0.5 - e["H"]))

    cur = obj(ev)
    if not math.isfinite(cur):
        return None, None
    best_u, best_ev, best_v = dict(u), ev, cur
    for step in range(steps):
        temp = 0.25 * (1 - step / steps) + 0.002
        sigma = 1.0 * (1 - step / steps) + 0.05
        cand = dict(u)
        z = rng.random()
        if z < 0.5:
            a = rng.choice(list(cand))
            cand[a] = max(cand[a], 1e-8) * math.exp(rng.gauss(0, sigma))
        elif z < 0.65:
            a = min(cand, key=cand.get)
            cand[a] = cand[a] * math.exp(rng.gauss(0, 2 * sigma))
        elif z < 0.82:
            a = rng.choice(list(cand))
            b = a ^ (1 << rng.randrange(n))
            w = cand.pop(a)
            cand[b] = cand.get(b, 0.0) + w
        elif z < 0.94 and len(cand) < k_max:
            a = rng.choice(list(cand))
            b = a ^ (1 << rng.randrange(n))
            if rng.random() < 0.5:
                b ^= 1 << rng.randrange(n)
            if b not in cand:
                cand[b] = math.exp(rng.uniform(-6, 0))
        elif len(cand) > 6:
            cand.pop(min(cand, key=cand.get))
        cand = clamp_u(cand)
        ev2 = agg_eval(cand, n, lam)
        v2 = obj(ev2)
        if v2 < cur or rng.random() < math.exp(-(v2 - cur) / temp):
            u, cur, ev = cand, v2, ev2
            if cur < best_v:
                best_u, best_ev, best_v = dict(u), ev, cur
    return best_u, best_ev


def load_round1_seed(n: int):
    f = DATA / f"aggprobe_partC_n{n}.json"
    if not f.exists():
        return None
    d = json.loads(f.read_text())
    b = d.get("best")
    if not b or not b.get("u"):
        return None
    return {int(k, 2): v for k, v in b["u"].items()}, b["lambda"], b["seed"]


def part_p3() -> dict:
    endpoints = []
    best = None
    for n in (24, 32):
        rng = random.Random(151000 + n)
        seeds = []
        r1 = load_round1_seed(n)
        if r1:
            seeds.append((f"round1_best({r1[2]})", r1[0], r1[1]))
        seeds.append(("MU_full_2.0", mu_ladder(n, min(n - 6, 18)), 2.0))
        seeds.append(("MU_full_3.5", mu_ladder(n, min(n - 6, 18)), 3.5))
        for name, u0, lam in seeds:
            t0 = time.monotonic()
            u_end, ev = anneal(u0, n, lam, rng)
            if ev is None:
                continue
            rec = {"n": n, "seed": name, "lambda": lam, "agg": ev["agg"],
                   "max_marginal": ev["max_marginal"],
                   "per_i_min": ev["per_i_min"], "argmin_i": ev["argmin_i"],
                   "H": ev["H"], "wsum": ev["wsum"],
                   "k_atoms": ev["k_atoms"],
                   "secs": round(time.monotonic() - t0, 1)}
            endpoints.append(rec)
            log(f"[P3] n={n} {name} lam={lam}: agg {ev['agg']:+.5f}, "
                f"mm {ev['max_marginal']:.4f}, per-i min "
                f"{ev['per_i_min']:+.4f} ({rec['secs']}s)")
            in_reg = ev["max_marginal"] < RECORD
            if in_reg and (best is None or ev["agg"] < best["agg"]):
                best = rec | {"u": {format(a, "b").zfill(n): w
                                    for a, w in u_end.items()}}
    nviol = sum(1 for e in endpoints
                if e["agg"] < -1e-9 and e["max_marginal"] < RECORD)
    log(f"[P3] deep climbs: {len(endpoints)} endpoints, in-regime "
        f"violations {nviol}, best in-regime "
        f"{best['agg'] if best else None}")
    save("aggprobe2_partP3.json", {"endpoints": endpoints, "best": best,
                                   "n_violations": nviol})
    return {"endpoints": endpoints, "best": best}


# ----------------------------------------------------------------------
# P4 -- lambda profiles of the extreme instances
# ----------------------------------------------------------------------

def part_p4(p2, p3) -> dict:
    grid = (0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0)
    out = {}
    jobs = []
    th2 = p2["best_by_lam"].get(2.0)
    if th2:
        jobs.append(("theta2_n96_ladder", None, 96, th2))
    if p3["best"]:
        b = p3["best"]
        jobs.append((f"p3_best_n{b['n']}",
                     {int(k, 2): v for k, v in b["u"].items()}, b["n"], None))
    for name, u_fixed, n, th in jobs:
        prof = []
        for lam in grid:
            if th is not None:
                th_n = tune_dil(n, n - 6, lam, th)
                u = mu_ladder_theta(n, n - 6, th_n)
            else:
                u = u_fixed
            ev = agg_eval(u, n, lam)
            prof.append({"lambda": lam, "agg": None if ev is None
                         else ev["agg"],
                         "max_marginal": None if ev is None
                         else ev["max_marginal"]})
        vals = [(p["lambda"], p["agg"]) for p in prof if p["agg"] is not None]
        log(f"[P4] lambda profile {name}: "
            + " ".join(f"{l}:{a:+.3f}" for l, a in vals))
        out[name] = prof
    save("aggprobe2_partP4.json", out)
    return out


# ----------------------------------------------------------------------
# P5 -- exact certification of the decisive points
# ----------------------------------------------------------------------

def tidy_u(u: dict[int, float], den: int = 10 ** 4) -> dict[int, Fraction]:
    """limit_denominator can round a tiny positive float to exactly 0; a
    zero-weight atom left in the support is a phantom atom that breaks the
    degeneracy dichotomy's support bookkeeping, so drop zeros AFTER
    rationalizing (caught by the dich_ok=False flag on the first P5 run's
    P3-best certificate; the kill certificates never had weights small
    enough to be affected)."""
    out = {a: Fraction(w).limit_denominator(den) for a, w in u.items()
           if w > 0}
    return {a: w for a, w in out.items() if w > 0}


def part_p5(p1, p2, p3) -> dict:
    t_by_lam = {2.0: Fraction(4), 3.5: Fraction(181, 16),
                5.0: Fraction(32)}
    jobs = []
    # the float-minimum ladder points at the largest n we can afford exactly
    lam_min_row = min(p1["rows"], key=lambda r: r["agg"])
    jobs.append((f"P1_min_MU({lam_min_row['n']})_lam{lam_min_row['lambda']}",
                 mu_ladder_theta(lam_min_row["n"], lam_min_row["r"],
                                 tune_dil(lam_min_row["n"], lam_min_row["r"],
                                          lam_min_row["lambda"])),
                 lam_min_row["n"], t_by_lam[lam_min_row["lambda"]]))
    tr = [r for r in p2["rows"] if r["stage"] == "transfer"]
    if tr:
        r_min = min(tr, key=lambda r: r["agg"])
        jobs.append((f"P2_min_theta_n{r_min['n']}_lam{r_min['lambda']}",
                     mu_ladder_theta(r_min["n"], r_min["n"] - 6,
                                     r_min["theta"]),
                     r_min["n"], t_by_lam[r_min["lambda"]]))
    if p3["best"]:
        b = p3["best"]
        jobs.append((f"P3_best_n{b['n']}_lam{b['lambda']}",
                     {int(k, 2): v for k, v in b["u"].items()},
                     b["n"], t_by_lam.get(b["lambda"],
                                          Fraction(2.0 ** b["lambda"]))))
    rows = []
    for name, u_f, n, t in jobs:
        u = tidy_u(u_f if isinstance(next(iter(u_f.values())), float)
                   else u_f)
        t0 = time.monotonic()
        r = exact_aggregate({a: w for a, w in u.items()}, n, t)
        verdict = ("CERTIFIED NEGATIVE (KILL)" if r["agg_hi"] < 0 else
                   "CERTIFIED POSITIVE" if r["agg_lo"] > 0 else
                   "STRADDLES 0")
        log(f"[P5] exact tidy-rational {name}: aggregate in "
            f"[{float(r['agg_lo']):+.9f}, {float(r['agg_hi']):+.9f}] "
            f"{verdict}; in-regime {r['in_regime']} (max marg "
            f"{float(r['max_marg']):.5f}); dichotomy {r['dich_ok']}; "
            f"{time.monotonic() - t0:.0f}s")
        rows.append({"name": name, "n": n, "t": str(t),
                     "agg_lo": dec_bound(r["agg_lo"], up=False),
                     "agg_hi": dec_bound(r["agg_hi"], up=True),
                     "agg_lo_float": float(r["agg_lo"]),
                     "agg_hi_float": float(r["agg_hi"]),
                     "in_regime": r["in_regime"],
                     "max_marg_float": float(r["max_marg"]),
                     "dich_ok": r["dich_ok"], "verdict": verdict,
                     "u_tidy": {format(a, "b").zfill(n): str(w)
                                for a, w in u.items()}})
    save("aggprobe2_partP5.json", {"rows": rows})
    return {"rows": rows}


# ----------------------------------------------------------------------
# P6 (--sweep) -- adversarial dilution-weight sweep: the tuned-family trend
# of P1/part B picks SOME feasible dilution weight per n, but the adversary
# owns mu, so the honest per-n ladder margin is the MINIMUM over the
# dilution weight subject to in-regime marginals.
# ----------------------------------------------------------------------

def part_p6() -> dict:
    grid = [0.5 * (2000 / 0.5) ** (k / 15) for k in range(16)]
    rows = []
    for lam in (2.0, 3.5):
        for n in (10, 16, 22, 28, 32, 40, 48, 64, 80, 96, 128):
            r = n - 6
            best = None
            th = dict(THETA_WITNESS)
            for w in grid:
                th["dil"] = w
                ev = agg_eval(mu_ladder_theta(n, r, th), n, lam)
                if ev is None or ev["max_marginal"] >= RECORD:
                    continue
                if best is None or ev["agg"] < best["agg"]:
                    best = {"wdil": w, "agg": ev["agg"],
                            "max_marginal": ev["max_marginal"],
                            "per_i_min": ev["per_i_min"]}
            rows.append({"n": n, "r": r, "lambda": lam} | (best or {}))
            log(f"[P6] MU({n},{r}) lam={lam}: min-over-wdil agg "
                f"{best['agg']:+.4f} at wdil {best['wdil']:.1f} "
                f"(mm {best['max_marginal']:.3f})" if best else
                f"[P6] MU({n},{r}) lam={lam}: no in-regime point")
        tr = [r for r in rows if r["lambda"] == lam and "agg" in r]
        log(f"[P6] lam={lam} adversarial-dil trend: "
            + " ".join(f"{r['n']}:{r['agg']:+.4f}" for r in tr))
    save("aggprobe2_partP6.json", {"rows": rows})
    return {"rows": rows}


# ----------------------------------------------------------------------
# P7 (--kill) -- drill into a float-negative aggregate found by P2: n-scan
# of the optimized theta, perturbation stability, and exact rational
# certification at exactly-rational tilt (t = 4, lambda = 2 exactly; both
# the float-exact dyadic weights and a tidy limit_denominator version).
# ----------------------------------------------------------------------

def part_p7() -> dict:
    d = json.loads((DATA / "aggprobe2_partP2.json").read_text())
    cands = [r for r in d["rows"] if r["lambda"] == 2.0]
    th = min(cands, key=lambda r: r["agg"])["theta"]
    lam = 2.0
    out: dict = {"theta": th}

    # (a) n-scan: where does the optimized ladder cross zero?
    scan = []
    for n in (48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 144, 160):
        th_n = tune_dil(n, n - 6, lam, th)
        ev = agg_eval(mu_ladder_theta(n, n - 6, th_n), n, lam)
        scan.append({"n": n, "agg": ev["agg"], "wdil": th_n["dil"],
                     "max_marginal": ev["max_marginal"],
                     "per_i_min": ev["per_i_min"], "H": ev["H"],
                     "wsum": ev["wsum"]})
        log(f"[P7a] n={n}: agg {ev['agg']:+.6f}, mm "
            f"{ev['max_marginal']:.3f}, H {ev['H']:.2f}")
    out["n_scan"] = scan
    neg = [s for s in scan if s["agg"] < -1e-9 and s["max_marginal"] < RECORD]
    tgt = min(neg, key=lambda s: s["agg"]) if neg else \
        min(scan, key=lambda s: s["agg"])
    n_t = tgt["n"]
    th_t = tune_dil(n_t, n_t - 6, lam, th)
    u_t = mu_ladder_theta(n_t, n_t - 6, th_t)
    log(f"[P7] drill target n={n_t}: float agg {tgt['agg']:+.6f}")

    # (b) perturbation stability: 20 x 3% multiplicative noise on ALL weights
    rng = random.Random(160001)
    pert = []
    for _ in range(20):
        up = {a: w * math.exp(rng.gauss(0, 0.03)) for a, w in u_t.items()}
        ev = agg_eval(up, n_t, lam)
        pert.append({"agg": ev["agg"], "max_marginal": ev["max_marginal"]})
    nneg = sum(1 for p in pert if p["agg"] < -1e-9
               and p["max_marginal"] < RECORD)
    log(f"[P7b] perturbation: {nneg}/20 remain in-regime negative "
        f"(range {min(p['agg'] for p in pert):+.6f} .. "
        f"{max(p['agg'] for p in pert):+.6f})")
    out["perturbation"] = {"rows": pert, "n_negative": nneg}

    # (c) cell floor at the deficit coordinates (underflow guard)
    pi = ORA.tilt_of_potential(u_t, lam)
    ev_t = agg_eval(u_t, n_t, lam, want_cens=True)
    out["cell_floor_at_argmin"] = ev_t["cell_floor_at_argmin"]
    log(f"[P7c] cell floor at argmin i={ev_t['argmin_i']}: "
        f"{ev_t['cell_floor_at_argmin']:.2e}")

    # (d) exact certification at t = 4 (lambda = 2 exactly)
    rows = []
    for tag, u_ex in (("float_exact_dyadic", {a: Fraction(w)
                                              for a, w in u_t.items()}),
                      ("tidy_1e4", tidy_u(u_t))):
        t0 = time.monotonic()
        r = exact_aggregate(u_ex, n_t, Fraction(4))
        verdict = ("CERTIFIED NEGATIVE (KILL)" if r["agg_hi"] < 0 else
                   "CERTIFIED POSITIVE" if r["agg_lo"] > 0 else
                   "STRADDLES 0")
        log(f"[P7d] exact {tag} n={n_t} t=4: aggregate in "
            f"[{float(r['agg_lo']):+.9f}, {float(r['agg_hi']):+.9f}] "
            f"{verdict}; in-regime {r['in_regime']} (max marg "
            f"{float(r['max_marg']):.5f}); dichotomy {r['dich_ok']}; "
            f"{time.monotonic() - t0:.0f}s")
        rows.append({"tag": tag, "n": n_t, "t": "4",
                     "agg_lo": dec_bound(r["agg_lo"], up=False),
                     "agg_hi": dec_bound(r["agg_hi"], up=True),
                     "agg_lo_float": float(r["agg_lo"]),
                     "agg_hi_float": float(r["agg_hi"]),
                     "in_regime": r["in_regime"],
                     "max_marg_float": float(r["max_marg"]),
                     "dich_ok": r["dich_ok"], "verdict": verdict})
    out["exact"] = rows
    out["u_target"] = {format(a, "b").zfill(n_t): w for a, w in u_t.items()}
    save("aggprobe2_partP7.json", out)
    return out


# ----------------------------------------------------------------------
# P8 (--cert96) -- exact certification of the theta*(2.0) ladder at the
# crossing point n = 96, float-exact dyadic weights, t = 4.  (First run as
# an inline snippet; kept here so every certificate has a reproducible
# command.)  Writes aggprobe2_n96cert.json.
# ----------------------------------------------------------------------

def part_p8() -> dict:
    d = json.loads((DATA / "aggprobe2_partP2.json").read_text())
    row = [r for r in d["rows"] if r["stage"] == "transfer"
           and r["lambda"] == 2.0 and r["n"] == 96][0]
    u = mu_ladder_theta(96, 90, row["theta"])
    t0 = time.monotonic()
    r = exact_aggregate({a: Fraction(w) for a, w in u.items()}, 96,
                        Fraction(4))
    verdict = ("CERTIFIED NEGATIVE (KILL)" if r["agg_hi"] < 0 else
               "CERTIFIED POSITIVE" if r["agg_lo"] > 0 else "STRADDLES 0")
    out = {"name": "theta2_n96_float_exact_dyadic_t4", "n": 96, "t": "4",
           "agg_lo": dec_bound(r["agg_lo"], up=False),
           "agg_hi": dec_bound(r["agg_hi"], up=True),
           "agg_lo_float": float(r["agg_lo"]),
           "agg_hi_float": float(r["agg_hi"]),
           "in_regime": r["in_regime"],
           "max_marg_float": float(r["max_marg"]),
           "dich_ok": r["dich_ok"], "verdict": verdict,
           "secs": round(time.monotonic() - t0, 1)}
    log(f"[P8] exact theta*(2.0) n=96 t=4: aggregate in "
        f"[{out['agg_lo_float']:+.9f}, {out['agg_hi_float']:+.9f}] "
        f"{verdict}; in-regime {out['in_regime']} "
        f"(max marg {out['max_marg_float']:.5f}); dichotomy "
        f"{out['dich_ok']}; {out['secs']}s")
    save("aggprobe2_n96cert.json", out)
    return out


# ----------------------------------------------------------------------
# P9 (--certp3) -- re-certification of the P3-best n = 32 endpoint with
# zero-dropped tidy weights (the first P5 pass left a phantom zero-weight
# atom after limit_denominator, tripping the exact dichotomy flag).
# Writes aggprobe2_p3best_cert.json.
# ----------------------------------------------------------------------

def part_p9() -> dict:
    d = json.loads((DATA / "aggprobe2_partP3.json").read_text())
    b = d["best"]
    u = tidy_u({int(k, 2): v for k, v in b["u"].items()})
    r = exact_aggregate(dict(u), b["n"], Fraction(4))
    verdict = ("CERTIFIED NEGATIVE (KILL)" if r["agg_hi"] < 0 else
               "CERTIFIED POSITIVE" if r["agg_lo"] > 0 else "STRADDLES 0")
    out = {"name": f"P3_best_n{b['n']}_lam2.0_tidy_zero_dropped",
           "n": b["n"], "t": "4",
           "agg_lo": dec_bound(r["agg_lo"], up=False),
           "agg_hi": dec_bound(r["agg_hi"], up=True),
           "agg_lo_float": float(r["agg_lo"]),
           "agg_hi_float": float(r["agg_hi"]),
           "in_regime": r["in_regime"],
           "max_marg_float": float(r["max_marg"]),
           "dich_ok": r["dich_ok"], "verdict": verdict}
    log(f"[P9] exact P3-best n={b['n']} t=4 (zero-dropped tidy): "
        f"[{out['agg_lo_float']:+.9f}, {out['agg_hi_float']:+.9f}] "
        f"{verdict}; dichotomy {out['dich_ok']}")
    save("aggprobe2_p3best_cert.json", out)
    return out


def main() -> None:
    if "--sweep" in sys.argv:
        part_p6()
        log(f"done in {time.monotonic() - T0:.0f}s")
        return
    if "--cert96" in sys.argv:
        part_p8()
        log(f"done in {time.monotonic() - T0:.0f}s")
        return
    if "--certp3" in sys.argv:
        part_p9()
        log(f"done in {time.monotonic() - T0:.0f}s")
        return
    if "--kill" in sys.argv:
        part_p7()
        log(f"done in {time.monotonic() - T0:.0f}s")
        return
    p1 = part_p1()
    p2 = part_p2()
    p3 = part_p3()
    part_p4(p2, p3)
    part_p5(p1, p2, p3)
    log(f"done in {time.monotonic() - T0:.0f}s")


if __name__ == "__main__":
    main()
