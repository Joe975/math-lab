#!/usr/bin/env python3
"""Attempt 014: large-n probe of the i-AGGREGATED odds-ratio control.

Object under test (007 lead 1 / 013 claim 6 -- the restated Gap 1 after the
per-i averaged form M_i >= lambda died float-free in 007/013):

    A(mu, lambda) = sum_{i<n, N_i nonempty} w_i (M_i - lambda) / sum w_i  >= 0
                    for all mu, all lambda > 0,

with M_i the well-posed history-averaged log2 odds ratio of 007 section 1
(mass-weighted mean of log2 OR_i(a,b) over the nondegenerate history block
N_i x N_i of the tilt coupling pi_lambda(A,B) prop u(A)u(B)2^{lambda|A^B|}),
and w_i = defined history mass at coordinate i.  The i = n term is excluded
(M_n = lambda identically, 005 Prop 1, so it contributes 0 anyway).

KILL condition (exactly as certified in 013 part C, sign reversed): an
instance with A(mu, lambda) < 0 certified in exact rational arithmetic
(rational tilt t, conjecture reads sum_i sum_rows m*log2(OR/t) < 0), with
all elementwise marginals < 38271/100000 (the record regime), H(mu) > 0,
positive defined mass.  Certified positives to date live at n <= 7 only
(+1.844669 on the 10-atom witness, 013 part C); 012's budget-census reversal
past n ~ 22 is the standing warning against small-n optimism -- hence this
probe BEFORE any proof effort.

Parts:
  A  anchors: witness aggregate at n = 7 (target +1.8446); MU(n,1) exact
     invariance at n = 20, 32 (padding coordinates carry w_i = 0);
     crash/mirror embedded at n = 8..32 (<= 4-atom supports, aggregate >= 0
     is a THEOREM there via 007 Theorem B -- pure engine check); dichotomy
     flags everywhere.
  B  constructive multi-unit witness ladder MU(n, r): r crash units
     (per-unit light slice {1,c_j,n} + block-b responder {c_j,n-1}) sharing
     the witness's two-block geometry, dilution auto-tuned into the record
     regime; aggregate and per-i profile vs n at lambda in {2, 3.5, 5} --
     the margin-vs-n trend for the recorded killer genre.
  C  seeded adversarial hill-climbs per n in {8,...,32}: seeds = MU ladders,
     witness-spread, sparsified Sawin geometric mixtures, delta_0+Bern
     genre, random sparse supports; support mutations (bit toggles, atom
     birth/death) + multiplicative weight noise; objective = aggregate +
     penalties (marginal cap 0.375, defined mass, entropy); dynamic-range
     clamp (1e-9 rel) and cell-floor diagnostics per 007's underflow guard.
     Checkpointed per n.
  D  permutation attack on the best endpoints (coordinate order is part of
     mu; the aggregate is NOT permutation-invariant): random + structured
     coordinate permutations, hunting an order that flips the sign.
  E  exact rational certification (engine imported from 013's
     uc_or_avg_skeptic: Fraction kernel + certified log2 enclosures) of the
     extreme points: witness anchor at t = 181/16 vs 013's +1.844669, the
     search minimum, and the largest-n ladder instances.

Standard library only.  Deterministic (fixed seeds, fixed step counts); a
wall-clock emergency abort exists but its firing is recorded in the
checkpoint.  Runtime ~45-75 min single process.

Usage:
    python problems/union-closed/explore/uc_or_agg_probe.py [--fast]
Checkpoints: problems/union-closed/data/aggprobe_part[ABCDE]*.json
             problems/union-closed/data/aggprobe_run.log (tee by hand)
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

import uc_or_avg as ORA                    # 007's engine (float census)
import uc_or_avg_skeptic as ORS            # 013's engine (exact rational)

RECORD = 0.38271
MARG_TARGET = 0.375
CLAMP_REL = 1e-9
T0 = time.monotonic()
FAST = "--fast" in sys.argv
# emergency abort only -- steps are fixed so the run is deterministic unless
# this fires (and then the checkpoint says so)
WALL_LIMIT = 300.0 if FAST else 4800.0


def save(name: str, obj) -> None:
    DATA.mkdir(exist_ok=True)
    (DATA / name).write_text(json.dumps(obj, indent=1, sort_keys=True,
                                        default=str) + "\n", encoding="utf-8")


def log(msg: str) -> None:
    print(msg, flush=True)


# ----------------------------------------------------------------------
# aggregate evaluation (float path; reuses 007's census verbatim)
# ----------------------------------------------------------------------

def agg_eval(u: dict[int, float], n: int, lam: float, want_cens=False):
    """Returns dict with aggregate over i<n, max marginal, entropy, defined
    mass sum, per-i minimum of M_i - lambda, worst cell floor; None if no
    defined coordinate."""
    pi = ORA.tilt_of_potential(u, lam)
    cens = ORA.census(pi, n)
    rows = [(c["defined_mass"], c["M_i"], c["i"]) for c in cens
            if c["M_i"] is not None and c["i"] < n]
    if not rows:
        return None
    wsum = sum(w for w, _v, _i in rows)
    agg = sum(w * (v - lam) for w, v, _i in rows) / wsum
    mu = ORA.mu_of(pi)
    mm = max(ORA.marginal_elements(mu, n))
    per_i_min, argmin_i = min(((v - lam, i) for w, v, i in rows))
    dich = all(c["dichotomy_ok"] for c in cens)
    out = {"agg": agg, "max_marginal": mm, "H": ORA.entropy_of(mu),
           "wsum": wsum, "per_i_min": per_i_min, "argmin_i": argmin_i,
           "dichotomy_ok": dich,
           "n_defined": len(rows), "k_atoms": len(u)}
    if want_cens:
        out["cens"] = cens
        out["cell_floor_at_argmin"] = ORA.cell_floor(pi, n, argmin_i)
    return out


def clamp_u(u: dict[int, float]) -> dict[int, float]:
    """dynamic-range clamp: drop atoms below CLAMP_REL of the max weight
    (honest support change, handled by the N_i bookkeeping)."""
    mx = max(u.values())
    return {a: w for a, w in u.items() if w >= mx * CLAMP_REL}


# ----------------------------------------------------------------------
# instance builders
# ----------------------------------------------------------------------

WIT = ORA.WITNESS_U          # 007's 10-atom witness, n = 7, lambda = 3.5


def mu_ladder(n: int, r: int, wdil: float = 20.0) -> dict[int, float]:
    """MU(n, r): multi-unit witness.  Coord 1 = block marker, coords 2-4 =
    dilution singletons, response coords 5..4+r, shared futures n-1, n.
    MU(7, 1) is exactly 007's witness.  Requires 1 <= r <= n - 6."""
    assert 1 <= r <= n - 6
    b_hi, b_lo = 1 << (n - 1), 1 << (n - 2)     # coords n, n-1
    u: dict[int, float] = {
        0b1 | b_lo: 3.802705744653898,          # {1, n-1}
        0b1 | b_hi: 8.094727225702503,          # {1, n}
        0b0: 0.2526660836353666,                # {}
        b_lo: 27.28641199464732,                # {n-1}
        b_hi: 12.529780963546393,               # {n}
    }
    for j in range(r):
        c = 1 << (4 + j)                        # response coord 5+j
        u[0b1 | c | b_hi] = 0.004792968160237178   # light {1, c_j, n}
        u[c | b_lo] = 0.3609234435898425           # {c_j, n-1}
    for d in range(3):
        u[1 << (1 + d)] = wdil                  # dilution {2},{3},{4}
    return u


def mu_ladder_tuned(n: int, r: int, lam: float):
    """Bisect the shared dilution weight so max marginal <= MARG_TARGET
    (dilution coords' own marginals rise as others fall; a feasible window
    usually exists).  Returns (u, info)."""
    lo, hi = 0.5, 400.0
    best = None
    for _ in range(22):
        w = math.sqrt(lo * hi)
        u = mu_ladder(n, r, w)
        ev = agg_eval(u, n, lam)
        mm = ev["max_marginal"]
        pi = ORA.tilt_of_potential(u, lam)
        mu = ORA.mu_of(pi)
        marg = ORA.marginal_elements(mu, n)
        dil_m = max(marg[1:4])
        oth_m = max(m for j, m in enumerate(marg) if j not in (1, 2, 3))
        best = (u, {"wdil": w, "max_marginal": mm, "dil_marg": dil_m,
                    "other_marg": oth_m} | {k: ev[k] for k in
                                            ("agg", "per_i_min", "argmin_i",
                                             "H", "wsum")})
        if oth_m > MARG_TARGET and dil_m < MARG_TARGET:
            lo = w
        elif dil_m > MARG_TARGET:
            hi = w
        else:
            break
    return best


def witness_spread(n: int) -> dict[int, float]:
    """007's witness embedded in [n] with its response coordinate moved to
    the middle and futures at the top: coords (1..7) -> (1, 2, 3, 4, m,
    n-1, n) with m = n//2.  Padding coords start empty (w_i = 0); climb
    mutations populate them."""
    m = n // 2
    perm = {0: 0, 1: 1, 2: 2, 3: 3, 4: m - 1, 5: n - 2, 6: n - 1}
    u = {}
    for a, w in WIT.items():
        b = 0
        for j in range(7):
            if a >> j & 1:
                b |= 1 << perm[j]
        u[b] = w
    return u


def sawin_sparse(n: int, k: int, rng) -> dict[int, float]:
    """Sparsified Sawin-genre geometric mixture: sample k atoms from
    sum_j (1-q) q^j Bern(p_j), p_j in a geometric ladder up to 0.92, and use
    the (unnormalized) mixture density as the potential u (fixed-u trick:
    any positive u is its own marginal's potential)."""
    q = 0.55
    ps = [0.92 * (0.55 ** j) + 0.02 for j in range(6)]
    u: dict[int, float] = {}
    for _ in range(k * 3):
        j = min(int(math.log(max(rng.random(), 1e-12)) / math.log(q)), 5)
        p = ps[j]
        a = 0
        for c in range(n):
            if rng.random() < p:
                a |= 1 << c
        dens = sum((1 - q) * q ** jj *
                   (ps[jj] ** ORA.popcount(a)) *
                   ((1 - ps[jj]) ** (n - ORA.popcount(a)))
                   for jj in range(6))
        u[a] = dens
        if len(u) >= k:
            break
    mx = max(u.values())
    return {a: w / mx for a, w in u.items()}


def d0bern_sparse(n: int, k: int, rng) -> dict[int, float]:
    """delta_0 + Bern(1/2 + eps) genre, sparsified: heavy empty-set atom
    plus k sampled Bernoulli atoms."""
    eps = 0.04
    u: dict[int, float] = {0: 1.0}
    while len(u) < k + 1:
        a = 0
        for c in range(n):
            if rng.random() < 0.5 + eps:
                a |= 1 << c
        u[a] = 0.35 / k
    return u


def rand_sparse(n: int, k: int, rng) -> dict[int, float]:
    atoms = rng.sample(range(1 << n), k) if n <= 20 else \
        [rng.getrandbits(n) for _ in range(k)]
    return {a: math.exp(rng.uniform(-1.5, 1.5)) for a in set(atoms)}


# ----------------------------------------------------------------------
# hill-climb toward an aggregate violation
# ----------------------------------------------------------------------

def objective(ev):
    if ev is None:
        return float("inf")
    return (ev["agg"]
            + 25.0 * max(0.0, ev["max_marginal"] - MARG_TARGET)
            + 10.0 * max(0.0, 0.05 - ev["wsum"])
            + 2.0 * max(0.0, 0.5 - ev["H"]))


def climb(u0: dict[int, float], n: int, lam: float, rng, steps: int,
          k_max: int = 64):
    u = clamp_u(dict(u0))
    ev = agg_eval(u, n, lam)
    cur = objective(ev)
    if not math.isfinite(cur):
        return None, None
    for step in range(steps):
        if time.monotonic() - T0 > WALL_LIMIT:
            break
        sigma = 1.0 * (1 - step / steps) + 0.05
        cand = dict(u)
        z = rng.random()
        if z < 0.55:                               # weight jiggle
            a = rng.choice(list(cand))
            cand[a] = max(cand[a], 1e-8) * math.exp(rng.gauss(0, sigma))
        elif z < 0.70:                             # light-slice push
            a = min(cand, key=cand.get)
            cand[a] = cand[a] * math.exp(rng.gauss(0, 2 * sigma))
        elif z < 0.85:                             # bit toggle (support move)
            a = rng.choice(list(cand))
            b = a ^ (1 << rng.randrange(n))
            w = cand.pop(a)
            cand[b] = cand.get(b, 0.0) + w
        elif z < 0.95 and len(cand) < k_max:       # atom birth
            a = rng.choice(list(cand))
            b = a ^ (1 << rng.randrange(n))
            if rng.random() < 0.5:
                b ^= 1 << rng.randrange(n)
            if b not in cand:
                cand[b] = math.exp(rng.uniform(-6, 0))
        elif len(cand) > 6:                        # atom death
            a = min(cand, key=cand.get)
            cand.pop(a)
        cand = clamp_u(cand)
        ev2 = agg_eval(cand, n, lam)
        v2 = objective(ev2)
        if v2 < cur:
            u, ev, cur = cand, ev2, v2
    return u, ev


# ----------------------------------------------------------------------
# part A -- anchors
# ----------------------------------------------------------------------

def part_a() -> dict:
    out: dict = {}
    ev = agg_eval(WIT, 7, 3.5, want_cens=True)
    log(f"[A] witness n=7 lam=3.5: aggregate {ev['agg']:+.6f} "
        f"(007/013: +1.8446 / +1.844669), per-i min {ev['per_i_min']:+.6f} "
        f"at i={ev['argmin_i']}, max marginal {ev['max_marginal']:.4f}, "
        f"dichotomy {ev['dichotomy_ok']}")
    out["witness_n7"] = {k: v for k, v in ev.items() if k != "cens"}

    inv = {}
    for n in (20, 32):
        e2 = agg_eval(mu_ladder(n, 1), n, 3.5)
        inv[n] = e2["agg"]
        log(f"[A] MU({n},1) padding invariance: aggregate {e2['agg']:+.6f} "
            f"(must equal witness {ev['agg']:+.6f}; diff "
            f"{abs(e2['agg'] - ev['agg']):.2e})")
    out["padding_invariance"] = {"witness": ev["agg"], "per_n": inv}

    thm_b = []
    for n in (8, 16, 24, 32):
        for name, mu in (("crash", ORA.crash_mu(n)), ("mirror", ORA.mirror_mu(n))):
            for lam in (0.5, 3.5):
                u, _r = ORA.sinkhorn(mu, lam, tol=1e-12)
                # sinkhorn returns pi; evaluate aggregate directly from pi
                cens = ORA.census(u, n)
                rows = [(c["defined_mass"], c["M_i"], c["i"]) for c in cens
                        if c["M_i"] is not None and c["i"] < n]
                agg = (sum(w * (v - lam) for w, v, _i in rows)
                       / sum(w for w, _v, _i in rows))
                thm_b.append({"n": n, "family": name, "lambda": lam,
                              "agg": agg})
    worst = min(r["agg"] for r in thm_b)
    log(f"[A] <=4-atom families (crash/mirror, n up to 32, Theorem-B-safe): "
        f"{len(thm_b)} instances, min aggregate {worst:+.3e} (>= 0 expected)")
    out["theorem_b_anchors"] = {"rows": thm_b, "min_agg": worst}
    save("aggprobe_partA.json", out)
    return out


# ----------------------------------------------------------------------
# part B -- the ladder trend
# ----------------------------------------------------------------------

def part_b() -> dict:
    rows = []
    ns = (7, 10, 13, 16, 19, 22, 25, 28, 32) if not FAST else (7, 10, 13)
    for lam in (2.0, 3.5, 5.0):
        for n in ns:
            for r_mode, r in (("full", n - 6), ("fixed6", min(6, n - 6))):
                if r < 1 or (r_mode == "fixed6" and n - 6 <= 6):
                    continue
                u, info = mu_ladder_tuned(n, r, lam)
                rows.append({"n": n, "r": r, "r_mode": r_mode, "lambda": lam}
                            | info)
        for r_mode in ("full", "fixed6"):
            tr = [r for r in rows if r["lambda"] == lam
                  and r["r_mode"] == r_mode]
            if tr:
                log(f"[B] MU ladder lam={lam} {r_mode}: agg vs n: "
                    + " ".join(f"{r['n']}:{r['agg']:+.4f}" for r in tr))
                log(f"[B]   per-i min vs n: "
                    + " ".join(f"{r['n']}:{r['per_i_min']:+.4f}" for r in tr))
                log(f"[B]   max marginal vs n: "
                    + " ".join(f"{r['n']}:{r['max_marginal']:.3f}" for r in tr))
    save("aggprobe_partB.json", {"rows": rows})
    return {"rows": rows}


# ----------------------------------------------------------------------
# part C -- seeded adversarial search per n
# ----------------------------------------------------------------------

def part_c() -> dict:
    ns = (8, 12, 16, 20, 24, 28, 32) if not FAST else (8, 12)
    summary = []
    best_overall = None
    for n in ns:
        rng = random.Random(140000 + n)
        steps = 500 if n <= 16 else (400 if n <= 24 else 300)
        if FAST:
            steps = 60
        seeds = []
        r_full = min(n - 6, 14)
        seeds.append(("MU_full_3.5", mu_ladder(n, r_full), 3.5))
        seeds.append(("MU_full_2.0", mu_ladder(n, r_full), 2.0))
        if n - 6 >= 3:
            seeds.append(("MU_r3_3.5", mu_ladder(n, 3), 3.5))
        seeds.append(("witness_spread_3.5", witness_spread(n), 3.5))
        seeds.append(("sawin_sparse_2.0", sawin_sparse(n, 18, rng), 2.0))
        seeds.append(("d0bern_3.5", d0bern_sparse(n, 15, rng), 3.5))
        seeds.append(("rand_2.0", rand_sparse(n, 12, rng), 2.0))
        seeds.append(("rand_5.0", rand_sparse(n, 10, rng), 5.0))
        endpoints = []
        best = None
        for name, u0, lam in seeds:
            t_start = time.monotonic()
            u_end, ev = climb(u0, n, lam, rng, steps)
            if ev is None:
                endpoints.append({"seed": name, "lambda": lam, "agg": None})
                continue
            rec = {"seed": name, "lambda": lam, "agg": ev["agg"],
                   "max_marginal": ev["max_marginal"], "H": ev["H"],
                   "wsum": ev["wsum"], "per_i_min": ev["per_i_min"],
                   "argmin_i": ev["argmin_i"], "k_atoms": ev["k_atoms"],
                   "secs": round(time.monotonic() - t_start, 1),
                   "steps": steps}
            endpoints.append(rec)
            in_reg = ev["max_marginal"] < RECORD
            if best is None or (in_reg and ev["agg"] < best["agg"]):
                best = rec | {"u": {format(a, "b").zfill(n): w
                                    for a, w in u_end.items()}}
        nviol = sum(1 for e in endpoints if e["agg"] is not None
                    and e["agg"] < -1e-9 and e["max_marginal"] < RECORD)
        lows = sorted(e["agg"] for e in endpoints if e["agg"] is not None)
        log(f"[C] n={n}: {len(endpoints)} climbs x {steps} steps, in-regime "
            f"aggregate violations: {nviol}; three lowest endpoints: "
            + " ".join(f"{v:+.5f}" for v in lows[:3])
            + f"  (deepest per-i min {min(e['per_i_min'] for e in endpoints if e['agg'] is not None):+.4f})")
        chk = {"n": n, "steps": steps, "endpoints": endpoints,
               "n_in_regime_violations": nviol, "best": best,
               "wall_abort": time.monotonic() - T0 > WALL_LIMIT}
        save(f"aggprobe_partC_n{n}.json", chk)
        summary.append({"n": n, "lowest": lows[0] if lows else None,
                        "n_violations": nviol})
        if best is not None and (best_overall is None
                                 or best["agg"] < best_overall["agg"]):
            best_overall = best | {"n": n}
        if time.monotonic() - T0 > WALL_LIMIT:
            log(f"[C] WALL LIMIT hit after n={n}; remaining n skipped")
            break
    save("aggprobe_partC_summary.json", {"summary": summary,
                                         "best_overall": best_overall})
    return {"summary": summary, "best_overall": best_overall}


# ----------------------------------------------------------------------
# part D -- permutation attack
# ----------------------------------------------------------------------

def apply_perm(u: dict[int, float], perm: list[int], n: int):
    out = {}
    for a, w in u.items():
        b = 0
        for j in range(n):
            if a >> j & 1:
                b |= 1 << perm[j]
        out[b] = out.get(b, 0.0) + w
    return out


def part_d(best_overall) -> dict:
    rng = random.Random(140999)
    targets = [("witness_n7", WIT, 7, 3.5)]
    for n in ((13, 22, 32) if not FAST else (13,)):
        if n - 6 >= 1:
            u, _ = mu_ladder_tuned(n, min(n - 6, 14), 3.5)
            targets.append((f"MU({n})", u, n, 3.5))
    if best_overall is not None and best_overall.get("u"):
        n = best_overall["n"]
        u = {int(k, 2): v for k, v in best_overall["u"].items()}
        targets.append((f"search_best_n{n}", u, n, best_overall["lambda"]))
    out_rows = []
    for name, u, n, lam in targets:
        base = agg_eval(u, n, lam)
        vals = []
        perms = []
        for _ in range(40 if not FAST else 6):
            p = list(range(n))
            rng.shuffle(p)
            perms.append(p)
        perms.append(list(reversed(range(n))))          # full reversal
        for p in perms:
            ev = agg_eval(apply_perm(u, p, n), n, lam)
            if ev is not None:
                vals.append((ev["agg"], ev["max_marginal"]))
        worst = min(v for v, _m in vals)
        nviol = sum(1 for v, m in vals if v < -1e-9 and m < RECORD)
        log(f"[D] permutation attack {name}: base {base['agg']:+.4f}, "
            f"{len(vals)} orders, min {worst:+.4f}, in-regime violations "
            f"{nviol}")
        out_rows.append({"target": name, "n": n, "lambda": lam,
                         "base_agg": base["agg"], "min_over_perms": worst,
                         "n_violations": nviol, "n_perms": len(vals)})
    save("aggprobe_partD.json", {"rows": out_rows})
    return {"rows": out_rows}


# ----------------------------------------------------------------------
# part E -- exact rational certification of the extreme points
# ----------------------------------------------------------------------

def exact_aggregate(u_float: dict[int, float], n: int, t: Fraction):
    """Certified enclosure of the aggregate at rational tilt t for the exact
    rational measure whose potential is the float-rationalization of u.
    Reuses 013's exact_census / log2 enclosures verbatim."""
    u = {a: Fraction(w) for a, w in u_float.items()}
    cens = ORS.exact_census(u, t, n)
    lo_sum = hi_sum = Fraction(0)
    wtot = Fraction(0)
    dich = True
    for c in cens:
        if c["i"] >= n or not c["rows"]:
            continue
        dich = dich and c["dich_ok"]
        for (_a, _b, m, orr) in c["rows"]:
            l, h = ORS.log2_enclosure(orr / t)
            lo_sum += m * l
            hi_sum += m * h
        wtot += c["def_mass"]
    _mu, marg = ORS.exact_marginals(u, t, n)
    return {"agg_lo": lo_sum / wtot, "agg_hi": hi_sum / wtot,
            "max_marg": max(marg), "in_regime": max(marg) < ORS.RECORD,
            "dich_ok": dich}


def dec_bound(x: Fraction, digits: int = 21, up: bool = False) -> str:
    """Directed-rounded decimal string for a certified bound (avoids
    serializing astronomically long exact fractions)."""
    s = 10 ** digits
    num, den = x.numerator * s, x.denominator
    q = -((-num) // den) if up else num // den
    sign = "-" if q < 0 else "+"
    q = abs(q)
    return f"{sign}{q // s}.{q % s:0{digits}d}"


def part_e(best_overall) -> dict:
    t_wit = Fraction(181, 16)          # lambda = log2(181/16) ~ 3.50022
    out: dict = {}
    jobs = [("witness_n7_t181/16", WIT, 7, t_wit)]
    for n in ((22, 32) if not FAST else (13,)):
        if n - 6 >= 1:
            u, _ = mu_ladder_tuned(n, min(n - 6, 14), 3.5)
            jobs.append((f"MU({n})_t181/16", u, n, t_wit))
    jobs.append(("MU(32)_t4", mu_ladder_tuned(32, 14, 2.0)[0], 32,
                 Fraction(4)) if not FAST else
                ("MU(13)_t4", mu_ladder_tuned(13, 7, 2.0)[0], 13,
                 Fraction(4)))
    if best_overall is not None and best_overall.get("u"):
        n = best_overall["n"]
        u = {int(k, 2): v for k, v in best_overall["u"].items()}
        lam = best_overall["lambda"]
        t = {2.0: Fraction(4), 3.5: t_wit, 5.0: Fraction(32),
             0.5: Fraction(181, 128)}.get(lam, Fraction(2.0 ** lam))
        jobs.append((f"search_best_n{n}_t{t}", u, n, t))
    rows = []
    for name, u, n, t in jobs:
        t_start = time.monotonic()
        r = exact_aggregate(u, n, t)
        lo, hi = float(r["agg_lo"]), float(r["agg_hi"])
        verdict = ("CERTIFIED NEGATIVE (KILL)" if r["agg_hi"] < 0 else
                   "CERTIFIED POSITIVE" if r["agg_lo"] > 0 else "STRADDLES 0")
        log(f"[E] exact {name}: aggregate in [{lo:+.9f}, {hi:+.9f}] "
            f"{verdict}; in-regime {r['in_regime']} "
            f"(max marg {float(r['max_marg']):.5f}); dichotomy {r['dich_ok']}; "
            f"{time.monotonic() - t_start:.0f}s")
        rows.append({"name": name, "n": n, "t": str(t),
                     "agg_lo": dec_bound(r["agg_lo"], up=False),
                     "agg_hi": dec_bound(r["agg_hi"], up=True),
                     "agg_lo_float": lo, "agg_hi_float": hi,
                     "in_regime": r["in_regime"],
                     "max_marg_float": float(r["max_marg"]),
                     "dich_ok": r["dich_ok"], "verdict": verdict})
    out["rows"] = rows
    save("aggprobe_partE.json", out)
    return out


def main() -> None:
    part_a()
    part_b()
    c = part_c()
    part_d(c["best_overall"])
    part_e(c["best_overall"])
    log(f"done in {time.monotonic() - T0:.0f}s; checkpoints in {DATA}")


if __name__ == "__main__":
    main()
