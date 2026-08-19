#!/usr/bin/env python3
"""Attempt 020: deep probes of the two surviving Gap-1 candidates, with
attack directions drafted by a dual-family LLM swarm (see
swarm020_attack_brief.md / swarm020_attack_values.txt; returns in
$MATHLAB_OUT/swarm/020-attack-*).

Part A -- queue 1(a) / 018 lead 1: the unsigned |sigma|-weighted control
  MM_abs at free supports n in {10, 12, 14, 16}, deep anneal (2500+ steps)
  with atom-count mutations, seeded from 018's tightest endpoint
  (gap1c_partC.json best_free, +0.0101 at n = 8) plus structural seeds
  suggested by the swarm's kill-geometry analysis (two-block nested /
  "frustrated mirrored blocks" genres) plus random restarts.
  Definite outcome per 018: an in-regime MM_abs < 0 witness kills the last
  margin-modulated reading; a floor surviving deep search licenses
  partial-proof effort.

Part B -- queue 1(b) / 018 lead 4, sharpened by the swarm's convergent
  suggestion (both families independently proposed prefix-signed
  cancellation): hunt a1-degenerate families for the lambda-window
  candidate.  Theorem C (007) gives the per-i first-order coefficient as
  || sum_a nu(a) D_a ||^2 -- a SIGNED VECTOR SUM that can cancel across
  prefixes without mu being a product.  We search small block measures
  (3-4 coordinates) for exact cancellation (a1 = 0) with the marginal cap,
  then measure the aggregate A(mu, lam) at small lambda to read a2's sign.
  The tilt kernel factorizes over disjoint coordinate blocks, so a
  product-of-blocks family has A additive across blocks (checked
  numerically in B3): a single block with a1 = 0, a2 < 0 kills the window
  candidate at EVERY n via replication.  a2 > 0 across all cancellation
  genres is structural survival evidence localizing the window question.

Standard library only.  Deterministic (fixed seeds).
Usage:
    python problems/union-closed/explore/uc_gap1_deep.py [--fast] [--partA|--partB]
Checkpoints: problems/union-closed/data/gap1deep_part*.json
"""

from __future__ import annotations

import json
import math
import random
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
sys.path.insert(0, str(HERE))

import uc_or_avg as ORA
from uc_gap1_candidates import census_mm, h2, z_plackett, dh_dlam, lamwin

FAST = "--fast" in sys.argv
ONLY_A = "--partA" in sys.argv
ONLY_B = "--partB" in sys.argv
T0 = time.monotonic()
MARG = 0.38271


def save(name: str, obj) -> None:
    DATA.mkdir(exist_ok=True)
    (DATA / name).write_text(json.dumps(obj, indent=1, sort_keys=True,
                                        default=str) + "\n", encoding="utf-8")


def log(msg: str) -> None:
    print(f"[{time.monotonic() - T0:7.0f}s] {msg}", flush=True)


# ======================================================================
# Part A -- MM_abs free-support deep attack
# ======================================================================

def best_free_seed() -> dict[int, float]:
    d = json.loads((DATA / "gap1c_partC.json").read_text())
    return {int(s, 2): w for s, w in d["best_free"]["u"].items()}


def frustrated_seed(n: int, rng: random.Random) -> dict[int, float]:
    """Swarm kill-geometry genre: selector bit 1; two future blocks whose
    intra-block alignment flips with the selector; empty-ish base atoms.
    Atoms only -- the anneal owns the weights afterwards."""
    half = (n - 1) // 2
    b1 = [1 + j for j in range(1, half + 1)]        # coords 2..half+1
    b2 = [1 + half + j for j in range(1, n - half)] # the rest
    def mask(bits):
        m = 0
        for j in bits:
            m |= 1 << (j - 1)
        return m
    atoms = {0: 30.0}
    sel = 1  # selector = coord 1
    # selector OFF: blocks aligned (both together)
    atoms[mask(b1)] = 3.0 * rng.uniform(0.5, 2)
    atoms[mask(b2)] = 3.0 * rng.uniform(0.5, 2)
    atoms[mask(b1) | mask(b2)] = 2.0 * rng.uniform(0.5, 2)
    # selector ON: blocks anti-aligned (one or the other)
    atoms[sel | mask(b1)] = 2.0 * rng.uniform(0.5, 2)
    atoms[sel | mask(b2)] = 2.0 * rng.uniform(0.5, 2)
    atoms[sel] = 1.0
    # light singles to keep coordinates nondegenerate
    for j in range(2, n + 1):
        if rng.random() < 0.5:
            atoms[1 << (j - 1)] = rng.uniform(2, 8)
    return atoms


def anneal_mm(u0: dict[int, float], n: int, lam: float, steps: int,
              rng: random.Random, atom_cap: int = 18):
    pen = lambda m: m["mm_abs"] + 300.0 * max(0.0, m["max_marginal"] - MARG)
    u = {a: w for a, w in u0.items() if w > 0}
    mm = census_mm(u, n, lam)
    if mm is None:
        return None, None
    cur = pen(mm)
    best = (dict(u), mm) if mm["max_marginal"] <= MARG else (None, None)
    for step in range(steps):
        sigma = 0.7 * (1 - step / steps) + 0.05
        cand = dict(u)
        r = rng.random()
        if r < 0.60 or len(cand) <= 3:                     # weight move
            a = rng.choice(list(cand))
            cand[a] = max(cand[a] * math.exp(rng.gauss(0, sigma + 0.2)), 1e-9)
        elif r < 0.75 and len(cand) < atom_cap:            # add atom
            cand[rng.getrandbits(n)] = math.exp(rng.uniform(-3, 1))
        elif r < 0.85 and len(cand) > 4:                   # delete atom
            del cand[rng.choice(list(cand))]
        else:                                              # bit-flip clone
            a = rng.choice(list(cand))
            b = a ^ (1 << rng.randrange(n))
            cand[b] = cand.get(b, 0.0) + cand[a] * rng.uniform(0.2, 1.0)
        mm2 = census_mm(cand, n, lam)
        if mm2 is None:
            continue
        v2 = pen(mm2)
        if v2 < cur or rng.random() < 0.02 * sigma:        # mild reheating
            if v2 < cur:
                u, cur, mm = cand, v2, mm2
            else:
                u, mm = cand, mm2
                cur = pen(mm2)
        if (mm["max_marginal"] <= MARG
                and (best[1] is None or mm["mm_abs"] < best[1]["mm_abs"])):
            best = (dict(u), dict(mm))
    return best


def part_a() -> dict:
    rows = []
    steps = 2500 if not FAST else 60
    ns = (10, 12, 14, 16) if not FAST else (10,)
    global_best = None
    seed_bf = best_free_seed()
    for n in ns:
        for lam in ((1.0, 2.0, 3.5) if not FAST else (2.0,)):
            for genre, mk in (
                ("bestfree", lambda rng: dict(seed_bf)),
                ("frustrated", lambda rng: frustrated_seed(n, rng)),
                ("random", lambda rng: {a: math.exp(rng.uniform(-1.5, 1.5))
                                        for a in rng.sample(range(1 << n),
                                                            12)}),
            ):
                rng = random.Random(202600 + 97 * n + int(10 * lam)
                                    + hash(genre) % 1000)
                u0 = mk(rng)
                got = anneal_mm(u0, n, lam, steps, rng)
                if got is None or got[0] is None:
                    log(f"[A] n={n} lam={lam} {genre}: no in-regime point")
                    continue
                u_best, mm = got
                row = {"n": n, "lambda": lam, "genre": genre,
                       "mm_abs": mm["mm_abs"], "mm_sec": mm["mm_sec"],
                       "agg_plain": mm["agg_plain"],
                       "max_marginal": mm["max_marginal"],
                       "atoms": len(u_best),
                       "u": {format(a, f"0{n}b"): w
                             for a, w in u_best.items()}}
                rows.append(row)
                log(f"[A] n={n:2d} lam={lam} {genre:10s}: best MM_abs "
                    f"{mm['mm_abs']:+.6e} (mm {mm['max_marginal']:.3f}, "
                    f"{len(u_best)} atoms)")
                if (global_best is None
                        or mm["mm_abs"] < global_best["mm_abs"]):
                    global_best = row
    out = {"rows": rows, "best": global_best, "steps": steps}
    save("gap1deep_partA.json", out)
    if global_best:
        log(f"[A] GLOBAL BEST MM_abs {global_best['mm_abs']:+.6e} at "
            f"n={global_best['n']} lam={global_best['lambda']} "
            f"({global_best['genre']})")
    return out


# ======================================================================
# Part B -- a1-cancellation block hunt for the window candidate
# ======================================================================

def prefix_stats(mu: dict[int, float], n: int):
    """Per-i: nu-weighted prefix future-bias sums (Theorem C first-order
    structure), computed directly from mu.  Returns cancel penalty
    sum_i ||sum_a nu(a) D_a||^2 and the per-i norms."""
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    per_i = []
    pens = 0.0
    for i in range(1, n + 1):
        pmask, bit = (1 << (i - 1)) - 1, 1 << (i - 1)
        pref: dict[int, list[float]] = {}
        for a, w in mu.items():
            c = a & pmask
            e = pref.setdefault(c, [0.0, 0.0])
            e[a >> (i - 1) & 1] += w
        act = {c for c, (w0, w1) in pref.items() if w0 > 0 and w1 > 0}
        if not act:
            per_i.append(0.0)
            continue
        nfut = n - i
        vec = [0.0] * nfut
        nu_tot = sum(sum(pref[c]) for c in act)
        for c in act:
            w0, w1 = pref[c]
            nu = (w0 + w1) / nu_tot
            for j in range(nfut):
                fb = 1 << (i + j)
                e1 = sum(w for a, w in mu.items()
                         if a & pmask == c and a >> (i - 1) & 1
                         and a & fb) / w1
                e0 = sum(w for a, w in mu.items()
                         if a & pmask == c and not a >> (i - 1) & 1
                         and a & fb) / w0
                vec[j] += nu * (e1 - e0)
        nrm = sum(v * v for v in vec)
        per_i.append(nrm)
        pens += nrm
    return pens, per_i


def census_pi(pi: dict[tuple[int, int], float], n: int, lam: float) -> dict:
    """Plain aggregate + MM values from a supplied coupling pi (for
    Sinkhorn couplings of a FIXED mu, unlike census_mm's potential
    convention)."""
    rho_t = 2.0 ** lam
    supp = {a for (a, _b) in pi}
    num_plain = num_abs = den_mass = den_absw = 0.0
    for i in range(1, n):
        pmask, bit = (1 << (i - 1)) - 1, 1 << (i - 1)
        tables: dict[tuple[int, int], list[float]] = {}
        for (A, B), w in pi.items():
            key = (A & pmask, B & pmask)
            tab = tables.setdefault(key, [0.0] * 4)
            tab[(2 if A & bit else 0) + (1 if B & bit else 0)] += w
        for (p00, p01, p10, p11) in tables.values():
            if min(p00, p01, p10, p11) <= 0.0:
                continue
            mass = p00 + p01 + p10 + p11
            l2or = (math.log2(p11) + math.log2(p00)
                    - math.log2(p10) - math.log2(p01))
            x, y = (p00 + p01) / mass, (p00 + p10) / mass
            d = l2or - lam
            w_der = dh_dlam(x, y, rho_t)
            num_plain += mass * d
            num_abs += mass * abs(w_der) * d
            den_mass += mass
            den_absw += mass * abs(w_der)
    if den_mass == 0.0:
        return None
    mu = ORA.mu_of(pi)
    return {"agg": num_plain / den_mass,
            "mm_abs": num_abs / (den_absw if den_absw > 0 else 1.0),
            "max_marginal": max(ORA.marginal_elements(mu, n))}


def agg_of_mu(mu: dict[int, float], n: int, lam: float):
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    pi, resid = ORA.sinkhorn(mu, lam)
    ev = census_pi(pi, n, lam)
    if ev is not None:
        ev["sinkhorn_resid"] = resid
    return ev


def cancel_residuals(mu: dict[int, float], n: int) -> list[float]:
    """The signed vector sums themselves (before squaring), concatenated
    over i -- the Gauss-Newton target."""
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    res = []
    for i in range(1, n + 1):
        pmask = (1 << (i - 1)) - 1
        pref: dict[int, list[float]] = {}
        for a, w in mu.items():
            c = a & pmask
            e = pref.setdefault(c, [0.0, 0.0])
            e[a >> (i - 1) & 1] += w
        act = {c for c, (w0, w1) in pref.items() if w0 > 0 and w1 > 0}
        if not act:
            continue
        nfut = n - i
        nu_tot = sum(sum(pref[c]) for c in act)
        for j in range(nfut):
            fb = 1 << (i + j)
            v = 0.0
            for c in act:
                w0, w1 = pref[c]
                nu = (w0 + w1) / nu_tot
                e1 = sum(w for a, w in mu.items()
                         if a & pmask == c and a >> (i - 1) & 1
                         and a & fb) / w1
                e0 = sum(w for a, w in mu.items()
                         if a & pmask == c and not a >> (i - 1) & 1
                         and a & fb) / w0
                v += nu * (e1 - e0)
            res.append(v)
    return res


def gauss_newton_polish(params: list[float], mk, n: int,
                        iters: int = 60) -> list[float]:
    """Drive cancel_residuals(mk(params)) to ~0 by damped Gauss-Newton with
    finite differences; free params only (zeros stay zero)."""
    free = [k for k, p in enumerate(params) if p > 0]
    ps = list(params)
    lam_lm = 1e-6

    def solve_step(A, g, m):
        A = [row[:] for row in A]
        g = g[:]
        for p in range(m):
            piv = max(range(p, m), key=lambda q: abs(A[q][p]))
            A[p], A[piv] = A[piv], A[p]
            g[p], g[piv] = g[piv], g[p]
            if abs(A[p][p]) < 1e-30:
                return None
            for q in range(p + 1, m):
                f = A[q][p] / A[p][p]
                for t in range(p, m):
                    A[q][t] -= f * A[p][t]
                g[q] -= f * g[p]
        dx = [0.0] * m
        for p in range(m - 1, -1, -1):
            dx[p] = (g[p] - sum(A[p][t] * dx[t]
                                for t in range(p + 1, m))) / A[p][p]
        return dx

    for _ in range(iters):
        mu = mk(ps)
        if len(mu) < 4:
            return ps
        r0 = cancel_residuals(mu, n)
        s0 = sum(v * v for v in r0)
        if s0 < 1e-26:
            break
        J = []
        for k in free:
            h = max(abs(ps[k]), 1e-4) * 1e-6
            q = list(ps)
            q[k] += h
            r1 = cancel_residuals(mk(q), n)
            J.append([(a - b) / h for a, b in zip(r1, r0)])
        m, d = len(free), len(r0)
        A0 = [[sum(J[p][t] * J[q][t] for t in range(d)) for q in range(m)]
              for p in range(m)]
        g = [sum(J[p][t] * r0[t] for t in range(d)) for p in range(m)]
        scale = 1 + max(A0[p][p] for p in range(m))
        improved = False
        for _try in range(25):
            A = [row[:] for row in A0]
            for p in range(m):
                A[p][p] += lam_lm * scale
            dx = solve_step(A, g, m)
            if dx is None:
                lam_lm *= 10
                continue
            cand = list(ps)
            ok = True
            for kk, k in enumerate(free):
                cand[k] = ps[k] - dx[kk]
                if cand[k] <= 0:
                    ok = False
                    break
            if ok:
                r2 = cancel_residuals(mk(cand), n)
                if sum(v * v for v in r2) < s0:
                    ps = cand
                    lam_lm = max(lam_lm / 3, 1e-12)
                    improved = True
                    break
            lam_lm *= 10
        if not improved:
            break
    return ps


def block_family(params: list[float]) -> dict[int, float]:
    """Parametric 3-coordinate block genre (swarm task-3 shape): coord 1 a
    selector; coords 2-3 aligned when selector off, anti-aligned when on.
    8 free weights over the 8 states; the anneal owns them."""
    w = [max(p, 0.0) for p in params]
    atoms = {0b000: w[0], 0b110: w[1], 0b011: w[2], 0b101: w[3],
             0b001: w[4], 0b010: w[5], 0b100: w[6], 0b111: w[7]}
    return {a: v for a, v in atoms.items() if v > 0}


def free_family(n: int):
    def mk(params: list[float]) -> dict[int, float]:
        return {a: p for a, p in enumerate(params) if p > 0}
    return mk


def nonproduct_score(mu: dict[int, float], n: int) -> float:
    """Max single-prefix future-bias norm: 0 iff every D_a = 0 (which holds
    at products) -- distinguishes genuine cancellation from product
    collapse."""
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    worst = 0.0
    for i in range(1, n + 1):
        pmask = (1 << (i - 1)) - 1
        pref: dict[int, list[float]] = {}
        for a, w in mu.items():
            c = a & pmask
            e = pref.setdefault(c, [0.0, 0.0])
            e[a >> (i - 1) & 1] += w
        for c, (w0, w1) in pref.items():
            if w0 <= 0 or w1 <= 0:
                continue
            for j in range(n - i):
                fb = 1 << (i + j)
                e1 = sum(w for a, w in mu.items()
                         if a & pmask == c and a >> (i - 1) & 1
                         and a & fb) / w1
                e0 = sum(w for a, w in mu.items()
                         if a & pmask == c and not a >> (i - 1) & 1
                         and a & fb) / w0
                worst = max(worst, abs(e1 - e0))
    return worst


def cancel_point(params, mk, n):
    """Anneal + polish to the cancellation manifold; returns
    (params, cancel_pen, marginals) or None."""
    def score(ps):
        mu = mk(ps)
        if len(mu) < 4:
            return None
        tot = sum(mu.values())
        marg = [sum(w for a, w in mu.items() if a >> j & 1) / tot
                for j in range(n)]
        if min(1 - m for m in marg) <= 0:
            return None
        pen_m = 40.0 * sum(max(0.0, m - MARG) for m in marg)
        r = cancel_residuals(mk(ps), n)
        cpen = sum(v * v for v in r)
        return cpen + pen_m, cpen, marg

    rng = random.Random(int(sum(params) * 1e6) % (1 << 30))
    ps = list(params)
    got = score(ps)
    if got is None:
        return None
    cur = got[0]
    for step in range(400 if not FAST else 150):
        sigma = 0.5 * (1 - step / 400) + 0.03
        cand = list(ps)
        k = rng.randrange(len(cand))
        if cand[k] == 0.0 and rng.random() < 0.9:
            continue
        cand[k] = max(cand[k] * math.exp(rng.gauss(0, sigma)), 0.0)
        got2 = score(cand)
        if got2 is not None and got2[0] < cur:
            ps, cur = cand, got2[0]
    ps = gauss_newton_polish(ps, mk, n)
    got = score(ps)
    if got is None:
        return None
    cur, cpen, marg = got
    if cpen > 1e-10 or max(marg) > MARG:
        return None
    return ps, cpen, marg


def a2_estimate(mu: dict[int, float], n: int, lam: float = 0.06):
    ev = agg_of_mu(mu, n, lam)
    if ev is None:
        return None
    return ev["agg"] / lam ** 2, ev


def part_b() -> dict:
    out = {"cancel_hunts": [], "descents": [], "lam_scans": [],
           "replication": []}

    genres = [("blk3", 3, block_family, 8),
              ("free4", 4, free_family(4), 16)]
    trials = 6 if not FAST else 1
    found = []            # (genre, n, params, mk)
    for gname, nb, mk, np_ in genres:
        for trial in range(trials):
            rng = random.Random(303000 + trial * 17 + hash(gname) % 999)
            params = [rng.uniform(0.02, 1.0) if rng.random() < 0.75 else 0.0
                      for _ in range(np_)]
            params[0] = max(params[0], 0.5)
            got = cancel_point(params, mk, nb)
            if got is None:
                continue
            ps, cpen, marg = got
            npscore = nonproduct_score(mk(ps), nb)
            a2, _ = a2_estimate(mk(ps), nb) or (None, None)
            rec = {"genre": gname, "trial": trial, "cancel_pen": cpen,
                   "marginals": marg, "nonproduct": npscore, "a2_est": a2,
                   "params": ps}
            out["cancel_hunts"].append(rec)
            log(f"[B1] {gname} trial {trial}: cancel {cpen:.1e} "
                f"nonprod {npscore:.3f} a2 {a2 if a2 is None else format(a2, '+.4e')} "
                f"margs {['%.3f' % m for m in marg]}")
            if npscore > 1e-3 and a2 is not None:
                found.append((gname, nb, ps, mk, a2))

    # B2: descend a2 WITHIN the cancellation manifold (mutate, re-polish,
    # accept if a2 drops; the real hunt -- can a2 go negative at a1 = 0?)
    dsteps = 120 if not FAST else 6
    for gname, nb, ps0, mk, a2_0 in found:
        rng = random.Random(404000 + hash(gname) % 999)
        ps, a2 = list(ps0), a2_0
        for step in range(dsteps):
            cand = list(ps)
            k = rng.randrange(len(cand))
            if cand[k] == 0.0:
                if rng.random() < 0.85:
                    continue
                cand[k] = 0.05
            cand[k] = max(cand[k] * math.exp(rng.gauss(0, 0.35)), 0.0)
            got = cancel_point(cand, mk, nb)
            if got is None:
                continue
            ps2, cpen2, marg2 = got
            if nonproduct_score(mk(ps2), nb) < 1e-3:
                continue
            est = a2_estimate(mk(ps2), nb)
            if est is None:
                continue
            a2_new, _ = est
            if a2_new < a2:
                ps, a2 = ps2, a2_new
                if a2 < 0:
                    log(f"[B2] {gname}: NEGATIVE a2 {a2:+.4e} at step {step}")
        rec = {"genre": gname, "a2_start": a2_0, "a2_end": a2,
               "params": ps,
               "mu": {format(a, f"0{nb}b"): w
                      for a, w in mk(ps).items()}}
        out["descents"].append(rec)
        log(f"[B2] {gname}: a2 descent {a2_0:+.4e} -> {a2:+.4e}")

        # lambda scan at the descent endpoint
        lams = (0.005, 0.01, 0.02, 0.04, 0.08, 0.16) if not FAST \
            else (0.02, 0.08)
        scan = {"genre": gname, "points": []}
        mu = mk(ps)
        for lam in lams:
            ev = agg_of_mu(mu, nb, lam)
            if ev is None:
                continue
            scan["points"].append({"lambda": lam, "agg": ev["agg"],
                                   "agg_over_lam2": ev["agg"] / lam ** 2,
                                   "mm_abs": ev["mm_abs"],
                                   "resid": ev["sinkhorn_resid"]})
        log(f"[B2] {gname} scan: A/lam^2 = "
            + " ".join(f"{q['agg_over_lam2']:+.3e}" for q in scan["points"]))
        out["lam_scans"].append(scan)

    # B3: replication -- evaluate the best endpoint tensored m times at
    # lam_win(n) (mean bookkeeping means this needs direct evaluation)
    if out["descents"]:
        best = min(out["descents"], key=lambda r: r["a2_end"])
        nb = len(next(iter(best["mu"])))
        mk_mu = {int(s, 2): w for s, w in best["mu"].items()}
        for m in ((2, 3) if not FAST else (2,)):
            n_tot = nb * m
            mu_m = {0: 1.0}
            for k in range(m):
                mu_m = {a | (b << (nb * k)): wa * wb
                        for a, wa in mu_m.items()
                        for b, wb in mk_mu.items()}
            if len(mu_m) > 400:
                break
            lw = lamwin(n_tot) if n_tot > 3 else 0.5
            for lam in (lw, 0.06):
                ev = agg_of_mu(mu_m, n_tot, lam)
                if ev:
                    row = {"m": m, "n": n_tot, "lambda": lam,
                           "agg": ev["agg"],
                           "max_marginal": ev["max_marginal"]}
                    out["replication"].append(row)
                    log(f"[B3] m={m} n={n_tot} lam={lam:.4f}: "
                        f"A {ev['agg']:+.5e} (mm {ev['max_marginal']:.3f})")

    save("gap1deep_partB.json", out)
    return out


# ======================================================================
# Part C -- kernel-guided first-order designer for MM_abs (swarm task 2's
# convergent reformulation: the first-order coefficient of MM_abs is a
# |sigma|-kernel-weighted Gram form, and the kernel is NOT PSD -- verified
# by direct scan; worst 2x2 violation near margins (0.2, 0.68))
# ======================================================================

def model_c1(mu: dict[int, float], n: int, lam: float):
    """First-order model of MM_abs's numerator: per coordinate i,
    sum_{a,b in N_i^2} nu(a)nu(b) K(x_a, x_b) D_a . D_b, with
    K = |dh_dlam| at rho = 2^lam and D_a the future-bias vector.
    Exact to O(lam) (verified against the census on toys)."""
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    rho = 2.0 ** lam
    c1 = 0.0
    for i in range(1, n + 1):
        pmask = (1 << (i - 1)) - 1
        pref: dict[int, list[float]] = {}
        for a, w in mu.items():
            c = a & pmask
            e = pref.setdefault(c, [0.0, 0.0])
            e[a >> (i - 1) & 1] += w
        act = [c for c, (w0, w1) in pref.items() if w0 > 0 and w1 > 0]
        if not act:
            continue
        nu_tot = sum(sum(pref[c]) for c in act)
        nfut = n - i
        info = {}
        for c in act:
            w0, w1 = pref[c]
            xz = w0 / (w0 + w1)
            if not (1e-6 < xz < 1 - 1e-6):
                continue
            Dv = []
            for j in range(nfut):
                fb = 1 << (i + j)
                e1 = sum(w for a, w in mu.items()
                         if a & pmask == c and a >> (i - 1) & 1
                         and a & fb) / w1
                e0 = sum(w for a, w in mu.items()
                         if a & pmask == c and not a >> (i - 1) & 1
                         and a & fb) / w0
                Dv.append(e1 - e0)
            info[c] = ((w0 + w1) / nu_tot, xz, Dv)
        for a in info:
            nua, xa, Da = info[a]
            for b in info:
                nub, xb, Db = info[b]
                K = abs(dh_dlam(xa, xb, rho))
                c1 += nua * nub * K * sum(u * v for u, v in zip(Da, Db))
    return c1


def part_c() -> dict:
    out = {"designs": [], "verified": []}
    lam_design = 0.05
    trials = 10 if not FAST else 2
    steps = 3000 if not FAST else 100
    best_rows = []
    for n in ((4, 5) if not FAST else (4,)):
        for trial in range(trials):
            rng = random.Random(505000 + 31 * n + trial)
            # seed: selector + response + futures, or free random
            mu = {a: math.exp(rng.uniform(-2, 1))
                  for a in rng.sample(range(1 << n), min(10, 1 << n))}

            def ev_obj(m):
                tot = sum(m.values())
                marg = [sum(w for a, w in m.items() if a >> j & 1) / tot
                        for j in range(n)]
                pen = 10.0 * sum(max(0.0, mg - MARG) for mg in marg)
                return model_c1(m, n, lam_design) + pen, max(marg)

            cur, mm = ev_obj(mu)
            for step in range(steps):
                sigma = 0.6 * (1 - step / steps) + 0.05
                cand = dict(mu)
                r = rng.random()
                if r < 0.7 or len(cand) <= 4:
                    a = rng.choice(list(cand))
                    cand[a] = max(cand[a] * math.exp(rng.gauss(0, sigma)),
                                  1e-6)
                elif r < 0.85 and len(cand) < 14:
                    cand[rng.getrandbits(n)] = math.exp(rng.uniform(-2, 0))
                else:
                    del cand[rng.choice(list(cand))]
                v2, mm2 = ev_obj(cand)
                if v2 < cur:
                    mu, cur, mm = cand, v2, mm2
            c1 = model_c1(mu, n, lam_design)
            row = {"n": n, "trial": trial, "model_c1": c1,
                   "max_marginal": mm,
                   "mu": {format(a, f"0{n}b"): w for a, w in mu.items()}}
            best_rows.append(row)
            log(f"[C] design n={n} trial {trial}: model c1 {c1:+.6e} "
                f"(mm {mm:.3f})")
    best_rows.sort(key=lambda r: r["model_c1"])
    out["designs"] = best_rows

    # verify the best negative designs against the real coupling
    for row in [r for r in best_rows if r["model_c1"] < 0][:6]:
        n = row["n"]
        mu = {int(s, 2): w for s, w in row["mu"].items()}
        ver = {"n": n, "model_c1": row["model_c1"], "points": []}
        for lam in (0.02, 0.05, 0.1, 0.2, 0.5):
            ev = agg_of_mu(mu, n, lam)
            if ev is None:
                continue
            ver["points"].append({"lambda": lam, "mm_abs": ev["mm_abs"],
                                  "agg": ev["agg"],
                                  "max_marginal": ev["max_marginal"],
                                  "resid": ev["sinkhorn_resid"]})
            log(f"[C] verify n={n} lam={lam}: MM_abs {ev['mm_abs']:+.6e} "
                f"agg {ev['agg']:+.6e} (mm {ev['max_marginal']:.3f})")
        out["verified"].append(ver)
    save("gap1deep_partC.json", out)
    return out


# ======================================================================
# Part E -- exact rational certification of the MM_abs kill
# ======================================================================

from fractions import Fraction
import uc_or_avg_skeptic as ORS
from uc_gap1_candidates import z_target_enclosure

F0 = Fraction(0)


def _imul(alo, ahi, blo, bhi):
    c = (alo * blo, alo * bhi, ahi * blo, ahi * bhi)
    return min(c), max(c)


def absw_enclosure(x: Fraction, y: Fraction, t: Fraction):
    """Certified [lo, hi] for |h2'(z_t) * dz/drho| at the Plackett root
    z_t(x, y) -- the MM_abs weight up to the global positive constant
    t*ln2, which cancels in the ratio and cannot affect the sign."""
    zlo, zhi = z_target_enclosure(x, y, t)
    # h2'(z) = log2((1-z)/z), strictly decreasing in z
    if zhi >= 1 or zlo <= 0:
        raise ValueError("z out of range")
    dlo_l, dlo_h = ORS.log2_enclosure((1 - zhi) / zhi)
    dhi_l, dhi_h = ORS.log2_enclosure((1 - zlo) / zlo)
    hp_lo, hp_hi = dlo_l, dhi_h            # h2' interval
    # dz/drho = (x-z)(y-z) / (1-x-y+2z + t(x+y-2z));  numerator decreasing
    # in z on (.., min(x,y)); denominator linear in z with slope 2-2t
    n_lo = (x - zhi) * (y - zhi)
    n_hi = (x - zlo) * (y - zlo)
    d_at = lambda z: 1 - x - y + 2 * z + t * (x + y - 2 * z)
    d1, d2 = d_at(zlo), d_at(zhi)
    d_lo, d_hi = min(d1, d2), max(d1, d2)
    if d_lo <= 0 or n_lo < 0:
        raise ValueError("enclosure degenerate")
    dz_lo, dz_hi = n_lo / d_hi, n_hi / d_lo
    s_lo, s_hi = _imul(hp_lo, hp_hi, dz_lo, dz_hi)   # h2' * dz/drho
    # |.|
    if s_lo >= 0:
        return s_lo, s_hi
    if s_hi <= 0:
        return -s_hi, -s_lo
    return F0, max(-s_lo, s_hi)


def exact_mm_abs(u_float: dict[int, float], n: int, t: Fraction,
                 rationalize=None):
    """Certified enclosure of the MM_abs numerator (up to the positive
    constant t*ln2) plus the plain aggregate, at rational tilt t, for the
    exact measure whose potential is the rationalization of u_float."""
    if rationalize is None:
        u = {a: Fraction(w) for a, w in u_float.items()}
    else:
        u = {a: rationalize(w) for a, w in u_float.items()}
        u = {a: w for a, w in u.items() if w > 0}
    cens = ORS.exact_census(u, t, n)
    wcache = {}
    num_lo = num_hi = F0
    den_lo = den_hi = F0
    plain_lo = plain_hi = F0
    wtot = F0
    dich = True
    n_rows = 0
    for c in cens:
        if c["i"] >= n:
            continue
        dich = dich and c["dich_ok"]
        for (a, b), tab in c["tables"].items():
            p00, p01, p10, p11 = tab
            if min(p00, p01, p10, p11) <= 0:
                continue
            mass = p00 + p01 + p10 + p11
            x, y = (p00 + p01) / mass, (p00 + p10) / mass
            key = (x, y)
            if key not in wcache:
                wcache[key] = absw_enclosure(x, y, t)
            wlo, whi = wcache[key]
            orr = (p11 * p00) / (p10 * p01)
            llo, lhi = ORS.log2_enclosure(orr / t)
            plo, phi = _imul(wlo, whi, llo, lhi)
            num_lo += mass * plo
            num_hi += mass * phi
            den_lo += mass * wlo
            den_hi += mass * whi
            plain_lo += mass * llo
            plain_hi += mass * lhi
            wtot += mass
            n_rows += 1
    _mu, marg = ORS.exact_marginals(u, t, n)
    out = {"num_lo": num_lo, "num_hi": num_hi,
           "den_lo": den_lo, "den_hi": den_hi,
           "agg_lo": plain_lo / wtot, "agg_hi": plain_hi / wtot,
           "max_marg": max(marg), "in_regime": max(marg) < ORS.RECORD,
           "dich_ok": dich, "n_rows": n_rows,
           "n_distinct_margins": len(wcache)}
    if den_lo > 0:
        out["mm_abs_lo"] = num_lo / den_hi if num_lo < 0 else num_lo / den_hi
        out["mm_abs_hi"] = num_hi / den_lo if num_hi > 0 else num_hi / den_lo
    return out


def part_e() -> dict:
    """Certify the two best part-C witnesses at rational tilts."""
    d = json.loads((DATA / "gap1deep_partC.json").read_text())
    allv = sorted(d["verified"], key=lambda v: min(p["mm_abs"]
                                                   for p in v["points"]))
    out = {"certs": []}
    for ver in allv[:3]:
        src = [r for r in d["designs"]
               if abs(r["model_c1"] - ver["model_c1"]) < 1e-12][0]
        n = src["n"]
        mu = {int(s, 2): w for s, w in src["mu"].items()}
        tot = sum(mu.values())
        mu = {a: w / tot for a, w in mu.items()}
        for t in (Fraction(3, 2), Fraction(6, 5)):
            lam = math.log2(t)
            pi, resid = ORA.sinkhorn(mu, lam)
            u = {A: math.sqrt(pi[(A, A)] / 2.0 ** (lam * bin(A).count("1")))
                 for A in mu}
            rec = {"n": n, "t": str(t), "lambda": lam,
                   "sinkhorn_resid": resid}
            enc = exact_mm_abs(
                u, n, t,
                rationalize=lambda w: Fraction(w).limit_denominator(10 ** 6))
            rec.update({k: (float(v) if isinstance(v, Fraction) else v)
                        for k, v in enc.items()})
            out["certs"].append(rec)
            log(f"[E] n={n} t={t}: MM_abs num in "
                f"[{float(enc['num_lo']):+.6e}, {float(enc['num_hi']):+.6e}] "
                f"agg in [{float(enc['agg_lo']):+.4e}, "
                f"{float(enc['agg_hi']):+.4e}] "
                f"max_marg {float(enc['max_marg']):.6f} "
                f"in_regime {enc['in_regime']} dich {enc['dich_ok']}")
    save("gap1deep_partE.json", out)
    return out


# ======================================================================
# Part F -- R+-constrained no-go witness (swarm task 5's reduction: any
# assembly-compatible weight equals sigma on R+ = {z_target < 1/2}, so a
# sigma-weighted violation with every history in R+ rules out EVERY
# per-history weighting at once)
# ======================================================================

def rplus_pen(mu: dict[int, float], n: int, lam: float) -> float:
    """Penalty for histories with z_target >= 1/2 (outside R+)."""
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    rho = 2.0 ** lam
    pen = 0.0
    for i in range(1, n + 1):
        pmask = (1 << (i - 1)) - 1
        pref: dict[int, list[float]] = {}
        for a, w in mu.items():
            c = a & pmask
            e = pref.setdefault(c, [0.0, 0.0])
            e[a >> (i - 1) & 1] += w
        xs = [w0 / (w0 + w1) for c, (w0, w1) in pref.items()
              if w0 > 0 and w1 > 0]
        for xa in xs:
            for xb in xs:
                zt = z_plackett(xa, xb, rho)
                if zt >= 0.5 - 0.01:
                    pen += 10.0 * (zt - (0.5 - 0.01)) + 0.1
    return pen


def part_f() -> dict:
    out = {"designs": [], "verified": [], "certs": []}
    lam_design = 0.5
    trials = 8 if not FAST else 2
    steps = 3000 if not FAST else 100
    rows = []
    for n in ((4, 5) if not FAST else (4,)):
        for trial in range(trials):
            rng = random.Random(606000 + 31 * n + trial)
            mu = {a: math.exp(rng.uniform(-2, 1))
                  for a in rng.sample(range(1 << n), min(10, 1 << n))}

            def ev_obj(m):
                tot = sum(m.values())
                marg = [sum(w for a, w in m.items() if a >> j & 1) / tot
                        for j in range(n)]
                pen = 10.0 * sum(max(0.0, mg - MARG) for mg in marg)
                pen += rplus_pen(m, n, lam_design)
                return model_c1(m, n, lam_design) + pen, max(marg)

            cur, mm = ev_obj(mu)
            for step in range(steps):
                sigma = 0.6 * (1 - step / steps) + 0.05
                cand = dict(mu)
                r = rng.random()
                if r < 0.7 or len(cand) <= 4:
                    a = rng.choice(list(cand))
                    cand[a] = max(cand[a] * math.exp(rng.gauss(0, sigma)),
                                  1e-6)
                elif r < 0.85 and len(cand) < 14:
                    cand[rng.getrandbits(n)] = math.exp(rng.uniform(-2, 0))
                else:
                    del cand[rng.choice(list(cand))]
                v2, mm2 = ev_obj(cand)
                if v2 < cur:
                    mu, cur, mm = cand, v2, mm2
            c1 = model_c1(mu, n, lam_design)
            rp = rplus_pen(mu, n, lam_design)
            row = {"n": n, "trial": trial, "model_c1": c1, "rplus_pen": rp,
                   "max_marginal": mm,
                   "mu": {format(a, f"0{n}b"): w for a, w in mu.items()}}
            rows.append(row)
            log(f"[F] design n={n} trial {trial}: model c1 {c1:+.6e} "
                f"rp_pen {rp:.2e} (mm {mm:.3f})")
    rows.sort(key=lambda r: r["model_c1"])
    out["designs"] = rows

    for row in [r for r in rows
                if r["model_c1"] < 0 and r["rplus_pen"] == 0.0][:4]:
        n = row["n"]
        mu = {int(s, 2): w for s, w in row["mu"].items()}
        tot = sum(mu.values())
        mu = {a: w / tot for a, w in mu.items()}
        ver = {"n": n, "model_c1": row["model_c1"], "points": []}
        for lam in (0.2, 0.5, 1.0):
            ev = agg_of_mu(mu, n, lam)
            if ev is None:
                continue
            ver["points"].append({"lambda": lam, "mm_abs": ev["mm_abs"],
                                  "agg": ev["agg"],
                                  "max_marginal": ev["max_marginal"]})
            log(f"[F] verify n={n} lam={lam}: MM_abs {ev['mm_abs']:+.6e} "
                f"(mm {ev['max_marginal']:.3f})")
        out["verified"].append(ver)
        # exact certification incl. certified z < 1/2 per history
        for t in (Fraction(3, 2), Fraction(2, 1)):
            lam = math.log2(t)
            pi, resid = ORA.sinkhorn(mu, lam)
            u = {A: math.sqrt(pi[(A, A)] / 2.0 ** (lam * bin(A).count("1")))
                 for A in mu}
            enc = exact_mm_abs(
                u, n, t,
                rationalize=lambda w: Fraction(w).limit_denominator(10 ** 6))
            # re-walk the census to certify z_target < 1/2 on every row
            uq = {a: Fraction(w).limit_denominator(10 ** 6)
                  for a, w in u.items()}
            cens = ORS.exact_census(uq, t, n)
            all_rplus = True
            for c in cens:
                if c["i"] >= n:
                    continue
                for (_ab, tab) in c["tables"].items():
                    p00, p01, p10, p11 = tab
                    if min(p00, p01, p10, p11) <= 0:
                        continue
                    mass = p00 + p01 + p10 + p11
                    x, y = (p00 + p01) / mass, (p00 + p10) / mass
                    zlo, zhi = z_target_enclosure(x, y, t)
                    if not zhi < Fraction(1, 2):
                        all_rplus = False
            rec = {"n": n, "t": str(t),
                   "num_lo": float(enc["num_lo"]),
                   "num_hi": float(enc["num_hi"]),
                   "max_marg": float(enc["max_marg"]),
                   "in_regime": enc["in_regime"],
                   "dich_ok": enc["dich_ok"],
                   "all_histories_Rplus_certified": all_rplus}
            out["certs"].append(rec)
            log(f"[F] cert n={n} t={t}: num [{rec['num_lo']:+.4e}, "
                f"{rec['num_hi']:+.4e}] R+ {all_rplus} "
                f"in_regime {rec['in_regime']}")
    save("gap1deep_partF.json", out)
    return out


ONLY_C = "--partC" in sys.argv
ONLY_E = "--partE" in sys.argv
ONLY_F = "--partF" in sys.argv


def main() -> None:
    log(f"start (FAST={FAST})")
    if ONLY_C:
        part_c()
        log("done")
        return
    if ONLY_E:
        part_e()
        log("done")
        return
    if ONLY_F:
        part_f()
        log("done")
        return
    if not ONLY_B:
        part_a()
    if not ONLY_A:
        part_b()
    log("done")


if __name__ == "__main__":
    main()
