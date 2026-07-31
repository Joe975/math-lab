#!/usr/bin/env python3
"""Attempt 007: skeptic review of 006 (design hunt past 135).  Default REFUTE.

Own verification stack, written for this review:

  * Exact unfolding by COMPOSED AFFINE MAPS with 2x2 rational matrices
    (z -> M z + v), a third implementation path: 006's exact layer is
    skeptic_orbit.py (002; reflected-vertex chains, complex-free tuples) and
    004's is complex composed maps.  Corridor = intersection of gate-endpoint
    projections onto n = (-tau_y, tau_x); same normalization as both, so
    exact widths must agree DIGIT-FOR-DIGIT with recorded certificates.
  * My own exact rational billiard simulator (own ray/segment solve with
    cross products, strict interior bounce check, closure of position and
    direction, struck-sequence check).
  * Certified gamma brackets via deathlaw_skeptic's interval stack (004's
    own Machin pi + Taylor enclosures -- NOT unfold.py, NOT 006's path).
  * Own word enumeration + canonicalization (independent recount of 006's
    universes), own samplers (seeded), own census-grid replica for the 001
    pinch-gap mechanism.

Nothing from design_*.py is imported in any verdict path.

Subcommands (all deterministic, seed 20260731):
  selftest    Fagnano + cross-implementation corridor agreement
  certs       re-verify the 5 exact certificates of 006 from data files
  birthdepth  W(4,3) birth vs sampler depth (is 135.0000076 a floor?)
  censusgap   replicate 001's census.py family sampler; window x-extent
  birthlaw    out-of-sample birth-law tests W(7,3), W(6,4); W(5,2) edges
  counts      independent recount of the screened universes; N1/N2 non-W
  deadsample  re-test a random subsample of "dead" words, own sampler/arcs
  coverage    formula-member aliveness at my own arc list
  newmembers  identities + case-tree search for (5,2), (6,2) via my ring
              driver (deathlaw_skeptic ring); Lemma C via its lemmac
"""

import argparse
import json
import math
import os
import random
import subprocess
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

SIDE_ENDS = {0: (1, 2), 1: (2, 0), 2: (0, 1)}
DATA = os.path.join(HERE, "..", "data")
SEED = 20260731


# ------------------------------------------------------------ affine maps
# isometry T = (M, v): z -> M z + v, M a 2x2 tuple ((a,b),(c,d)), z tuple.

def mat_apply(M, z):
    return (M[0][0] * z[0] + M[0][1] * z[1],
            M[1][0] * z[0] + M[1][1] * z[1])


def t_apply(T, z):
    M, v = T
    w = mat_apply(M, z)
    return (w[0] + v[0], w[1] + v[1])


def mat_mul(A, B):
    return ((A[0][0] * B[0][0] + A[0][1] * B[1][0],
             A[0][0] * B[0][1] + A[0][1] * B[1][1]),
            (A[1][0] * B[0][0] + A[1][1] * B[1][0],
             A[1][0] * B[0][1] + A[1][1] * B[1][1]))


def t_compose(T, S):
    """T o S: z -> T(S(z))."""
    M, v = T
    N, w = S
    Mw = mat_apply(M, w)
    return (mat_mul(M, N), (Mw[0] + v[0], Mw[1] + v[1]))


def reflection(P, Q):
    """Isometry reflecting across the line through P, Q."""
    dx, dy = Q[0] - P[0], Q[1] - P[1]
    n2 = dx * dx + dy * dy
    S = (((dx * dx - dy * dy) / n2, 2 * dx * dy / n2),
         (2 * dx * dy / n2, (dy * dy - dx * dx) / n2))
    SP = mat_apply(S, P)
    return (S, (P[0] - SP[0], P[1] - SP[1]))


def unfold(word, apex, exact):
    """(gates, tau) or (None, reason). gates = [(u, v)] images of base sides.
    exact=True: Fraction arithmetic, exact identity test."""
    one, zero = (F(1), F(0)) if exact else (1.0, 0.0)
    V = ((zero, zero), (one, zero), apex)
    T = (((one, zero), (zero, one)), (zero, zero))
    gates = []
    for s in word:
        i, j = SIDE_ENDS[s]
        gates.append((t_apply(T, V[i]), t_apply(T, V[j])))
        T = t_compose(T, reflection(V[i], V[j]))
    M, v = T
    if exact:
        ident = M == ((F(1), F(0)), (F(0), F(1)))
    else:
        ident = all(abs(M[r][c] - (1.0 if r == c else 0.0)) < 1e-9
                    for r in range(2) for c in range(2))
    if not ident:
        return None, "linear part is not the identity"
    if v == (zero, zero) or (not exact and abs(v[0]) + abs(v[1]) < 1e-14):
        return None, "zero translation"
    return gates, v


def corridor(word, apex, exact=True):
    """(lo, hi, tau, gates) with projections onto n = (-tau_y, tau_x)."""
    gates, tau = unfold(word, apex, exact)
    if gates is None:
        return None
    nx, ny = -tau[1], tau[0]
    lo = hi = None
    for u, v in gates:
        a = u[0] * nx + u[1] * ny
        b = v[0] * nx + v[1] * ny
        if a > b:
            a, b = b, a
        lo = a if lo is None else (a if a > lo else lo)
        hi = b if hi is None else (b if b < hi else hi)
    return lo, hi, tau, gates


def fwidth(word, x, y):
    c = corridor(word, (x, y), exact=False)
    return -math.inf if c is None else c[1] - c[0]


# ------------------------------------------------------------- simulator

def reflect_vec(d, e):
    """Reflect direction d across direction e."""
    n2 = e[0] * e[0] + e[1] * e[1]
    t = (d[0] * e[0] + d[1] * e[1]) / n2
    return (2 * t * e[0] - d[0], 2 * t * e[1] - d[1])


def my_verify_orbit(word, apex):
    """Own exact simulation. (ok, msg)."""
    c = corridor(word, apex, exact=True)
    if c is None:
        return False, "not a translation word"
    lo, hi, tau, gates = c
    if not lo < hi:
        return False, "corridor not positive"
    m = (lo + hi) / 2
    nx, ny = -tau[1], tau[0]
    V = ((F(0), F(0)), (F(1), F(0)), (F(apex[0]), F(apex[1])))
    i, j = SIDE_ENDS[word[0]]
    P, Q = V[i], V[j]
    den = (Q[0] - P[0]) * nx + (Q[1] - P[1]) * ny
    if den == 0:
        return False, "first gate parallel to axis"
    s = (m - (P[0] * nx + P[1] * ny)) / den
    if not (0 < s < 1):
        return False, "midline misses first gate interior"
    p0 = (P[0] + s * (Q[0] - P[0]), P[1] + s * (Q[1] - P[1]))
    d0 = reflect_vec(tau, (Q[0] - P[0], Q[1] - P[1]))
    p, d, cur = p0, d0, word[0]
    struck = []
    for _ in range(len(word)):
        best = None
        for sd in (0, 1, 2):
            if sd == cur:
                continue
            a, b = SIDE_ENDS[sd]
            U, W = V[a], V[b]
            ex, ey = W[0] - U[0], W[1] - U[1]
            det = d[0] * ey - d[1] * ex
            if det == 0:
                continue
            rx, ry = U[0] - p[0], U[1] - p[1]
            t = (rx * ey - ry * ex) / det
            r = (rx * d[1] - ry * d[0]) / det
            if t <= 0:
                continue
            if not (0 < r < 1):
                continue
            if best is None or t < best[0]:
                best = (t, sd, (p[0] + t * d[0], p[1] + t * d[1]))
        if best is None:
            return False, "escaped or vertex hit"
        _, cur, p = best
        a, b = SIDE_ENDS[cur]
        d = reflect_vec(d, (V[b][0] - V[a][0], V[b][1] - V[a][1]))
        struck.append(cur)
    if struck != list(word[1:]) + [word[0]]:
        return False, "wrong bounce sequence"
    if p != p0 or d != d0:
        return False, "did not close"
    return True, f"own exact simulation: closed periodic orbit, {len(word)} bounces"


# ----------------------------------------------------------- 135-arc facts

def on_arc_135(x, y):
    """gamma(x, y) == 135 deg exactly. Fraction. Hand-derived equivalence:
    circle (2x-1)^2+(2y+1)^2=2 <=> dot = x^2-x+y^2 = -y; then
    ca2 = x - y, cb2 = 1 - x - y, ca2*cb2 = 2 y^2 = 2 dot^2, dot < 0."""
    if not y > 0:
        return False
    dot = x * x - x + y * y
    if not dot < 0:
        return False
    return 2 * dot * dot == (x * x + y * y) * ((1 - x) ** 2 + y * y)


# ---------------------------------------------------------------- sampling

def apex_from_angles(al, be):
    ta = math.tan(math.radians(al))
    tb = math.tan(math.radians(be))
    return tb / (ta + tb), ta * tb / (ta + tb)


def arc_probe_alphas(g, uniform=60, jmax=14, depth=24, rng=None, nrand=0,
                     extra_edges=()):
    th = 180.0 - g
    als = [th * (k + 0.5) / uniform for k in range(uniform)]
    edges = [90.0 / j for j in range(1, jmax + 1)] + list(extra_edges)
    for e in edges:
        for d in range(1, depth + 1):
            eps = 0.5 * 2.0 ** -d
            for al in (e - eps, e + eps, th - e - eps, th - e + eps):
                if 0 < al < th:
                    als.append(al)
    if rng and nrand:
        als += [rng.uniform(1e-6, th - 1e-6) for _ in range(nrand)]
    return als


def best_width_on_arc(word, g, **kw):
    best, arg = -math.inf, None
    for al in arc_probe_alphas(g, **kw):
        be = (180.0 - g) - al
        if al <= 0 or be <= 0:
            continue
        x, y = apex_from_angles(al, be)
        v = fwidth(word, x, y)
        if v > best:
            best, arg = v, (al, be)
    return best, arg


TOL = 1e-12


def family_word(a, b):
    return ((0,) + (1, 2) * a + (0, 2) * b) * 2


# ------------------------------------------------------------ enumeration

def canonical(word):
    n = len(word)
    best = None
    for w in (word, word[::-1]):
        for sw in (False, True):
            ww = tuple((1 - c if c < 2 else 2) for c in w) if sw else w
            for r in range(n):
                cand = ww[r:] + ww[:r]
                if best is None or cand < best:
                    best = cand
    return best


def is_proper_power(u):
    n = len(u)
    for p in range(1, n):
        if n % p == 0 and u == u[:p] * (n // p):
            return True
    return False


def universe_all(L):
    """Canonical doubled words u^2, |u| = L odd; own recursive generator."""
    out = set()

    def rec(u):
        if len(u) == L:
            if u[-1] != u[0] and not is_proper_power(u):
                out.add(canonical(u * 2))
            return
        for c in (0, 1, 2):
            if c != u[-1]:
                rec(u + (c,))
    for c in (0, 1, 2):
        rec((c,))
    return out


def universe_s1(k):
    out = set()
    for bits in range(2 ** (k + 1)):
        y = [(bits >> i) & 1 for i in range(k + 1)]
        if y[0] == y[1]:
            continue
        u = (y[0], y[1], 2)
        for i in range(2, k + 1):
            u += (y[i], 2)
        if not is_proper_power(u):
            out.add(canonical(u * 2))
    return out


# ------------------------------------------------------------ subcommands

def cmd_selftest(_a):
    ok = True
    # Fagnano at an acute apex
    good, msg = my_verify_orbit((0, 1, 2, 0, 1, 2), (F(1, 2), F(4, 5)))
    print("[fagnano]", good, msg)
    ok &= good
    # obtuse fagnano must fail
    good, _ = my_verify_orbit((0, 1, 2, 0, 1, 2), (F(1, 2), F(1, 10)))
    print("[fagnano obtuse rejected]", not good)
    ok &= not good
    # cross-implementation exact corridor agreement (vs 002's and 004's)
    import skeptic_orbit
    import deathlaw_skeptic as DLS
    rng = random.Random(SEED)
    words = [family_word(2, 1), family_word(3, 3), family_word(4, 2),
             (0, 2, 1, 2, 1, 2, 0, 2, 0, 2, 1) * 2]
    agree = 0
    tot = 0
    for w in words:
        for _ in range(25):
            x = F(rng.randint(1, 2047), 4096)
            y = F(rng.randint(1, 1300), 4096)
            c1 = corridor(w, (x, y), exact=True)
            c2 = skeptic_orbit.corridor(w, (x, y))
            c3 = DLS.corridor_apex_exact(w, x, y)
            w1 = None if c1 is None else c1[1] - c1[0]
            w2 = None if c2 is None else c2[1] - c2[0]
            w3 = None if c3 is None else c3[1] - c3[0]
            tot += 1
            if w1 == w2 == w3:
                agree += 1
    print(f"[exact corridor 3-way agreement {agree}/{tot}]")
    ok &= agree == tot
    # W(3,3) death near 135 with own float corridor
    b, _ = best_width_on_arc(family_word(3, 3), 134.99, uniform=200, depth=20)
    d, _ = best_width_on_arc(family_word(3, 3), 135.01, uniform=200, depth=20)
    print(f"[W(3,3) alive 134.99: {b:.3g} dead 135.01: {d:.3g}]")
    ok &= b > TOL and d < TOL
    print("SELFTEST", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


def cmd_certs(args):
    import deathlaw_skeptic as DLS
    rows = []
    allok = True
    specs = [
        ("design_cert_W43_pinchgap.json", "bracket", family_word(4, 3)),
        ("design_cert_W43_gapdeep.json", "bracket", family_word(4, 3)),
        ("design_cert_W54_touch144.json", "bracket", family_word(5, 4)),
        ("design_exact135_W52.json", "exact135", family_word(5, 2)),
        ("design_exact135_N1.json", "exact135",
         ((0,) + (1, 2) * 3 + (0, 2) * 3 + (1, 2) * 2 + (0, 2) * 3) * 2),
    ]
    for fname, kind, expect_word in specs:
        d = json.load(open(os.path.join(DATA, fname)))
        w = tuple(int(c) for c in d["word"])
        x, y = F(d["apex"][0]), F(d["apex"][1])
        row = {"file": fname, "len": len(w),
               "word_matches_claimed_family": w == expect_word}
        c = corridor(w, (x, y), exact=True)
        alive = c is not None and c[1] > c[0]
        row["my_exact_corridor_positive"] = alive
        row["width_matches_record_digit_for_digit"] = (
            alive and str(c[1] - c[0]) == d["exact_corridor_width"])
        good, msg = my_verify_orbit(w, (x, y))
        row["my_own_exact_simulation"] = good
        row["sim_msg"] = msg
        if kind == "bracket":
            glo = F(d["certified_gamma_bracket_deg"][0])
            ghi = F(d["certified_gamma_bracket_deg"][1])
            lo_ok, hi_ok = DLS.certify_gamma_bounds(x, y, glo, ghi)
            row["my_gamma_gt"] = [str(glo), bool(lo_ok)]
            row["my_gamma_lt"] = [str(ghi), bool(hi_ok)]
            row["ok"] = all([alive, good, lo_ok, hi_ok,
                             row["word_matches_claimed_family"],
                             row["width_matches_record_digit_for_digit"]])
        else:
            on = on_arc_135(x, y)
            row["my_on_arc_135_exact"] = on
            circ = (2 * x - 1) ** 2 + (2 * y + 1) ** 2 == 2
            row["circle_identity_exact"] = circ
            tan_alpha = y / x
            row["tan_alpha"] = str(tan_alpha)
            row["niven_applicable"] = tan_alpha not in (0, 1, -1)
            row["ok"] = all([alive, good, on, circ,
                             row["word_matches_claimed_family"],
                             row["width_matches_record_digit_for_digit"]])
        allok &= row["ok"]
        rows.append(row)
        print(fname, "OK" if row["ok"] else row)
    out = {"command": "certs", "rows": rows, "all_ok": allok}
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    print("CERTS", "ALL CONFIRMED" if allok else "FAILURE")
    sys.exit(0 if allok else 1)


def birth_bisect(word, g_dead, g_alive, edges, iters=46, uniform=120,
                 depth=30):
    lo, hi = g_dead, g_alive
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        v, _ = best_width_on_arc(word, mid, uniform=uniform, depth=depth,
                                 extra_edges=edges)
        if v > TOL:
            hi = mid
        else:
            lo = mid
    return lo, hi


def cmd_birthdepth(args):
    """W(4,3) birth measured at increasing sampler depths."""
    w = family_word(4, 3)
    rows = []
    for depth in (16, 24, 32, 40):
        lo, hi = birth_bisect(w, 134.9, 135.2, edges=(22.5,), iters=52,
                              uniform=60, depth=depth)
        floor = 0.5 * 2.0 ** -depth
        rows.append({"depth": depth, "sampler_floor_deg": floor,
                     "birth_between": [lo, hi],
                     "birth_minus_135": 0.5 * (lo + hi) - 135.0})
        print(rows[-1])
    out = {"command": "birthdepth", "word": "W(4,3)", "rows": rows}
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)


def y_for_gamma(x, g):
    lo, hi = 1e-12, 0.5
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        ax, ay = -x, -mid
        bx, by = 1.0 - x, -mid
        gg = math.degrees(math.atan2(abs(ax * by - ay * bx),
                                     ax * bx + ay * by))
        if gg > g:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def census_grid_alive(word, g, points=400, x_lo=0.02):
    """Replica of census.py cmd_family's alive_max design (my corridor)."""
    best = -math.inf
    for k in range(points):
        x = x_lo + (0.5 - x_lo) * (k + 0.5) / points
        v = fwidth(word, x, y_for_gamma(x, g))
        if v > best:
            best = v
    return best


def cmd_censusgap(args):
    w = family_word(4, 3)
    out = {"command": "censusgap"}
    # 1. replicate 001's birth bisection with its own sampler design
    alive0 = None
    g = 135.25
    while g < 141.0:
        if census_grid_alive(w, g) > 0:
            alive0 = g
            break
        g += 0.5
    lo, hi = alive0 - 0.5, alive0
    for _ in range(16):
        mid = 0.5 * (lo + hi)
        if census_grid_alive(w, mid) > 0:
            hi = mid
        else:
            lo = mid
    out["replicated_001_birth_between"] = [lo, hi]
    print("replicated 001-style birth:", lo, hi, "(001 recorded 135.0486)")
    # 2. window x-extent at sample gammas vs the census grid
    grid_max_x = 0.02 + 0.48 * (400 - 0.5) / 400
    gap = 0.48 / 400
    out["census_grid"] = {"max_x": grid_max_x, "spacing": gap}
    rows = []
    for g in (135.005, 135.01, 135.02, 135.03, 135.04, 135.0486, 135.06):
        # find window in x by dense scan near 0.5 then edge bisection
        xs = [0.5 - 0.002 * (k + 0.5) / 4000 for k in range(4000)]
        alive_xs = [x for x in xs if fwidth(w, x, y_for_gamma(x, g)) > 0]
        if not alive_xs:
            rows.append({"gamma": g, "window": None})
            print(g, "no alive x in [0.498, 0.5) at scan resolution")
            continue
        xin, xax = min(alive_xs), max(alive_xs)
        # refine edges
        lo1, hi1 = xin - 0.002 / 4000, xin
        for _ in range(40):
            m = 0.5 * (lo1 + hi1)
            if fwidth(w, m, y_for_gamma(m, g)) > 0:
                hi1 = m
            else:
                lo1 = m
        lo2, hi2 = xax, min(xax + 0.002 / 4000, 0.5 - 1e-15)
        for _ in range(40):
            m = 0.5 * (lo2 + hi2)
            if fwidth(w, m, y_for_gamma(m, g)) > 0:
                lo2 = m
            else:
                hi2 = m
        row = {"gamma": g, "x_window": [hi1, lo2],
               "x_window_width": lo2 - hi1,
               "below_grid_max_x": lo2 < grid_max_x,
               "contains_grid_point": any(
                   hi1 <= 0.02 + 0.48 * (k + 0.5) / 400 <= lo2
                   for k in range(400))}
        rows.append(row)
        print(row)
    out["window_rows"] = rows
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)


def cmd_birthlaw(args):
    """Out-of-sample members: measure birth with own sampler, both halves."""
    rows = []
    for a, b in ((7, 3), (6, 4)):
        pred_birth = 180.0 - 90.0 * (a + b + 1) / (a * (b + 1))
        pred_death = 180.0 - 90.0 * (a + b) / (a * (b + 1))
        w = family_word(a, b)
        edges = (90.0 / a, 90.0 / (b + 1))
        mid = 0.5 * (pred_birth + pred_death)
        v_mid, arg = best_width_on_arc(w, mid, uniform=120, depth=30,
                                       extra_edges=edges)
        lo, hi = birth_bisect(w, pred_birth - 0.3, mid, edges=edges,
                              iters=50, uniform=120, depth=34)
        v_below, _ = best_width_on_arc(w, pred_birth - 1e-3, uniform=400,
                                       depth=34, extra_edges=edges)
        row = {"a": a, "b": b, "len": len(w),
               "predicted_birth": pred_birth, "predicted_death": pred_death,
               "alive_at_window_mid": v_mid > TOL, "mid_width": v_mid,
               "measured_birth_between": [lo, hi],
               "birth_minus_predicted": 0.5 * (lo + hi) - pred_birth,
               "dead_sampled_at_birth_minus_1e-3": v_below < TOL,
               "best_width_below": v_below}
        rows.append(row)
        print(row)
    # W(5,2): window edges per law = [132, 138]
    w = family_word(5, 2)
    edges = (90.0 / 5, 90.0 / 3)
    v_in, _ = best_width_on_arc(w, 132.001, uniform=120, depth=34,
                                extra_edges=edges)
    v_out, _ = best_width_on_arc(w, 131.999, uniform=400, depth=34,
                                 extra_edges=edges)
    v135, _ = best_width_on_arc(w, 135.0, uniform=120, depth=20,
                                extra_edges=edges)
    rows.append({"a": 5, "b": 2, "alive_132.001": v_in > TOL,
                 "dead_sampled_131.999": v_out < TOL,
                 "alive_135": v135 > TOL, "width_135": v135})
    print(rows[-1])
    out = {"command": "birthlaw", "rows": rows}
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)


def cmd_counts(args):
    out = {"command": "counts"}
    tot13 = 0
    per = {}
    for L in (3, 5, 7, 9, 11, 13):
        n = len(universe_all(L))
        per[L] = n
        tot13 += n
    out["all_odd_by_half_len"] = per
    out["total_half_le_13"] = tot13
    out["claim_245"] = tot13 == 245
    n15 = len(universe_all(15))
    out["half_15"] = n15
    out["claim_811"] = n15 == 811
    ns1 = sum(len(universe_s1(k)) for k in range(8, 13))
    out["s1_k8_12"] = ns1
    out["claim_2008"] = ns1 == 2008
    # N1, N2 genuinely non-W?
    N1 = ((0,) + (1, 2) * 3 + (0, 2) * 3 + (1, 2) * 2 + (0, 2) * 3) * 2
    N2 = ((0,) + (1, 2) * 3 + (0, 2) * 3 + (1, 2) * 2 + (0, 2) * 4) * 2
    res = {}
    for name, N in (("N1", N1), ("N2", N2)):
        cN = canonical(N)
        s = (len(N) - 2) // 4
        clash = [(a, s - a) for a in range(1, s)
                 if canonical(family_word(a, s - a)) == cN]
        res[name] = {"len": len(N), "W_members_same_len_checked": s - 1,
                     "canonical_clash_with_W": clash}
        print(name, res[name])
    out["nonW"] = res
    print({k: v for k, v in out.items() if k != "nonW"})
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)


def cmd_deadsample(args):
    rng = random.Random(SEED)
    uni = sorted(universe_all(3) | universe_all(5) | universe_all(7)
                 | universe_all(9) | universe_all(11) | universe_all(13))
    sample = rng.sample(uni, args.n)
    arcs = [135.0, 135.02, 136.5, 138.0, 140.0, 140.7, 142.5, 144.0,
            145.0, 148.0, 151.0, 155.0, 159.0,
            137.3, 141.9, 146.2, 149.4, 152.7, 156.5]  # 6 arcs of my own
    hits = []
    best_overall = -math.inf
    for idx, w in enumerate(sample):
        for g in arcs:
            v, arg = best_width_on_arc(w, g, uniform=90, depth=26,
                                       rng=rng, nrand=40)
            best_overall = max(best_overall, v)
            if v > TOL:
                hits.append({"word": "".join(map(str, w)), "gamma": g,
                             "width": v, "alpha": arg[0]})
        if (idx + 1) % 20 == 0:
            print(f"# {idx + 1}/{len(sample)}", flush=True)
    # attack the recorded near-misses too
    near = json.load(open(os.path.join(DATA,
                                       "design_screen_all13.json")))["near_misses"]
    near_hits = []
    for nm in near:
        w = tuple(int(c) for c in nm["word"])
        for g in (nm["gamma"], nm["gamma"] - 0.35, nm["gamma"] + 0.35):
            v, arg = best_width_on_arc(w, g, uniform=400, depth=34,
                                       rng=rng, nrand=200)
            if v > TOL:
                near_hits.append({"word": nm["word"], "gamma": g, "width": v})
    out = {"command": "deadsample", "n_sampled": len(sample),
           "arcs": arcs, "hits": hits,
           "best_width_seen": best_overall,
           "near_miss_words_attacked": len(near),
           "near_miss_hits": near_hits}
    print(f"deadsample: {len(hits)} hits in sample, "
          f"{len(near_hits)} near-miss hits, best width {best_overall:.3g}")
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    sys.exit(0)


def pick_member(gamma):
    """Shortest W(a,b) with gamma strictly inside the formula window."""
    t = (F(180) - gamma) / 90
    best = None
    for a in range(1, 100):
        for b in range(1, 100):
            v = t * a * (b + 1)
            if a + b < v < a + b + 1:
                if best is None or a + b < best[0] + best[1]:
                    best = (a, b)
    return best


def cmd_coverage(args):
    arcs = [F("90.7"), F("93.3"), F("97.9"), F("101.1"), F("107.4"),
            F("112.5"), F("118.7"), F("120"), F("124.9"), F("127.3"),
            F("130"), F("133.9"), F("135"), F("135.02"), F("135.0486"),
            F("138.6"), F("141.3"), F("144"), F("147.5"), F("150"),
            F(1080, 7), F("157.1"), F("160.9"), F("162.5"), F("164.8")]
    rows = []
    fails = 0
    for g in arcs:
        ab = pick_member(g)
        if ab is None:
            rows.append({"gamma": str(g), "member": None})
            fails += 1
            continue
        a, b = ab
        w = family_word(a, b)
        v, arg = best_width_on_arc(w, float(g), uniform=90, depth=40,
                                   extra_edges=(90.0 / a, 90.0 / (b + 1)))
        row = {"gamma": str(g), "a": a, "b": b, "len": len(w),
               "alive": v > TOL, "width": v}
        if not row["alive"]:
            fails += 1
        rows.append(row)
        print(row)
    out = {"command": "coverage", "n_arcs": len(arcs), "failures": fails,
           "rows": rows}
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    print(f"coverage: {len(arcs) - fails}/{len(arcs)} alive")
    sys.exit(0 if fails == 0 else 1)


def cmd_newmembers(args):
    import deathlaw_skeptic as DLS
    rng = random.Random(SEED)
    out = {"command": "newmembers", "members": []}
    allok = True
    for a, b in ((5, 2), (6, 2)):
        scope = a >= b >= 1 and a >= 2 and a <= 2 * b + 3
        res = DLS.identity_residuals(DLS.Ring, a, b)
        ring_ok = {k: DLS.Ring.is_zero(v) for k, v in res.items()}
        spots_ok = True
        for _ in range(10):
            u0 = F(rng.randint(2, 40), rng.randint(1, 7))
            v0 = F(rng.randint(2, 40), rng.randint(1, 7))
            if u0 == v0:
                v0 += 1
            pres = DLS.identity_residuals(DLS.PairOps(u0, v0), a, b)
            if not all(p.is_zero() for p in pres.values()):
                spots_ok = False
                break
        # case-tree adversarial sign search at theta <= theta_d
        td = float(DLS.theta_d(a, b))
        ct_hits = 0
        total = 0
        for _ in range(args.ctsamples):
            r = rng.random()
            theta = td * r if rng.random() < 0.5 else \
                td * (1 - 10 ** -rng.uniform(0, 12))
            alpha = (90.0 / a - 10 ** -rng.uniform(0, 12)
                     if rng.random() < 0.4 else theta * rng.random())
            beta = theta - alpha
            if alpha <= 0 or beta <= 0:
                continue
            total += 1
            N1, N2, N3 = DLS.signs_from_identities(a, b, alpha, beta)
            if (N1 > 1e-9 and N2 > 1e-9 and N3 > 1e-9) or \
               (N1 < -1e-9 and N2 < -1e-9 and N3 < -1e-9):
                ct_hits += 1
        ok = scope and all(ring_ok.values()) and spots_ok and ct_hits == 0
        allok &= ok
        entry = {"a": a, "b": b, "scope_precondition_a_le_2b+3": scope,
                 "ring_identities": ring_ok, "pair_spots_ok": spots_ok,
                 "casetree_points": total, "casetree_counterexamples": ct_hits,
                 "ok": ok}
        out["members"].append(entry)
        print(entry)
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    print("NEWMEMBERS", "CONFIRMED" if allok else "FAILURE")
    sys.exit(0 if allok else 1)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest").set_defaults(fn=cmd_selftest)
    for name, fn in (("certs", cmd_certs), ("birthdepth", cmd_birthdepth),
                     ("censusgap", cmd_censusgap), ("birthlaw", cmd_birthlaw),
                     ("counts", cmd_counts), ("coverage", cmd_coverage)):
        c = sub.add_parser(name)
        c.add_argument("--out")
        c.set_defaults(fn=fn)
    c = sub.add_parser("deadsample")
    c.add_argument("--n", type=int, default=60)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_deadsample)
    c = sub.add_parser("newmembers")
    c.add_argument("--ctsamples", type=int, default=300000)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_newmembers)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
