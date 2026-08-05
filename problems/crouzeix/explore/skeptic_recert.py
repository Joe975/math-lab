#!/usr/bin/env python3
"""Independent re-certification of 001's endpoint enclosures (skeptic pass).

Written for the review of attempt 001 from the PROBLEM.md verification
contract, deliberately sharing NO algorithm with the tier-0 harness:

  quantity              harness route (ratio.py)      this route
  --------------------  ----------------------------  ---------------------------
  lambda_max upper      bisection on an exact LDL     bisection on Sylvester's
                        definiteness test             criterion (leading principal
                                                      minors by cofactor expansion)
  lambda_max lower      same bisection                exact Rayleigh quotient at a
                                                      rationalized power-iteration
                                                      vector (any vector is valid)
  rational sqrt bounds  bisection                     integer math.isqrt scaling
  W(A) outer polygon    32 Gaussian-integer coprime   40 rounded-trig rational
                        directions, |a|,|b| <= 3      directions (denominator 64),
                                                      support values rounded UP to
                                                      denominator 10^8 (safe side)
  max|p|^2 on polygon   interval Horner on the edges  Bernstein-coefficient bounds
                        of an ordered convex hull     on ALL pairwise vertex
                                                      segments (no hull-ordering
                                                      code at all: every hull edge
                                                      is a pairwise segment, and
                                                      every pairwise segment lies
                                                      in the hull)

Everything on the certifying path is exact rational arithmetic.  Floats are
used only to CHOOSE Rayleigh vectors and inner points; soundness never
depends on those choices.

Soundness directions (the contract's crux):
  * outer2 >= (max_{W(A)} |p|)^2  because the polygon contains W(A)
    (halfplanes from upper eigenvalue bounds, offsets rounded UP, vertex
    filter has enlarging slack) and |p| is subharmonic so its max over the
    hull is attained on the boundary, which the pairwise segments cover.
  * n2_lo <= ||p(A)||^2 by Rayleigh.  Hence ratio_lo <= true ratio.
  * n2_hi >= ||p(A)||^2 by Sylvester bisection; inner2 <= (max_W |p|)^2
    because inner points x*Ax/x*x lie in W(A).  Hence ratio_hi >= true ratio.
  * refutation test in Q: n2_lo > 4*outer2 (same decision as the contract).

Usage:
  skeptic_recert.py --selftest
  skeptic_recert.py --anchors          re-certify the two calibration anchors
  skeptic_recert.py --sample [--out F] stratified re-certification of 001's
                                       endpoints (reads data/classified.jsonl
                                       and data/certified.jsonl)

Standard library only.  Deterministic (fixed RNG seed).
"""

import argparse
import cmath
import json
import math
import os
import random
import sys
import time
from fractions import Fraction as Q

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

# ----------------------------------------------------------------------
# exact complex rationals as (re, im) Fraction pairs
# ----------------------------------------------------------------------

def C(re=0, im=0):
    return (Q(re), Q(im))


def cadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def csub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def cconj(a):
    return (a[0], -a[1])


def cnorm2(a):
    return a[0] * a[0] + a[1] * a[1]


def cdiv(a, b):
    n = cnorm2(b)
    if n == 0:
        raise ZeroDivisionError
    return ((a[0] * b[0] + a[1] * b[1]) / n, (a[1] * b[0] - a[0] * b[1]) / n)


ZERO, ONE = C(0), C(1)


def mmulq(a, b):
    n = len(a)
    return [[sum_c(cmul(a[i][k], b[k][j]) for k in range(n)) for j in range(n)]
            for i in range(n)]


def sum_c(it):
    re, im = Q(0), Q(0)
    for x in it:
        re += x[0]
        im += x[1]
    return (re, im)


def madjq(a):
    n = len(a)
    return [[cconj(a[j][i]) for j in range(n)] for i in range(n)]


def maddq(a, b):
    return [[cadd(x, y) for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]


def mscaleq(s, a):
    return [[cmul(s, x) for x in row] for row in a]


def eyeq(n):
    return [[C(1) if i == j else C(0) for j in range(n)] for i in range(n)]


def herm_partq(a):
    aa = madjq(a)
    return [[((x[0] + y[0]) / 2, (x[1] + y[1]) / 2)
             for x, y in zip(ra, rb)] for ra, rb in zip(a, aa)]


def poly_matq(coeffs, a):
    n = len(a)
    out = [[C(0)] * n for _ in range(n)]
    power = eyeq(n)
    for c in coeffs:
        out = maddq(out, mscaleq(c, power))
        power = mmulq(power, a)
    return out


# ----------------------------------------------------------------------
# integer-sqrt rational bounds (no bisection)
# ----------------------------------------------------------------------

def sqrt_upper(f):
    """Rational >= sqrt(f) for f >= 0, via math.isqrt."""
    if f < 0:
        raise ValueError
    a, b = f.numerator, f.denominator
    s = math.isqrt(a * b)
    if s * s < a * b:
        s += 1
    return Q(s, b)


def sqrt_lower(f, scale=10 ** 12):
    """Rational <= sqrt(f), floor at 1/scale resolution."""
    if f <= 0:
        return Q(0)
    a, b = f.numerator, f.denominator
    return Q(math.isqrt(a * scale * scale // b), scale)


def sqrt_upper_scaled(f, scale=10 ** 12):
    if f < 0:
        raise ValueError
    a, b = f.numerator, f.denominator
    num = a * scale * scale
    q = -(-num // b)           # ceil
    s = math.isqrt(q)
    if s * s < q:
        s += 1
    return Q(s, scale)


# ----------------------------------------------------------------------
# eigenvalue bounds: Rayleigh from below, Sylvester bisection from above
# ----------------------------------------------------------------------

def det_cofactor(m):
    """Determinant of a small complex-rational matrix by cofactor expansion.
    Deliberately not elimination and not LDL."""
    n = len(m)
    if n == 1:
        return m[0][0]
    if n == 2:
        return csub(cmul(m[0][0], m[1][1]), cmul(m[0][1], m[1][0]))
    total = C(0)
    for j in range(n):
        minor = [[m[i][k] for k in range(n) if k != j] for i in range(1, n)]
        term = cmul(m[0][j], det_cofactor(minor))
        if j % 2:
            term = (-term[0], -term[1])
        total = cadd(total, term)
    return total


def sylvester_pd(h):
    """Positive definiteness via leading principal minors, all exact."""
    n = len(h)
    for k in range(1, n + 1):
        sub = [row[:k] for row in h[:k]]
        d = det_cofactor(sub)
        # Hermitian principal minors are real; the imaginary part must vanish
        if d[1] != 0:
            raise ArithmeticError("non-real principal minor of a Hermitian matrix")
        if d[0] <= 0:
            return False
    return True


def rowsum_bound(h):
    """lambda_max <= max_i sum_j |h_ij| for Hermitian h (infinity norm)."""
    best = Q(0)
    for row in h:
        s = Q(0)
        for x in row:
            s += sqrt_upper(cnorm2(x))
        best = max(best, s)
    return best


def rayleigh(h, x):
    """Exact (x* h x)/(x* x) for a rational complex vector x; real for
    Hermitian h.  A certified LOWER bound on lambda_max for any x != 0."""
    n = len(h)
    num = C(0)
    for i in range(n):
        for j in range(n):
            num = cadd(num, cmul(cmul(cconj(x[i]), h[i][j]), x[j]))
    den = sum(cnorm2(v) for v in x)
    if den == 0:
        return None
    return num[0] / den


def top_vec_float(hf, iters=200):
    """Float power iteration for the top eigenvector of a Hermitian matrix
    given as python-complex rows.  Shifted to make lambda_max dominant."""
    n = len(hf)
    shift = 1.0 + max(sum(abs(x) for x in row) for row in hf)
    m = [[hf[i][j] + (shift if i == j else 0.0) for j in range(n)]
         for i in range(n)]
    v = [complex(1.0, 0.3 * k + 0.1) for k in range(n)]
    for _ in range(iters):
        w = [sum(m[i][j] * v[j] for j in range(n)) for i in range(n)]
        nw = math.sqrt(sum(abs(x) ** 2 for x in w))
        if nw < 1e-300:
            return [1.0] + [0.0] * (n - 1)
        v = [x / nw for x in w]
    return v


def to_float_mat(h):
    return [[complex(float(x[0]), float(x[1])) for x in row] for row in h]


def rationalize_vec(vf, den=10 ** 6):
    return [C(Q(round(x.real * den), den), Q(round(x.imag * den), den))
            for x in vf]


def lam_max_bounds(h, tol=Q(1, 10 ** 7)):
    """Certified [lo, hi] around lambda_max(h): Rayleigh below, Sylvester
    bisection above."""
    x = rationalize_vec(top_vec_float(to_float_mat(h)))
    lo = rayleigh(h, x)
    if lo is None:
        lo = -rowsum_bound(h)
    hi = rowsum_bound(h)
    if hi < lo:
        hi = lo + tol
    n = len(h)
    while hi - lo > tol:
        mid = (lo + hi) / 2
        shifted = [[csub(((mid, Q(0)) if i == j else C(0)), h[i][j])
                    for j in range(n)] for i in range(n)]
        if sylvester_pd(shifted):
            hi = mid
        else:
            lo = mid
    return lo, hi


# ----------------------------------------------------------------------
# outer polygon of W(A) from rounded-trig rational directions
# ----------------------------------------------------------------------

def trig_directions(count=40, den=64):
    dirs = []
    seen = set()
    for k in range(count):
        th = 2 * math.pi * (k + 0.5) / count
        u = (Q(round(math.cos(th) * den), den), Q(round(math.sin(th) * den), den))
        if u == (Q(0), Q(0)) or u in seen:
            continue
        seen.add(u)
        dirs.append(u)
    return dirs


def support_upper(a, u, tol=Q(1, 10 ** 7), cap_den=10 ** 8):
    """Upper bound on sup Re(conj(u) z) over W(A): lambda_max(Herm(conj(u)A)),
    rounded UP to denominator cap_den (weakening keeps the halfplane valid)."""
    h = herm_partq(mscaleq(cconj(u), a))
    _, hi = lam_max_bounds(h, tol)
    return Q(-(-(hi * cap_den).numerator // (hi * cap_den).denominator), cap_den)


def outer_vertices(a, dirs, tol=Q(1, 10 ** 7), slack=Q(1, 10 ** 9)):
    planes = [(u, support_upper(a, u, tol)) for u in dirs]
    pts = set()
    m = len(planes)
    for i in range(m):
        u1, h1 = planes[i]
        for j in range(i + 1, m):
            u2, h2 = planes[j]
            det = u1[0] * u2[1] - u1[1] * u2[0]
            if det == 0:
                continue
            x = (h1 * u2[1] - h2 * u1[1]) / det
            y = (u1[0] * h2 - u2[0] * h1) / det
            if all(u[0] * x + u[1] * y <= h + slack for u, h in planes):
                pts.add((x, y))
    return sorted(pts), planes


# ----------------------------------------------------------------------
# Bernstein bound for |p|^2 on segments
# ----------------------------------------------------------------------

def binom(n, k):
    return math.comb(n, k)


def seg_abs2_poly(coeffs, z0, z1):
    """|p(z0 + s (z1 - z0))|^2 as real rational coefficients in s, ascending."""
    d = csub(z1, z0)
    cur = [C(1)]
    acc = [C(0)] * len(coeffs)
    for c in coeffs:
        for i, t in enumerate(cur):
            acc[i] = cadd(acc[i], cmul(c, t))
        nxt = [C(0)] * (len(cur) + 1)
        for i, t in enumerate(cur):
            nxt[i] = cadd(nxt[i], cmul(t, z0))
            nxt[i + 1] = cadd(nxt[i + 1], cmul(t, d))
        cur = nxt
    deg = len(acc)
    out = [Q(0)] * (2 * deg - 1)
    for i in range(deg):
        for j in range(deg):
            out[i + j] += acc[i][0] * acc[j][0] + acc[i][1] * acc[j][1]
    return out


def bernstein_max_01(q, pieces=8):
    """Upper bound for the rational polynomial q on [0,1]: subdivide, then on
    each piece the max Bernstein coefficient bounds the polynomial (Bernstein
    basis is a nonnegative partition of unity)."""
    n = len(q) - 1
    if n <= 0:
        return q[0] if q else Q(0)
    best = None
    w = Q(1, pieces)
    for k in range(pieces):
        c = Q(k, pieces)
        # shift: a'_j = sum_{m>=j} q_m C(m,j) c^(m-j) w^j
        ap = []
        for j in range(n + 1):
            s = Q(0)
            for m in range(j, n + 1):
                s += q[m] * binom(m, j) * c ** (m - j)
            ap.append(s * w ** j)
        for i in range(n + 1):
            b = Q(0)
            for j in range(i + 1):
                b += Q(binom(i, j), binom(n, j)) * ap[j]
            if best is None or b > best:
                best = b
    return best


def outer_max_abs2(coeffs, verts, pieces=8):
    """Upper bound for |p|^2 over the hull of verts: Bernstein bound over
    every pairwise segment (superset of the hull's edges, subset of the hull)."""
    if not verts:
        return Q(0)
    if len(verts) == 1:
        z = (verts[0][0], verts[0][1])
        return eval_abs2(coeffs, z)
    best = Q(0)
    m = len(verts)
    for i in range(m):
        z0 = (verts[i][0], verts[i][1])
        for j in range(i + 1, m):
            z1 = (verts[j][0], verts[j][1])
            q = seg_abs2_poly(coeffs, z0, z1)
            b = bernstein_max_01(q, pieces)
            if b > best:
                best = b
    return best


def eval_abs2(coeffs, z):
    val = C(0)
    power = C(1)
    for c in coeffs:
        val = cadd(val, cmul(c, power))
        power = cmul(power, z)
    return cnorm2(val)


# ----------------------------------------------------------------------
# inner points (for the upper end of the enclosure)
# ----------------------------------------------------------------------

def inner_points(a, extra_float_dirs=16):
    n = len(a)
    vecs = []
    for i in range(n):
        e = [C(0)] * n
        e[i] = C(1)
        vecs.append(e)
    for i in range(n):
        for j in range(i + 1, n):
            for s in (C(1), C(-1), C(0, 1), C(0, -1)):
                e = [C(0)] * n
                e[i] = C(1)
                e[j] = s
                vecs.append(e)
    af = to_float_mat(a)
    for k in range(extra_float_dirs):
        th = 2 * math.pi * k / extra_float_dirs
        e = cmath.exp(-1j * th)
        hf = [[(e * af[i][j] + (e * af[j][i]).conjugate()) / 2
               for j in range(n)] for i in range(n)]
        vecs.append(rationalize_vec(top_vec_float(hf), den=10 ** 4))
    pts = []
    for x in vecs:
        den = sum(cnorm2(v) for v in x)
        if den == 0:
            continue
        num = C(0)
        for i in range(n):
            for j in range(n):
                num = cadd(num, cmul(cmul(cconj(x[i]), a[i][j]), x[j]))
        pts.append((num[0] / den, num[1] / den))
    return pts


# ----------------------------------------------------------------------
# the enclosure
# ----------------------------------------------------------------------

def certify(a, coeffs, ndirs=40, dir_den=64, tol=Q(1, 10 ** 7), pieces=8):
    m = poly_matq(coeffs, a)
    gram = mmulq(madjq(m), m)
    n2_lo, n2_hi = lam_max_bounds(gram, tol)
    n2_lo = max(n2_lo, Q(0))

    verts, planes = outer_vertices(a, trig_directions(ndirs, dir_den), tol)
    outer2 = outer_max_abs2(coeffs, verts, pieces)

    inner2 = Q(0)
    for z in inner_points(a):
        # sanity: inner points must satisfy every halfplane (they are in W(A))
        v = eval_abs2(coeffs, z)
        if v > inner2:
            inner2 = v

    if outer2 <= 0 or inner2 <= 0:
        return {"degenerate": True, "outer2": outer2, "inner2": inner2}
    return {
        "degenerate": False,
        "norm2": (n2_lo, n2_hi),
        "outer2": outer2,
        "inner2": inner2,
        "ratio_lo": sqrt_lower(n2_lo / outer2),
        "ratio_hi": sqrt_upper_scaled(n2_hi / inner2),
        "refutes": n2_lo > 4 * outer2,
        "nverts": len(verts),
    }


# ----------------------------------------------------------------------
# endpoint rationalization -- must match 001's convention exactly
# (den-256 entrywise rounding; p normalized to unit float norm first)
# ----------------------------------------------------------------------

def endpoint_to_rational(x, den=256):
    a = [[C(0)] * 3 for _ in range(3)]
    idx = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    for k, (i, j) in enumerate(idx):
        a[i][j] = C(Q(round(x[2 * k] * den), den), Q(round(x[2 * k + 1] * den), den))
    p = [complex(x[12 + 2 * k], x[13 + 2 * k]) for k in range(4)]
    nrm = math.sqrt(sum(abs(c) ** 2 for c in p))
    coeffs = [C(Q(round(c.real / nrm * den), den), Q(round(c.imag / nrm * den), den))
              for c in p]
    return a, coeffs


# ----------------------------------------------------------------------
# reference cases
# ----------------------------------------------------------------------

def jordanq(n):
    m = [[C(0)] * n for _ in range(n)]
    for i in range(n - 1):
        m[i][i + 1] = C(1)
    return m


def diagq(vals):
    n = len(vals)
    m = [[C(0)] * n for _ in range(n)]
    for i, v in enumerate(vals):
        m[i][i] = v
    return m


def selftest():
    ok = True

    def check(label, good, detail=""):
        nonlocal ok
        ok &= good
        print(f"[{label}] {detail} {'OK' if good else 'FAIL'}")

    r = certify(jordanq(2), [C(0), C(1)])
    check("J2, p=z: encloses 2", r["ratio_lo"] <= 2 <= r["ratio_hi"],
          f"[{float(r['ratio_lo']):.9f}, {float(r['ratio_hi']):.9f}]")
    check("J2, p=z: no false refutation", not r["refutes"])
    check("J2: ||A||^2 encloses 1", r["norm2"][0] <= 1 <= r["norm2"][1])

    r3 = certify(jordanq(3), [C(0), C(1)])
    rt2 = Q(14142135623, 10 ** 10)
    check("J3, p=z: encloses sqrt2", r3["ratio_lo"] <= rt2 <= r3["ratio_hi"],
          f"[{float(r3['ratio_lo']):.9f}, {float(r3['ratio_hi']):.9f}]")

    # normal control (attack b): for diagonal A the exact denominator is
    # known: max |p| over hull of the spectrum.  Here A = diag(1, i, -1),
    # p = z: true max = 1.  The polygon must OVER-estimate (outer2 >= 1)
    # and the inner points must UNDER-estimate (inner2 <= 1), so the ratio
    # enclosure brackets the true ratio 1 from the safe sides.
    dn = certify(diagq([C(1), C(0, 1), C(-1)]), [C(0), C(1)])
    check("normal control: outer2 >= 1 >= inner2",
          dn["outer2"] >= 1 >= dn["inner2"],
          f"outer2={float(dn['outer2']):.9f} inner2={float(dn['inner2']):.9f}")
    check("normal control: ratio encloses 1",
          dn["ratio_lo"] <= 1 <= dn["ratio_hi"],
          f"[{float(dn['ratio_lo']):.9f}, {float(dn['ratio_hi']):.9f}]")
    check("normal control: no refutation", not dn["refutes"])

    # a wrong claim must be rejectable: ratio_lo of the normal case must not
    # reach 2 (the enclosure can actually discriminate)
    check("normal control: lower end well below 2", dn["ratio_lo"] < Q(3, 2))

    print("SELFTEST", "PASSED" if ok else "FAILED")
    return ok


def anchors():
    """Re-certify 001's two calibration anchors with ~80 directions."""
    out = {}
    print("[anchor 1: J3, p = z^2]  001 claims [1.980557335, 2.25], true value 2")
    r = certify(jordanq(3), [C(0), C(0), C(1)], ndirs=80, dir_den=128)
    print(f"  independent enclosure [{float(r['ratio_lo']):.9f}, "
          f"{float(r['ratio_hi']):.9f}] refutes={r['refutes']}")
    ok1 = r["ratio_lo"] <= 2 <= r["ratio_hi"] and not r["refutes"] \
        and abs(r["ratio_lo"] - Q(1980557335, 10 ** 9)) < Q(2, 100)
    print(f"  contains 2 and lower end within 0.02 of theirs: {'OK' if ok1 else 'FAIL'}")
    out["J3_z2"] = (str(r["ratio_lo"]), str(r["ratio_hi"]), ok1)

    print("[anchor 2: J2 (+) [1/4], p = z]  001 claims [1.990266652, 2.0], true value 2")
    a = jordanq(3)
    a[0][1] = C(1)
    a[1][2] = C(0)
    a[2][2] = C(Q(1, 4))
    r = certify(a, [C(0), C(1)], ndirs=80, dir_den=128)
    print(f"  independent enclosure [{float(r['ratio_lo']):.9f}, "
          f"{float(r['ratio_hi']):.9f}] refutes={r['refutes']}")
    ok2 = r["ratio_lo"] <= 2 <= r["ratio_hi"] and not r["refutes"] \
        and abs(r["ratio_lo"] - Q(1990266652, 10 ** 9)) < Q(2, 100)
    print(f"  contains 2 and lower end within 0.02 of theirs: {'OK' if ok2 else 'FAIL'}")
    out["J2c_z"] = (str(r["ratio_lo"]), str(r["ratio_hi"]), ok2)
    return ok1 and ok2


def pick_sample(classified, certified):
    """Stratified sample: all 4 below2-max, 5 one-exact, 10 near-2 including
    the certified-lower-bound record holder, the rest random.  Deterministic."""
    bycl = {}
    for r in classified:
        bycl.setdefault(r["class"], []).append(r)
    cert_by_seed = {c["seed"]: c for c in certified}
    rng = random.Random(20260804)
    sample = []
    sample += bycl.get("below2-max", [])
    one = sorted(bycl.get("one-exact", []), key=lambda r: r["seed"])
    sample += rng.sample(one, 5)
    near2 = [r for r in classified if r["class"].startswith("near2")]
    near2s = sorted(near2, key=lambda r: -float(cert_by_seed[r["seed"]]["ratio_lo"]))
    picked = near2s[:5]                      # top certified incl. record holder
    picked += [r for r in near2 if r["class"] == "near2-triple"
               if r not in picked][:2]
    rest2 = [r for r in near2 if r not in picked]
    picked += rng.sample(rest2, 10 - len(picked))
    sample += picked
    remaining = [r for r in classified if r not in sample]
    sample += rng.sample(remaining, 11)
    return sample, cert_by_seed


def run_sample(outfile):
    classified = [json.loads(l) for l in
                  open(os.path.join(DATA, "classified.jsonl"))]
    certified = [json.loads(l) for l in
                 open(os.path.join(DATA, "certified.jsonl"))]
    sample, cert_by_seed = pick_sample(classified, certified)
    print(f"re-certifying {len(sample)} endpoints "
          f"({[r['class'] for r in sample].count('below2-max')} below2-max, "
          f"{[r['class'] for r in sample].count('one-exact')} one-exact)")
    n_bad = 0
    with open(outfile, "w") as out:
        for i, rec in enumerate(sample):
            theirs = cert_by_seed[rec["seed"]]
            a, coeffs = endpoint_to_rational(rec["x"])
            t0 = time.time()
            r = certify(a, coeffs)
            dt = time.time() - t0
            their_lo = Q(theirs["ratio_lo"])
            their_hi = Q(theirs["ratio_hi_f"]).limit_denominator(10 ** 12)
            checks = {
                # their certified lower bound must not exceed my certified
                # upper bound (both bound the same true ratio)
                "their_lo_le_my_hi": their_lo <= r["ratio_hi"],
                # my certified lower bound must not exceed their upper end
                "my_lo_le_their_hi": r["ratio_lo"] <= their_hi + Q(1, 10 ** 6),
                "refutes_agree": bool(r["refutes"]) == bool(theirs["refutes"]),
                "my_lo_below_2": r["ratio_lo"] < 2,
            }
            bad = not all(checks.values())
            n_bad += bad
            row = {
                "seed": rec["seed"], "family": rec["family"],
                "class": rec["class"], "float_ratio": rec["diag"]["ratio"],
                "their_lo": theirs["ratio_lo"], "their_hi_f": theirs["ratio_hi_f"],
                "my_lo": f"{float(r['ratio_lo']):.12f}",
                "my_hi": f"{float(r['ratio_hi']):.12f}",
                "my_refutes": r["refutes"], "their_refutes": theirs["refutes"],
                "checks": checks, "consistent": not bad, "sec": round(dt, 2),
            }
            out.write(json.dumps(row) + "\n")
            out.flush()
            print(f"  [{i+1}/{len(sample)}] seed={rec['seed']} {rec['class']:12s} "
                  f"theirs=[{theirs['ratio_lo'][:8]},{theirs['ratio_hi_f']:.4f}] "
                  f"mine=[{float(r['ratio_lo']):.6f},{float(r['ratio_hi']):.6f}] "
                  f"{'OK' if not bad else '*** INCONSISTENT ***'} ({dt:.1f}s)")
    print(f"done: {len(sample)} endpoints, {n_bad} inconsistent")
    return n_bad == 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--anchors", action="store_true")
    ap.add_argument("--sample", action="store_true")
    ap.add_argument("--out", default=os.path.join(DATA, "skeptic_recert.jsonl"))
    args = ap.parse_args()
    if args.selftest:
        sys.exit(0 if selftest() else 1)
    if args.anchors:
        sys.exit(0 if anchors() else 1)
    if args.sample:
        sys.exit(0 if run_sample(args.out) else 1)
    ap.print_help()


if __name__ == "__main__":
    main()
