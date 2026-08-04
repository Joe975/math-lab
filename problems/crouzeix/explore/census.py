#!/usr/bin/env python3
"""Local-maxima census for the Crouzeix ratio, n = 3, deg p <= 3.

Float-only EXPLORATION tool (per the verification contract, nothing here is a
claim).  It maximizes the Crouzeix ratio

    f(A, p) = ||p(A)|| / max_{z in W(A)} |p(z)|

over 3x3 upper-triangular complex A (Schur form is WLOG: the ratio is
unitary-similarity invariant) and polynomials p of degree <= 3, by
Nelder-Mead restarts plus a compass-search polish, from several *structured*
start families.  Endpoints are re-evaluated at high boundary resolution and
classified by invariant features.  Certification of representative endpoints
is done separately by the tier-0 harness (see certify_endpoints.py).

Pure standard library.  Usage:
  census.py --pilot            a few starts, sanity output
  census.py --run N --seed S --out FILE     N starts per family, JSON lines
"""

import argparse
import cmath
import json
import math
import random
import sys
import time

# ----------------------------------------------------------------------
# 3x3 complex linear algebra (lists of lists of python complex)
# ----------------------------------------------------------------------

def mat_mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)]
            for i in range(3)]


def mat_add(a, b):
    return [[a[i][j] + b[i][j] for j in range(3)] for i in range(3)]


def mat_scale(s, a):
    return [[s * a[i][j] for j in range(3)] for i in range(3)]


def adjoint(a):
    return [[a[j][i].conjugate() for j in range(3)] for i in range(3)]


def eye():
    return [[1.0 + 0j if i == j else 0j for j in range(3)] for i in range(3)]


def fro(a):
    return math.sqrt(sum(abs(a[i][j]) ** 2 for i in range(3) for j in range(3)))


def herm_part(a):
    aa = adjoint(a)
    return [[(a[i][j] + aa[i][j]) / 2 for j in range(3)] for i in range(3)]


def poly_mat(coeffs, a):
    out = [[0j] * 3 for _ in range(3)]
    power = eye()
    for c in coeffs:
        out = mat_add(out, mat_scale(c, power))
        power = mat_mul(power, a)
    return out


# ----------------------------------------------------------------------
# Hermitian 3x3 eigenvalues (trig cubic) and top eigenvector
# ----------------------------------------------------------------------

def herm_eigs(h):
    """Real eigenvalues of Hermitian 3x3, ascending."""
    # charpoly: l^3 - c2 l^2 + c1 l - c0
    c2 = (h[0][0] + h[1][1] + h[2][2]).real
    c1 = 0.0
    for (i, j) in ((0, 1), (0, 2), (1, 2)):
        c1 += (h[i][i] * h[j][j] - h[i][j] * h[j][i]).real
    c0 = det3(h).real
    # shift l = m + c2/3 -> m^3 + p m + q
    p = c1 - c2 * c2 / 3.0
    q = -2.0 * c2 ** 3 / 27.0 + c2 * c1 / 3.0 - c0
    # m^3 + p m + q = 0, all roots real for Hermitian
    p = min(p, 0.0)
    if p > -1e-300:
        r = -q ** (1.0 / 3.0) if q >= 0 else (-q) ** (1.0 / 3.0)
        roots = [r, r, r] if q >= 0 else [r, r, r]
        # p ~ 0: triple root of shifted cubic
        m = (-q) ** (1.0 / 3.0) if q <= 0 else -(q ** (1.0 / 3.0))
        roots = [m, m, m]
    else:
        a = math.sqrt(-p / 3.0)
        arg = 3.0 * q / (2.0 * p * a)
        arg = max(-1.0, min(1.0, arg))
        phi = math.acos(arg) / 3.0
        roots = [2 * a * math.cos(phi - 2 * math.pi * k / 3.0)
                 for k in range(3)]
    roots = sorted(r + c2 / 3.0 for r in roots)
    return roots


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def cross(u, v):
    """Formal (unconjugated) cross product; lies in ker M when rows of M
    span rank <= 2 and u, v are rows of M."""
    return [u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0]]


def top_eigvec(h, lam):
    """Eigenvector of Hermitian 3x3 h for eigenvalue lam (float tolerant)."""
    m = [[h[i][j] - (lam if i == j else 0.0) for j in range(3)]
         for i in range(3)]
    best, bn = None, -1.0
    rows = m
    for (i, j) in ((0, 1), (0, 2), (1, 2)):
        v = cross(rows[i], rows[j])
        n = sum(abs(x) ** 2 for x in v)
        if n > bn:
            best, bn = v, n
    if bn > 1e-24:
        n = math.sqrt(bn)
        return [x / n for x in best]
    # rank <= 1 (multiplicity 2): any vector plain-orthogonal to the
    # largest row is in the kernel
    r = max(rows, key=lambda row: sum(abs(x) ** 2 for x in row))
    rn = sum(abs(x) ** 2 for x in r)
    if rn < 1e-24:
        return [1.0, 0.0, 0.0]
    # two independent candidates, pick the longer
    cands = [[-r[1], r[0], 0j], [-r[2], 0j, r[0]], [0j, -r[2], r[1]]]
    v = max(cands, key=lambda w: sum(abs(x) ** 2 for x in w))
    n = math.sqrt(sum(abs(x) ** 2 for x in v))
    return [x / n for x in v]


def quad_form(a, v):
    """v* A v for unit v."""
    out = 0j
    for i in range(3):
        for j in range(3):
            out += v[i].conjugate() * a[i][j] * v[j]
    return out


# ----------------------------------------------------------------------
# the ratio (float, boundary-sampled)
# ----------------------------------------------------------------------

def spec_norm(m):
    g = mat_mul(adjoint(m), m)
    return math.sqrt(max(0.0, herm_eigs(g)[2]))


def peval(coeffs, z):
    out = 0j
    for c in reversed(coeffs):
        out = out * z + c
    return out


def boundary_points(a, ntheta):
    pts = []
    for k in range(ntheta):
        th = 2 * math.pi * k / ntheta
        e = cmath.exp(-1j * th)
        h = herm_part(mat_scale(e, a))
        lam = herm_eigs(h)[2]
        v = top_eigvec(h, lam)
        pts.append(quad_form(a, v))
    return pts


def crouzeix_ratio_float(a, coeffs, ntheta=96, pts=None):
    num = spec_norm(poly_mat(coeffs, a))
    if pts is None:
        pts = boundary_points(a, ntheta)
    den = max(abs(peval(coeffs, z)) for z in pts)
    if den < 1e-14:
        return 0.0
    return num / den


# ----------------------------------------------------------------------
# parametrization: x in R^20 -> (A upper triangular, p)
# ----------------------------------------------------------------------

def unpack(x):
    a = [[0j] * 3 for _ in range(3)]
    a[0][0] = complex(x[0], x[1])
    a[1][1] = complex(x[2], x[3])
    a[2][2] = complex(x[4], x[5])
    a[0][1] = complex(x[6], x[7])
    a[0][2] = complex(x[8], x[9])
    a[1][2] = complex(x[10], x[11])
    p = [complex(x[12 + 2 * k], x[13 + 2 * k]) for k in range(4)]
    return a, p


def pack(a, p):
    x = [a[0][0].real, a[0][0].imag, a[1][1].real, a[1][1].imag,
         a[2][2].real, a[2][2].imag, a[0][1].real, a[0][1].imag,
         a[0][2].real, a[0][2].imag, a[1][2].real, a[1][2].imag]
    for c in p:
        x += [c.real, c.imag]
    return x


def objective(x, ntheta=96):
    a, p = unpack(x)
    na = fro(a)
    if na < 1e-9 or na > 50.0:
        return -1.0
    np_ = math.sqrt(sum(abs(c) ** 2 for c in p))
    if np_ < 1e-9:
        return -1.0
    p = [c / np_ for c in p]          # pure gauge: ratio is p-scale invariant
    return crouzeix_ratio_float(a, p, ntheta)


# ----------------------------------------------------------------------
# optimizers: Nelder-Mead (maximize) + compass polish
# ----------------------------------------------------------------------

def nelder_mead(f, x0, step=0.25, maxit=2500, ftol=1e-10):
    n = len(x0)
    simplex = [list(x0)]
    for i in range(n):
        y = list(x0)
        y[i] += step
        simplex.append(y)
    vals = [f(s) for s in simplex]
    for _ in range(maxit):
        order = sorted(range(n + 1), key=lambda i: -vals[i])
        simplex = [simplex[i] for i in order]
        vals = [vals[i] for i in order]
        if vals[0] - vals[-1] < ftol:
            break
        cen = [sum(s[i] for s in simplex[:-1]) / n for i in range(n)]
        worst = simplex[-1]
        refl = [2 * cen[i] - worst[i] for i in range(n)]
        fr = f(refl)
        if fr > vals[0]:
            exp = [3 * cen[i] - 2 * worst[i] for i in range(n)]
            fe = f(exp)
            if fe > fr:
                simplex[-1], vals[-1] = exp, fe
            else:
                simplex[-1], vals[-1] = refl, fr
        elif fr > vals[-2]:
            simplex[-1], vals[-1] = refl, fr
        else:
            con = [(cen[i] + worst[i]) / 2 for i in range(n)]
            fc = f(con)
            if fc > vals[-1]:
                simplex[-1], vals[-1] = con, fc
            else:
                best = simplex[0]
                for i in range(1, n + 1):
                    simplex[i] = [(simplex[i][j] + best[j]) / 2
                                  for j in range(n)]
                    vals[i] = f(simplex[i])
    order = sorted(range(n + 1), key=lambda i: -vals[i])
    return simplex[order[0]], vals[order[0]]


def compass(f, x0, v0, step=0.02, minstep=1e-6, maxit=4000):
    x, v = list(x0), v0
    n = len(x)
    it = 0
    while step > minstep and it < maxit:
        improved = False
        for i in range(n):
            for s in (step, -step):
                it += 1
                y = list(x)
                y[i] += s
                fy = f(y)
                if fy > v + 1e-13:
                    x, v = y, fy
                    improved = True
                    break
        if not improved:
            step /= 2.0
    return x, v


# ----------------------------------------------------------------------
# start families
# ----------------------------------------------------------------------

def rnd_c(rng, s=1.0):
    return complex(rng.gauss(0, s), rng.gauss(0, s))


def schur_of(a):
    """No exact Schur needed: our families are built upper triangular."""
    return a


def start(family, rng):
    if family == "random":
        # baseline: dense-equivalent random upper triangular + random p
        a = [[rnd_c(rng) if j >= i else 0j for j in range(3)] for i in range(3)]
        p = [rnd_c(rng) for _ in range(4)]
    elif family == "jordan3":
        # perturbed nilpotent Jordan block, p near z^2 (known 3x3 extremal)
        eps = 0.3
        a = [[eps * rnd_c(rng) if j >= i else 0j for j in range(3)]
             for i in range(3)]
        a[0][1] += 1.0
        a[1][2] += 1.0
        p = [eps * rnd_c(rng) for _ in range(4)]
        p[2] += 1.0
    elif family == "jordan2":
        # perturbed J2 (+) [c], p near z (the embedded 2x2 extremal)
        eps = 0.3
        a = [[eps * rnd_c(rng) if j >= i else 0j for j in range(3)]
             for i in range(3)]
        a[0][1] += 1.0
        a[2][2] += 0.3 * rnd_c(rng)
        p = [eps * rnd_c(rng) for _ in range(4)]
        p[1] += 1.0
    elif family == "nearnormal":
        # normal (diagonal) + small strictly-upper perturbation; p a cubic
        # with roots near the spectrum (peaks |p| at several boundary points)
        d = [rnd_c(rng) for _ in range(3)]
        eps = 0.15
        a = [[0j] * 3 for _ in range(3)]
        for i in range(3):
            a[i][i] = d[i]
        a[0][1], a[0][2], a[1][2] = (eps * rnd_c(rng), eps * rnd_c(rng),
                                     eps * rnd_c(rng))
        # p with roots at perturbed midpoints of the spectrum
        r1 = (d[0] + d[1]) / 2 + 0.2 * rnd_c(rng)
        r2 = (d[1] + d[2]) / 2 + 0.2 * rnd_c(rng)
        r3 = (d[0] + d[2]) / 2 + 0.2 * rnd_c(rng)
        p = poly_from_roots([r1, r2, r3])
    elif family == "companion":
        # companion matrix of a random cubic (upper Hessenberg made upper
        # triangular is wrong; instead use its transpose-free triangular
        # surrogate: random triangular with the same spectrum), p lacunary z^3
        roots = [rnd_c(rng) for _ in range(3)]
        a = [[0j] * 3 for _ in range(3)]
        for i in range(3):
            a[i][i] = roots[i]
        # strong coupling, as in a companion form
        a[0][1], a[1][2], a[0][2] = 1.0 + 0j, 1.0 + 0j, rnd_c(rng, 0.5)
        p = [0j, 0j, 0j, 1.0 + 0j]
        p = [p[k] + 0.1 * rnd_c(rng) for k in range(4)]
    elif family == "peaks":
        # elliptic-ish W(A); p = product of linear factors just outside the
        # boundary at spread angles -> multiple near-coincident |p| peaks
        a = [[0j] * 3 for _ in range(3)]
        a[0][0], a[1][1], a[2][2] = -1.0 + 0j, 0j, 1.0 + 0j
        a[0][1] = 0.8 + 0.4 * rng.random()
        a[1][2] = 0.8 + 0.4 * rng.random()
        a[0][2] = 0.3 * rnd_c(rng)
        pts = boundary_points(a, 24)
        zs = rng.sample(pts, 3)
        # push each root outward by 10-40%
        c = sum(pts, 0j) / len(pts)
        roots = [c + (z - c) * (1.1 + 0.3 * rng.random()) for z in zs]
        p = poly_from_roots(roots)
    elif family == "real":
        # real upper triangular, real polynomial
        a = [[complex(rng.gauss(0, 1), 0) if j >= i else 0j
              for j in range(3)] for i in range(3)]
        p = [complex(rng.gauss(0, 1), 0) for _ in range(4)]
    else:
        raise ValueError(family)
    return pack(a, p)


def poly_from_roots(roots):
    p = [1.0 + 0j]
    for r in roots:
        new = [0j] * (len(p) + 1)
        for i, c in enumerate(p):
            new[i + 1] += c
            new[i] -= r * c
        p = new
    return p


# ----------------------------------------------------------------------
# endpoint diagnostics (invariant features)
# ----------------------------------------------------------------------

def diagnostics(x):
    a, p = unpack(x)
    np_ = math.sqrt(sum(abs(c) ** 2 for c in p))
    p = [c / np_ for c in p]
    pts = boundary_points(a, 512)
    r = crouzeix_ratio_float(a, p, pts=pts)
    na = fro(a)
    # non-normality (unitary-similarity invariant)
    comm = mat_add(mat_mul(adjoint(a), a), mat_scale(-1, mat_mul(a, adjoint(a))))
    nonnormal = fro(comm) / (na * na)
    # eigenvalues of A: diagonal (A is triangular throughout the search)
    evs = [a[0][0], a[1][1], a[2][2]]
    scale = max(abs(z1 - z2) for z1 in pts for z2 in pts)  # diam W(A)
    gaps = sorted(abs(evs[i] - evs[j]) / max(scale, 1e-12)
                  for i in range(3) for j in range(i + 1, 3))
    # disc-ness of W(A): fit support function h(th) ~ Re(e^{-i th} c) + R
    n = len(pts)
    hs, cs, ss = [], [], []
    for k in range(n):
        th = 2 * math.pi * k / n
        e = cmath.exp(-1j * th)
        h = herm_part(mat_scale(e, a))
        hs.append(herm_eigs(h)[2])
        cs.append(math.cos(th))
        ss.append(math.sin(th))
    R = sum(hs) / n
    cr = 2 * sum(h * c for h, c in zip(hs, cs)) / n
    ci = 2 * sum(h * s for h, s in zip(hs, ss)) / n
    resid = math.sqrt(sum((h - (R + cr * c + ci * s)) ** 2
                          for h, c, s in zip(hs, cs, ss)) / n)
    discness = resid / max(R - math.hypot(cr, ci) + 1e-12, 1e-12)
    # |p| peak structure on the boundary
    vals = [abs(peval(p, z)) for z in pts]
    vmax = max(vals)
    peaks = 0
    for k in range(n):
        if vals[k] >= 0.99 * vmax and vals[k] >= vals[k - 1] \
                and vals[k] >= vals[(k + 1) % n]:
            peaks += 1
    # degree of p and root positions relative to W(A) centre/diam
    deg = max((k for k, c in enumerate(p) if abs(c) > 1e-7), default=0)
    return {
        "ratio": r,
        "nonnormal": nonnormal,
        "eig_gaps": gaps,
        "discness": discness,
        "peaks": peaks,
        "deg": deg,
        "diamW": scale,
    }


# ----------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------

FAMILIES = ["random", "jordan3", "jordan2", "nearnormal", "companion",
            "peaks", "real"]


def one_start(family, seed):
    rng = random.Random(seed)
    x0 = start(family, rng)
    f96 = lambda x: objective(x, 96)
    x1, v1 = nelder_mead(f96, x0, step=0.2, maxit=1500)
    x2, v2 = nelder_mead(f96, x1, step=0.05, maxit=1000)
    f192 = lambda x: objective(x, 192)
    x3, v3 = compass(f192, x2, f192(x2), step=0.01, maxit=1500)
    d = diagnostics(x3)
    return {"family": family, "seed": seed, "x": x3,
            "ratio96": v2, "ratio_final": d["ratio"], "diag": d}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--run", type=int, default=0, metavar="N",
                    help="starts per family")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--families", type=str, default=",".join(FAMILIES))
    ap.add_argument("--out", type=str, default="")
    args = ap.parse_args()

    fams = args.families.split(",")
    if args.pilot:
        for fam in fams:
            t0 = time.time()
            rec = one_start(fam, 1000 + hash(fam) % 1000)
            dt = time.time() - t0
            print(f"{fam:11s} ratio={rec['ratio_final']:.6f} "
                  f"({dt:.1f}s) diag={ {k: (round(v,4) if isinstance(v,float) else v) for k,v in rec['diag'].items()} }")
        return

    out = open(args.out, "w") if args.out else sys.stdout
    # interleave families so an early cutoff stays unbiased across them
    for i in range(args.run):
        for fam in fams:
            rec = one_start(fam, args.seed * 100000 + FAMILIES.index(fam) * 1000 + i)
            out.write(json.dumps(rec) + "\n")
            out.flush()
    if args.out:
        out.close()


if __name__ == "__main__":
    main()
