#!/usr/bin/env python3
"""Float exploration for an explicit ABK-type witness configuration.

Candidate machinery ONLY (contract: floats locate, never certify).  Searches
the two-parameter family

    unit charges at e1, e2, e3;  charges q at (1/3,1/3,1/3) +- t*(1,1,1)

which is the arXiv:2607.27197 configuration mapped onto the rational
equilateral embedding (their frame has circumradius 1 and axial height eps;
ours has circumradius sqrt(6)/3, so eps = 3t/sqrt(2), and charges are
scale-invariant).  Their charge law: q_eps = (3/4)eps^3 - (5/32)eps^5.

For a grid of (t, lambda) with q = lambda * q_eps(3t/sqrt(2)), runs a global
Newton census of equilibria (dense seeds near the centroid cluster + coarse
global seeds in the localization ball), dedups converged zeros, and reports
count, minimum pairwise separation, minimum distance to a charge, worst
Jacobian conditioning, and the float index sum (should be 1-n = -4 if the
census is complete and everything is nondegenerate).

Usage:
  abk_witness_scan.py --scan                 coarse (t, lambda) table
  abk_witness_scan.py --probe T LAMBDA       heavy census at one point
  abk_witness_scan.py --probe-q T Q          heavy census at explicit q
"""

import argparse
import math
import random
import sys

C = (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)


def make_charges(t: float, q: float):
    return [
        (1.0, (1.0, 0.0, 0.0)),
        (1.0, (0.0, 1.0, 0.0)),
        (1.0, (0.0, 0.0, 1.0)),
        (q, (C[0] + t, C[1] + t, C[2] + t)),
        (q, (C[0] - t, C[1] - t, C[2] - t)),
    ]


def q_law(t: float) -> float:
    eps = 3.0 * t / math.sqrt(2.0)
    return 0.75 * eps ** 3 - (5.0 / 32.0) * eps ** 5


def field(ch, x):
    out = [0.0, 0.0, 0.0]
    for q, p in ch:
        d = [x[k] - p[k] for k in range(3)]
        r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
        if r2 < 1e-24:
            return None
        r3 = r2 ** 1.5
        for k in range(3):
            out[k] += q * d[k] / r3
    return out


def jac(ch, x):
    J = [[0.0] * 3 for _ in range(3)]
    for q, p in ch:
        d = [x[k] - p[k] for k in range(3)]
        r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
        if r2 < 1e-24:
            return None
        r3 = r2 ** 1.5
        r5 = r3 * r2
        for k in range(3):
            for l in range(3):
                J[k][l] += q * ((1.0 / r3 if k == l else 0.0) - 3.0 * d[k] * d[l] / r5)
    return J


def solve3(J, b):
    # Gaussian elimination with partial pivoting; returns None if singular.
    A = [row[:] + [b[i]] for i, row in enumerate(J)]
    for col in range(3):
        piv = max(range(col, 3), key=lambda r: abs(A[r][col]))
        if abs(A[piv][col]) < 1e-300:
            return None
        A[col], A[piv] = A[piv], A[col]
        for r in range(col + 1, 3):
            m = A[r][col] / A[col][col]
            for cc in range(col, 4):
                A[r][cc] -= m * A[col][cc]
    x = [0.0] * 3
    for r in range(2, -1, -1):
        s = A[r][3] - sum(A[r][cc] * x[cc] for cc in range(r + 1, 3))
        x[r] = s / A[r][r]
    return x


def det3(J):
    a, b, c = J[0]
    d, e, f = J[1]
    g, h, i = J[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def newton(ch, x0, iters=80, tol=1e-13):
    x = list(x0)
    for _ in range(iters):
        F = field(ch, x)
        if F is None:
            return None
        J = jac(ch, x)
        if J is None:
            return None
        step = solve3(J, F)
        if step is None:
            return None
        # damping for large steps
        s = max(abs(v) for v in step)
        lam = 1.0 if s < 0.1 else 0.1 / s
        for k in range(3):
            x[k] -= lam * step[k]
        if s < tol:
            F = field(ch, x)
            if F is not None and max(abs(v) for v in F) < 1e-10:
                return tuple(x)
            return None
    return None


def census(t, q, n_cluster=6000, n_global=1500, seed=20260804):
    ch = make_charges(t, q)
    rng = random.Random(seed)
    zeros = []

    def add(z):
        if z is None:
            return
        if sum(v * v for v in z) > 1.0:  # localization ball |x| < R0 = 1
            return
        for q2, p in ch:
            if sum((z[k] - p[k]) ** 2 for k in range(3)) < 1e-12:
                return
        for w in zeros:
            if sum((z[k] - w[k]) ** 2 for k in range(3)) < 1e-16:
                return
        zeros.append(z)

    # dense seeds around the centroid cluster (radius ~ 6t)
    rad = max(6.0 * t, 0.05)
    for _ in range(n_cluster):
        u = [rng.gauss(0, 1) for _ in range(3)]
        nrm = math.sqrt(sum(v * v for v in u))
        r = rad * rng.random() ** (1.0 / 3.0)
        add(newton(ch, [C[k] + r * u[k] / nrm for k in range(3)]))

    # coarse global seeds in the ball
    for _ in range(n_global):
        u = [rng.gauss(0, 1) for _ in range(3)]
        nrm = math.sqrt(sum(v * v for v in u))
        r = rng.random() ** (1.0 / 3.0)
        add(newton(ch, [r * u[k] / nrm for k in range(3)]))

    # diagnostics
    minsep = min((math.dist(a, b) for i, a in enumerate(zeros)
                  for b in zeros[i + 1:]), default=float("inf"))
    mindch = min((math.dist(z, p) for z in zeros for _, p in ch), default=float("inf"))
    dets = [det3(jac(ch, z)) for z in zeros]
    idxsum = sum(1 if d > 0 else -1 for d in dets)
    mindet = min((abs(d) for d in dets), default=0.0)
    return {
        "count": len(zeros), "minsep": minsep, "min_dist_charge": mindch,
        "index_sum": idxsum, "min_abs_det": mindet, "zeros": sorted(zeros),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", action="store_true")
    ap.add_argument("--probe", nargs=2, type=float, metavar=("T", "LAMBDA"))
    ap.add_argument("--probe-q", nargs=2, type=float, metavar=("T", "Q"))
    ap.add_argument("--ts", default="0.25,0.2,0.15,0.125,0.1,0.075")
    ap.add_argument("--lams", default="0.6,0.8,0.9,1.0,1.1,1.2,1.4")
    args = ap.parse_args()

    if args.scan:
        ts = [float(v) for v in args.ts.split(",")]
        lams = [float(v) for v in args.lams.split(",")]
        print(f"{'t':>7} {'lam':>5} {'q':>12} {'count':>5} {'minsep':>9} "
              f"{'mind(chg)':>9} {'idxsum':>6} {'min|det|':>10}")
        for t in ts:
            for lam in lams:
                q = lam * q_law(t)
                r = census(t, q, n_cluster=2500, n_global=600)
                print(f"{t:7.4f} {lam:5.2f} {q:12.6e} {r['count']:5d} "
                      f"{r['minsep']:9.5f} {r['min_dist_charge']:9.5f} "
                      f"{r['index_sum']:6d} {r['min_abs_det']:10.3e}")
            print()
        return 0

    if args.probe or args.probe_q:
        if args.probe:
            t, lam = args.probe
            q = lam * q_law(t)
        else:
            t, q = args.probe_q
        r = census(t, q, n_cluster=20000, n_global=5000)
        print(f"t={t} q={q!r} count={r['count']} minsep={r['minsep']:.6f} "
              f"min_dist_charge={r['min_dist_charge']:.6f} "
              f"index_sum={r['index_sum']} min|det|={r['min_abs_det']:.3e}")
        for z in r["zeros"]:
            d = math.dist(z, C)
            print(f"  ({z[0]:+.9f}, {z[1]:+.9f}, {z[2]:+.9f})  |z-c|={d:.6f}")
        return 0

    ap.error("need --scan, --probe or --probe-q")


if __name__ == "__main__":
    sys.exit(main())
