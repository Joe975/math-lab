#!/usr/bin/env python3
"""How the closed forms in tp_curvature were FOUND: exact rational-function
reconstruction from sampled values.

beta, gamma and c are rational functions of (x, y, r) = (k2/k1, k3/k1,
sqrt(m2)) -- every step of their construction is field arithmetic over Q -- but
tp_curvature only computes them pointwise.  This module recovers the symbolic
forms from those points, which is what turned "c is exactly rational at every
rational instance" into the printed formulas.

Method, one variable at a time:

 1. Sample beta(r) at enough rational r for a fixed sigma and solve the
    homogeneous linear system  P(r_i) - beta_i Q(r_i) = 0  over Q for the
    smallest total degree with a one-dimensional null space; validate the
    candidate on every sampled point, not just the ones used to fit.
 2. Factor the denominator by testing the conjectured roots 1, -p, -p*q with
    p = (1+x)/(y-x) and q = (y-1)/(x-1) -- read off from three triples and then
    confirmed by exact deflation everywhere.
 3. Read the remaining shape numbers (K, n1) off the deflated fit, and fit
    THOSE against k2 and k3 the same way.

Nothing here is load-bearing: the resulting formulas are checked as identities
against the jet expansion by `tp_curvature.py --closed-form` (2610 exact
instances at the time of writing), so a reconstruction bug would show up there
rather than silently propagating.

Usage:
    python tp_ratrecon.py --derive
    python tp_ratrecon.py --selftest
Standard library only.
"""

from __future__ import annotations

import argparse
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tp_curvature as TC  # noqa: E402

RS = [Fr(a, b) for a, b in
      [(1, 9), (1, 7), (1, 6), (1, 5), (2, 7), (1, 3), (3, 8), (2, 5), (3, 7),
       (1, 2), (4, 7), (3, 5), (5, 8), (2, 3), (5, 7), (3, 4), (4, 5), (5, 6),
       (6, 7), (7, 8), (8, 9), (9, 10)]]


def nullspace(M):
    """Exact null-space basis of M (rows of Fractions)."""
    A = [row[:] for row in M]
    rows, cols = len(A), len(A[0]) if A else 0
    piv, r = [], 0
    for c in range(cols):
        pr = next((i for i in range(r, rows) if A[i][c] != 0), None)
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == rows:
            break
    out = []
    for f in (c for c in range(cols) if c not in piv):
        v = [Fr(0)] * cols
        v[f] = Fr(1)
        for i, c in enumerate(piv):
            v[c] = -A[i][f]
        out.append(v)
    return out


def ev(co, x):
    return sum(c * x ** j for j, c in enumerate(co))


def fit_auto(xs, ys, maxtot=10):
    """Smallest-total-degree rational fit validated on ALL sampled points."""
    for tot in range(maxtot + 1):
        for p in range(tot + 1):
            q = tot - p
            if p + q + 2 > len(xs):
                continue
            M = [[x ** j for j in range(p + 1)] + [-y * x ** j for j in range(q + 1)]
                 for x, y in zip(xs[:p + q + 2], ys[:p + q + 2])]
            ns = nullspace(M)
            if len(ns) != 1:
                continue
            P, Q = ns[0][:p + 1], ns[0][p + 1:]
            if not any(Q):
                continue
            if all(ev(P, x) == y * ev(Q, x) for x, y in zip(xs, ys)):
                return P, Q
    return None


def deflate(co, root):
    """Divide an ascending-coefficient polynomial by (t - root)."""
    n = len(co) - 1
    out = [Fr(0)] * n
    acc = Fr(0)
    for i in range(n, 0, -1):
        acc = co[i] if i == n else co[i] + acc * root
        out[i - 1] = acc
    return out, co[0] + acc * root


def sample(sig, key):
    xs, ys = [], []
    for r in RS:
        xs.append(r)
        ys.append(TC.expand(sig, r, deg=2)[key])
    return xs, ys


def shape(trip, key):
    """Fit key(r) at a fixed triple and deflate the conjectured roots."""
    k1, k2, k3 = trip
    sig = tuple(Fr(x) for x in trip)
    P, Q = fit_auto(*sample(sig, key))
    lead = Q[-1]
    P, Q = [x / lead for x in P], [x / lead for x in Q]
    p = Fr(k1 + k2) / (k3 - k2)
    q = Fr(k3 - k1) / (k2 - k1)
    roots = [Fr(1), -p, -p] + ([-p * q] if key == "beta" else [])
    rest = Q
    for root in roots:
        rest, rem = deflate(rest, root)
        assert rem == 0, (trip, key, root, "conjectured root failed")
    return P, rest, p, q


def derive():
    trips = [(1, 2, 3), (1, 2, 5), (1, 2, 7), (1, 3, 7), (1, 4, 9), (2, 5, 11)]
    print("beta = K (r + n1) / [(1-r)(r+p)^2 (r+pq)]\n")
    print(f"{'sigma':>12} {'p':>8} {'q':>8} {'K':>18} {'n1':>12} "
          f"{'K pred':>18} {'n1 pred':>12}")
    for trip in trips:
        k1, k2, k3 = trip
        P, rest, p, q = shape(trip, "beta")
        assert rest == [Fr(1)], rest
        K, n1 = -P[1], P[0] / P[1]
        x, y = Fr(k2, k1), Fr(k3, k1)
        Kp = k1 * p ** 3 * (y - 1) ** 2 * (y + 1) / 4
        n1p = (1 + x) * (y - 1) / ((x - 1) * (y + 1))
        print(f"{str(trip):>12} {str(p):>8} {str(q):>8} {str(K):>18} "
              f"{str(n1):>12} {str(Kp):>18} {str(n1p):>12}")
        assert (K, n1) == (Kp, n1p), trip

    print("\ngamma = K r (A - r) / [(1-r)(r+p)^2 q2(r)],  q2 = C2 + B2r r - r^2\n")
    print(f"{'sigma':>12} {'A':>10} {'A pred':>10} {'C2':>16} {'B2r':>16}")
    for trip in trips:
        k1, k2, k3 = trip
        P, rest, p, q = shape(trip, "gamma")
        assert P[0] == 0 and len(P) == 3, P
        A = -P[1] / P[2]
        assert -P[2] == k1 * p ** 3 * (Fr(k3, k1) - 1) ** 2 * (Fr(k3, k1) + 1) / 4, trip
        C2, B2r = -rest[0], -rest[1]
        cf = TC.closed_forms(tuple(Fr(t) for t in trip), Fr(1, 2))
        Ap = p * (Fr(k3, k1) - 1)
        print(f"{str(trip):>12} {str(A):>10} {str(Ap):>10} {str(C2):>16} "
              f"{str(B2r):>16}")
        assert A == Ap, (trip, A, Ap)
        assert cf["q2"] == C2 + B2r * Fr(1, 2) - Fr(1, 4), trip
    print("\nAll deflations and shape identifications hold exactly; the assembled "
          "formulas are the ones in tp_curvature.closed_forms.")


def selftest():
    derive()
    print("selftest: the reconstruction reproduces tp_curvature's closed forms")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--derive", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.derive:
        derive()
    if a.selftest:
        selftest()
