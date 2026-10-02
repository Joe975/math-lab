#!/usr/bin/env python3
"""Integer-relation search by LLL (exact integer arithmetic, stdlib only).

Used to test whether a float value x (known to ~eps) satisfies a small
integer polynomial. find_poly(x, deg, digits) returns the shortest LLL
vector of the lattice [I | round(10^digits * x^i)] and its residual.
A hit is only a CANDIDATE: a relation found at the precision limit can be
an artefact. Accept one only if it (a) has small coefficients relative to
10^(digits/(deg+1)) and (b) still holds when the value is recomputed to more
digits by an independent route.
"""
from __future__ import annotations

from fractions import Fraction


def lll(B, delta=Fraction(3, 4)):
    B = [list(map(int, b)) for b in B]
    n = len(B)

    def dot(u, v):
        return sum(a * b for a, b in zip(u, v))

    def gso(B):
        Bs = []; mu = [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            v = [Fraction(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = dot(B[i], Bs[j]) / dot(Bs[j], Bs[j])
                v = [a - mu[i][j] * b for a, b in zip(v, Bs[j])]
            Bs.append(v)
        return Bs, mu

    Bs, mu = gso(B); k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                B[k] = [a - q * b for a, b in zip(B[k], B[j])]
                Bs, mu = gso(B)
        if dot(Bs[k], Bs[k]) >= (delta - mu[k][k - 1] ** 2) * dot(Bs[k - 1], Bs[k - 1]):
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            Bs, mu = gso(B); k = max(k - 1, 1)
    return B


def find_poly(x: float, deg: int, digits: int = 13):
    scale = 10 ** digits
    rows = []
    for i in range(deg + 1):
        r = [0] * (deg + 1); r[i] = 1
        rows.append(r + [round(scale * x ** i)])
    best = min(lll(rows), key=lambda v: sum(c * c for c in v[:-1]))
    coeffs = best[:-1]
    resid = sum(c * x ** i for i, c in enumerate(coeffs))
    return coeffs, resid
