#!/usr/bin/env python3
"""Exact eigenvalue signature of DF at the centroid for the abk family.

Configuration family (attempt 001): unit charges at e1, e2, e3, and two
charges q > 0 at c +- t*(1,1,1), where c = (1/3,1/3,1/3) is the centroid.

Facts derived here (proofs in attempt 002; every load-bearing inequality is
instantiated below in exact rational arithmetic):

  1. F(c) = 0 for every (t, q): the three triangle terms have equal radii
     r_t = |c - e_i| and directions summing to 3c - (e1+e2+e3) = 0; the two
     axial terms are opposite with equal radii r_a = t*sqrt(3).

  2. DF(c) = A*I + B*J with J the all-ones matrix (D3 symmetry + the
     rank-one structure of each charge term q*(I/r^3 - 3 d d^T/r^5)):
        sum_tri d_i d_i^T = I - J/3   (exact, verified below),
        axial d d^T = t^2 J.
     Since trace DF = -laplacian V = 0, B = -A; equivalently the identity
     3/r_t^3 - 3/r_t^5 = -(3/2)/r_t^3 holds exactly because r_t^2 = 2/3
     (verified below). Eigenvalues: lambda_ax = A + 3B = -2A on (1,1,1),
     lambda_perp = A (multiplicity 2).

  3. lambda_ax = 3/r_t^3 - 4q/r_a^3. Both terms positive, so
        sign(lambda_ax) = sign(9*r_a^6 - 16 q^2 r_t^6)
                        = sign(6561 t^6 - 128 q^2)      (exact rational),
     using r_a^6 = 27 t^6 and r_t^6 = 8/27.
     6561 t^6 > 128 q^2: signature (+,-,-), det DF > 0 (nondegenerate saddle)
     6561 t^6 < 128 q^2: signature (-,+,+), det DF < 0
     equality: DEGENERATE, at q*(t) = 81*sqrt(2)*t^3/16.

Usage: skeptic_centroid_hessian.py  [--t 17/200] [--q 4367/1000000 ...]
Prints the exact comparison for each q, and q*(t) exactly and as a decimal.
Standard library only, exact except the labelled float cross-check.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import isqrt


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--t", default="17/200")
    ap.add_argument("--q", nargs="*",
                    default=["4360/1000000", "4367/1000000", "4400/1000000"])
    args = ap.parse_args()
    t = Fraction(args.t)

    # -- fact 1: F(c) = 0, exact pieces -------------------------------------
    c = (Fraction(1, 3),) * 3
    tri = [(Fraction(1), Fraction(0), Fraction(0)),
           (Fraction(0), Fraction(1), Fraction(0)),
           (Fraction(0), Fraction(0), Fraction(1))]
    dirs = [tuple(c[k] - p[k] for k in range(3)) for p in tri]
    assert all(sum(d[k] for d in dirs) == 0 for k in range(3)), "triangle dirs"
    r2s = {sum(d[k] ** 2 for k in range(3)) for d in dirs}
    assert r2s == {Fraction(2, 3)}, "equal triangle radii, r_t^2 = 2/3"
    print("F(c) = 0: triangle directions sum to 0 at equal radii (r_t^2 = 2/3),")
    print("          axial terms opposite at equal radii (r_a^2 = 3 t^2)  [exact]")

    # -- fact 2: sum_tri d d^T = I - J/3, and the traceless identity --------
    S = [[sum(d[a] * d[b] for d in dirs) for b in range(3)] for a in range(3)]
    expect = [[Fraction(2, 3) if a == b else Fraction(-1, 3) for b in range(3)]
              for a in range(3)]
    assert S == expect, "sum_tri d d^T != I - J/3"
    # identity 3/r^3 - 3/r^5 = -(3/2)/r^3  <=>  r^2 = 2/3 (with r = r_t):
    assert Fraction(2, 3) * (3 + Fraction(3, 2)) == 3, "traceless identity"
    print("DF(c) = A*I + B*J with B = -A: rank-one sums and the r_t^2 = 2/3")
    print("          traceless identity verified exactly")

    # -- fact 3: the sign of lambda_ax --------------------------------------
    lhs = 6561 * t ** 6
    print(f"\nt = {t}:   6561 t^6 = {lhs}  ~ {float(lhs):.6e}")
    for qs in args.q:
        q = Fraction(qs)
        rhs = 128 * q ** 2
        if lhs > rhs:
            verdict = "lambda_ax > 0: signature (+,-,-), det DF(c) > 0 (index +1)"
        elif lhs < rhs:
            verdict = "lambda_ax < 0: signature (-,+,+), det DF(c) < 0 (index -1)"
        else:
            verdict = "DEGENERATE (lambda_ax = 0)"
        print(f"q = {qs}: 128 q^2 = {rhs} ~ {float(rhs):.6e}\n"
              f"    -> {verdict}   [exact rational comparison]")

    # -- the degeneracy locus ------------------------------------------------
    print(f"\nq*(t) = 81*sqrt(2)/16 * t^3  (exact algebraic form), t^3 = {t**3}")
    # certified decimal enclosure of q* via integer sqrt at 2^-120
    v = 2 * (Fraction(81, 16) * t ** 3) ** 2      # q*^2 = 2*(81 t^3/16)^2
    s = 1 << 120
    r = isqrt(v.numerator * v.denominator * s * s)
    lo = Fraction(r, v.denominator * s)
    hi = Fraction(r + 1, v.denominator * s)
    assert lo ** 2 <= v <= hi ** 2
    print(f"q*({t}) in [{float(lo):.12e}, {float(hi):.12e}]  (certified)")

    # -- float cross-check (EVIDENCE only, never certifies) ------------------
    tf, r_t3, = float(t), (2.0 / 3.0) ** 1.5
    r_a3 = (3.0 * tf * tf) ** 1.5
    for qs in args.q:
        qf = float(Fraction(qs))
        lam_ax = 3.0 / r_t3 - 4.0 * qf / r_a3
        print(f"float check q={qs}: lambda_ax ~ {lam_ax:+.6e}, "
              f"lambda_perp ~ {-lam_ax / 2:+.6e} (x2)")


if __name__ == "__main__":
    main()
