#!/usr/bin/env python3
"""Exact (rational-coefficient) algebra for the overlap pattern of the floor.

Hypothesis from relations_check.py: at the four-basis floor every defective
pair has rows that are permutations of
    p1 = 4c^2/3,  p2 = (2c-1)^2,  p3 = c(3-4c)/3  (x4),   c = cos^2(theta).
This script works with polynomials in c over Q (Fractions), exactly:
  * checks p1 + p2 + 4 p3 = 1 identically (rows of a unistochastic matrix);
  * forms f(c) = 6[(p1-1/6)^2 + (p2-1/6)^2 + 4(p3-1/6)^2]  (one pair's defect)
    and L(c) = 3 f(c);
  * compares f'(c) with the cubic 112c^3 - 144c^2 + 63c - 9 found numerically;
  * evaluates L at the cubic's relevant root by bisection in Fractions.
Stdlib only.
"""
from fractions import Fraction as F


def padd(a, b):
    n = max(len(a), len(b)); a = a + [0] * (n - len(a)); b = b + [0] * (n - len(b))
    return [x + y for x, y in zip(a, b)]


def pmul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def pscale(a, k): return [x * k for x in a]
def pder(a): return [i * a[i] for i in range(1, len(a))]
def pev(a, x): return sum(c * x ** i for i, c in enumerate(a))
def trim(a):
    while len(a) > 1 and a[-1] == 0: a = a[:-1]
    return a


p1 = [F(0), F(0), F(4, 3)]
p2 = [F(1), F(-4), F(4)]
p3 = [F(0), F(1), F(-4, 3)]
print("p1 + p2 + 4 p3 =", trim(padd(padd(p1, p2), pscale(p3, 4))))

sixth = [F(-1, 6)]
def sq(p): q = padd(p, sixth); return pmul(q, q)
f = pscale(padd(padd(sq(p1), sq(p2)), pscale(sq(p3), 4)), 6)
L = pscale(f, 3)
print("L(c) =", " + ".join(f"({c})c^{i}" for i, c in enumerate(trim(L))))
dL = trim(pder(L))
cubic = [F(-9), F(63), F(-144), F(112)]
ratio = dL[-1] / cubic[-1]
print("L'(c) =", dL, "; ratio to cubic:", ratio,
      "; proportional:", all(a == ratio * b for a, b in zip(dL, cubic)))

# root of the cubic in [0.25, 0.35] (c = 1 - 0.69455...)
lo, hi = F(1, 4), F(7, 20)
assert pev(cubic, lo) * pev(cubic, hi) < 0
for _ in range(200):
    mid = (lo + hi) / 2
    if pev(cubic, lo) * pev(cubic, mid) <= 0: hi = mid
    else: lo = mid
c0 = (lo + hi) / 2
print(f"c* = {float(c0):.17g}  (s* = 1 - c* = {float(1 - c0):.17g})")
print(f"L(c*) = {float(pev(L, c0)):.17g}")
print(f"L''(c*) = {float(pev(pder(dL), c0)):.6g}  (> 0: minimum along the family)")
# all real roots of L' in [0,1]
pts = [F(i, 1000) for i in range(1001)]
sgn = [pev(dL, x) for x in pts]
print("sign changes of L' on [0,1] at c ~", [float(pts[i]) for i in range(1000) if sgn[i] * sgn[i + 1] < 0])
for cc in (F(0), F(1, 4), F(3, 4), F(1)):
    print(f"  L({cc}) = {float(pev(L, cc)):.6g}   overlaps {[float(pev(p, cc)) for p in (p1, p2, p3)]}")
