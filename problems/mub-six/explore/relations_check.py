#!/usr/bin/env python3
"""Exact-form checks on the 50-digit polished floor (see data/hiprec_d6k4.txt).

Runs integer-relation searches (LLL, exact integers) at 45 digits and
evaluates every candidate at the full 50 digits. Stdlib only.
"""
from decimal import Decimal as D, getcontext
from intrel import lll

getcontext().prec = 70
p1 = D("0.12439794566570321461037700050123282717085012795674")
p2 = D("0.15140197747747829881064201194373798309990083907013")
p3 = D("0.18105001921420462164474524688875729743231225824328")
L = D("0.051249218996283827782882716219040904801850493547422")
c4 = (1 + 7 * L / 3) / 12
s2 = 1 - c4.sqrt()          # RLE's sin^2(theta_opt), via their closed form


def relation(xs, digits=45):
    """Shortest integer vector c with sum c_i xs_i ~ 0."""
    scale = D(10) ** digits
    rows = []
    for i, x in enumerate(xs):
        r = [0] * len(xs); r[i] = 1
        rows.append(r + [int((scale * x).to_integral_value())])
    best = min(lll(rows), key=lambda v: sum(c * c for c in v[:-1]))
    c = best[:-1]
    return c, sum(D(ci) * x for ci, x in zip(c, xs))


def ev(co, x):
    return sum(D(c) * x ** i for i, c in enumerate(co))


print("cubic residuals at 50 digits:")
for name, co, x in [("p1", [-12, 153, -552, 784], p1), ("p2", [-1, -15, 24, 784], p2),
                    ("p3", [-3, 45, -228, 392], p3), ("sin2", [-22, 111, -192, 112], s2)]:
    print(f"  {name}: {ev(co, x):.2e}")

print("minimal polynomial search for L (deg <= 6, 45 digits):")
for deg in (2, 3, 4, 6):
    c, r = relation([L ** i for i in range(deg + 1)])
    print(f"  deg {deg}: {c}  residual {r:.2e}")

print("express p1, p2, p3, L in the basis 1, s, s^2 (s = sin^2 theta_opt):")
for name, x in [("p1", p1), ("p2", p2), ("p3", p3), ("L", L)]:
    c, r = relation([x, D(1), s2, s2 * s2])
    print(f"  {name}: {c}  residual {r:.2e}")

# Raynal-Lu-Englert Eq. 21 radical form [T]: r = (21 sqrt3 - 36)^(1/3),
# sin^2 theta_opt = (3 + 16 r - r^2) / (28 r). Compared against s2 above,
# which came from our 50-digit polish (via their ASD formula) -- an
# independent published route to the same number.
r = (21 * D(3).sqrt() - 36) ** (D(1) / 3)
s_rle = (3 + 16 * r - r * r) / (28 * r)
print(f"RLE radical sin^2 = {s_rle:.50f}")
print(f"ours            = {s2:.50f}")
print(f"difference      = {s_rle - s2:.2e}")
c = 1 - s_rle
print(f"overlaps from RLE radical: p1 = {4*c*c/3:.40f}")
print(f"                           p2 = {(2*c-1)**2:.40f}")
print(f"                           p3 = {c*(3-4*c)/3:.40f}")
print(f"  vs polished p1 - 4c^2/3 = {p1 - 4*c*c/3:.2e}, p2 = {p2 - (2*c-1)**2:.2e}, p3 = {p3 - c*(3-4*c)/3:.2e}")
Lq = 448*c**4 - 768*c**3 + 504*c**2 - 144*c + 15
print(f"L(c) quartic at RLE radical c: {Lq:.50f}; vs polished L: {Lq - L:.2e}")
