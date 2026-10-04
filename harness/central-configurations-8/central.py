"""Planar unit-mass central configurations; exact local certificates.

Gauge: lambda=1, y_0=0; omit the y_0 force equation, require x_0 != 0.
Pairwise torque cancellation recovers the omitted equation. Summing the full
force equations recovers center of mass zero. Certificates assert one root in
one gauge box, NOT a global configuration count. Verification is stdlib only.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from math import isqrt
from pathlib import Path


class Interval:
    def __init__(self, lo=0, hi=None):
        self.lo = Q(lo)
        self.hi = self.lo if hi is None else Q(hi)
        if self.lo > self.hi:
            raise ValueError('reversed interval')

    @staticmethod
    def cast(x):
        return x if isinstance(x, Interval) else Interval(x)

    def __add__(self, other):
        b = self.cast(other)
        return Interval(self.lo+b.lo, self.hi+b.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.cast(other)

    def __rsub__(self, other):
        return self.cast(other) + -self

    def __mul__(self, other):
        b = self.cast(other)
        v = [self.lo*b.lo, self.lo*b.hi, self.hi*b.lo, self.hi*b.hi]
        return Interval(min(v), max(v))

    __rmul__ = __mul__

    def __truediv__(self, other):
        b = self.cast(other)
        if b.lo <= 0 <= b.hi:
            raise ValueError('division across zero')
        return self * Interval(1/b.hi, 1/b.lo)

    def square(self):
        lower = 0 if self.lo <= 0 <= self.hi else min(self.lo**2, self.hi**2)
        return Interval(lower, max(self.lo**2, self.hi**2))

    def sqrt(self, bits=100):
        if self.lo < 0:
            raise ValueError('negative square root')
        scale = 1 << bits
        a = isqrt((self.lo.numerator*scale*scale)//self.lo.denominator)
        b = isqrt((self.hi.numerator*scale*scale)//self.hi.denominator)
        return Interval(Q(a, scale), Q(b+1, scale))

    def magnitude(self):
        return max(abs(self.lo), abs(self.hi))


def unpack(z):
    if len(z) < 3 or len(z) % 2 != 1:
        raise ValueError('expected 2n-1 coordinates')
    return [[z[0], z[0]*0]] + [[z[i], z[i+1]] for i in range(1, len(z), 2)]


def interval_system(z):
    """Direct inverse-distance residual and its analytic interval Jacobian."""
    p = unpack(z)
    n = len(p)
    f = [[a, b] for a, b in p]
    j = [[Interval(int(a == b)) for b in range(2*n)] for a in range(2*n)]
    for a in range(n):
        for b in range(a):
            d = [p[a][k]-p[b][k] for k in range(2)]
            r2 = d[0].square()+d[1].square()
            if r2.lo <= 0:
                raise ValueError('box may contain collision')
            inv3 = Interval(1)/(r2*r2.sqrt())
            for k in range(2):
                force = d[k]*inv3
                f[a][k] = f[a][k]-force
                f[b][k] = f[b][k]+force
                for ell in range(2):
                    product = d[k].square() if k == ell else d[k]*d[ell]
                    h = inv3*(int(k == ell)-3*product/r2)
                    j[2*a+k][2*a+ell] = j[2*a+k][2*a+ell]-h
                    j[2*b+k][2*b+ell] = j[2*b+k][2*b+ell]-h
                    j[2*a+k][2*b+ell] = j[2*a+k][2*b+ell]+h
                    j[2*b+k][2*a+ell] = j[2*b+k][2*a+ell]+h
    keep = [i for i in range(2*n) if i != 1]
    flat = [x for row in f for x in row]
    return [flat[i] for i in keep], [[j[i][k] for k in keep] for i in keep]


def verify(cert, details=False):
    """Banach self-map and contraction proof in exact rational arithmetic."""
    try:
        if cert['version'] != 1 or cert['normalization'] != 'unit-masses-lambda-1-y0-0':
            return False
        n = cert['n']
        if type(n) is not int or not 2 <= n <= 20:
            return False
        x = [Q(v) for v in cert['center']]
        r = Q(cert['radius'])
        c = [[Q(v) for v in row] for row in cert['inverse']]
        d = 2*n-1
        if len(x) != d or len(c) != d or any(len(row) != d for row in c):
            return False
        if r <= 0 or abs(x[0]) <= r:
            return False
        f0, _ = interval_system([Interval(a) for a in x])
        _, jac = interval_system([Interval(a-r, a+r) for a in x])
        norm = Q(0)
        movement = Q(0)
        for i in range(d):
            rownorm = sum((Interval(int(i == j))-sum(c[i][k]*jac[k][j]
                          for k in range(d))).magnitude() for j in range(d))
            shift = sum(c[i][k]*f0[k] for k in range(d)).magnitude()
            norm = max(norm, rownorm)
            movement = max(movement, shift+rownorm*r)
        valid = norm < 1 and movement < r
        if details:
            return {'valid': valid, 'contraction_bound': str(norm),
                    'self_map_ratio': str(movement/r)}
        return valid
    except (KeyError, TypeError, ValueError, ZeroDivisionError, OverflowError):
        return False


def certify(points, radius='1/100000000'):
    """Untrusted numeric inverse proposal, followed by exact verification."""
    import numpy as np
    p = np.asarray(points, dtype=float)
    if p.ndim != 2 or p.shape[1] != 2 or p[0, 1] != 0:
        raise ValueError('expected points in y0=0 gauge')
    x = [str(float(v)) for v in np.delete(p.ravel(), 1)]
    _, j = interval_system([Interval(v) for v in x])
    inv = np.linalg.inv([[float((v.lo+v.hi)/2) for v in row] for row in j])
    cert = {'version': 1, 'normalization': 'unit-masses-lambda-1-y0-0',
            'n': len(p), 'center': x, 'radius': radius,
            'inverse': [[str(float(v)) for v in row] for row in inv]}
    if not verify(cert):
        raise ValueError('failed to certify root')
    return cert


def distance_intervals(cert):
    """Squared distances enclosing the unique root certified by cert."""
    x, r = [Q(v) for v in cert['center']], Q(cert['radius'])
    p = unpack([Interval(v-r, v+r) for v in x])
    return sorted(((p[i][0]-p[j][0]).square()+(p[i][1]-p[j][1]).square()
                   for i in range(len(p)) for j in range(i)), key=lambda v: v.lo)


def distinct_by_distances(a, b):
    """One-way proof of inequivalence, never a proof of equivalence.

    Sorting lower and upper endpoint lists separately encloses each order
    statistic, even if individual distance intervals overlap.
    """
    da, db = distance_intervals(a), distance_intervals(b)
    if len(da) != len(db):
        return True
    alo, ahi = sorted(v.lo for v in da), sorted(v.hi for v in da)
    blo, bhi = sorted(v.lo for v in db), sorted(v.hi for v in db)
    return any(h < l or hh < ll for ll, h, l, hh in zip(alo, ahi, blo, bhi))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', type=Path)
    args = parser.parse_args()
    result = verify(json.loads(args.certificate.read_text()), details=True)
    print(json.dumps(result))
    return 0 if result and result['valid'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
