"""Exact rational contraction checker for equal unit-mass central configurations.

The lifted variables are (x0,x1,y1,...), followed by positive inverse
distances s_ij.  We set y0=0 and lambda=1.  The equations are
q_i-sum_j (q_i-q_j)s_ij**3=0 (except its y0 component), and
s_ij**2 |q_i-q_j|**2=1.  Positivity selects physical inverse distances.

The omitted equation follows from the torque identity sum q_i cross f_i=0:
all other components vanish, leaving x0*f0y=0.  Summing full force equations
then establishes center of mass zero.  A nonzero x0 makes the rotational
slice transverse; uniqueness is local in this slice, not modulo all global
rotations or permutations.  No enumeration completeness is asserted.

Verification uses only stdlib Fraction, including interval Jacobian bounds.
Floating point is used solely by certify() to propose a preconditioner.
"""
from fractions import Fraction as F
from itertools import combinations


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def neg(a):
    return -a[1], -a[0]


def mul(a, b):
    products = [x*y for x in a for y in b]
    return min(products), max(products)


def scale(a, k):
    return mul(a, (k, k))


def sq(a):
    return (F(0) if a[0] <= 0 <= a[1] else min(a[0]**2, a[1]**2),
            max(a[0]**2, a[1]**2))


ZERO = (F(0), F(0))
ONE = (F(1), F(1))


def system(v, n):
    """Interval residual and sparse interval Jacobian, independently derived."""
    coords = [(i, axis) for i in range(n) for axis in range(2)
              if (i, axis) != (0, 1)]
    index = {key: k for k, key in enumerate(coords)}
    pairs = list(combinations(range(n), 2))
    size = len(coords) + len(pairs)
    out = list(v[:len(coords)]) + [ZERO]*len(pairs)
    jac = [{k: ONE} for k in range(len(coords))] + [{} for _ in pairs]

    def put(row, col, value):
        if col is not None:
            jac[row][col] = add(jac[row].get(col, ZERO), value)

    for p, (i, j) in enumerate(pairs):
        col = len(coords) + p
        s = v[col]
        s2, s3 = sq(s), mul(sq(s), s)
        distance = ZERO
        for axis in range(2):
            a, b = index.get((i, axis)), index.get((j, axis))
            delta = add(v[a] if a is not None else ZERO,
                        neg(v[b] if b is not None else ZERO))
            distance = add(distance, sq(delta))
            for row, sign in ((a, -1), (b, 1)):
                if row is not None:
                    out[row] = add(out[row], scale(mul(delta, s3), sign))
                    put(row, a, scale(s3, sign))
                    put(row, b, scale(s3, -sign))
                    put(row, col, scale(mul(delta, s2), 3*sign))
            put(col, a, scale(mul(s2, delta), 2))
            put(col, b, scale(mul(s2, delta), -2))
        out[col] = add(mul(s2, distance), neg(ONE))
        put(col, col, scale(mul(s, distance), 2))
    assert len(out) == size
    return out, jac


def verify(cert):
    """Return True only after strict inclusion and contraction are proved."""
    try:
        if cert['version'] != 1 or type(cert['n']) is not int:
            return False
        n = cert['n']
        if not 2 <= n <= 32 or len(cert['masses']) != n:
            return False
        if any(F(m) != 1 for m in cert['masses']):
            return False
        d = 2*n-1+n*(n-1)//2
        x = [F(t) for t in cert['center']]
        radius = F(cert['radius'])
        c = [[F(t) for t in row] for row in cert['inverse']]
        if len(x) != d or len(c) != d or any(len(row) != d for row in c):
            return False
        if radius <= 0 or abs(x[0]) <= radius:
            return False
        if any(s <= radius for s in x[2*n-1:]):
            return False
        f, _ = system([(t, t) for t in x], n)
        _, jac = system([(t-radius, t+radius) for t in x], n)
        for i in range(d):
            product = [ZERO]*d
            for k in range(d):
                for j, entry in jac[k].items():
                    product[j] = add(product[j], scale(entry, c[i][k]))
            norm = F(0)
            for j in range(d):
                entry = add((F(i == j), F(i == j)), neg(product[j]))
                norm += max(abs(entry[0]), abs(entry[1]))
            correction = abs(sum(c[i][k]*f[k][0] for k in range(d)))
            if norm >= 1 or correction + radius*norm >= radius:
                return False
        return True
    except (KeyError, TypeError, ValueError, ZeroDivisionError, OverflowError):
        return False


def certify(points, radius='1/100000000'):
    """Propose and check a JSON-serializable certificate; numpy is proposer-only."""
    import math
    import numpy as np
    n = len(points)
    if n < 2 or any(len(p) != 2 for p in points) or points[0][1] != 0:
        raise ValueError('Expected points in gauge y0=0')
    x = [F(str(points[i][a])) for i in range(n) for a in range(2)
         if (i, a) != (0, 1)]
    for i, j in combinations(range(n), 2):
        distance = math.dist(points[i], points[j])
        if distance == 0:
            raise ValueError('Collision')
        x.append(F(str(1/distance)))
    _, jac = system([(t, t) for t in x], n)
    matrix = np.zeros((len(x), len(x)))
    for i, row in enumerate(jac):
        for j, value in row.items():
            matrix[i, j] = float(value[0])
    inverse = np.linalg.inv(matrix)
    cert = dict(version=1, n=n, masses=['1']*n,
                center=list(map(str, x)), radius=str(radius),
                inverse=[[str(F(str(t))) for t in row] for row in inverse])
    if not verify(cert):
        raise ValueError('Proposed enclosure failed exact contraction verification')
    return cert
