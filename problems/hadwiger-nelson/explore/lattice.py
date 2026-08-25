"""Triangular-lattice patches, exact rotations, and unit-circle intersections.

Everything here is generic construction machinery for unit-distance graphs over
an exact field: it enumerates and rotates point sets and computes compass
intersections, and decides nothing by floating point.

Conventions
-----------
The triangular lattice is the ring of Eisenstein integers E = Z + Z*zeta with
zeta = exp(i*pi/3) = 1/2 + i*sqrt3/2.  The lattice point (a, b) is a + b*zeta,
with Cartesian coordinates ((2a+b)/2, b*sqrt3/2) and squared length the integer
norm N(a, b) = a^2 + a*b + b^2.  Nearest neighbours are the six units, so the
unit-distance graph on E is the usual triangular grid.

For a positive integer n, ``rotation_for_norm(n)`` is the rotation about the
origin through the angle theta with cos(theta) = 1 - 1/(2n).  It is exactly the
angle that moves a point at distance sqrt(n) from the origin to a point at
distance 1 from where it started, since |p - R p|^2 = 2n(1 - cos theta) = 1.
Its matrix entries lie in Q(sqrt(4n-1)), so the union of a lattice patch with
its image lives in Q(sqrt3, sqrt(4n-1)).
"""

from __future__ import annotations

from fractions import Fraction

from exact_field import MQField, field_for, squarefree_part
from unit_distance import Point, eisenstein_norm, eisenstein_point


def eisenstein_patch(max_norm):
    """All (a, b) with N(a, b) <= max_norm, in a deterministic order."""
    out = []
    bound = int(2 * (max_norm ** 0.5)) + 3
    for a in range(-bound, bound + 1):
        for b in range(-bound, bound + 1):
            if eisenstein_norm(a, b) <= max_norm:
                out.append((a, b))
    out.sort(key=lambda ab: (eisenstein_norm(*ab), ab))
    return out


def norms_up_to(limit):
    """Integers that are norms of Eisenstein integers, 1 <= n <= limit."""
    seen = set()
    bound = int(limit ** 0.5) + 2
    for a in range(-2 * bound, 2 * bound + 1):
        for b in range(-2 * bound, 2 * bound + 1):
            n = eisenstein_norm(a, b)
            if 1 <= n <= limit:
                seen.add(n)
    return sorted(seen)


class Rotation:
    """Exact plane rotation with entries in a multiquadratic field."""

    __slots__ = ("field", "c", "s", "label")

    def __init__(self, field, c, s, label=""):
        self.field = field
        self.c = c
        self.s = s
        self.label = label

    def __call__(self, p: Point) -> Point:
        return Point(self.c * p.x - self.s * p.y, self.s * p.x + self.c * p.y)

    def is_orthogonal(self):
        """c^2 + s^2 == 1, checked exactly."""
        return self.c * self.c + self.s * self.s == self.field.one()

    def __repr__(self):
        return f"Rotation({self.label or ''}, cos={self.c}, sin={self.s})"


def rotation_for_norm(n, field=None):
    """The rotation with cos(theta) = 1 - 1/(2n); see the module docstring."""
    n = int(n)
    if n < 1:
        raise ValueError("n must be a positive integer")
    if field is None:
        field = field_for(3, 4 * n - 1)
    c = field.rational(Fraction(2 * n - 1, 2 * n))
    s = field.sqrt_int(4 * n - 1) * field.rational(Fraction(1, 2 * n))
    return Rotation(field, c, s, label=f"n={n}")


def rotation_field(n):
    """The field Q(sqrt3, sqrt(4n-1)) in coprime-generator form."""
    return field_for(3, 4 * n - 1)


def patch_points(field, patch):
    return [eisenstein_point(field, a, b) for a, b in patch]


def rotated_union(field, patch, rot, copies=2):
    """Points of ``patch``, together with its images under rot, rot^2, ...

    Returns (points, tags) where tags[i] = (copy_index, (a, b)).  Duplicate
    points across copies are removed by exact equality, so the vertex set is a
    genuine set of distinct plane points.
    """
    base = patch_points(field, patch)
    points, tags, seen = [], [], {}
    cur = base
    for k in range(copies):
        if k > 0:
            cur = [rot(p) for p in cur]
        for p, ab in zip(cur, patch):
            if p in seen:
                continue
            seen[p] = len(points)
            points.append(p)
            tags.append((k, ab))
    return points, tags


class RadicandEscape(Exception):
    """Raised when a compass construction needs a square root outside the field."""

    def __init__(self, radicand):
        super().__init__(f"sqrt({radicand}) is not in the working field")
        self.radicand = radicand


def unit_circle_intersections(p: Point, q: Point):
    """The two points at distance 1 from both p and q, exactly.

    Returns [] when the circles do not meet in two real points (squared
    distance 0 or >= 4).  Raises :class:`RadicandEscape` when the intersection
    is real but its coordinates need a square root the field does not contain --
    which is the interesting case, so the radicand is carried on the exception.

    Implemented for rational squared distance, which is what the lattice
    supplies; a non-rational squared distance raises NotImplementedError rather
    than guessing.
    """
    field = p.field
    d = Point(q.x - p.x, q.y - p.y)
    d2 = d.x * d.x + d.y * d.y
    if d2 == field.zero():
        return []
    if not d2.is_rational():
        raise NotImplementedError("non-rational squared distance")
    r = d2.as_fraction()
    if r >= 4:
        return []
    # t = sqrt((4 - r) / (4r)); the intersections are (p+q)/2 +- t * perp(q - p).
    frac = Fraction(4 - r, 4 * r)
    radicand = frac.numerator * frac.denominator
    try:
        root = field.sqrt_int(radicand)
    except ValueError:
        raise RadicandEscape(squarefree_part(radicand)[1])
    t = root * field.rational(Fraction(1, frac.denominator))
    mx = (p.x + q.x) * field.rational(Fraction(1, 2))
    my = (p.y + q.y) * field.rational(Fraction(1, 2))
    px, py = -d.y, d.x  # perp(q - p)
    return [
        Point(mx + t * px, my + t * py),
        Point(mx - t * px, my - t * py),
    ]
