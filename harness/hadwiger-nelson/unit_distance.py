"""Unit-distance graphs over an exact real field, with a two-sided certifier.

A unit-distance graph is a finite point set in the plane with an edge between
two points exactly when their distance is 1.  The only thing that can go wrong
computationally is a distance that is nearly 1 being treated as 1, so every
decision here is an exact equality test in a multiquadratic field
(``exact_field``): squared distance == 1 or it is not.

The certifier checks both directions, because only checking the edges would let
a graph declare fewer edges than the point set actually has -- which can turn a
non-3-colourable point set into a 3-colourable graph.
"""

from __future__ import annotations

from exact_field import MQ, MQField


class Point:
    """A plane point with both coordinates in one :class:`MQField`."""

    __slots__ = ("x", "y")

    def __init__(self, x: MQ, y: MQ):
        if x.field != y.field:
            raise ValueError("coordinates live in different fields")
        self.x = x
        self.y = y

    @property
    def field(self):
        return self.x.field

    def __eq__(self, other):
        return isinstance(other, Point) and self.x == other.x and self.y == other.y

    def __hash__(self):
        return hash((self.x, self.y))

    def __repr__(self):
        return f"Point({self.x}, {self.y})"

    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y)

    def to_floats(self):
        return (self.x.to_float(), self.y.to_float())


def sq_dist(p: Point, q: Point) -> MQ:
    dx = p.x - q.x
    dy = p.y - q.y
    return dx * dx + dy * dy


def distinct_points(points):
    """True iff all points are pairwise distinct as exact field elements."""
    return len(set(points)) == len(points)


def build_graph(points):
    """Exhaustive exact construction: every pair tested, both directions decided.

    Returns (edges, stats) where edges is a sorted list of (i, j), i < j.
    Every pair of indices is tested; there is no numeric prefilter, so the
    non-edges are certified by construction.
    """
    n = len(points)
    one = points[0].field.one() if n else None
    edges = []
    tested = 0
    for i in range(n):
        pi = points[i]
        for j in range(i + 1, n):
            tested += 1
            if sq_dist(pi, points[j]) == one:
                edges.append((i, j))
    return edges, {"vertices": n, "pairs_tested": tested, "edges": len(edges)}


def certify(points, edges):
    """Independently re-check a claimed unit-distance graph.

    Recomputes every pair from the coordinates and compares against the claimed
    edge set.  Returns a report dict; ``report['ok']`` is the verdict.
    """
    n = len(points)
    claimed = set()
    for a, b in edges:
        if a == b or not (0 <= a < n and 0 <= b < n):
            return {"ok": False, "reason": f"bad edge index pair {(a, b)}"}
        claimed.add((min(a, b), max(a, b)))
    if not distinct_points(points):
        return {"ok": False, "reason": "points are not pairwise distinct"}
    one = points[0].field.one()
    bad_edges, bad_nonedges = [], []
    for i in range(n):
        for j in range(i + 1, n):
            is_unit = sq_dist(points[i], points[j]) == one
            if is_unit and (i, j) not in claimed:
                bad_nonedges.append((i, j))
            if (not is_unit) and (i, j) in claimed:
                bad_edges.append((i, j))
    return {
        "ok": not bad_edges and not bad_nonedges,
        "vertices": n,
        "edges_claimed": len(claimed),
        "edges_wrong": bad_edges,
        "unit_pairs_missing": bad_nonedges,
        "pairs_checked": n * (n - 1) // 2,
    }


def numeric_separation(points, edges, prec=120):
    """Report how far the non-edges are from being edges, at high precision.

    Reporting only -- the exact test above is what decides.  This exists so a
    reader can see that the exact verdict is not a knife-edge: a genuine
    non-edge should have |d^2 - 1| enormously larger than the working epsilon,
    and any pair where it is not deserves a second look.
    """
    from decimal import Decimal, getcontext

    getcontext().prec = prec
    claimed = {(min(a, b), max(a, b)) for a, b in edges}
    worst_edge = None
    worst_nonedge = None
    n = len(points)
    for i in range(n):
        for j in range(i + 1, n):
            d2 = sq_dist(points[i], points[j]).to_decimal(prec)
            gap = abs(d2 - Decimal(1))
            if (i, j) in claimed:
                if worst_edge is None or gap > worst_edge[0]:
                    worst_edge = (gap, (i, j))
            else:
                if worst_nonedge is None or gap < worst_nonedge[0]:
                    worst_nonedge = (gap, (i, j))
    return {
        "prec": prec,
        "max_edge_residual": worst_edge,
        "min_nonedge_gap": worst_nonedge,
    }


def check_colouring(n, edges, colouring):
    """Linear check that ``colouring`` is proper.  Written to be trivially auditable."""
    if len(colouring) != n:
        return False
    for a, b in edges:
        if colouring[a] == colouring[b]:
            return False
    return True


def adjacency(n, edges):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    return adj


def eisenstein_point(field: MQField, a: int, b: int) -> Point:
    """The triangular-lattice point a + b*zeta, zeta = exp(i*pi/3).

    x = a + b/2, y = b*sqrt3/2.  Requires sqrt3 in ``field``.
    """
    from fractions import Fraction as _F

    s3 = field.sqrt_int(3)
    x = field.rational(_F(2 * a + b, 2))
    y = s3 * field.rational(_F(b, 2))
    return Point(x, y)


def eisenstein_norm(a: int, b: int) -> int:
    """|a + b*zeta|^2 for zeta = exp(i*pi/3): the integer a^2 + a*b + b^2."""
    return a * a + a * b + b * b
