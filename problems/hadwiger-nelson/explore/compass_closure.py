"""What can a unit compass alone build, starting from a unit segment?

Operation: given two existing points at distance strictly between 0 and 2, add
the two points at distance exactly 1 from both (intersection of the two unit
circles about them).  Iterate to a fixpoint inside a bounded disc.

The script checks three things exactly:

  1. the operation never needs a square root outside Q(sqrt3) (no escape);
  2. every point produced is an Eisenstein lattice point;
  3. the resulting unit-distance graph is 3-chromatic, with an explicit
     colouring checked by the independent linear checker.

It also checks the closure lemma directly: for every ordered pair of lattice
points at squared distance < 4 inside a patch, both intersection points are
again lattice points.

Run:  python explore/compass_closure.py [radius_squared_bound]
"""

from __future__ import annotations

import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))


def _harness_dir():
    """Locate harness/hadwiger-nelson by walking up from this file."""
    d = HERE
    while True:
        cand = os.path.join(d, "harness", "hadwiger-nelson")
        if os.path.isdir(cand):
            return cand
        parent = os.path.dirname(d)
        if parent == d:
            raise RuntimeError("could not locate harness/hadwiger-nelson")
        d = parent


sys.path.insert(0, _harness_dir())
sys.path.insert(0, HERE)

import colouring as C  # noqa: E402
from exact_field import field_for  # noqa: E402
from lattice import (  # noqa: E402
    RadicandEscape,
    eisenstein_patch,
    unit_circle_intersections,
)
from unit_distance import (  # noqa: E402
    Point,
    build_graph,
    certify,
    eisenstein_norm,
    eisenstein_point,
    sq_dist,
)

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))


def as_eisenstein(p: Point):
    """Return (a, b) if p = a + b*zeta with integers a, b, else None."""
    xr, xs = p.x.coords()          # basis (1, sqrt3)
    yr, ys = p.y.coords()
    if xs != 0 or yr != 0:
        return None
    b2 = ys * 2                    # y = (b/2) sqrt3  =>  b = 2*ys
    if b2.denominator != 1:
        return None
    b = int(b2)
    a2 = xr * 2 - b                # x = (2a+b)/2
    if a2.denominator != 1 or int(a2) % 2 != 0:
        return None
    return int(a2) // 2, b


def closure(field, bound_sq, max_rounds=12):
    """Fixpoint of the unit-compass operation inside |p|^2 <= bound_sq."""
    pts = [
        Point(field.rational(0), field.rational(0)),
        Point(field.rational(1), field.rational(0)),
    ]
    seen = set(pts)
    radicands = set()
    escapes = []
    rounds = 0
    while rounds < max_rounds:
        rounds += 1
        new = []
        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                d2 = sq_dist(pts[i], pts[j])
                if not d2.is_rational():
                    raise AssertionError("non-rational squared distance in closure")
                r = d2.as_fraction()
                if r == 0 or r >= 4:
                    continue
                radicands.add(r)
                try:
                    cands = unit_circle_intersections(pts[i], pts[j])
                except RadicandEscape as exc:
                    escapes.append((r, exc.radicand))
                    continue
                for q in cands:
                    if q in seen:
                        continue
                    if (q.x * q.x + q.y * q.y).as_fraction() > bound_sq:
                        continue
                    seen.add(q)
                    new.append(q)
        if not new:
            break
        pts.extend(new)
    return pts, sorted(radicands), escapes, rounds


def lemma_check(field, patch_norm):
    """Directly verify closure of the lattice under the operation."""
    patch = eisenstein_patch(patch_norm)
    pts = {ab: eisenstein_point(field, *ab) for ab in patch}
    distances = set()
    bad = []
    checked = 0
    for i, ab in enumerate(patch):
        for cd in patch[i + 1:]:
            d2 = sq_dist(pts[ab], pts[cd]).as_fraction()
            distances.add(d2)
            if d2 == 0 or d2 >= 4:
                continue
            checked += 1
            for q in unit_circle_intersections(pts[ab], pts[cd]):
                if as_eisenstein(q) is None:
                    bad.append((ab, cd, repr(q)))
    return {
        "patch_norm": patch_norm,
        "patch_size": len(patch),
        "pairs_with_intersections": checked,
        "squared_distances_below_4": sorted(d for d in distances if 0 < d < 4),
        "non_lattice_results": bad,
    }


def main():
    bound_sq = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    field = field_for(3)
    print(f"field                 : {field}")

    pts, radicands, escapes, rounds = closure(field, bound_sq)
    print(f"closure within |p|^2<={bound_sq}: {len(pts)} points after {rounds} rounds")
    print(f"squared distances used: {[str(r) for r in radicands]}")
    print(f"field escapes needed  : {len(escapes)}")

    tags = [as_eisenstein(p) for p in pts]
    non_lattice = [repr(p) for p, t in zip(pts, tags) if t is None]
    print(f"points off the lattice: {len(non_lattice)}")
    assert not non_lattice, non_lattice[:3]

    edges, stats = build_graph(pts)
    rep = certify(pts, edges)
    print(f"unit-distance graph   : {stats}, certificate ok={rep['ok']}")
    assert rep["ok"]

    chi = C.chromatic_number(len(pts), edges)
    explicit = [(a - b) % 3 for (a, b) in tags]
    ok_explicit = C.verify_colouring(len(pts), edges, explicit)
    print(f"chi(closure graph)    : {chi}")
    print(f"explicit (a-b) mod 3  : proper = {ok_explicit}")
    adj = C.adjacency(len(pts), edges)
    has_triangle = any(adj[u] & adj[v] for u, v in edges)
    print(f"triangle present      : {has_triangle} (forces chi >= 3)")

    lem = lemma_check(field, 21)
    print(f"lemma patch norm<=21  : {lem['patch_size']} points, "
          f"{lem['pairs_with_intersections']} intersecting pairs, "
          f"distances^2 < 4 seen = {[str(d) for d in lem['squared_distances_below_4']]}, "
          f"non-lattice results = {len(lem['non_lattice_results'])}")
    assert not lem["non_lattice_results"]

    # Which norms below 4 exist at all?
    norms_below_4 = sorted({eisenstein_norm(a, b)
                            for a in range(-8, 9) for b in range(-8, 9)
                            if 0 < eisenstein_norm(a, b) < 4})
    print(f"lattice norms in (0,4): {norms_below_4}")

    os.makedirs(OUT, exist_ok=True)
    rec = {
        "bound_sq": bound_sq,
        "closure_points": len(pts),
        "rounds_to_fixpoint": rounds,
        "squared_distances_used": [str(r) for r in radicands],
        "field_escapes": escapes,
        "all_points_are_lattice_points": True,
        "unit_distance_graph": stats,
        "chi": chi,
        "explicit_colouring_rule": "(a - b) mod 3 on the lattice coordinates",
        "explicit_colouring_proper": ok_explicit,
        "lemma": lem,
        "lattice_norms_strictly_between_0_and_4": norms_below_4,
    }
    path = os.path.join(OUT, "compass_closure.json")
    with open(path, "w") as fh:
        json.dump(rec, fh, indent=2, default=str)
    print(f"wrote                 : {path}")


if __name__ == "__main__":
    main()
