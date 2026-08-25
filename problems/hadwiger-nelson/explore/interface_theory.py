"""Why the patch/rotated-patch interface is bounded, and an integer re-derivation.

Write lattice points as p = a1 + b1*zeta, q = a2 + b2*zeta, and let

    N1 = a1^2 + a1 b1 + b1^2,          N2 = a2^2 + a2 b2 + b2^2,
    2D = 2 a1 a2 + a1 b2 + a2 b1 + 2 b1 b2      (twice the dot product),
    K  = a2 b1 - a1 b2                          (twice the cross product / sqrt3).

Then, expanding |p - rho_n q|^2 with cos = (2n-1)/(2n) and sin = sqrt(4n-1)/(2n),

    |p - rho_n q|^2 = N1 + N2 - ((2n-1)/n) D - (sqrt(3(4n-1)) / (2n)) K.

When 3(4n-1) is not a perfect square -- equivalently 4n-1 != 3m^2, the
"non-commensurable" case -- the irrational part must vanish on its own, so

    |p - rho_n q| = 1   <=>   K = 0  and  2n(N1 + N2) - (2n-1)(2D) = 2n,

a pair of conditions in integers only.  This module implements that test and
compares its edge count against the count the exact-field geometry produces.
Two independently written routes to the same number: one does plane geometry in
a quadratic field, the other never leaves Z.

The bound.  K = 0 says p and q are parallel, so q = mu p for a rational mu, and
the second condition becomes N1 (1 + mu^2 - 2 c mu) = 1 with c = (2n-1)/(2n).
Since 1 + mu^2 - 2 c mu >= 1 - c^2 = (4n-1)/(4n^2),

    N1 <= 4n^2 / (4n - 1)  <  n + 1.

So every cross edge has both endpoints at squared distance at most n from the
origin, whatever the patch is.  The interface is bounded by a function of n
alone and cannot grow with the patch.

Run:  python problems/hadwiger-nelson/explore/interface_theory.py
"""

from __future__ import annotations

import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))


def _harness_dir():
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

from lattice import (  # noqa: E402
    eisenstein_patch,
    norms_up_to,
    rotated_union,
    rotation_field,
    rotation_for_norm,
)
from unit_distance import build_graph, eisenstein_norm  # noqa: E402

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))


def commensurable(n):
    v = 4 * n - 1
    if v % 3:
        return False
    q = v // 3
    r = int(q ** 0.5)
    return any(c >= 0 and c * c == q for c in (r - 1, r, r + 1))


def cross_edges_integer(n, patch):
    """Cross-copy unit edges, counted with integer arithmetic only.

    Only valid in the non-commensurable case, where the irrational part of the
    squared distance is forced to vanish.
    """
    if commensurable(n):
        raise ValueError("integer criterion does not apply to commensurable n")
    count = 0
    witnesses = []
    for a1, b1 in patch:
        N1 = eisenstein_norm(a1, b1)
        for a2, b2 in patch:
            K = a2 * b1 - a1 * b2
            if K != 0:
                continue
            N2 = eisenstein_norm(a2, b2)
            D2 = 2 * a1 * a2 + a1 * b2 + a2 * b1 + 2 * b1 * b2
            if 2 * n * (N1 + N2) - (2 * n - 1) * D2 == 2 * n:
                # (a1,b1) is in the base copy, (a2,b2) is rotated.  The origin
                # is shared between the copies, so an edge whose rotated
                # endpoint is the origin is not a cross edge.
                if (a2, b2) == (0, 0):
                    continue
                count += 1
                witnesses.append(((a1, b1), (a2, b2), N1, N2))
    return count, witnesses


def bound_on_N1(n):
    return Fraction(4 * n * n, 4 * n - 1)


def main():
    patches = [3, 7, 13, 21, 28, 48, 79]
    ns = [n for n in norms_up_to(49) if not commensurable(n)]
    rows = []
    hdr = (f"{'n':>4} {'patch':>6} {'|patch|':>8} {'integer':>8} {'geometry':>9} "
           f"{'agree':>6} {'N1 bound':>9} {'max N1':>7} {'r(n)':>5}")
    print(hdr)
    print("-" * len(hdr))
    for n in ns:
        field = rotation_field(n)
        rot = rotation_for_norm(n, field)
        for m in patches:
            patch = eisenstein_patch(m)
            pred, wit = cross_edges_integer(n, patch)
            points, tags = rotated_union(field, patch, rot, copies=2)
            edges, _ = build_graph(points)
            meas = sum(1 for a, b in edges if tags[a][0] != tags[b][0])
            maxN1 = max((w[2] for w in wit), default=0)
            rn = sum(1 for ab in patch if eisenstein_norm(*ab) == n)
            row = {
                "n": n, "patch_norm": m, "patch_size": len(patch),
                "integer_criterion": pred, "exact_geometry": meas,
                "agree": pred == meas,
                "N1_bound": str(bound_on_N1(n)),
                "max_N1_seen": maxN1,
                "lattice_points_of_norm_n_in_patch": rn,
            }
            rows.append(row)
            print(f"{n:>4} {m:>6} {len(patch):>8} {pred:>8} {meas:>9} "
                  f"{str(pred == meas):>6} {str(bound_on_N1(n)):>9} {maxN1:>7} {rn:>5}")
    bad = [r for r in rows if not r["agree"]]
    over = [r for r in rows if r["max_N1_seen"] > bound_on_N1(r["n"])]
    print(f"\ndisagreements: {len(bad)}   N1-bound violations: {len(over)}")
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "interface_theory.json")
    with open(path, "w") as fh:
        json.dump(rows, fh, indent=2)
    print(f"wrote {path}")
    return 1 if (bad or over) else 0


if __name__ == "__main__":
    raise SystemExit(main())
