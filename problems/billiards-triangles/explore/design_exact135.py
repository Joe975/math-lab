#!/usr/bin/env python3
"""Attempt 006: exact certificate of a stable periodic orbit at a triangle
whose obtuse angle is EXACTLY 135 degrees (gamma = 3 pi / 4).

The apexes seeing AB = [(0,0),(1,0)] at angle 135 deg form the open arc of
the circle (2x-1)^2 + (2y+1)^2 = 2 with y > 0 (inscribed-angle theorem;
proven here directly from the dot product, no trigonometry: for C = (x,y),
cos(gamma) = -1/sqrt(2)  <=>  CA.CB < 0 and 2 (CA.CB)^2 = |CA|^2 |CB|^2).
The circle has the rational point (x,y) = (1,0) (vertex B), hence a full
rational parametrization by the slope t of chords through it:

    s = 2(1-t)/(1+t^2),  x = 1 - s/2,  y = t s / 2.

So rational apexes are DENSE on the 135-arc.  For a word alive on that arc
we pick t rational near the float argmax, and verify EXACTLY (Fractions):

  (E1) the apex lies on the arc: y > 0 and CA.CB < 0 and
       2 (CA.CB)^2 == |CA|^2 |CB|^2   -- so gamma = 135 deg exactly;
  (E2) the word's corridor at that apex has positive width
       (skeptic_orbit.corridor, 002's independent exact implementation);
  (E3) the periodic orbit is re-derived by exact rational billiard
       simulation (skeptic_orbit.verify_orbit -- no unfolding).

Positive corridor width is an open condition, so the orbit is stable under
perturbation of the triangle; and a rational apex on this arc generically
has alpha, beta irrational multiples of pi (not checked, not claimed).

Usage:
  python3 .../design_exact135.py --word <digits> [--denom-bits B] [--out F]

Stdlib only, deterministic.
"""

import argparse
import json
import math
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from design_neighbours import reduced_width  # noqa: E402
from skeptic_family import width as sk_width  # noqa: E402
import skeptic_orbit  # noqa: E402


def apex_from_t(t):
    """Rational apex on the gamma = 135 deg arc, from rational chord slope
    t in (0, 1) (t -> 0+ degenerates to vertex A, t -> 1- to vertex B)."""
    s = 2 * (1 - t) / (1 + t * t)
    return 1 - s / 2, t * s / 2


def on_arc_exact(x, y):
    """(E1): gamma = 135 deg exactly at apex (x, y). Pure Fraction check."""
    if not y > 0:
        return False
    dot = (-x) * (1 - x) + y * y
    if not dot < 0:
        return False
    ca2 = x * x + y * y
    cb2 = (1 - x) ** 2 + y * y
    return 2 * dot * dot == ca2 * cb2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--word", required=True)
    ap.add_argument("--denom-bits", type=int, default=40,
                    help="t is snapped to denominator 2^bits")
    ap.add_argument("--out")
    args = ap.parse_args()
    w = tuple(int(c) for c in args.word)

    # float hunt for the widest apex ON the arc, parametrized by t
    doubled = (len(w) % 2 == 0 and (len(w) // 2) % 2 == 1
               and w == w[:len(w) // 2] * 2)
    def fwidth(x, y):
        return reduced_width(w, x, y)[0] if doubled else sk_width(w, x, y)
    best_t, best_v = None, -math.inf
    N = 4000
    for k in range(1, N):
        t = k / N
        x, y = apex_from_t(t)
        v = fwidth(x, y)
        if v > best_v:
            best_v, best_t = v, t
    # golden-section-ish refinement by trisection
    lo, hi = best_t - 1.0 / N, best_t + 1.0 / N
    for _ in range(60):
        t1, t2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
        v1 = fwidth(*apex_from_t(t1))
        v2 = fwidth(*apex_from_t(t2))
        if v1 < v2:
            lo = t1
        else:
            hi = t2
    best_t = 0.5 * (lo + hi)
    if not best_v > 0:
        raise SystemExit("word has no positive float corridor on the arc")

    tq = F(round(best_t * 2 ** args.denom_bits), 2 ** args.denom_bits)
    xq, yq = apex_from_t(tq)
    assert on_arc_exact(xq, yq), "exact on-arc check failed"
    c = skeptic_orbit.corridor(w, (xq, yq))
    if c is None or not c[1] > c[0]:
        raise SystemExit("exact corridor not positive at the snapped apex")
    ok, msg = skeptic_orbit.verify_orbit(w, (xq, yq))
    res = {
        "word": args.word, "len": len(w),
        "t": str(tq),
        "apex": [str(xq), str(yq)],
        "gamma_exactly_135_deg": True,
        "on_arc_identity": "2*(CA.CB)^2 == |CA|^2*|CB|^2 and CA.CB < 0, exact",
        "exact_corridor_width_float": float(c[1] - c[0]),
        "exact_corridor_width": str(c[1] - c[0]),
        "orbit_verified_by_independent_simulation": bool(ok),
        "orbit_msg": str(msg),
        "conclusion": (f"stable periodic orbit (word length {len(w)}) "
                       f"certified at a rational-apex triangle with obtuse "
                       f"angle EXACTLY 135 deg" if ok else "INCOMPLETE"),
    }
    if args.out:
        json.dump(res, open(args.out, "w"), indent=1)
    print(json.dumps(res, indent=1))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
