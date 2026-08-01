#!/usr/bin/env python3
"""Attempt 006: exact alive certificate for an arbitrary word at a target
gamma, with a certified rational gamma bracket.

Same certification pattern as deathlaw_exact.py (003), generalized from the
W(a,b) family to any even bounce word and any target arc:

  1. float search for the widest apex on the arc gamma = --gamma, using
     design_neighbours' whole-arc sampler PLUS deep geometric accumulation
     at the caller-supplied corner --corner-alpha (the windows this attempt
     hunts pinch onto corners like alpha = beta = 22.5);
  2. snap the apex to denominator 2^--bits; corridor positivity is then
     checked in EXACT Fraction arithmetic (skeptic_orbit.corridor, 002's
     independent implementation);
  3. the periodic orbit is re-derived by exact rational simulation
     (skeptic_orbit.verify_orbit -- independent of any unfolding);
  4. gamma is certified strictly inside the rational bracket
     (--glo, --ghi) via exact rational cos^2 against interval enclosures
     (deathlaw_exact.gamma_bounds_certified; single non-iterated interval
     evaluations, tier-0 unfold.py interval cosine + rational pi).

The result "word alive at a triangle with gamma in (glo, ghi), exactly" is
what this prints and writes.

Usage example (the pinch-gap certificate of attempt 006):
  python3 problems/billiards-triangles/explore/design_certify.py \
      --word 012121212020202012121212020202 \
      --gamma 135.02 --corner-alpha 22.5 \
      --glo 135000001/1000000 --ghi 135048/1000 \
      --out problems/billiards-triangles/data/design_cert_W43_gap.json

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
from design_neighbours import (apex_from_angles, reduced_width,  # noqa: E402
                               arc_alphas)
from deathlaw_exact import gamma_bounds_certified, snap  # noqa: E402
import skeptic_orbit  # noqa: E402


def best_apex(word, g, corner_alpha, uniform=400, depth=48):
    th = 180.0 - g
    als = arc_alphas(g, uniform=uniform, jmax=12, depth=14)
    if corner_alpha is not None:
        for d in range(1, depth + 1):
            eps = 0.5 * 2.0 ** -d
            for al in (corner_alpha - eps, corner_alpha + eps,
                       th - corner_alpha - eps, th - corner_alpha + eps):
                if 0 < al < th:
                    als.append(al)
    best, arg = -math.inf, None
    doubled = (len(word) % 2 == 0 and (len(word) // 2) % 2 == 1
               and word == word[:len(word) // 2] * 2)
    for al in als:
        be = th - al
        x, y = apex_from_angles(al, be)
        if doubled:
            v, _ = reduced_width(word, x, y)
        else:
            from skeptic_family import width as sk_w
            v = sk_w(word, x, y)
        if v > best:
            best, arg = v, (al, be)
    return best, arg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--word", required=True)
    ap.add_argument("--gamma", type=float, required=True)
    ap.add_argument("--corner-alpha", type=float, default=None)
    ap.add_argument("--glo", required=True,
                    help="rational lower gamma bound to certify (deg)")
    ap.add_argument("--ghi", required=True,
                    help="rational upper gamma bound to certify (deg)")
    ap.add_argument("--bits", type=int, default=48)
    ap.add_argument("--out")
    args = ap.parse_args()
    w = tuple(int(c) for c in args.word)
    glo, ghi = F(args.glo), F(args.ghi)

    bw, arg = best_apex(w, args.gamma, args.corner_alpha)
    if not bw > 0:
        raise SystemExit(f"no alive float apex on arc gamma={args.gamma}")
    x, y = apex_from_angles(*arg)
    xq, yq = snap(x, args.bits), snap(y, args.bits)
    c = skeptic_orbit.corridor(w, (xq, yq))
    if c is None or not c[1] > c[0]:
        raise SystemExit("exact corridor not positive after snap; "
                         "raise --bits or move --gamma")
    ok, msg = skeptic_orbit.verify_orbit(w, (xq, yq))
    lo_ok, hi_ok = gamma_bounds_certified(xq, yq, glo, ghi)
    res = {
        "word": args.word, "len": len(w),
        "apex": [str(xq), str(yq)],
        "float_alpha_beta_deg": [arg[0], arg[1]],
        "exact_corridor_width": str(c[1] - c[0]),
        "exact_corridor_width_float": float(c[1] - c[0]),
        "orbit_verified_by_independent_simulation": bool(ok),
        "orbit_msg": str(msg),
        "certified_gamma_bracket_deg": [str(glo), str(ghi)],
        "gamma_gt_glo_certified": bool(lo_ok),
        "gamma_lt_ghi_certified": bool(hi_ok),
        "conclusion": (f"word alive (exact positive corridor + independent "
                       f"orbit simulation) at a triangle with gamma in "
                       f"({glo}, {ghi}) deg, certified"
                       if (ok and lo_ok and hi_ok) else "INCOMPLETE"),
    }
    if args.out:
        json.dump(res, open(args.out, "w"), indent=1)
    print(json.dumps(res, indent=1))
    sys.exit(0 if (ok and lo_ok and hi_ok) else 1)


if __name__ == "__main__":
    main()
