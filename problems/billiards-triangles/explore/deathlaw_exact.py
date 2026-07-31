#!/usr/bin/env python3
"""Attempt 003: exact-arithmetic brackets for W(a,b) death angles.

Two subcommands, both stdlib-only and deterministic:

  alive --a A --b B --dg DG
      Exhibit a RATIONAL apex whose exact (Fraction) corridor for W(a,b) is
      positive and which realizes an exact periodic orbit (skeptic_orbit's
      independent simulator), together with a CERTIFIED lower bound on its
      obtuse angle gamma.  This proves: death(W(a,b)) >= gamma_lo exactly.
      The apex is found near the float argmax at gamma = D_pred - DG and
      snapped to denominator 2^40; positivity is then exact, not float.

      The gamma bound is certified as follows.  With dot = CA.CB < 0 exact
      and r = dot^2/(|CA|^2 |CB|^2) exact rational, gamma > g iff
      r > cos^2(g), since cos is decreasing and both cosines are negative.
      cos^2(g) is enclosed with unfold.py's interval cosine and rational pi
      enclosure; the comparison uses the enclosure's upper end.  These are
      single (non-iterated) interval evaluations, so plain intervals are the
      right tool here (the affine-form lesson applies to iterated maps).

  deadscan --a A --b B --dg DG [--uniform N] [--geo K]
      Sample rational apexes with CERTIFIED gamma in [D_pred + DG/2,
      D_pred + 2*DG] along the arc, concentrated geometrically at the
      measured last-life window edge alpha = 90/a deg (plus a uniform grid
      and the beta = 90/(b+1) edge), and check every exact corridor is
      empty.  This is sample-bounded EVIDENCE that the word is dead above
      the predicted angle, never a proof.

Usage examples (from the repo root):
  python3 problems/billiards-triangles/explore/deathlaw_exact.py alive \
      --a 10 --b 9 --dg 1e-8 --out .../deathlaw_exact_alive_W109.json
  python3 problems/billiards-triangles/explore/deathlaw_exact.py deadscan \
      --a 10 --b 9 --dg 1e-4 --out .../deathlaw_exact_dead_W109.json
"""

import argparse
import json
import math
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..",
                                "harness", "billiards-triangles"))
from skeptic_family import family_word  # noqa: E402
import skeptic_orbit  # noqa: E402
import unfold  # noqa: E402


def predicted_death_exact(a, b):
    """General law of attempt 003: gamma_d = 180 - 90(a+b)/(a(b+1)).
    Reduces to 001's laws at b = a and b = a-1."""
    return F(180) - F(90 * (a + b), a * (b + 1))


def apex_from_angles(alpha_deg, beta_deg):
    ta = math.tan(math.radians(alpha_deg))
    tb = math.tan(math.radians(beta_deg))
    return tb / (ta + tb), ta * tb / (ta + tb)


def snap(v, bits=40):
    return F(round(v * 2 ** bits), 2 ** bits)


def cos2_encl(g_deg):
    """Interval enclosure of cos^2(g_deg * pi / 180) for rational g_deg."""
    x = unfold.Iv(g_deg) * unfold.Iv(unfold.PI_LO, unfold.PI_HI) / 180
    return unfold.square(unfold.cos_iv(x))


def gamma_bounds_certified(x, y, g_lo, g_hi):
    """(gamma > g_lo certified, gamma < g_hi certified) for the exact
    rational apex (x, y); both g in degrees, rational, in (90, 180)."""
    dot = (-x) * (1 - x) + y * y           # CA.CB with CA=(-x,-y), CB=(1-x,-y)
    if dot >= 0:
        return False, False
    r = dot * dot / ((x * x + y * y) * ((1 - x) ** 2 + y * y))
    return r > cos2_encl(g_lo).hi, r < cos2_encl(g_hi).lo


def cmd_alive(args):
    a, b = args.a, args.b
    w = family_word(a, b)
    d_pred = predicted_death_exact(a, b)
    g_target = float(d_pred) - args.dg
    theta = 180.0 - g_target
    # float search for the widest apex on the target arc, then snap.
    # near death the argmax is near alpha:beta = a:b, window edge alpha=90/a.
    from deathlaw_measure import best_alive
    bw, (al, be) = best_alive(w, a, b, g_target, full=(b < a - 1))
    if not bw > 0:
        raise SystemExit(f"no alive float apex found at gamma={g_target}")
    x, y = apex_from_angles(al, be)
    xq, yq = snap(x, args.bits), snap(y, args.bits)
    c = skeptic_orbit.corridor(w, (xq, yq))
    if c is None or not c[1] > c[0]:
        raise SystemExit("exact corridor not positive after snap; "
                         "reduce --dg or raise --bits")
    ok, msg = skeptic_orbit.verify_orbit(w, (xq, yq))
    # certified gamma lower bound: try successively tighter rational bounds
    g_lo_claim = None
    for k in (2, 3, 4, 6, 8, 10, 12):
        cand = d_pred - F(10) ** -k if F(10) ** -k > F(args.dg / 4) else None
        if cand is None:
            break
        lo_ok, _ = gamma_bounds_certified(xq, yq, cand, F(180))
        if lo_ok:
            g_lo_claim = cand
        else:
            break
    _, hi_ok = gamma_bounds_certified(xq, yq, F(90), d_pred)
    res = {
        "a": a, "b": b, "len": len(w),
        "apex": [str(xq), str(yq)],
        "exact_corridor_width": str(c[1] - c[0]),
        "exact_corridor_width_float": float(c[1] - c[0]),
        "orbit_verified_by_simulation": bool(ok),
        "orbit_msg": msg,
        "predicted_death": str(d_pred),
        "certified_gamma_lower_bound_deg": str(g_lo_claim),
        "certified_gamma_below_prediction": bool(hi_ok),
        "conclusion": (f"death(W({a},{b})) >= {g_lo_claim} deg, exactly; "
                       f"apex gamma certified < predicted death: {hi_ok}"),
    }
    if args.out:
        json.dump(res, open(args.out, "w"), indent=1)
    print(json.dumps(res, indent=1))
    if not ok:
        sys.exit(1)


def cmd_deadscan(args):
    a, b = args.a, args.b
    w = family_word(a, b)
    d_pred = predicted_death_exact(a, b)
    band_lo = d_pred + F(args.dg) / 2
    band_hi = d_pred + 2 * F(args.dg)
    g_mid = float(d_pred) + args.dg
    theta = 180.0 - g_mid
    edge_a, edge_b = 90.0 / a, 90.0 / (b + 1)
    lo_cut = 0.0 if b < a - 1 else theta / 2
    alphas = [lo_cut + (theta - lo_cut) * (k + 0.5) / args.uniform
              for k in range(args.uniform)]
    for k in range(args.geo):
        d = 0.5 * 2.0 ** -k
        if lo_cut <= edge_a - d < theta:
            alphas.append(edge_a - d)
        if lo_cut <= theta - (edge_b - d) < theta:
            alphas.append(theta - (edge_b - d))
    n_in_band = n_alive = n_skipped = 0
    worst = None
    for al in alphas:
        x, y = apex_from_angles(al, theta - al)
        xq, yq = snap(x, args.bits), snap(y, args.bits)
        lo_ok, hi_ok = gamma_bounds_certified(xq, yq, band_lo, band_hi)
        if not (lo_ok and hi_ok):
            n_skipped += 1
            continue
        n_in_band += 1
        c = skeptic_orbit.corridor(w, (xq, yq))
        wd = c[1] - c[0] if c else None
        if wd is not None and wd > 0:
            n_alive += 1
            print(f"ALIVE apex above prediction: {xq},{yq} width={wd}")
        if wd is not None and (worst is None or wd > worst):
            worst = wd
    res = {
        "a": a, "b": b, "len": len(w),
        "band_deg": [str(band_lo), str(band_hi)],
        "n_apexes_certified_in_band": n_in_band,
        "n_skipped_gamma_uncertified": n_skipped,
        "n_alive": n_alive,
        "max_exact_width_float": float(worst) if worst is not None else None,
        "conclusion": ("REFUTES law upward" if n_alive else
                       f"all {n_in_band} sampled apexes with certified gamma "
                       f"in band are exactly dead (sample-bounded evidence)"),
    }
    if args.out:
        json.dump(res, open(args.out, "w"), indent=1)
    print(json.dumps(res, indent=1))
    sys.exit(1 if n_alive else 0)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("alive")
    c.add_argument("--a", type=int, required=True)
    c.add_argument("--b", type=int, required=True)
    c.add_argument("--dg", type=float, default=1e-8)
    c.add_argument("--bits", type=int, default=40)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_alive)
    c = sub.add_parser("deadscan")
    c.add_argument("--a", type=int, required=True)
    c.add_argument("--b", type=int, required=True)
    c.add_argument("--dg", type=float, default=1e-4)
    c.add_argument("--uniform", type=int, default=120)
    c.add_argument("--geo", type=int, default=40)
    c.add_argument("--bits", type=int, default=40)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_deadscan)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
