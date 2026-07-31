#!/usr/bin/env python3
"""Attempt 003: high-precision float measurement of W(a,b) death angles.

W(a,b) = (0 (12)^a (02)^b)^2, the family of attempt 001.  The death angle is
the supremum of the obtuse angle gamma over apexes at which the word's
corridor is non-empty.  Attempts 001/002 measured it by bisection over a
fixed 400-point x-scan per constant-gamma arc, which under-estimates the
death when the alive x-window is narrower than the scan spacing (their
measured deaths sat 0.001-0.006 deg below the conjectured law).

Measured mechanism this tool exploits (and re-verifies): near death the alive
window in the angle alpha at vertex A is an interval whose upper edge is
alpha = 90/a deg, pinching onto it (a != b) or staying wide (a == b).  The
alive test therefore samples, on each arc:
  * a uniform alpha-grid (fallback, catches any other window), and
  * geometric grids accumulating at alpha = 90/a and at beta = 90/(b+1),
    the two conjectured window edges, so a window of width 2^-k deg is found
    with ~k samples instead of 2^k.
A positive verdict is a genuinely positive float corridor width at a concrete
apex (exact-checkable downstream); a DEAD verdict is bounded by this sampling
design.  Float only - nothing here is a certificate.  Exact statements are
made by deathlaw_exact.py and deathlaw_symbolic.py.

Corridor evaluation is skeptic_family.corridor (002's independent
implementation), reused deliberately so this is not census.py's code path.

Usage:
  deathlaw_measure.py death --a A --b B [--out FILE]   bisect the death angle
  deathlaw_measure.py profile --a A --b B --gamma G    max width on one arc

Deterministic.  Runtime: W(12,12) (word length 98) death takes ~2 min.
"""

import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from skeptic_family import width, family_word  # noqa: E402


def apex_from_angles(alpha_deg, beta_deg):
    """Apex (x, y) of the triangle with angle alpha at A=(0,0), beta at
    B=(1,0).  Requires 0 < alpha, beta and alpha + beta < 180."""
    ta = math.tan(math.radians(alpha_deg))
    tb = math.tan(math.radians(beta_deg))
    return tb / (ta + tb), ta * tb / (ta + tb)


def best_alive(w, a, b, g, uniform=800, geo=64, full=False):
    """(best width, best (alpha, beta)) over the sampled apexes of the arc
    gamma = g.  By default restricted to alpha >= beta (x <= 1/2, the half
    of parameter space 001/002 measured; the mirror half holds the
    0<->1-swapped word's congruent tile).  full=True scans the whole arc:
    the b <= a-2 members (e.g. W(5,3), W(4,2)) have their alive window on
    the mirror half, pinching onto the corner alpha = 90/a from below.
    Alive == best width > ALIVE_TOL, which sits far above float noise
    (~1e-15 for these word lengths) and far below true widths at the gamma
    resolution this tool reports (width ~ 0.1 * dg near death)."""
    theta = 180.0 - g                       # alpha + beta on this arc
    edge_a = 90.0 / a                       # conjectured upper edge in alpha
    edge_b = 90.0 / (b + 1)                 # conjectured upper edge in beta
    lo_cut = 0.0 if full else theta / 2
    alphas = []
    for k in range(uniform):
        alphas.append(lo_cut + (theta - lo_cut) * (k + 0.5) / uniform)
    for k in range(geo):
        d = 0.5 * 2.0 ** -k
        if lo_cut <= edge_a - d < theta:
            alphas.append(edge_a - d)
        if lo_cut <= theta - (edge_b - d) < theta:
            alphas.append(theta - (edge_b - d))
    best, arg = -math.inf, None
    for al in alphas:
        be = theta - al
        if al <= 0 or be <= 0:
            continue
        x, y = apex_from_angles(al, be)
        v = width(w, x, y)
        if v > best:
            best, arg = v, (al, be)
    return best, arg


ALIVE_TOL = 1e-12


def predicted_death(a, b):
    """The general law of attempt 003: gamma_d = 180 - 90(a+b)/(a(b+1)).
    Reduces to 001's two laws at b = a and b = a-1."""
    return 180.0 - 90.0 * (a + b) / (a * (b + 1))


def cmd_death(args):
    w = family_word(args.a, args.b)
    pred = predicted_death(args.a, args.b)
    g_lo = g_hi = None
    if pred is not None:
        wl, _ = best_alive(w, args.a, args.b, pred - 0.5, args.uniform,
                           full=args.full)
        wh, _ = best_alive(w, args.a, args.b, pred + 0.5, args.uniform,
                           full=args.full)
        if wl > ALIVE_TOL and wh <= ALIVE_TOL:
            g_lo, g_hi = pred - 0.5, pred + 0.5
        else:
            print(f"# prediction bracket failed (w({pred - 0.5})={wl:.3g}, "
                  f"w({pred + 0.5})={wh:.3g}); falling back to coarse walk")
    if g_lo is None:
        alive = []
        g = 90.25
        while g < 179.5:
            if best_alive(w, args.a, args.b, g, args.uniform,
                          full=args.full)[0] > ALIVE_TOL:
                alive.append(g)
            g += 0.5
        if not alive:
            print(json.dumps({"a": args.a, "b": args.b, "alive": None}))
            return
        g_lo, g_hi = alive[-1], alive[-1] + 0.5
    last_arg = None
    for _ in range(args.iters):
        mid = 0.5 * (g_lo + g_hi)
        mw, ar = best_alive(w, args.a, args.b, mid, args.uniform,
                            full=args.full)
        if mw > ALIVE_TOL:
            g_lo, last_arg = mid, ar
        else:
            g_hi = mid
    mw, ar = best_alive(w, args.a, args.b, g_lo, args.uniform,
                        full=args.full)
    res = {
        "a": args.a, "b": args.b, "len": len(w),
        "death_between": [g_lo, g_hi],
        "uniform_points": args.uniform, "geo_points": 64,
        "full_arc": args.full,
        "bisect_iters": args.iters,
        "alive_at_death_lo": {"alpha_deg": ar[0], "beta_deg": ar[1],
                              "alpha_beta_ratio": ar[0] / ar[1],
                              "width": mw},
        "predicted_death": predicted_death(args.a, args.b),
    }
    if res["predicted_death"] is not None:
        res["measured_minus_predicted"] = (
            0.5 * (g_lo + g_hi) - res["predicted_death"])
    if args.out:
        json.dump(res, open(args.out, "w"), indent=1)
    print(json.dumps(res))


def cmd_profile(args):
    w = family_word(args.a, args.b)
    mw, ar = best_alive(w, args.a, args.b, args.gamma, args.uniform,
                        full=args.full)
    x, y = apex_from_angles(*ar)
    print(json.dumps({
        "a": args.a, "b": args.b, "gamma": args.gamma,
        "best_width": mw, "alpha_deg": ar[0], "beta_deg": ar[1],
        "x": x, "y": y}))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("death")
    c.add_argument("--a", type=int, required=True)
    c.add_argument("--b", type=int, required=True)
    c.add_argument("--uniform", type=int, default=800)
    c.add_argument("--iters", type=int, default=44)
    c.add_argument("--full", action="store_true",
                   help="scan the whole arc, not just alpha >= beta")
    c.add_argument("--out")
    c.set_defaults(fn=cmd_death)
    c = sub.add_parser("profile")
    c.add_argument("--a", type=int, required=True)
    c.add_argument("--b", type=int, required=True)
    c.add_argument("--gamma", type=float, required=True)
    c.add_argument("--uniform", type=int, default=800)
    c.add_argument("--full", action="store_true")
    c.set_defaults(fn=cmd_profile)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
