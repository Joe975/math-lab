#!/usr/bin/env python3
"""SKEPTIC float-side re-measurements: family death angles and sliver probes.

Own implementation throughout (all-three-vertex reflection unfolding, as in
skeptic_orbit.py but in floats; census.py uses the foot-parameter shortcut).
Float here measures frontiers only, exactly as the record's own float layer
does; nothing here is a certificate.

Subcommands:
  death  --a A --b B      bisect the death angle of W(a,b) = (0(12)^a(02)^b)^2
  sliver --gamma G --words FILE --xpoints N
                          dense x-scan of every word in FILE along the
                          constant-gamma arc: max corridor width found
"""

import argparse
import json
import math

ENDS = ((1, 2), (2, 0), (0, 1))


def corridor(word, ax, ay):
    """(lo, hi) corridor projection interval, or None if tau degenerate."""
    tri = ((0.0, 0.0), (1.0, 0.0), (ax, ay))
    gates = []
    for s in word:
        i, j = ENDS[s]
        u, v = tri[i], tri[j]
        gates.append((u, v))
        dx, dy = v[0] - u[0], v[1] - u[1]
        n2 = dx * dx + dy * dy
        new = []
        for p in tri:
            t = ((p[0] - u[0]) * dx + (p[1] - u[1]) * dy) / n2
            new.append((2 * (u[0] + t * dx) - p[0],
                        2 * (u[1] + t * dy) - p[1]))
        tri = tuple(new)
    taux, tauy = tri[0][0], tri[0][1]
    if taux * taux + tauy * tauy < 1e-16:
        return None
    nx, ny = -tauy, taux
    lo, hi = -math.inf, math.inf
    for u, v in gates:
        a = u[0] * nx + u[1] * ny
        b = v[0] * nx + v[1] * ny
        if a > b:
            a, b = b, a
        if a > lo:
            lo = a
        if b < hi:
            hi = b
    return lo, hi


def width(word, ax, ay):
    c = corridor(word, ax, ay)
    return -math.inf if c is None else c[1] - c[0]


def gamma_deg(x, y):
    ax, ay = -x, -y
    bx, by = 1.0 - x, -y
    return math.degrees(math.atan2(abs(ax * by - ay * bx),
                                   ax * bx + ay * by))


def y_for_gamma(x, g):
    lo, hi = 1e-12, 0.5
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if gamma_deg(x, mid) > g:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def family_word(a, b):
    return ((0,) + (1, 2) * a + (0, 2) * b) * 2


def alive_max(word, g, K, x_lo=0.02):
    best = -math.inf
    for k in range(K):
        x = x_lo + (0.5 - x_lo) * (k + 0.5) / K
        w = width(word, x, y_for_gamma(x, g))
        if w > best:
            best = w
    return best


def cmd_death(args):
    w = family_word(args.a, args.b)
    # coarse scan for the alive interval
    alive = []
    g = 90.25
    while g < 179.5:
        if alive_max(w, g, args.points) > 0:
            alive.append(g)
        g += args.step
    if not alive:
        print(json.dumps({"a": args.a, "b": args.b, "alive": None}))
        return
    lo, hi = alive[-1], alive[-1] + args.step
    for _ in range(args.iters):
        mid = 0.5 * (lo + hi)
        if alive_max(w, mid, args.points) > 0:
            lo = mid
        else:
            hi = mid
    res = {"a": args.a, "b": args.b, "len": len(w),
           "alive_coarse": [alive[0], alive[-1]],
           "death_between": [lo, hi],
           "points": args.points}
    if args.b == args.a:
        res["predicted_death"] = 180.0 * args.a / (args.a + 1)
    elif args.b == args.a - 1:
        res["predicted_death"] = 180.0 * (1 - (2 * args.a - 1)
                                          / (2.0 * args.a * args.a))
    if args.out:
        json.dump(res, open(args.out, "w"), indent=1)
    print(json.dumps(res))


def cmd_sliver(args):
    words = [tuple(int(c) for c in line.strip())
             for line in open(args.words) if line.strip()]
    best = (-math.inf, None, None)
    npos = 0
    for k in range(args.xpoints):
        x = 0.002 + (0.5 - 0.002) * (k + 0.5) / args.xpoints
        y = y_for_gamma(x, args.gamma)
        for w in words:
            wd = width(w, x, y)
            if wd > 0:
                npos += 1
            if wd > best[0]:
                best = (wd, "".join(map(str, w)), x)
    print(json.dumps({"gamma": args.gamma, "n_words": len(words),
                      "xpoints": args.xpoints,
                      "n_positive_samples": npos,
                      "max_width": best[0], "at_word": best[1],
                      "at_x": best[2]}))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("death")
    c.add_argument("--a", type=int, required=True)
    c.add_argument("--b", type=int, required=True)
    c.add_argument("--points", type=int, default=400)
    c.add_argument("--step", type=float, default=0.5)
    c.add_argument("--iters", type=int, default=26)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_death)
    c = sub.add_parser("sliver")
    c.add_argument("--gamma", type=float, required=True)
    c.add_argument("--words", required=True)
    c.add_argument("--xpoints", type=int, default=2000)
    c.set_defaults(fn=cmd_sliver)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
