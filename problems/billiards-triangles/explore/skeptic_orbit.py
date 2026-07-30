#!/usr/bin/env python3
"""SKEPTIC re-implementation: exact-rational periodic-orbit checker.

Written from scratch for the adversarial review of the census attempt.
Deliberately different implementation choices from both tier-0 tools and
census.py:

  * unfolding reflects ALL THREE vertices of the current copy across the
    gate line (unfold.py moves only the opposite vertex via a precomputed
    foot parameter; census.py's float layer copies that trick);
  * the translation property is checked GEOMETRICALLY: after unfolding,
    all three vertex displacements must be equal (unfold.py decides it by
    integer bookkeeping on the word alone);
  * the corridor is intersected from scratch, and the orbit is then
    reconstructed and verified by my own billiard simulator (own
    ray/segment solve, own direction reflection), checking the struck
    sequence, exact closure of position and direction, and that every
    bounce is strictly interior to its side.

Everything is Fraction arithmetic: a PASS/FAIL at a rational apex is exact.

Only shared convention (not code): side s is opposite vertex s, i.e.
side 0 = BC, side 1 = CA, side 2 = AB, with A=(0,0), B=(1,0), C the apex.

Subcommands:
  selfcheck                    calibration on known cases (Fagnano etc.)
  check-anchors FILES          verify sampled certificates (centre + 4
                               near-corner apexes of each certified box)
  invariance                   exact rotation/reversal corridor invariance
  pointsweep-recheck           exact corridor sign for every word of one
                               length at apexes taken from a pointsweep file
"""

import argparse
import json
import random
import sys
from fractions import Fraction as F

SIDE_ENDS = {0: (1, 2), 1: (2, 0), 2: (0, 1)}


# ---------------------------------------------------------------- geometry

def reflect_point(p, u, v):
    """Reflect point p across the line through u and v. Exact."""
    dx, dy = v[0] - u[0], v[1] - u[1]
    n2 = dx * dx + dy * dy
    wx, wy = p[0] - u[0], p[1] - u[1]
    t = (wx * dx + wy * dy) / n2
    fx, fy = u[0] + t * dx, u[1] + t * dy   # foot of perpendicular
    return (2 * fx - p[0], 2 * fy - p[1])


def reflect_dir(d, u, v):
    """Reflect direction d across the direction of the line (u,v)."""
    ex, ey = v[0] - u[0], v[1] - u[1]
    n2 = ex * ex + ey * ey
    t = (d[0] * ex + d[1] * ey) / n2
    return (2 * t * ex - d[0], 2 * t * ey - d[1])


def unfold_chain(word, apex):
    """Unfold the triangle along the word, reflecting all three vertices
    each step. Returns (gates, tau) or (None, reason).

    tau is accepted only if the three vertex displacements after the full
    word are identical (geometric translation test) and nonzero."""
    tri = ((F(0), F(0)), (F(1), F(0)), (F(apex[0]), F(apex[1])))
    first = tri
    gates = []
    for s in word:
        i, j = SIDE_ENDS[s]
        u, v = tri[i], tri[j]
        gates.append((u, v))
        tri = tuple(reflect_point(p, u, v) for p in tri)
    shifts = [(tri[k][0] - first[k][0], tri[k][1] - first[k][1])
              for k in range(3)]
    if not (shifts[0] == shifts[1] == shifts[2]):
        return None, "unfolding is not a translation (geometric test)"
    tau = shifts[0]
    if tau == (0, 0):
        return None, "translation vector is zero"
    return gates, tau


def corridor(word, apex):
    """Exact corridor (lo, hi) of the word at a rational apex, in the
    projection onto the normal of tau; None if not a translation word."""
    gates_tau = unfold_chain(word, apex)
    if gates_tau[0] is None:
        return None
    gates, tau = gates_tau
    nx, ny = -tau[1], tau[0]
    lo, hi = None, None
    for u, v in gates:
        a = u[0] * nx + u[1] * ny
        b = v[0] * nx + v[1] * ny
        if a > b:
            a, b = b, a
        lo = a if lo is None else max(lo, a)
        hi = b if hi is None else min(hi, b)
    return lo, hi, tau, gates


def is_bounce_word(word):
    return (len(word) >= 2 and all(s in (0, 1, 2) for s in word)
            and all(word[i] != word[(i + 1) % len(word)]
                    for i in range(len(word))))


# ---------------------------------------------------------------- simulator

def bounce(tri, p, d, skip):
    """First side struck from p in direction d, ignoring side `skip`.
    Returns (side, point, along) with `along` strictly in (0,1), or None."""
    best = None
    for s in (0, 1, 2):
        if s == skip:
            continue
        i, j = SIDE_ENDS[s]
        u, v = tri[i], tri[j]
        ex, ey = v[0] - u[0], v[1] - u[1]
        den = d[0] * ey - d[1] * ex
        if den == 0:
            continue
        rx, ry = u[0] - p[0], u[1] - p[1]
        t = (rx * ey - ry * ex) / den
        a = (rx * d[1] - ry * d[0]) / den
        if t <= 0:
            continue
        if not (F(0) < a < F(1)):
            continue        # vertex hit (a in {0,1}) or off the segment
        if best is None or t < best[0]:
            best = (t, s, (p[0] + t * d[0], p[1] + t * d[1]), a)
    if best is None:
        return None
    return best[1], best[2], best[3]


def verify_orbit(word, apex):
    """Full exact verification of a stable periodic orbit with this word at
    a rational apex. Returns (True, info) or (False, reason)."""
    if not is_bounce_word(word) or len(word) % 2:
        return False, "not an even bounce word"
    if not (apex[1] > 0):
        return False, "degenerate apex"
    c = corridor(word, apex)
    if c is None:
        return False, "not a translation word here"
    lo, hi, tau, gates = c
    if not lo < hi:
        return False, f"corridor empty (width {float(hi - lo):.3g})"
    m = (lo + hi) / 2
    nx, ny = -tau[1], tau[0]
    u, v = gates[0]
    ex, ey = v[0] - u[0], v[1] - u[1]
    den = ex * nx + ey * ny
    if den == 0:
        return False, "first gate parallel to corridor"
    a = (m - (u[0] * nx + u[1] * ny)) / den
    if not (F(0) < a < F(1)):
        return False, "corridor midline misses the first gate interior"
    start = (u[0] + a * ex, u[1] + a * ey)
    tri = ((F(0), F(0)), (F(1), F(0)), (F(apex[0]), F(apex[1])))
    i0, j0 = SIDE_ENDS[word[0]]
    d = reflect_dir(tau, tri[i0], tri[j0])

    p, dcur, skip = start, d, word[0]
    struck = []
    for _ in range(len(word)):
        hit = bounce(tri, p, dcur, skip)
        if hit is None:
            return False, "trajectory escaped or hit a vertex"
        s, p, _ = hit
        dcur = reflect_dir(dcur, tri[SIDE_ENDS[s][0]], tri[SIDE_ENDS[s][1]])
        struck.append(s)
        skip = s
    want = list(word[1:]) + [word[0]]
    if struck != want:
        return False, f"struck {struck} != word rotation {want}"
    if p != start or dcur != d:
        return False, "did not close (position or direction)"
    return True, f"exact periodic orbit confirmed, {len(word)} bounces"


def box_samples(centre, half):
    """Centre plus four corner apexes pulled in to 3/4 of the half-width."""
    cx, cy = F(centre[0]), F(centre[1])
    h = F(half) * 3 / 4
    return [(cx, cy), (cx - h, cy - h), (cx - h, cy + h),
            (cx + h, cy - h), (cx + h, cy + h)]


# ---------------------------------------------------------------- commands

def cmd_selfcheck(_args):
    ok = True
    fag = (0, 1, 2, 0, 1, 2)
    good, msg = verify_orbit(fag, (F(1, 2), F(4, 5)))
    print("fagnano acute:", good, msg); ok &= good
    good, msg = verify_orbit(fag, (F(1, 2), F(1, 10)))
    print("fagnano obtuse rejected:", not good, msg); ok &= not good
    good, msg = verify_orbit((0, 1, 0, 2, 1, 2), (F(1, 2), F(4, 5)))
    print("junk word rejected:", not good, msg); ok &= not good
    # translation test must fail on a single traversal
    c = corridor((0, 1, 2), (F(1, 2), F(4, 5)))
    print("single traversal not a translation:", c is None); ok &= c is None
    print("SELFCHECK", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


def iter_anchor_records(paths):
    for path in paths:
        if path.endswith(".json"):
            r = json.load(open(path))
            if r.get("anchor"):
                yield {"src": path, "mid": r["anchor_apex"],
                       "word": r["anchor"]["word"],
                       "half_width": r["anchor"]["half_width"]}
            continue
        for line in open(path):
            if not line.strip():
                continue
            r = json.loads(line)
            if r.get("status") == "ANCHORED":
                yield {"src": path, "mid": r["mid"], "word": r["word"],
                       "half_width": r["half_width"]}
            elif r.get("anchor"):
                yield {"src": path, "mid": r["anchor_at"],
                       "word": r["anchor"]["word"],
                       "half_width": r["anchor"]["half_width"]}


def cmd_check_anchors(args):
    recs = list(iter_anchor_records(args.src))
    rng = random.Random(args.seed)
    if args.sample and args.sample < len(recs):
        recs = rng.sample(recs, args.sample)
    n_ok = n_fail = 0
    out = open(args.out, "w") if args.out else None
    for r in recs:
        word = tuple(int(c) for c in r["word"])
        h = F(r["half_width"])
        entry = {"src": r["src"].split("/")[-1], "mid": r["mid"],
                 "len": len(word), "half_width": r["half_width"]}
        results = []
        for apex in box_samples(r["mid"], h):
            good, msg = verify_orbit(word, apex)
            results.append(good)
            if not good:
                entry["fail_msg"] = msg
                entry["fail_apex"] = [str(apex[0]), str(apex[1])]
                break
        entry["ok"] = all(results)
        entry["n_apexes"] = len(results)
        n_ok += entry["ok"]
        n_fail += not entry["ok"]
        if out:
            out.write(json.dumps(entry) + "\n")
        print(("OK  " if entry["ok"] else "FAIL"),
              entry["src"], r["mid"], "len", len(word))
    if out:
        out.close()
    print(f"check-anchors: {n_ok} confirmed, {n_fail} failed "
          f"of {len(recs)} sampled certificates x 5 apexes")
    sys.exit(1 if n_fail else 0)


def cmd_invariance(args):
    """Exact corridor width must be identical for every rotation and every
    rotation of the reversal, at random rational apexes."""
    words = [tuple(int(c) for c in line.strip())
             for line in open(args.words) if line.strip()]
    rng = random.Random(args.seed)
    words = rng.sample(words, min(args.nwords, len(words)))
    apexes = [(F(rng.randint(1, 2047), 4096), F(rng.randint(1, 1300), 4096))
              for _ in range(args.napex)]
    bad = 0
    for w in words:
        n = len(w)
        variants = [w[r:] + w[:r] for r in range(n)]
        rw = w[::-1]
        variants += [rw[r:] + rw[:r] for r in range(n)]
        for apex in apexes:
            widths = set()
            for v in variants:
                c = corridor(v, apex)
                if c is None:
                    widths.add(None)
                else:
                    widths.add(c[1] - c[0])
            if len(widths) != 1:
                bad += 1
                print("INVARIANCE FAIL", "".join(map(str, w)), apex, widths)
    print(f"invariance: {len(words)} words x {len(apexes)} apexes x "
          f"all rotations+reversals: {bad} violations")
    sys.exit(1 if bad else 0)


def cmd_pointsweep_recheck(args):
    """Exact corridor sign of EVERY word in a length file at apexes drawn
    from an existing pointsweep file (positions re-used, verdicts mine)."""
    words = [tuple(int(c) for c in line.strip())
             for line in open(args.words) if line.strip()]
    src = json.load(open(args.src))
    pts = src["points"]
    rng = random.Random(args.seed)
    if args.sample and args.sample < len(pts):
        pts = rng.sample(pts, args.sample)
    n_pos = 0
    pos_detail = []
    for pt in pts:
        x, y = F(pt["x"]), F(pt["y"])
        for w in words:
            c = corridor(w, (x, y))
            if c is not None and c[1] > c[0]:
                n_pos += 1
                pos_detail.append(("".join(map(str, w)), pt["x"], pt["y"]))
    print(f"pointsweep-recheck {args.src}: {len(pts)} apexes x "
          f"{len(words)} words -> {n_pos} positive exact corridors")
    for d in pos_detail[:10]:
        print("  POSITIVE", d)
    theirs = set()
    for pt in src["points"]:
        for tw in pt.get("true_words", []):
            theirs.add(tw)
    print("  (their file reports n_true_total =", src["n_true_total"], ")")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selfcheck").set_defaults(fn=cmd_selfcheck)
    c = sub.add_parser("check-anchors")
    c.add_argument("--src", nargs="+", required=True)
    c.add_argument("--sample", type=int, default=0)
    c.add_argument("--seed", type=int, default=20260730)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_check_anchors)
    c = sub.add_parser("invariance")
    c.add_argument("--words", required=True)
    c.add_argument("--nwords", type=int, default=6)
    c.add_argument("--napex", type=int, default=3)
    c.add_argument("--seed", type=int, default=20260730)
    c.set_defaults(fn=cmd_invariance)
    c = sub.add_parser("pointsweep-recheck")
    c.add_argument("--words", required=True)
    c.add_argument("--src", required=True)
    c.add_argument("--sample", type=int, default=0)
    c.add_argument("--seed", type=int, default=20260730)
    c.set_defaults(fn=cmd_pointsweep_recheck)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
