#!/usr/bin/env python3
"""Word census by length for periodic billiard orbits in obtuse triangles.

Standard-library only.  Everything here is a *candidate pipeline* feeding the
tier-0 harness: the only object ever claimed is a (word, box) pair for which
`harness/billiards-triangles/unfold.py` returns TRUE in exact arithmetic, and
those claims are then attacked with `verify_cover.py` by direct simulation.
The float code in this file is a prefilter and never appears in a certificate.

SCOPE: an open-region certificate only sees orbits stable under perturbation
of the triangle.  Nothing here says anything about unstable/isolated orbits.

Parameter space.  Every obtuse triangle is similar to exactly one triangle
with the obtuse vertex at C, A = (0,0), B = (1,0), apex C = (x, y) in the
open half-disk (x-1/2)^2 + y^2 < 1/4, y > 0, and (up to the reflection that
swaps A and B, i.e. congruence) x <= 1/2.  The census region is therefore

    R = { (x, y) : 0 < x <= 1/2, y > 0, (x-1/2)^2 + y^2 < 1/4 },

of area pi/16.  The largest angle gamma is the angle at C; gamma -> 90 deg at
the arc, gamma -> 180 deg as y -> 0.

Word universe and pruning (the exact rules; see `words` subcommand):
  U(L) = all tuples w over {0,1,2} with
    (W1) no two cyclically-adjacent letters equal (bounce words);
    (W2) |w| even (odd-length unfoldings are never translations);
    (W3) the composed unfolding isometry is a translation for every triangle
         at once: linear part is a rotation by p*alpha+q*beta+r*pi with
         p = q = 0 and r even (integer bookkeeping identical to
         unfold.linear_part, re-implemented incrementally and cross-checked
         against the harness by `selftest`);
    (W4) canonical under rotation and reversal: w == min over all cyclic
         rotations of w and of reversed(w).  Sound because (a) the corridor
         of a rotated word is the same set: the unfolded chain of rot(w)
         differs by replacing gate 1 with gate n+1 = gate 1 + tau, and tau
         projects to 0 on the corridor normal; (b) the unfolded chain of
         reversed(w) is congruent to that of w (run the chain backwards and
         reflect), so its corridor is isometric.
    (W5) primitive-root pruning: w = u^k with k >= 2 is dropped when u itself
         satisfies (W1)-(W3), because the gate-projection multiset of u^k
         equals that of u, so the corridors are identical.  (u^k with u NOT a
         translation word is kept: e.g. Fagnano 012012 = (012)^2.)
  The 0<->1 letter swap is NOT quotiented: it maps the certified region to
  its mirror x -> 1-x, and the census region is already the x <= 1/2 half, so
  swapped words act differently there and both members are kept.

Two-layer census (and why).  The harness's affine enclosure of the unfolding
translation vector grows rapidly with box width: measured on the length-14
word 01212020121202 at apex ~(0.44, 0.38), the tau-enclosure has width ~470
on a box of size 1/64, ~26 at 1/256, and the verdict only becomes TRUE at box
size ~1/4096 -- while the *point* corridor there has width ~1.04.  Direct
whole-cell certification at map resolution is therefore hopeless for the word
lengths that matter, and the census is built in two layers:

  layer 1 (map)     -- float corridor sweep at the midpoint of every grid
                       cell, all canonical words by ascending length: the
                       minimal length with a positive corridor ("candidate
                       length") per cell.  Never claimed as a certificate.
  layer 2 (anchor)  -- for every cell with a candidate, an exact harness
                       certificate on a small centred dyadic box, trying a
                       ladder of half-widths (1/32768 up to 1/128) and
                       keeping the largest TRUE.  These small boxes are the
                       exactly-certified regions; a cell is "anchored" when
                       its midpoint neighbourhood is certified.

Subcommands:
  selftest              validation suite (cross-checks against the harness)
  words --max-len L     enumerate the universe, JSONL per length + counts
  map ...               layer-1 candidate sweep over the cell grid
  anchor ...            layer-2 exact certificates at candidate cells
  growth ...            minimal certificate length along a line of increasing
                        gamma (the kill-condition measurement)
  simcheck ...          re-verify claimed certificates by direct simulation
  summarize ...         aggregate everything into the coverage map
"""

import argparse
import json
import math
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
# tier-0 harness (the certificate path; this file never certifies on its own)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..",
                                "harness", "billiards-triangles"))
import unfold  # noqa: E402
import verify_cover  # noqa: E402
from unfold import Iv, Q, TRUE  # noqa: E402


# ----------------------------------------------------------------------
# word enumeration
# ----------------------------------------------------------------------

# Incremental copy of unfold.INITIAL_DIRECTIONS / linear_part bookkeeping.
INITIAL_DIRECTIONS = ((0, -1, 1), (1, 0, 0), (0, 0, 0))


def canonical(w):
    """Lexicographic minimum over all rotations of w and of reversed(w)."""
    n = len(w)
    best = w
    for base in (w, w[::-1]):
        for r in range(n):
            cand = base[r:] + base[:r]
            if cand < best:
                best = cand
    return best


def _is_translation_state(ang):
    return ang[0] == 0 and ang[1] == 0 and ang[2] % 2 == 0


def _primitive_root_is_translation(w):
    """True when w = u^k, k >= 2, with u a bounce word satisfying (W1)-(W3)."""
    n = len(w)
    for p in range(2, n):
        if n % p:
            continue
        if w != w[:p] * (n // p):
            continue
        u = w[:p]
        if p % 2:
            continue  # odd u is never a translation word
        if u[0] == u[-1]:
            continue  # not cyclically valid
        if unfold.is_translation_word(u):
            return True
    return False


def enumerate_translation_words(max_len, emit):
    """DFS over bounce words, calling emit(word) for every canonical
    translation word of even length in [4, max_len].  (Length 2 words 01, 02,
    12 are never translations: check r parity -- the DFS just finds none.)

    The state carried down the tree is exactly the (ang, dirs) recursion of
    unfold.linear_part; kind alternates with depth so even depth = rotation.
    """
    sys.setrecursionlimit(max_len + 100)

    def rec(word, ang, dirs, depth):
        if depth >= 2 and depth % 2 == 0 and word[0] != word[-1] \
                and _is_translation_state(ang):
            w = tuple(word)
            if canonical(w) == w and not _primitive_root_is_translation(w):
                emit(w)
        if depth == max_len:
            return
        prev = word[-1] if word else -1
        for s in (0, 1, 2):
            if s == prev:
                continue
            t = dirs[s]
            nang = (2 * t[0] - ang[0], 2 * t[1] - ang[1], 2 * t[2] - ang[2])
            ndirs = tuple((2 * t[0] - d[0], 2 * t[1] - d[1], 2 * t[2] - d[2])
                          for d in dirs)
            word.append(s)
            rec(word, nang, ndirs, depth + 1)
            word.pop()

    rec([], (0, 0, 0), INITIAL_DIRECTIONS, 0)


def cmd_words(args):
    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)
    by_len = {}

    def emit(w):
        by_len.setdefault(len(w), []).append(w)

    enumerate_translation_words(args.max_len, emit)
    counts = {}
    for n in sorted(by_len):
        words = sorted(by_len[n])
        counts[n] = len(words)
        path = os.path.join(out_dir, f"words_L{n:02d}.jsonl")
        with open(path + ".partial", "w") as f:
            for w in words:
                f.write("".join(map(str, w)) + "\n")
        os.replace(path + ".partial", path)
    with open(os.path.join(out_dir, "counts.json"), "w") as f:
        json.dump({"max_len": args.max_len,
                   "pruning": "W1-W5 as in census.py docstring",
                   "counts_by_length": counts,
                   "total": sum(counts.values())}, f, indent=1)
    print(json.dumps(counts))
    print("total canonical translation words:", sum(counts.values()))


def load_words(words_dir, max_len):
    """Canonical words by length, ascending, from the words stage."""
    out = {}
    for n in range(4, max_len + 1, 2):
        path = os.path.join(words_dir, f"words_L{n:02d}.jsonl")
        if not os.path.exists(path):
            continue
        with open(path) as f:
            out[n] = [tuple(int(c) for c in line.strip())
                      for line in f if line.strip()]
    return out


# ----------------------------------------------------------------------
# float prefilter (never part of any certificate)
# ----------------------------------------------------------------------

def float_corridor(word, ax, ay):
    """Corridor width of the unfolding at one apex, in floats.

    Mirrors unfold.unfold/certify at a point: foot parameters once, reflect
    only the vertex opposite the mirror, project gates on the normal of tau.
    Returns the corridor width (positive = candidate) or None (degenerate).
    """
    A = (0.0, 0.0)
    B = (1.0, 0.0)
    C = (ax, ay)
    tri = [A, B, C]
    ends = ((1, 2), (2, 0), (0, 1))
    feet = []
    for s in (0, 1, 2):
        i, j = ends[s]
        u, v = tri[i], tri[j]
        dx, dy = v[0] - u[0], v[1] - u[1]
        p = tri[s]
        feet.append(((p[0] - u[0]) * dx + (p[1] - u[1]) * dy)
                    / (dx * dx + dy * dy))
    t = tri
    gates = []
    for s in word:
        i, j = ends[s]
        u, v = t[i], t[j]
        gates.append((u, v))
        c2 = 2.0 * feet[s]
        p = t[s]
        np_ = (2.0 * u[0] + c2 * (v[0] - u[0]) - p[0],
               2.0 * u[1] + c2 * (v[1] - u[1]) - p[1])
        t2 = list(t)
        t2[s] = np_
        t = t2
    tau = (t[0][0] - 0.0, t[0][1] - 0.0)  # translation word: A's shift
    if tau[0] * tau[0] + tau[1] * tau[1] < 1e-18:
        return None
    nx, ny = -tau[1], tau[0]
    lo = -math.inf
    hi = math.inf
    for u, v in gates:
        a = u[0] * nx + u[1] * ny
        b = v[0] * nx + v[1] * ny
        if a > b:
            a, b = b, a
        if a > lo:
            lo = a
        if b < hi:
            hi = b
        if lo >= hi:
            # keep going is pointless; width is already non-positive
            return lo, hi, hi - lo
    return lo, hi, hi - lo


def float_width(word, ax, ay):
    r = float_corridor(word, ax, ay)
    if r is None:
        return -math.inf
    # normalize by |tau|?  No: only the sign is used downstream.
    return r[2]


# ----------------------------------------------------------------------
# region geometry
# ----------------------------------------------------------------------

def box_in_region(x0, x1, y0, y1):
    """Whole closed box inside the open census region R?  Disk is convex, so
    corner checks suffice; exact rational arithmetic."""
    if y0 <= 0 or x1 > Q(1, 2):
        return False
    for x in (x0, x1):
        for y in (y0, y1):
            if (x - Q(1, 2)) ** 2 + y ** 2 >= Q(1, 4):
                return False
    return True


def box_meets_region(x0, x1, y0, y1):
    """Does the open box meet R at all?  Nearest point of box to the disk
    centre (1/2, 0); exact."""
    if y1 <= 0 or x0 >= Q(1, 2):
        # x0 >= 1/2 boxes belong to the mirrored half (x1 > 1/2 handled by
        # caller clipping at 1/2)
        return False
    cx, cy = Q(1, 2), Q(0)
    nx = min(max(cx, x0), x1)
    ny = min(max(cy, y0), y1)
    return (nx - cx) ** 2 + (ny - cy) ** 2 < Q(1, 4)


def gamma_deg(x, y):
    """Angle at C in degrees (float, for reporting only)."""
    ax, ay = -float(x), -float(y)
    bx, by = 1.0 - float(x), -float(y)
    d = ax * bx + ay * by
    c = ax * by - ay * bx
    return math.degrees(math.atan2(abs(c), d))


def box_gamma_range(x0, x1, y0, y1):
    gs = [gamma_deg(x, y) for x in (x0, x1) for y in (y0, y1)]
    return min(gs), max(gs)


# ----------------------------------------------------------------------
# layer 1: float candidate scan at a point
# ----------------------------------------------------------------------

def scan_point(words_by_len, fmx, fmy, top=4):
    """Float sweep of the whole universe at one apex.  Returns
    {length: (n_positive, [(width, word) descending, at most `top`])}
    for lengths with at least one positive-corridor word, plus the minimal
    such length (or None)."""
    out = {}
    cand_len = None
    for n in sorted(words_by_len):
        pos = []
        for w in words_by_len[n]:
            width = float_width(w, fmx, fmy)
            if width > 0.0:
                pos.append((width, w))
        if pos:
            pos.sort(key=lambda t: (-t[0], t[1]))
            out[n] = (len(pos), pos[:top])
            if cand_len is None:
                cand_len = n
    return cand_len, out


# ----------------------------------------------------------------------
# layer 2: exact anchor certification
# ----------------------------------------------------------------------

def variants(w):
    """The word, its half-rotation, its reversal, and the reversal's
    half-rotation.  The true corridor is the same set for all rotations and
    reversals, but the affine *enclosure* is order-dependent (gates late in
    the chain carry wider enclosures), so when one variant is UNKNOWN on a
    box another can still be TRUE.  Four variants is a fixed, recorded cap.
    """
    n = len(w)
    k = n // 2
    out = []
    for v in (w, w[k:] + w[:k], w[::-1], (w[::-1])[k:] + (w[::-1])[:k]):
        if v not in out:
            out.append(v)
    return out


def anchor_at(words_by_len, mx, my, scan=None, exact_cap=4):
    """Minimal-length exact certificate on a small box centred at (mx, my).

    For each length with float candidates, in ascending order: take the top
    `exact_cap` candidates by float corridor width, and for each try harness
    certification on a centred dyadic box (in up to 4 rotation/reversal
    variants).  The starting half-width is word-length adaptive: measured on
    this harness, the box size at which the affine enclosure of a length-n
    word's unfolding becomes decidable shrinks by ~2 bits per letter
    (TRUE at half-width ~2^-13 for n=14, ~2^-19 for n=20, ~2^-25 for n=26),
    so the first attempt is at half-width 2^-(n+2), with one retry at
    2^-(n+6) for thin corridors.  On success, climb by doubling and keep the
    largest half-width still TRUE (a certificate provably survives
    shrinking, so the climb only ever adds area; cap 1/128 = half a map
    cell).  Returns a dict or None.
    """
    if scan is None:
        _, scan = scan_point(words_by_len, float(mx), float(my),
                             top=exact_cap)
    for n in sorted(scan):
        _, tops = scan[n]
        for _, w in tops[:exact_cap]:
            for v in variants(w):
                start = None
                for e in (n + 2, n + 6):
                    h0 = Q(1, 2 ** e)
                    if not box_in_region(mx - h0, mx + h0,
                                         my - h0, my + h0):
                        continue
                    if unfold.certify(
                            v, (Iv(mx - h0, mx + h0),
                                Iv(my - h0, my + h0)))[0] == TRUE:
                        start = e
                        break
                if start is None:
                    continue
                best = Q(1, 2 ** start)
                for e in range(start - 1, 6, -1):
                    h = Q(1, 2 ** e)
                    if not box_in_region(mx - h, mx + h, my - h, my + h):
                        break
                    if unfold.certify(
                            v, (Iv(mx - h, mx + h),
                                Iv(my - h, my + h)))[0] == TRUE:
                        best = h
                    else:
                        break
                return {"len": n, "word": "".join(map(str, v)),
                        "half_width": str(best)}
    return None


# ----------------------------------------------------------------------
# map subcommand (layer 1)
# ----------------------------------------------------------------------

def cmd_map(args):
    words = load_words(args.words, args.max_len)
    if not words:
        sys.exit("no word files found; run the words stage first")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    done = set()
    if os.path.exists(args.out):
        with open(args.out) as f:
            for line in f:
                if line.strip():
                    r = json.loads(line)
                    done.add((r["i"], r["j"]))
    g = args.grid  # cells of size 1/(2g) covering [0,1/2] x [0,1/2]
    h = Q(1, 2 * g)
    rows = range(args.row_lo, min(args.row_hi, g))
    out = open(args.out, "a")
    n_done = 0
    for j in rows:               # j indexes y
        y0, y1 = j * h, (j + 1) * h
        for i in range(g):       # i indexes x
            if (i, j) in done:
                continue
            x0, x1 = i * h, (i + 1) * h
            if not box_meets_region(x0, x1, y0, y1):
                continue
            mx, my = (x0 + x1) / 2, (y0 + y1) / 2
            rec = {"i": i, "j": j, "grid": g,
                   "box": [str(x0), str(x1), str(y0), str(y1)],
                   "mid": [str(mx), str(my)]}
            if box_in_region(x0, x1, y0, y1):
                cand_len, scan = scan_point(words, float(mx), float(my),
                                            top=args.top)
                rec["status"] = "CAND" if cand_len is not None else "NOCAND"
                rec["cand_len"] = cand_len
                rec["pos_counts"] = {str(n): scan[n][0] for n in scan}
                rec["top"] = {str(n): [[round(wd, 6),
                                        "".join(map(str, w))]
                                       for wd, w in scan[n][1]]
                              for n in scan}
                rec["gamma"] = [round(v, 3) for v in
                                box_gamma_range(x0, x1, y0, y1)]
                rec["gamma_mid"] = round(gamma_deg(mx, my), 3)
            else:
                rec["status"] = "BOUNDARY"
            out.write(json.dumps(rec) + "\n")
            out.flush()
            n_done += 1
    out.close()
    print(f"map: wrote {n_done} new cells to {args.out}")


# ----------------------------------------------------------------------
# anchor subcommand (layer 2)
# ----------------------------------------------------------------------

def cmd_anchor(args):
    words = load_words(args.words, args.max_len)
    cells = []
    with open(args.src) as f:
        for line in f:
            if line.strip():
                r = json.loads(line)
                if r.get("status") == "CAND":
                    cells.append(r)
    done = set()
    if os.path.exists(args.out):
        with open(args.out) as f:
            for line in f:
                if line.strip():
                    r = json.loads(line)
                    done.add((r["i"], r["j"]))
    out = open(args.out, "a")
    n_done = 0
    for r in cells:
        if (r["i"], r["j"]) in done:
            continue
        if not (args.row_lo <= r["j"] < args.row_hi):
            continue
        mx, my = Q(r["mid"][0]), Q(r["mid"][1])
        res = anchor_at(words, mx, my, exact_cap=args.exact_cap)
        rec = {"i": r["i"], "j": r["j"], "mid": r["mid"],
               "gamma_mid": r.get("gamma_mid"),
               "cand_len": r.get("cand_len")}
        if res is None:
            rec["status"] = "NO_ANCHOR"
        else:
            rec["status"] = "ANCHORED"
            rec.update(res)
        out.write(json.dumps(rec) + "\n")
        out.flush()
        n_done += 1
    out.close()
    print(f"anchor: wrote {n_done} cells to {args.out}")


# ----------------------------------------------------------------------
# growth subcommand (kill-condition measurement)
# ----------------------------------------------------------------------

def y_for_gamma(x, gdeg):
    """y with angle at C = gdeg on the vertical line at x (float; probe
    placement only, then snapped to a dyadic)."""
    lo, hi = 1e-9, 0.5
    for _ in range(80):
        mid = (lo + hi) / 2
        if gamma_deg(x, mid) > gdeg:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def cmd_growth(args):
    words = load_words(args.words, args.max_len)
    x = Q(args.x)
    done = set()
    if os.path.exists(args.out):
        with open(args.out) as f:
            for line in f:
                if line.strip():
                    done.add(json.loads(line)["gamma_target"])
    out = open(args.out, "a")
    g = args.gamma_lo
    while g <= args.gamma_hi + 1e-9:
        gt = round(g, 4)
        if gt in done:
            g += args.gamma_step
            continue
        yf = y_for_gamma(float(x), gt)
        # snap centre to a dyadic just off the target curve
        yc = Q(round(yf * args.denom), args.denom)
        rec = {"gamma_target": gt, "x": str(x), "yc": str(yc)}
        if not box_in_region(x - Q(1, 32768), x + Q(1, 32768),
                             yc - Q(1, 32768), yc + Q(1, 32768)):
            rec["status"] = "OUTSIDE"
        else:
            rec["gamma_mid"] = round(gamma_deg(x, yc), 4)
            cand_len, scan = scan_point(words, float(x), float(yc),
                                        top=args.exact_cap)
            rec["cand_len"] = cand_len
            rec["pos_counts"] = {str(n): scan[n][0] for n in scan}
            res = anchor_at(words, x, yc, scan=scan,
                            exact_cap=args.exact_cap)
            if res is None:
                rec["status"] = "CAND_ONLY" if cand_len else "NOCAND"
            else:
                rec["status"] = "ANCHORED"
                rec.update(res)
        out.write(json.dumps(rec) + "\n")
        out.flush()
        g += args.gamma_step
    out.close()
    print(f"growth: appended to {args.out}")


def cmd_arcsweep(args):
    """Minimal candidate/certified word length as a function of the obtuse
    angle gamma: sample --points apexes along each constant-gamma arc
    (equally spaced in x over (x_lo, 1/2)), record cand_len at each, and
    exactly anchor the best (minimal cand_len, ties by larger float width)
    point of the arc.  This is the kill-condition measurement: how fast does
    the needed word length grow with gamma?"""
    words = load_words(args.words, args.max_len)
    done = set()
    if os.path.exists(args.out):
        with open(args.out) as f:
            for line in f:
                if line.strip():
                    done.add(json.loads(line)["gamma_target"])
    out = open(args.out, "a")
    g = args.gamma_lo
    while g <= args.gamma_hi + 1e-9:
        gt = round(g, 4)
        if gt in done:
            g += args.gamma_step
            continue
        pts = []
        for k in range(args.points):
            xf = args.x_lo + (0.5 - args.x_lo) * (k + 0.5) / args.points
            yf = y_for_gamma(xf, gt)
            x = Q(round(xf * 32768), 32768)
            y = Q(round(yf * 32768), 32768)
            e = Q(1, 2 ** 20)
            if not box_in_region(x - e, x + e, y - e, y + e):
                continue
            cand_len, scan = scan_point(words, float(x), float(y),
                                        top=args.exact_cap)
            pts.append({"x": str(x), "y": str(y),
                        "gamma": round(gamma_deg(x, y), 4),
                        "cand_len": cand_len,
                        "best_width": (round(scan[cand_len][1][0][0], 6)
                                       if cand_len else None),
                        "_scan": scan, "_xy": (x, y)})
        rec = {"gamma_target": gt,
               "samples": [{k: p[k] for k in
                            ("x", "y", "gamma", "cand_len", "best_width")}
                           for p in pts]}
        with_c = [p for p in pts if p["cand_len"] is not None]
        rec["min_cand_len"] = (min(p["cand_len"] for p in with_c)
                               if with_c else None)
        rec["n_samples"] = len(pts)
        rec["n_with_cand"] = len(with_c)
        if with_c and args.anchor_best:
            best = min(with_c,
                       key=lambda p: (p["cand_len"], -p["best_width"]))
            res = anchor_at(words, *best["_xy"], scan=best["_scan"],
                            exact_cap=args.exact_cap)
            rec["anchor_at"] = [best["x"], best["y"]]
            rec["anchor"] = res
        out.write(json.dumps(rec) + "\n")
        out.flush()
        g += args.gamma_step
    out.close()
    print(f"arcsweep: appended to {args.out}")


def family_word(a, b):
    """W(a,b) = (0 (12)^a (02)^b)^2.  W(1,0) is Fagnano's word.  The census
    found W(2,1), W(2,2), W(3,2), W(3,3) to be the last-surviving words at
    lengths 14, 18, 22, 26; this builds the family beyond the enumerated
    universe."""
    prim = (0,) + (1, 2) * a + (0, 2) * b
    return prim * 2


def cmd_family(args):
    """Scan one family word W(a,b): confirm it is a translation word, map its
    alive gamma-interval on the arc grid (float), bisect its death angle,
    and exactly anchor + simulate it at its widest sample point."""
    w = family_word(args.a, args.b)
    rec = {"a": args.a, "b": args.b, "word": "".join(map(str, w)),
           "len": len(w)}
    rec["is_word"] = unfold.is_word(w)
    rec["is_translation_word"] = unfold.is_translation_word(w)
    if not (rec["is_word"] and rec["is_translation_word"]):
        print(json.dumps(rec, indent=1))
        return

    def alive_max(g, K):
        best = (-1.0, None)
        for k in range(K):
            x = args.x_lo + (0.5 - args.x_lo) * (k + 0.5) / K
            y = y_for_gamma(x, g)
            wd = float_width(w, x, y)
            if wd > best[0]:
                best = (wd, x)
        return best

    # coarse alive scan
    alive_g = []
    g = 90.25
    while g < 179.5:
        wd, _ = alive_max(g, args.points)
        if wd > 0:
            alive_g.append(round(g, 4))
        g += args.step
    rec["alive_gammas_coarse"] = [alive_g[0], alive_g[-1]] if alive_g \
        else None
    rec["n_alive_gammas"] = len(alive_g)
    if alive_g:
        # bisect death upward from the last alive gamma
        lo, hi = alive_g[-1], alive_g[-1] + args.step
        for _ in range(args.iters):
            mid = (lo + hi) / 2
            if alive_max(mid, args.points)[0] > 0:
                lo = mid
            else:
                hi = mid
        rec["death_between"] = [lo, hi]
        # bisect birth downward from the first alive gamma
        lo2, hi2 = alive_g[0] - args.step, alive_g[0]
        for _ in range(args.iters):
            mid = (lo2 + hi2) / 2
            if alive_max(mid, args.points)[0] > 0:
                hi2 = mid
            else:
                lo2 = mid
        rec["birth_between"] = [lo2, hi2]
        # anchor at the widest point of the widest arc
        best = (-1.0, None, None)
        for g in alive_g:
            wd, x = alive_max(g, args.points)
            if wd > best[0]:
                best = (wd, x, g)
        wd, xf, g = best
        x = Q(round(xf * 32768), 32768)
        y = Q(round(y_for_gamma(float(x), g) * 32768), 32768)
        scan = {len(w): (1, [(wd, w)])}
        res = anchor_at({len(w): [w]}, x, y, scan=scan)
        rec["anchor_apex"] = [str(x), str(y)]
        rec["anchor_gamma"] = round(gamma_deg(x, y), 4)
        rec["anchor"] = res
        if res:
            simok, msg = verify_cover.verify_at(
                tuple(int(c) for c in res["word"]), (x, y))
            rec["sim_ok"] = simok
            rec["sim_msg"] = msg
    with open(args.out, "w") as f:
        json.dump(rec, f, indent=1)
    print(json.dumps({k: v for k, v in rec.items()
                      if k not in ("alive_gammas_coarse",)}, indent=1))


def cmd_bisect_death(args):
    """Measure the gamma at which the last length-n candidate dies (float).

    Phase 1: full scan of all length-n canonical words at K points on the
    gamma-lo arc; words positive somewhere form the 'alive set'.
    Phase 2: bisect gamma on the alive set (max over K points and alive
    words of the corridor width; positive = alive).
    Phase 3: full scan of ALL length-n words at gamma_death + margin to
    confirm the alive-set restriction did not miss a late-blooming family.
    Float throughout -- this measures a frontier, it certifies nothing.
    """
    words = load_words(args.words, args.len)[args.len]

    def xs(K):
        return [args.x_lo + (0.5 - args.x_lo) * (k + 0.5) / K
                for k in range(K)]

    def alive_words(g, wl, K):
        out = []
        for w in wl:
            for x in xs(K):
                if float_width(w, x, y_for_gamma(x, g)) > 0:
                    out.append(w)
                    break
        return out

    alive = alive_words(args.gamma_lo, words, args.points)
    print(f"alive at gamma={args.gamma_lo}: {len(alive)}/{len(words)}")
    lo, hi = args.gamma_lo, args.gamma_hi
    for _ in range(args.iters):
        mid = (lo + hi) / 2
        if alive_words(mid, alive, args.points):
            lo = mid
        else:
            hi = mid
    confirm = alive_words(hi + args.margin, words, args.points)
    last = alive_words(lo, alive, args.points)
    res = {"len": args.len, "points_per_arc": args.points,
           "x_lo": args.x_lo,
           "alive_at_gamma_lo": len(alive),
           "death_between": [lo, hi],
           "last_alive_words": ["".join(map(str, w)) for w in last],
           "confirm_gamma": hi + args.margin,
           "confirm_alive_full_scan": len(confirm),
           "confirm_words": ["".join(map(str, w)) for w in confirm]}
    with open(args.out, "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps({k: res[k] for k in
                      ("death_between", "last_alive_words",
                       "confirm_alive_full_scan")}, indent=1))


def cmd_pointsweep(args):
    """Exact-arithmetic refutation sweep: certify EVERY canonical word of
    one length at exact rational apexes snapped near a constant-gamma arc.
    At a zero-width box the harness enclosure is exact, so each verdict is
    decisive, and by rotation/reversal invariance of the corridor at a point
    a non-TRUE verdict for the canonical word rules out its whole class.
    Used to confirm float-detected family frontiers in exact arithmetic."""
    words = load_words(args.words, args.len)[args.len]
    denom = 32768
    results = []
    n_true = 0
    for k in range(args.points):
        xf = args.x_lo + (0.5 - args.x_lo) * (k + 0.5) / args.points
        yf = y_for_gamma(xf, args.gamma)
        x = Q(round(xf * denom), denom)
        y = Q(round(yf * denom), denom)
        if not ((x - Q(1, 2)) ** 2 + y ** 2 < Q(1, 4) and 0 < x <= Q(1, 2)
                and y > 0):
            continue
        pt = (Iv(x), Iv(y))
        true_words = []
        for w in words:
            if unfold.certify(w, pt)[0] == TRUE:
                true_words.append("".join(map(str, w)))
        n_true += len(true_words)
        results.append({"x": str(x), "y": str(y),
                        "gamma": round(gamma_deg(x, y), 4),
                        "true_words": true_words})
    with open(args.out, "w") as f:
        json.dump({"gamma_target": args.gamma, "len": args.len,
                   "n_words": len(words), "n_points": len(results),
                   "n_true_total": n_true, "points": results}, f, indent=1)
    print(f"pointsweep gamma={args.gamma} len={args.len}: "
          f"{len(results)} exact apexes x {len(words)} words, "
          f"{n_true} TRUE verdicts")


# ----------------------------------------------------------------------
# simcheck subcommand: attack certified boxes with the simulator
# ----------------------------------------------------------------------

def cmd_simcheck(args):
    """For every ANCHORED record, re-verify by direct simulation at the
    anchor centre (independent code path in verify_cover.py: exact rational
    billiard simulation, no unfolding).  With --full N, additionally run the
    9-sample check over the whole certified anchor box on every N-th record.
    """
    recs = []
    for path in args.src:
        with open(path) as f:
            for line in f:
                if line.strip():
                    r = json.loads(line)
                    if r.get("status") == "ANCHORED":
                        recs.append(r)
                    elif r.get("anchor"):  # arcsweep shape
                        recs.append({"status": "ANCHORED",
                                     "mid": r["anchor_at"],
                                     "word": r["anchor"]["word"],
                                     "half_width":
                                         r["anchor"]["half_width"]})
    ok = fail = 0
    nfull = 0
    failures = []
    out = open(args.out, "w")
    for k, r in enumerate(recs):
        mx, my = Q(r["mid"][0]), Q(r["mid"][1])
        w = tuple(int(c) for c in r["word"])
        good, msg = verify_cover.verify_at(w, (mx, my))
        entry = {"mid": r["mid"], "word": r["word"], "mid_ok": good,
                 "msg": msg}
        if good and args.full and k % args.full == 0:
            h = Q(r["half_width"])
            box = (Iv(mx - h, mx + h), Iv(my - h, my + h))
            entry["full_ok"] = all(
                verify_cover.verify_at(w, a)[0]
                for a in verify_cover.samples(box, 3))
            nfull += 1
            good = good and entry["full_ok"]
        ok += good
        fail += not good
        if not good:
            failures.append(entry)
        out.write(json.dumps(entry) + "\n")
    out.close()
    print(f"simcheck: {ok} confirmed, {fail} contradicted "
          f"out of {len(recs)} certificates ({nfull} with full 9-sample box "
          f"check)")
    for f_ in failures[:20]:
        print("  FAIL", f_)
    sys.exit(1 if fail else 0)


# ----------------------------------------------------------------------
# summarize
# ----------------------------------------------------------------------

def cmd_summarize(args):
    cells = [json.loads(line) for line in open(args.map) if line.strip()]
    anchors = {}
    if args.anchor and os.path.exists(args.anchor):
        for line in open(args.anchor):
            if line.strip():
                r = json.loads(line)
                anchors[(r["i"], r["j"])] = r

    region_area = math.pi / 16

    def area(rec):
        x0, x1, y0, y1 = (Q(v) for v in rec["box"])
        return float((x1 - x0) * (y1 - y0))

    interior = [r for r in cells if r["status"] != "BOUNDARY"]
    interior_area = sum(area(r) for r in interior)
    boundary_area_est = region_area - interior_area

    # candidate coverage and anchored coverage by max word length,
    # area-weighted over interior cells
    lengths = sorted({r["cand_len"] for r in interior if r.get("cand_len")})
    rows = []
    for L in lengths:
        ca = sum(area(r) for r in interior
                 if r.get("cand_len") and r["cand_len"] <= L)
        aa = sum(area(r) for r in interior
                 if anchors.get((r["i"], r["j"]), {}).get("status")
                 == "ANCHORED"
                 and anchors[(r["i"], r["j"])]["len"] <= L)
        rows.append({"max_len": L,
                     "candidate_cells_area": ca,
                     "candidate_fraction_of_interior": ca / interior_area,
                     "anchored_cells_area": aa,
                     "anchored_fraction_of_interior": aa / interior_area})

    # exactly-certified area: the union of anchor boxes (disjoint: one per
    # cell, half-width <= 1/128 = half the cell size)
    cert_area = sum(float((2 * Q(a["half_width"])) ** 2)
                    for a in anchors.values()
                    if a.get("status") == "ANCHORED")

    bands = [(90, 91), (91, 92), (92, 95), (95, 100), (100, 105),
             (105, 112.3), (112.3, 120), (120, 135), (135, 150), (150, 180)]
    band_rows = []
    for glo, ghi in bands:
        tot = cand = anch = 0.0
        lens = []
        for r in interior:
            gmid = r.get("gamma_mid")
            if gmid is None or not (glo <= gmid < ghi):
                continue
            a = area(r)
            tot += a
            if r.get("cand_len"):
                cand += a
                lens.append(r["cand_len"])
            ar = anchors.get((r["i"], r["j"]))
            if ar and ar.get("status") == "ANCHORED":
                anch += a
        band_rows.append(
            {"gamma_band": [glo, ghi], "interior_area": tot,
             "candidate_fraction": (cand / tot) if tot else None,
             "anchored_fraction": (anch / tot) if tot else None,
             "cand_len_min": min(lens) if lens else None,
             "cand_len_max": max(lens) if lens else None})

    nocand = [{"box": r["box"], "gamma": r.get("gamma"),
               "gamma_mid": r.get("gamma_mid")}
              for r in interior if r["status"] == "NOCAND"]
    nocand.sort(key=lambda r: r["gamma_mid"] or 999)
    cand_unanchored = [
        {"box": r["box"], "gamma_mid": r.get("gamma_mid"),
         "cand_len": r.get("cand_len")}
        for r in interior
        if r["status"] == "CAND"
        and anchors.get((r["i"], r["j"]), {}).get("status") != "ANCHORED"]

    summary = {"region": "obtuse-at-C half-disk, x<=1/2, area pi/16",
               "grid": cells[0]["grid"] if cells else None,
               "interior_area": interior_area,
               "boundary_unresolved_area": boundary_area_est,
               "n_interior_cells": len(interior),
               "coverage_by_max_len": rows,
               "exactly_certified_area_anchors": cert_area,
               "gamma_bands": band_rows,
               "n_nocand_cells": len(nocand),
               "n_cand_unanchored_cells": len(cand_unanchored),
               "nocand_cells": nocand,
               "cand_unanchored_cells": cand_unanchored}
    with open(args.out, "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps({"coverage_by_max_len": rows[-1:],
                      "gamma_bands": band_rows}, indent=1))
    print("n_nocand:", len(nocand), " n_cand_unanchored:",
          len(cand_unanchored))


def cmd_ascii(args):
    """Human-readable rendering of the layer-1 map."""
    recs = [json.loads(l) for l in open(args.map) if l.strip()]
    anchors = {}
    if args.anchor and os.path.exists(args.anchor):
        for line in open(args.anchor):
            if line.strip():
                r = json.loads(line)
                anchors[(r["i"], r["j"])] = r
    grid = {(r["i"], r["j"]): r for r in recs}
    g = recs[0]["grid"]
    sym = {14: "d", 16: "e", 18: "f", 20: "g", 22: "h", 24: "i", 26: "j",
           6: "a", 10: "b", 12: "c"}
    lines = []
    for j in range(g - 1, -1, -1):
        row = []
        for i in range(g):
            r = grid.get((i, j))
            if r is None:
                row.append(" ")
            elif r["status"] == "BOUNDARY":
                row.append(".")
            elif r["cand_len"] is None:
                row.append("#")
            else:
                c = sym[r["cand_len"]]
                a = anchors.get((i, j), {})
                # anchored cells upper-case
                row.append(c.upper() if a.get("status") == "ANCHORED"
                           else c)
        lines.append(f"{j:2d} {''.join(row)}")
    lines.append("")
    lines.append("x: i/64..(i+1)/64, y: j/64..(j+1)/64;"
                 " apex parameter region x<=1/2 of the obtuse half-disk")
    lines.append("letter = minimal candidate word length at cell midpoint:"
                 " d=14 e=16 f=18 g=20 h=22 i=24 j=26")
    lines.append("UPPERCASE = midpoint neighbourhood carries an exact"
                 " harness certificate (anchor)")
    lines.append("# = no candidate word of length <= 26; . = cell straddles"
                 " the region boundary")
    text = "\n".join(lines)
    with open(args.out, "w") as f:
        f.write(text + "\n")
    print(text)


# ----------------------------------------------------------------------
# selftest
# ----------------------------------------------------------------------

def cmd_selftest(args):
    ok = True

    def check(label, got, want):
        nonlocal ok
        good = got == want
        ok &= good
        print(f"[{label}] {got} (expect {want}) {'OK' if good else 'FAIL'}")

    # 1. incremental translation bookkeeping == harness, all words len <= 10
    import itertools
    mism = 0
    mine = set()
    enumerate_translation_words(10, lambda w: mine.add(w))
    theirs = set()
    for n in range(2, 11):
        for w in itertools.product((0, 1, 2), repeat=n):
            if unfold.is_word(w) and unfold.is_translation_word(w):
                c = canonical(w)
                if not _primitive_root_is_translation(c):
                    theirs.add(c)
    check("canonical translation words agree with harness (len<=10)",
          mine == theirs, True)
    if mine != theirs:
        print("   only-mine:", sorted(mine - theirs)[:5],
              "only-harness:", sorted(theirs - mine)[:5])

    # 2. every canonical rep is itself a valid translation word per harness
    bad = [w for w in mine if not (unfold.is_word(w)
                                   and unfold.is_translation_word(w))]
    check("all reps valid per harness", bad, [])

    # 3. rotation/reversal corridor invariance.  The *true* corridor is the
    # same set for every rotation and reversal, so on a box where one variant
    # is TRUE no variant may be FALSE (enclosure slack can only demote a
    # verdict to UNKNOWN, never flip it).  At a zero-width box (a point) the
    # enclosure is exact, so all variants must agree outright.
    fag = (0, 1, 2, 0, 1, 2)
    acute = (Iv(Q(49, 100), Q(51, 100)), Iv(Q(79, 100), Q(81, 100)))
    verdicts = set()
    for base in (fag, fag[::-1]):
        for r in range(6):
            w = base[r:] + base[:r]
            verdicts.add(unfold.certify(w, acute)[0])
    check("no contradictory verdict among fagnano variants (box)",
          "FALSE" in verdicts or TRUE not in verdicts, False)
    pt = (Iv(Q(1, 2)), Iv(Q(4, 5)))
    verdicts_pt = set()
    for base in (fag, fag[::-1]):
        for r in range(6):
            w = base[r:] + base[:r]
            verdicts_pt.add(unfold.certify(w, pt)[0])
    check("fagnano variants all agree at a point", verdicts_pt, {TRUE})

    # 4. float corridor sign agrees with exact point certification on a
    # sample of words and apexes
    words = sorted(mine)
    apexes = [(Q(2, 5), Q(1, 5)), (Q(1, 4), Q(3, 10)), (Q(9, 20), Q(2, 5)),
              (Q(1, 2), Q(2, 5)), (Q(3, 10), Q(1, 10))]
    dis = 0
    tested = 0
    for w in words:
        for ax, ay in apexes:
            fw = float_width(w, float(ax), float(ay))
            if abs(fw) < 1e-9:
                continue  # too close to call in floats; screen is agnostic
            v, _ = unfold.certify(w, (Iv(ax), Iv(ay)))
            tested += 1
            if (fw > 0) != (v == TRUE):
                dis += 1
    check(f"float screen agrees with exact point verdicts ({tested} cases)",
          dis, 0)

    # 5. region predicates: exact circle tests
    check("box inside region", box_in_region(Q(1, 4), Q(5, 16), Q(1, 8),
                                             Q(3, 16)), True)
    check("box crossing arc not inside",
          box_in_region(Q(1, 4), Q(1, 2), Q(1, 4), Q(1, 2)), False)
    check("box outside disk does not meet region",
          box_meets_region(Q(0), Q(1, 16), Q(7, 16), Q(1, 2)), False)

    # 6. anchor pipeline end-to-end at one apex: finds an exact certificate,
    # the certificate survives shrinking (harness contract), and the direct
    # simulator confirms the orbit at the anchor centre
    upto14 = set()
    enumerate_translation_words(14, lambda w: upto14.add(w))
    wl = {n: sorted(w for w in upto14 if len(w) == n)
          for n in (6, 10, 12, 14)}
    wl = {n: v for n, v in wl.items() if v}
    mx, my = Q(7, 16) + Q(1, 128), Q(3, 8) + Q(1, 128)
    res = anchor_at(wl, mx, my)
    good = res is not None
    ok &= good
    print(f"[anchor at (57/128, 49/128)] {res} {'OK' if good else 'FAIL'}")
    if res:
        w = tuple(int(c) for c in res["word"])
        h = Q(res["half_width"]) / 2
        v = unfold.certify(w, (Iv(mx - h, mx + h), Iv(my - h, my + h)))[0]
        check("anchor survives halving the box", v, TRUE)
        simok, msg = verify_cover.verify_at(w, (mx, my))
        check("simulator confirms anchor orbit", simok, True)

    # 7. gamma sanity: right-angle arc gives 90, thin gives ~180
    g = gamma_deg(Q(1, 4), math.sqrt(3) / 4)
    ok &= abs(g - 90.0) < 1e-9
    print(f"[gamma on arc] {g:.9f} (expect 90) {'OK' if abs(g-90)<1e-9 else 'FAIL'}")
    g2 = gamma_deg(Q(1, 2), Q(1, 1000))
    ok &= g2 > 179.5
    print(f"[gamma thin] {g2:.3f} (expect >179.5) "
          f"{'OK' if g2 > 179.5 else 'FAIL'}")

    print("SELFTEST", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


# ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("selftest")
    s.set_defaults(func=cmd_selftest)

    s = sub.add_parser("words")
    s.add_argument("--max-len", type=int, required=True)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_words)

    s = sub.add_parser("map")
    s.add_argument("--words", required=True)
    s.add_argument("--max-len", type=int, required=True)
    s.add_argument("--grid", type=int, default=32,
                   help="cells per half-unit; cell size = 1/(2*grid)")
    s.add_argument("--row-lo", type=int, default=0)
    s.add_argument("--row-hi", type=int, default=10 ** 9)
    s.add_argument("--top", type=int, default=4)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_map)

    s = sub.add_parser("anchor")
    s.add_argument("--words", required=True)
    s.add_argument("--max-len", type=int, required=True)
    s.add_argument("--src", required=True)
    s.add_argument("--row-lo", type=int, default=0)
    s.add_argument("--row-hi", type=int, default=10 ** 9)
    s.add_argument("--exact-cap", type=int, default=4)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_anchor)

    s = sub.add_parser("growth")
    s.add_argument("--words", required=True)
    s.add_argument("--max-len", type=int, required=True)
    s.add_argument("--x", required=True, help="fixed x as a fraction")
    s.add_argument("--denom", type=int, default=1024,
                   help="box half-width = 1/denom, centre snapped")
    s.add_argument("--gamma-lo", type=float, required=True)
    s.add_argument("--gamma-hi", type=float, required=True)
    s.add_argument("--gamma-step", type=float, default=1.0)
    s.add_argument("--exact-cap", type=int, default=6)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_growth)

    s = sub.add_parser("arcsweep")
    s.add_argument("--words", required=True)
    s.add_argument("--max-len", type=int, required=True)
    s.add_argument("--points", type=int, default=16)
    s.add_argument("--x-lo", type=float, default=0.02)
    s.add_argument("--gamma-lo", type=float, required=True)
    s.add_argument("--gamma-hi", type=float, required=True)
    s.add_argument("--gamma-step", type=float, default=1.0)
    s.add_argument("--exact-cap", type=int, default=4)
    s.add_argument("--anchor-best", action="store_true")
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_arcsweep)

    s = sub.add_parser("family")
    s.add_argument("--a", type=int, required=True)
    s.add_argument("--b", type=int, required=True)
    s.add_argument("--points", type=int, default=400)
    s.add_argument("--step", type=float, default=0.5)
    s.add_argument("--iters", type=int, default=16)
    s.add_argument("--x-lo", type=float, default=0.02)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_family)

    s = sub.add_parser("bisect_death")
    s.add_argument("--words", required=True)
    s.add_argument("--len", type=int, required=True)
    s.add_argument("--gamma-lo", type=float, required=True)
    s.add_argument("--gamma-hi", type=float, required=True)
    s.add_argument("--points", type=int, default=400)
    s.add_argument("--iters", type=int, default=14)
    s.add_argument("--margin", type=float, default=0.25)
    s.add_argument("--x-lo", type=float, default=0.02)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_bisect_death)

    s = sub.add_parser("pointsweep")
    s.add_argument("--words", required=True)
    s.add_argument("--len", type=int, required=True)
    s.add_argument("--gamma", type=float, required=True)
    s.add_argument("--points", type=int, default=40)
    s.add_argument("--x-lo", type=float, default=0.02)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_pointsweep)

    s = sub.add_parser("simcheck")
    s.add_argument("--src", nargs="+", required=True)
    s.add_argument("--full", type=int, default=0,
                   help="also 9-sample-check every N-th certificate")
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_simcheck)

    s = sub.add_parser("ascii")
    s.add_argument("--map", required=True)
    s.add_argument("--anchor", default=None)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_ascii)

    s = sub.add_parser("summarize")
    s.add_argument("--map", required=True)
    s.add_argument("--anchor", default=None)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_summarize)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
