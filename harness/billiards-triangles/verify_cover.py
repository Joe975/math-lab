#!/usr/bin/env python3
"""Adversarial re-verification of unfolding certificates, by direct simulation.

`unfold.py` decides whether a word certifies a region by folding the triangle
out flat and intersecting projected gates.  That argument has a documented soft
spot: a line meeting every gate is the standard criterion, but the criterion is
about the *set* of gates, and it does not by itself confirm that a trajectory
crosses them in the word's order.  Re-running the same unfolding would not
notice.  So this tool does not unfold at all.

Instead it bounces the ball.  At a rational apex every quantity in a billiard
trajectory is rational -- the ray/segment intersections and the direction
reflections are linear solves over Q -- so a trajectory can be simulated with
no enclosures and no error terms at all, and the answer is not an
approximation.  For each sampled triangle inside a claimed region the tool
checks, exactly:

  * the sequence of sides struck is the claimed word, cyclically,
  * after that many bounces the trajectory is back at its starting point with
    its starting direction, so the orbit really closes,
  * no bounce lands on a vertex, where the billiard is undefined,
  * and the trajectory's own period is measured and reported.

That last one is reported rather than judged, for a reason worth recording: a
trajectory whose return map reverses orientation closes in half the letters of
its word, because the unfolding only closes after going round twice.  Fagnano's
orbit is exactly that -- three bounces, six letters -- so a shorter period is
ordinary, not a defect.

What this establishes and what it does not.  A region certificate is a claim
about infinitely many triangles, and finitely many exact samples cannot confirm
it -- that is the certificate's job.  What the samples do is independent: they
catch a certificate that is vacuous, mis-derived, or right about the gates and
wrong about their order, none of which the unfolding can catch about itself.  A
claim that passes here is a claim whose certificate has been attacked and
survived, not a claim that has been proved twice.

The one thing taken from `unfold.py` is the starting point: the corridor tells
this tool where to aim.  Everything that happens after that is computed here.

Usage:
  verify_cover.py --selftest                 run the validation suite
  verify_cover.py --check WORD X0,X1,Y0,Y1   re-verify a claimed certificate
  verify_cover.py --at WORD X,Y              simulate at one exact apex
"""

import argparse
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unfold
from unfold import Iv, Q, SIDE_ENDS, TRUE, is_word, parse_box


# ----------------------------------------------------------------------
# exact rational geometry, independent of the unfolding machinery
# ----------------------------------------------------------------------

def corners(apex):
    """The triangle at an exact rational apex: A = (0,0), B = (1,0), C."""
    return ((Q(0), Q(0)), (Q(1), Q(0)), (Q(apex[0]), Q(apex[1])))


def side_points(tri, s):
    i, j = SIDE_ENDS[s]
    return tri[i], tri[j]


def reflect_direction(d, tri, s):
    """Reflect a direction vector in the line of a side.  A direction has no
    position, so only the linear part of the reflection applies."""
    u, v = side_points(tri, s)
    e = (v[0] - u[0], v[1] - u[1])
    n2 = e[0] * e[0] + e[1] * e[1]
    f = 2 * (d[0] * e[0] + d[1] * e[1]) / n2
    return (f * e[0] - d[0], f * e[1] - d[1])


def next_hit(p, d, tri, skip):
    """The first side struck from p in direction d, ignoring the side just
    left.  Returns (side, point, position along that side) or None.

    Solving p + t*d = u + s*(v-u) is a 2x2 rational system, so "which side is
    struck" is decided exactly rather than by a tolerance.
    """
    best = None
    for s in (0, 1, 2):
        if s == skip:
            continue
        u, v = side_points(tri, s)
        e = (v[0] - u[0], v[1] - u[1])
        det = d[0] * (-e[1]) - d[1] * (-e[0])
        if det == 0:
            continue
        rx, ry = u[0] - p[0], u[1] - p[1]
        t = (rx * (-e[1]) - ry * (-e[0])) / det
        w = (d[0] * ry - d[1] * rx) / det
        if t <= 0 or not (0 < w < 1):
            continue
        if best is None or t < best[0]:
            best = (t, s, (p[0] + t * d[0], p[1] + t * d[1]), w)
    if best is None:
        return None
    _, s, point, w = best
    return s, point, w


def simulate(tri, p, d, steps, skip):
    """Bounce a ball, returning (sides struck, end point, end direction) or a
    failure string."""
    struck = []
    for _ in range(steps):
        hit = next_hit(p, d, tri, skip)
        if hit is None:
            return "trajectory left the triangle or struck a vertex"
        s, p, _ = hit
        d = reflect_direction(d, tri, s)
        struck.append(s)
        skip = s
    return struck, p, d


def launch(word, apex):
    """A starting point and outgoing direction for a word at an exact apex.

    Taken from the unfolding: the corridor says which lines cross every gate,
    and the trajectory is launched from where the middle such line meets the
    first gate.  Because the ball leaves from a point *on* the first gate, the
    sides it then strikes are the word rotated by one.
    """
    box = (Iv(apex[0]), Iv(apex[1]))
    verdict, detail = unfold.certify(word, box)
    if verdict != TRUE:
        return None, f"no certificate at this apex ({verdict}: {detail})"
    lower, upper = detail
    copies, gates = unfold.unfold(word, box)
    tau = (unfold.iv(copies[-1][0][0] - copies[0][0][0]).lo,
           unfold.iv(copies[-1][0][1] - copies[0][0][1]).lo)
    normal = (-tau[1], tau[0])
    u, v = gates[0]
    u = (unfold.iv(u[0]).lo, unfold.iv(u[1]).lo)
    v = (unfold.iv(v[0]).lo, unfold.iv(v[1]).lo)
    x = (lower + upper) / 2
    edge = (v[0] - u[0], v[1] - u[1])
    denom = edge[0] * normal[0] + edge[1] * normal[1]
    if denom == 0:
        return None, "first gate is parallel to the corridor"
    s = (x - (u[0] * normal[0] + u[1] * normal[1])) / denom
    if not (0 < s < 1):
        return None, "corridor does not meet the first gate strictly inside"
    start = (u[0] + s * edge[0], u[1] + s * edge[1])
    tri = corners(apex)
    return (start, reflect_direction(tau, tri, word[0])), None


def verify_at(word, apex):
    """Exact verification at one rational apex.  Returns (ok, message)."""
    if not is_word(word):
        return False, "not a bounce word"
    tri = corners(apex)
    if tri[2][1] <= 0:
        return False, "degenerate triangle"
    launched, err = launch(word, apex)
    if launched is None:
        return False, err
    start, direction = launched
    result = simulate(tri, start, direction, len(word), word[0])
    if isinstance(result, str):
        return False, result
    struck, end, end_dir = result

    expected = list(word[1:]) + [word[0]]
    if struck != expected:
        return False, f"struck {struck}, word predicts {expected}"
    if end != start:
        return False, "trajectory did not return to its starting point"
    if end_dir != direction:
        return False, "trajectory returned with a different direction"

    # The trajectory's own period can be shorter than the word, and that is
    # not an error.  An orbit whose return map reverses orientation closes as
    # a trajectory while its unfolding only closes after going round twice --
    # Fagnano's orbit is exactly this, a 3-bounce path with a 6-letter word.
    # So the period is reported, not judged.
    period = len(word)
    for k in range(1, len(word)):
        partial = simulate(tri, start, direction, k, word[0])
        if not isinstance(partial, str) and partial[1] == start \
                and partial[2] == direction:
            period = k
            break
    return True, (f"closed after {period} bounces "
                  f"(word length {len(word)}) at {tuple(map(str, start))}")


# ----------------------------------------------------------------------
# region checking
# ----------------------------------------------------------------------

def samples(box, n=3):
    """Rational sample apexes inside a box: the corners pulled in, and an
    interior grid.  Corners are pulled in because the certificate is about the
    closed box while an orbit at an exact corner may be degenerate."""
    x, y = box
    out = []
    for i in range(n):
        for j in range(n):
            fx = Q(i + 1, n + 1)
            fy = Q(j + 1, n + 1)
            out.append((x.lo + fx * (x.hi - x.lo), y.lo + fy * (y.hi - y.lo)))
    return out


def verify_region(word, box, n=3):
    """Re-verify a claimed certificate at sampled triangles inside it."""
    verdict, detail = unfold.certify(word, box)
    print(f"claim: word {''.join(map(str, word))} certifies "
          f"x in {box[0]}, y in {box[1]}")
    if isinstance(detail, tuple):
        detail = f"corridor ({float(detail[0]):.6g}, {float(detail[1]):.6g})"
    print(f"  unfolding says {verdict}: {detail}")
    if verdict != TRUE:
        print("  nothing to verify: this is not a certificate")
        return False
    ok = True
    for apex in samples(box, n):
        good, message = verify_at(word, apex)
        ok &= good
        mark = "OK  " if good else "FAIL"
        print(f"  {mark} apex ({apex[0]}, {apex[1]}): {message}")
    print(f"  {'CONFIRMED' if ok else 'CONTRADICTED'} at {n * n} sampled "
          f"triangles; the region claim itself rests on the certificate")
    return ok


# ----------------------------------------------------------------------

def selftest():
    ok = True

    def check(label, got, want):
        nonlocal ok
        good = got == want
        ok &= good
        print(f"[{label}] {got} (expect {want}) {'OK' if good else 'FAIL'}")

    acute = (Iv(Q(49, 100), Q(51, 100)), Iv(Q(79, 100), Q(81, 100)))
    print("[fagnano on the acute box]")
    ok &= verify_region(unfold.FAGNANO, acute)

    # The known orbit, checked against its known geometry: in an acute
    # triangle Fagnano's orbit is the orthic triangle, so each bounce point is
    # the foot of an altitude.  This is an answer from the literature, not
    # from either of our two implementations.
    apex = (Q(1, 2), Q(4, 5))
    tri = corners(apex)
    launched, err = launch(unfold.FAGNANO, apex)
    good = launched is not None
    ok &= good
    print(f"[launch at (1/2, 4/5)] {'OK' if good else 'FAIL: ' + str(err)}")
    if launched:
        start, direction = launched
        struck, end, end_dir = simulate(tri, start, direction, 3,
                                        unfold.FAGNANO[0])
        # after three bounces the orbit is at the antipodal phase, and the
        # first bounce point must be the foot of an altitude
        first = next_hit(start, direction, tri, unfold.FAGNANO[0])
        s, point, _ = first
        u, v = side_points(tri, s)
        edge = (v[0] - u[0], v[1] - u[1])
        opp = tri[s]
        to_opp = (opp[0] - point[0], opp[1] - point[1])
        perpendicular = edge[0] * to_opp[0] + edge[1] * to_opp[1] == 0
        ok &= perpendicular
        print(f"[orthic] first bounce at {tuple(map(str, point))} is the foot "
              f"of the altitude from the opposite vertex: {perpendicular} "
              f"{'OK' if perpendicular else 'FAIL'}")

    # Negative controls.  A verifier that cannot fail is not a verifier.
    check("wrong word rejected", verify_at((0, 1, 0, 2, 1, 2), apex)[0], False)
    check("word with a repeat rejected", verify_at((0, 0, 1), apex)[0], False)
    check("obtuse apex has no certificate",
          verify_at(unfold.FAGNANO, (Q(1, 2), Q(1, 10)))[0], False)

    # A verifier that cannot fail is not a verifier.  Aim the same word from
    # a point outside the corridor and the struck sequence must stop matching:
    # if it still matched, the check would be confirming the word rather than
    # the trajectory.
    launched, _ = launch(unfold.FAGNANO, apex)
    start, direction = launched
    u, v = side_points(tri, unfold.FAGNANO[0])
    off = (u[0] + Q(1, 50) * (v[0] - u[0]), u[1] + Q(1, 50) * (v[1] - u[1]))
    result = simulate(tri, off, direction, len(unfold.FAGNANO),
                      unfold.FAGNANO[0])
    expected = list(unfold.FAGNANO[1:]) + [unfold.FAGNANO[0]]
    diverged = isinstance(result, str) or result[0] != expected \
        or result[1] != off
    ok &= diverged
    print(f"[off-corridor launch diverges] "
          f"{result if isinstance(result, str) else result[0]} "
          f"{'OK' if diverged else 'FAIL'}")

    # A doubled word is a valid but redundant certificate: same trajectory,
    # twice the letters.  It must be accepted, and its true period reported.
    doubled = unfold.FAGNANO * 2
    good, message = verify_at(doubled, apex)
    reports_period = good and "closed after 3 bounces" in message
    ok &= reports_period
    print(f"[doubled word] {message} {'OK' if reports_period else 'FAIL'}")

    print("SELFTEST", "PASSED" if ok else "FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--check", nargs=2, metavar=("WORD", "X0,X1,Y0,Y1"))
    ap.add_argument("--at", nargs=2, metavar=("WORD", "X,Y"))
    ap.add_argument("--samples", type=int, default=3)
    args = ap.parse_args()

    if args.selftest:
        sys.exit(0 if selftest() else 1)
    if args.check:
        word = tuple(int(c) for c in args.check[0])
        sys.exit(0 if verify_region(word, parse_box(args.check[1]),
                                    args.samples) else 1)
    if args.at:
        word = tuple(int(c) for c in args.at[0])
        x, y = (Fraction(v) for v in args.at[1].split(","))
        good, message = verify_at(word, (x, y))
        print(f"{'CONFIRMED' if good else 'REJECTED'}: {message}")
        sys.exit(0 if good else 1)
    ap.print_help()


if __name__ == "__main__":
    main()
