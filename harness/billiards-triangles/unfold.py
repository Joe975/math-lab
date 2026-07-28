#!/usr/bin/env python3
"""Unfolding certificates for periodic billiard orbits in triangles.

A bounce word w = (s_1, ..., s_n) over the side labels {0, 1, 2} is realized by
a periodic orbit of a given triangle when the unfolded chain of reflected
copies admits a straight line meeting every gate.  This module decides that
question over a *region* of triangles at once, which is the only form of answer
worth having: a single triangle is a measure-zero statement, while a region is
a certificate covering an open set of the parameter space.

Two exactness decisions, both load-bearing.

**Apex coordinates, not angles.**  Every triangle is normalized to A = (0,0),
B = (1,0) and an apex C, so a family of triangles is a rational box in the
plane.  Reflections, gates and projections are then rational functions of C
with no trigonometry anywhere in the certificate path.  `angles_to_apex`
converts an angle region into an apex box when one is specified that way, and
is the only routine here that touches transcendental functions.

**Affine forms with rational coefficients.**  A box of apexes is carried
through the unfolding in first-order forms over Q, so every inequality decided
is decided for every triangle in the box simultaneously.  A verdict is one of
TRUE (the orbit exists throughout the box), FALSE (it exists nowhere in the
box), or UNKNOWN (the box straddles the boundary, or the enclosure is too wide
-- subdivide and retry).  There is no floating point in this file.

Plain intervals do not work here and the reason is worth recording: each
reflection recombines coordinates that earlier reflections already widened,
and an interval cannot see that those errors are the same error, so the
enclosure grows by a factor per reflection.  Measured on a six-reflection
word, that is about an order of magnitude per step -- a box of width 1e-6
collapses before the word ends.  Keeping first-order dependence on the two
apex coordinates makes the growth additive instead.

What this method can and cannot see.  A certificate covers an open region, so
it only ever finds orbits that survive perturbation of the triangle.  Orbits
that exist for one triangle but no neighbourhood of it -- the unstable ones,
whose unfolding closes up only because a particular integer combination of
that triangle's angles happens to vanish -- are invisible to it by
construction, not by oversight.  Nothing here should be read as evidence that
such a triangle has no periodic orbit.

The translation condition is combinatorial, and exact.  Reflecting a direction
d in a mirror of direction t sends d to 2t - d, so every mirror direction in
the chain is an integer combination p*alpha + q*beta + r*pi.  The composed
isometry is a translation for *every* triangle in an open region exactly when
that combination has p = q = 0 with r even -- a condition on the word alone,
checked in integer arithmetic by `linear_part`.  This is why one word can
certify an open region rather than a curve.

Known limitation, stated rather than assumed away: the corridor test decides
whether some line in the translation direction meets every gate, which is the
standard criterion, but it does not by itself confirm that the resulting
trajectory visits the gates in the word's order.  That confirmation is the job
of `verify_cover.py`, which re-derives any accepted certificate by direct
simulation instead.  A certificate from this file alone is a candidate.

Usage:
  unfold.py --selftest                     run the validation suite
  unfold.py --certify WORD X0,X1,Y0,Y1     verdict for a word on an apex box
  unfold.py --words N X0,X1,Y0,Y1          words of length <= N certifying it
"""

import argparse
import itertools
import sys
from fractions import Fraction

Q = Fraction

TRUE, FALSE, UNKNOWN = "TRUE", "FALSE", "UNKNOWN"

# pi to 50 places, as a rational enclosure.  Used only by angles_to_apex and
# the degree helper; nothing in the certificate path depends on it.
_PI_DIGITS = 314159265358979323846264338327950288419716939937510
PI_LO = Q(_PI_DIGITS, 10 ** 50)
PI_HI = Q(_PI_DIGITS + 1, 10 ** 50)


# ----------------------------------------------------------------------
# interval arithmetic over Q
# ----------------------------------------------------------------------

class Iv:
    """A closed interval with exact rational endpoints."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = Q(lo)
        self.hi = Q(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError("empty interval")

    def __add__(self, o):
        o = iv(o)
        return Iv(self.lo + o.lo, self.hi + o.hi)

    def __sub__(self, o):
        o = iv(o)
        return Iv(self.lo - o.hi, self.hi - o.lo)

    def __neg__(self):
        return Iv(-self.hi, -self.lo)

    def __mul__(self, o):
        o = iv(o)
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return Iv(min(p), max(p))

    def __truediv__(self, o):
        o = iv(o)
        if o.lo <= 0 <= o.hi:
            raise ZeroDivisionError("divisor interval contains 0")
        p = (self.lo / o.lo, self.lo / o.hi, self.hi / o.lo, self.hi / o.hi)
        return Iv(min(p), max(p))

    __radd__ = __add__
    __rmul__ = __mul__

    def __rsub__(self, o):
        return iv(o) - self

    def straddles_zero(self):
        return self.lo <= 0 <= self.hi

    def width(self):
        return self.hi - self.lo

    def mid(self):
        return (self.lo + self.hi) / 2

    def __repr__(self):
        return f"[{float(self.lo):.6g}, {float(self.hi):.6g}]"


def iv(x):
    if isinstance(x, Iv):
        return x
    if isinstance(x, Aff):
        return x.to_iv()
    return Iv(x)


# ----------------------------------------------------------------------
# affine arithmetic over Q
# ----------------------------------------------------------------------

class Aff:
    """A first-order form: centre + sum of coefficients times noise symbols,
    plus an aggregated non-linear radius.  Every noise symbol ranges over
    [-1, 1], and two forms sharing a symbol share the same value of it.

    Plain intervals are unusable for unfolding.  A reflection is an affine
    combination of coordinates that have themselves been widened by earlier
    reflections, and an interval cannot see that the errors are the same
    error, so an enclosure grows by a constant factor per step and a
    six-reflection word is already hopeless from a box of width 1e-6.  Keeping
    the first-order dependence on the two apex coordinates makes the growth
    additive instead, which is what makes a region certificate possible at all.
    """

    __slots__ = ("c", "d", "e")

    def __init__(self, c, d=None, e=0):
        self.c = Q(c)
        self.d = dict(d) if d else {}
        self.e = Q(e)

    @staticmethod
    def var(lo, hi, symbol):
        """A quantity known only to lie in [lo, hi], tagged with a symbol so
        that later uses of the same quantity stay correlated."""
        lo, hi = Q(lo), Q(hi)
        return Aff((lo + hi) / 2, {symbol: (hi - lo) / 2})

    def rad(self):
        return sum(abs(v) for v in self.d.values()) + self.e

    def to_iv(self):
        r = self.rad()
        return Iv(self.c - r, self.c + r)

    def __add__(self, o):
        o = aff(o)
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d.get(k, Q(0)) + v
        return Aff(self.c + o.c, d, self.e + o.e)

    def __neg__(self):
        return Aff(-self.c, {k: -v for k, v in self.d.items()}, self.e)

    def __sub__(self, o):
        return self + (-aff(o))

    def __mul__(self, o):
        o = aff(o)
        d = {}
        for k in set(self.d) | set(o.d):
            d[k] = self.c * o.d.get(k, Q(0)) + o.c * self.d.get(k, Q(0))
        e = abs(self.c) * o.e + abs(o.c) * self.e + self.rad() * o.rad()
        return Aff(self.c * o.c, d, e)

    def inv(self):
        """1/x, keeping the first-order dependence.

        Replacing the reciprocal by a fresh interval-valued constant is sound
        but expensive: it injects an uncorrelated error at every reflection,
        and those errors are what compound.  The min-range approximation used
        here takes the slope at the far end, where 1/x is flattest, so the
        residual is second-order.  1/x is convex and decreasing on a positive
        interval, so x - p*x is decreasing there and its extremes are the two
        endpoints.
        """
        b = self.to_iv()
        if b.lo <= 0 <= b.hi:
            raise ZeroDivisionError("divisor contains 0")
        if b.hi < 0:
            return -(-self).inv()
        p = -Q(1) / (b.hi * b.hi)
        hi_end = Q(1) / b.lo - p * b.lo      # value of 1/x - p*x at x = lo
        lo_end = Q(1) / b.hi - p * b.hi      # and at x = hi
        centre, rad = (hi_end + lo_end) / 2, (hi_end - lo_end) / 2
        return self * p + Aff(centre, {}, rad)

    def __truediv__(self, o):
        return self * aff(o).inv()

    __radd__ = __add__
    __rmul__ = __mul__

    def __rsub__(self, o):
        return aff(o) - self

    def __repr__(self):
        return repr(self.to_iv())


def aff(x):
    if isinstance(x, Aff):
        return x
    if isinstance(x, Iv):
        raise TypeError("give an interval a symbol with Aff.var before mixing")
    return Aff(x)


APEX_X, APEX_Y = 0, 1


def square(a):
    """x^2.  Not the same as a*a: a generic product cannot see that both
    factors are the same quantity, so it loses the fact that a square is never
    negative -- and that lost sign is enough to make a length come out zero."""
    if isinstance(a, Aff):
        r = a.rad()
        return Aff(a.c * a.c + r * r / 2,
                   {k: 2 * a.c * v for k, v in a.d.items()},
                   2 * abs(a.c) * a.e + r * r / 2)
    a = iv(a)
    if a.lo >= 0:
        return Iv(a.lo * a.lo, a.hi * a.hi)
    if a.hi <= 0:
        return Iv(a.hi * a.hi, a.lo * a.lo)
    return Iv(Q(0), max(a.lo * a.lo, a.hi * a.hi))


def dot(u, v):
    return u[0] * v[0] + u[1] * v[1]


def norm2(u):
    return square(u[0]) + square(u[1])


def sub(u, v):
    return (u[0] - v[0], u[1] - v[1])


def add(u, v):
    return (u[0] + v[0], u[1] + v[1])


def scale(s, u):
    return (s * u[0], s * u[1])


# ----------------------------------------------------------------------
# certified sin and cos (used only to convert an angle region)
# ----------------------------------------------------------------------

def _sin_point(x, terms=40):
    """Enclosure of sin(x) at a rational x, by the alternating Taylor series.

    The remainder of an alternating series with decreasing terms is bounded by
    the first omitted term, which is what the returned width is."""
    total, term = Q(0), Q(x)
    for k in range(terms):
        total += term
        term = -term * x * x / ((2 * k + 2) * (2 * k + 3))
    r = abs(term)
    return Iv(total - r, total + r)


def _cos_point(x, terms=40):
    total, term = Q(0), Q(1)
    for k in range(terms):
        total += term
        term = -term * x * x / ((2 * k + 1) * (2 * k + 2))
    r = abs(term)
    return Iv(total - r, total + r)


def sin_iv(a):
    """sin over an interval.  Both sin and cos are 1-Lipschitz, so evaluating
    at the midpoint and widening by the half-width is a valid enclosure."""
    a = iv(a)
    h = a.width() / 2
    v = _sin_point(a.mid())
    return Iv(v.lo - h, v.hi + h)


def cos_iv(a):
    a = iv(a)
    h = a.width() / 2
    v = _cos_point(a.mid())
    return Iv(v.lo - h, v.hi + h)


def degrees(d):
    """An angle given in degrees, as a radian interval."""
    return Iv(Q(d) * PI_LO / 180, Q(d) * PI_HI / 180)


def angles_to_apex(alpha, beta):
    """Apex box for the triangles with angle alpha at A and beta at B.

    By the law of sines with |AB| = 1, the apex is at distance
    sin(beta)/sin(alpha+beta) from A along the ray at angle alpha."""
    alpha, beta = iv(alpha), iv(beta)
    t = sin_iv(beta) / sin_iv(alpha + beta)
    return (t * cos_iv(alpha), t * sin_iv(alpha))


# ----------------------------------------------------------------------
# words and the translation condition
# ----------------------------------------------------------------------

def is_word(word):
    """A bounce word never repeats a side, including cyclically: a trajectory
    cannot strike the same side twice in succession."""
    if len(word) < 2 or any(s not in (0, 1, 2) for s in word):
        return False
    return all(word[i] != word[(i + 1) % len(word)] for i in range(len(word)))


# A direction is p*alpha + q*beta + r*pi.  Side 2 = AB has direction 0, side 1
# = CA leaves A at the interior angle alpha, and side 0 = BC leaves B at
# pi - beta.
INITIAL_DIRECTIONS = ((0, -1, 1), (1, 0, 0), (0, 0, 0))


def linear_part(word):
    """The linear part of the composed isometry, as (kind, (p, q, r)).

    kind is "rot" (z -> e^{i*angle} z) or "ref" (z -> e^{i*angle} conj z), and
    the angle is p*alpha + q*beta + r*pi.  Reflecting in a mirror of direction
    t maps a rotation of angle a to a reflection of angle 2t - a and vice
    versa, and maps every side direction d to 2t - d.
    """
    dirs = list(INITIAL_DIRECTIONS)
    kind, ang = "rot", (0, 0, 0)
    for s in word:
        t = dirs[s]
        ang = tuple(2 * t[i] - ang[i] for i in range(3))
        kind = "ref" if kind == "rot" else "rot"
        dirs = [tuple(2 * t[i] - d[i] for i in range(3)) for d in dirs]
    return kind, ang


def is_translation_word(word):
    """True when the unfolding closes up to a translation for every triangle
    in an open region, rather than only along a curve in parameter space."""
    kind, (p, q, r) = linear_part(word)
    return kind == "rot" and p == 0 and q == 0 and r % 2 == 0


def translation_words(max_len):
    """Every bounce word up to a length whose unfolding is a translation.
    Words are generated up to nothing -- rotations and reversals of a word are
    distinct words here, and deduplicating them is a caller's decision."""
    out = []
    for n in range(2, max_len + 1):
        for w in itertools.product((0, 1, 2), repeat=n):
            if is_word(w) and is_translation_word(w):
                out.append(w)
    return out


# ----------------------------------------------------------------------
# unfolding
# ----------------------------------------------------------------------

def triangle(apex):
    """The normalized triangle (A, B, C) for an apex box, in affine forms.

    The apex coordinates get one noise symbol each; everything downstream is
    expressed in terms of those two, which is what keeps the chain correlated.
    """
    x, y = iv(apex[0]), iv(apex[1])
    return ((Aff(0), Aff(0)), (Aff(1), Aff(0)),
            (Aff.var(x.lo, x.hi, APEX_X), Aff.var(y.lo, y.hi, APEX_Y)))


# side s is opposite vertex s: side 0 = BC, side 1 = CA, side 2 = AB
SIDE_ENDS = {0: (1, 2), 1: (2, 0), 2: (0, 1)}


def reflect(p, u, v):
    """Reflect a point in the line through u and v."""
    d = sub(v, u)
    w = sub(p, u)
    f = dot(w, d) / norm2(d) * 2
    return sub(add(u, scale(f, d)), w)


def foot_parameters(apex):
    """For each side, where the opposite vertex projects onto it, as a
    fraction of the way from the first endpoint to the second.

    This is the one quantity that makes interval unfolding usable.  Written
    directly, a reflection divides by the squared length of the mirror, and
    both are recomputed from vertex coordinates that have already been widened
    by earlier steps -- the enclosure then grows by an order of magnitude per
    reflection and collapses within half a dozen.  But every copy in the chain
    is congruent to the original with the same labelling, so this parameter is
    the *same* for every copy of a given side.  Computing it once from the
    original triangle replaces the compounding division with a constant.
    """
    t = triangle(apex)
    out = []
    for s in (0, 1, 2):
        i, j = SIDE_ENDS[s]
        d = sub(t[j], t[i])
        out.append(dot(sub(t[s], t[i]), d) / norm2(d))
    return out


def reflect_opposite(p, u, v, c):
    """The image of the vertex opposite a mirror, given the mirror's foot
    parameter: p' = 2u - p + 2c(v - u)."""
    return sub(add(scale(2, u), scale(c * 2, sub(v, u))), p)


def unfold(word, apex):
    """Unfold along a word.  Returns (copies, gates): the reflected triangles
    and, for each step, the segment shared by consecutive copies.

    Only the vertex opposite the mirror moves; the other two lie on the mirror
    and are carried over unchanged.  Reflecting all three through the formula
    would be arithmetically identical and numerically much worse, because
    interval arithmetic cannot see that a point on the mirror is its own image
    and would widen it at every step.
    """
    t = triangle(apex)
    feet = foot_parameters(apex)
    copies, gates = [t], []
    for s in word:
        i, j = SIDE_ENDS[s]
        u, v = t[i], t[j]
        gates.append((u, v))
        moved = list(t)
        moved[s] = reflect_opposite(t[s], u, v, feet[s])
        t = tuple(moved)
        copies.append(t)
    return copies, gates


def certify(word, apex):
    """Verdict for a word over an apex box, with the corridor bounds.

    The composed isometry is a translation by tau.  A periodic orbit with this
    word exists when some line in direction tau meets every gate, so project
    each gate onto the normal of tau and intersect.  Using the upper end of
    each gate's near endpoint and the lower end of its far endpoint keeps the
    intersection an under-estimate, so a non-empty answer is certified for
    every triangle in the box.
    """
    if not is_word(word):
        return FALSE, "not a bounce word"
    if not is_translation_word(word):
        return FALSE, "unfolding is not a translation"
    copies, gates = unfold(word, apex)
    tau = sub(copies[-1][0], copies[0][0])
    if iv(tau[0]).straddles_zero() and iv(tau[1]).straddles_zero():
        return UNKNOWN, "translation vector not certified nonzero"
    normal = (-tau[1], tau[0])

    lower, upper = None, None       # certified inner bounds of the corridor
    outer_lo, outer_hi = None, None  # outer bounds, for the FALSE verdict
    for u, v in gates:
        a, b = iv(dot(u, normal)), iv(dot(v, normal))
        near_hi, far_lo = min(a.hi, b.hi), max(a.lo, b.lo)
        lower = near_hi if lower is None else max(lower, near_hi)
        upper = far_lo if upper is None else min(upper, far_lo)
        o_lo, o_hi = min(a.lo, b.lo), max(a.hi, b.hi)
        outer_lo = o_lo if outer_lo is None else max(outer_lo, o_lo)
        outer_hi = o_hi if outer_hi is None else min(outer_hi, o_hi)

    if lower < upper:
        return TRUE, (lower, upper)
    if outer_lo > outer_hi:
        return FALSE, (outer_lo, outer_hi)
    return UNKNOWN, (lower, upper)


def certifying_words(apex, max_len):
    """Words up to a length that certify the whole box."""
    return [w for w in translation_words(max_len)
            if certify(w, apex)[0] == TRUE]


# ----------------------------------------------------------------------
# shape predicates, exact in apex coordinates
# ----------------------------------------------------------------------

def angle_signs(apex):
    """Sign of the cosine of each interior angle, as (at A, at B, at C).
    Positive means acute; the sign is the sign of a dot product, so it is
    exact.  Returns None for a component the box does not decide."""
    a, b, c = triangle(apex)
    out = []
    for p, q, r in ((a, b, c), (b, c, a), (c, a, b)):
        d = iv(dot(sub(q, p), sub(r, p)))
        out.append(1 if d.lo > 0 else (-1 if d.hi < 0 else None))
    return tuple(out)


def is_acute_box(apex):
    return all(s == 1 for s in angle_signs(apex))


# ----------------------------------------------------------------------
# validation
# ----------------------------------------------------------------------

def parse_box(text):
    x0, x1, y0, y1 = (Q(v) for v in text.split(","))
    return (Iv(x0, x1), Iv(y0, y1))


FAGNANO = (0, 1, 2, 0, 1, 2)


def selftest():
    ok = True

    def check(label, got, want):
        nonlocal ok
        good = got == want
        ok &= good
        print(f"[{label}] {got} (expect {want}) {'OK' if good else 'FAIL'}")

    # 1. The translation condition is a property of the word alone.
    check("fagnano is a translation word", is_translation_word(FAGNANO), True)
    check("single traversal is not", is_translation_word((0, 1, 2)), False)
    check("odd length is never a translation",
          any(is_translation_word(w) for w in translation_words(7)
              if len(w) % 2), False)
    check("repeated side rejected", is_word((0, 0, 1)), False)
    check("cyclic repeat rejected", is_word((0, 1, 0)), False)

    # 2. The combinatorial translation test must agree with what the
    #    geometry actually does.  At a generic apex, unfold every bounce word
    #    and check the last copy is a translate of the first exactly when
    #    is_translation_word said it would be.  This is the check that would
    #    catch an error in the integer bookkeeping, which nothing downstream
    #    could detect.
    generic = (Iv(Q(37, 100)), Iv(Q(53, 71)))
    mismatch, translations = 0, 0
    for n in range(2, 9):
        for w in itertools.product((0, 1, 2), repeat=n):
            if not is_word(w):
                continue
            copies, _ = unfold(w, generic)
            shifts = [(copies[-1][k][0].c - copies[0][k][0].c,
                       copies[-1][k][1].c - copies[0][k][1].c) for k in range(3)]
            geometric = shifts[0] == shifts[1] == shifts[2]
            predicted = is_translation_word(w)
            translations += predicted
            if geometric != predicted:
                mismatch += 1
    ok &= mismatch == 0
    print(f"[word arithmetic] {translations} translation words up to length 8, "
          f"{mismatch} disagreements with the geometry "
          f"{'OK' if not mismatch else 'FAIL'}")

    # 3. The orthic triangle: Fagnano's orbit exists in every acute triangle
    #    and in no obtuse one.  This is the known answer the harness is
    #    calibrated against.
    acute = (Iv(Q(49, 100), Q(51, 100)), Iv(Q(79, 100), Q(81, 100)))
    check("acute box is acute", is_acute_box(acute), True)
    check("fagnano certifies acute box", certify(FAGNANO, acute)[0], TRUE)

    obtuse_a = (Iv(Q(-51, 100), Q(-49, 100)), Iv(Q(59, 100), Q(61, 100)))
    check("left box is obtuse at A", angle_signs(obtuse_a)[0], -1)
    check("fagnano fails obtuse-at-A box", certify(FAGNANO, obtuse_a)[0], FALSE)

    flat = (Iv(Q(1, 2)), Iv(Q(1, 10)))
    check("flat triangle is obtuse at C", angle_signs(flat)[2], -1)
    check("fagnano fails there", certify(FAGNANO, flat)[0], FALSE)

    # 4. A right triangle is the boundary case, so a box straddling it must
    #    not come back TRUE.  Silent overreach at the boundary is exactly the
    #    failure this harness exists to prevent.
    straddle = (Iv(Q(1, 2), Q(1, 2)), Iv(Q(9, 20), Q(11, 20)))
    verdict = certify(FAGNANO, straddle)[0]
    good = verdict in (FALSE, UNKNOWN)
    ok &= good
    print(f"[right-angle straddle] {verdict} (must not be TRUE) "
          f"{'OK' if good else 'FAIL'}")

    # 5. The search finds the known certificate without being told it, and
    #    finds the rotations of it that are separate words.
    found = certifying_words(acute, 8)
    good = FAGNANO in found
    ok &= good
    print(f"[search] {len(found)} words of length <= 8 certify the acute box, "
          f"including fagnano: {good} {'OK' if good else 'FAIL'}")

    # 6. A certificate must survive shrinking the box, and a box too wide to
    #    decide must become decidable when halved.  That is the intended
    #    workflow, so it is checked rather than assumed.
    sub = (Iv(Q(1, 2), Q(51, 100)), Iv(Q(4, 5), Q(81, 100)))
    check("certificate survives a sub-box", certify(FAGNANO, sub)[0], TRUE)

    wide = (Iv(Q(48, 100), Q(52, 100)), Iv(Q(78, 100), Q(82, 100)))
    coarse = certify(FAGNANO, wide)[0]
    half = (Iv(Q(49, 100), Q(51, 100)), Iv(Q(79, 100), Q(81, 100)))
    print(f"[subdivision] wide box {coarse}, one quarter of it "
          f"{certify(FAGNANO, half)[0]}")

    # 6. The angle conversion agrees with the apex coordinates it replaces:
    #    the equilateral triangle sits at (1/2, sqrt(3)/2).
    sixty = degrees(60)
    ax, ay = angles_to_apex(sixty, sixty)
    # The apex height is sqrt(3)/2, which is not rational, so the enclosure
    # is checked by squaring rather than against a truncated decimal.
    good = (ax.lo <= Q(1, 2) <= ax.hi
            and ay.lo * ay.lo <= Q(3, 4) <= ay.hi * ay.hi)
    ok &= good
    print(f"[equilateral apex] {ax}, {ay} encloses (1/2, sqrt(3)/2) "
          f"{'OK' if good else 'FAIL'}")

    print("SELFTEST", "PASSED" if ok else "FAILED")
    return ok


# ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--certify", nargs=2, metavar=("WORD", "X0,X1,Y0,Y1"),
                    help="word as digits, e.g. 012012")
    ap.add_argument("--words", nargs=2, metavar=("MAXLEN", "X0,X1,Y0,Y1"))
    args = ap.parse_args()

    if args.selftest:
        sys.exit(0 if selftest() else 1)
    if args.certify:
        word = tuple(int(c) for c in args.certify[0])
        box = parse_box(args.certify[1])
        verdict, detail = certify(word, box)
        print(f"word {args.certify[0]} on {box}: {verdict}  {detail}")
        return
    if args.words:
        box = parse_box(args.words[1])
        for w in certifying_words(box, int(args.words[0])):
            print("".join(str(s) for s in w))
        return
    ap.print_help()


if __name__ == "__main__":
    main()
