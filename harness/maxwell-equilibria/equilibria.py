#!/usr/bin/env python3
"""Certified equilibrium counting for point-charge fields (reference tool).

Setting.  Distinct points x_1..x_n in R^3 carry nonzero rational charges
q_1..q_n.  The field of interest is

    F(x) = sum_i q_i * (x - x_i) / |x - x_i|^3        (= E, minus grad of
                                                       V(x) = sum q_i/|x-x_i|)

and an equilibrium is a zero of F away from the charges.  This tool takes a
configuration with rational data and produces a *certificate* for the number
of equilibria in a search region: a binary-subdivision tree of the region in
which every leaf is one of

    excluded-sign     interval evaluation of F over the (closed) leaf box
                      proves some component is bounded away from 0;
    excluded-newton   the Krawczyk image K(B) is disjoint from B, so B
                      contains no zero (any zero lies in K(B), by the
                      mean-value form);
    excluded-ball     the leaf lies inside a punctured ball around charge i
                      on which the domination lemma below proves F != 0;
    excluded-outside  the leaf lies outside the localization sphere of the
                      all-positive lemma below;
    isolated          the Krawczyk operator proves the leaf contains exactly
                      one zero of F, in its interior;
    unresolved        nothing above succeeded at the minimum width (reported,
                      never dropped: the count is then a lower bound).

Leaves are interior-disjoint by construction, so if no leaf is unresolved the
number of equilibria in the region is exactly the number of isolated leaves.

Arithmetic.  Outward-rounded interval arithmetic with dyadic rational
endpoints (fractions.Fraction; products/quotients/roots are rounded outward
to a fixed 2^-PREC_BITS grid, sums and differences are exact).  Square roots
use integer isqrt enclosures:  for f = p/q >= 0 and S = 2^PREC_BITS, with
s = isqrt(p*q*S^2),  (s/(qS))^2 <= f <= ((s+1)/(qS))^2.  Squaring is its own
operation so that squares of intervals straddling 0 stay nonnegative.

Krawczyk certificate (standard; see A. Neumaier, "Interval Methods for
Systems of Equations", CUP 1990, section 5.2, or R. Moore, R. Kearfott,
M. Cloud, "Introduction to Interval Analysis", SIAM 2009, chapter 8).  For a
box B with midpoint m, an arbitrary rational matrix Y, interval Jacobian
enclosure J(B) and interval evaluation F(m):

    K(B) = m - Y*F(m) + (I - Y*J(B)) * (B - m).

If K(B) is contained in the interior of B, then F has exactly one zero in B
(existence by Brouwer via the mean-value extension, uniqueness because the
contraction forces every matrix in J(B) to be nonsingular).  If K(B) and B
are disjoint, F has no zero in B.  Y is a floating-point approximate inverse
of the midpoint Jacobian, rationalized; the *certificate* does not depend on
Y being accurate, only the success rate does.

Localization lemma (all charges positive).  If q_i > 0 for all i then every
equilibrium satisfies |x| < R_0 := max_i |x_i|:  for |x| >= R_0 and x not a
charge point,

    x . F(x) = sum_i q_i * x.(x - x_i) / |x - x_i|^3 > 0,

because x.(x - x_i) >= |x| (|x| - |x_i|) >= 0 with equality in every term
simultaneously only if x equals every x_i.  So a complete census may take
the search region to be any box containing the closed ball of radius R_0.
(Mixed-sign configurations get no automatic region; the caller must supply
one, and completeness claims are then scoped to it.)

Domination lemma (any signs).  Let d_ij = |x_i - x_j| and pick rho_i with
0 < rho_i < min_j d_ij such that

    |q_i| / rho_i^2  >  sum_{j != i} |q_j| / (d_ij - rho_i)^2.

Then F has no zero x with 0 < |x - x_i| <= rho_i:  writing r = |x - x_i|,
the i-th term of F has magnitude |q_i|/r^2 >= |q_i|/rho_i^2 while the others
sum to at most the right-hand side (|x - x_j| >= d_ij - rho_i).  The tool
finds such rho_i by halving from min_j d_ij / 2, instantiating the
inequality in exact rational arithmetic (with certified lower bounds for
d_ij).  Without this lemma, subdivision near a charge would never terminate:
interval evaluations blow up as the singularity is approached.

Index-sum identity (all charges positive; classical Poincare-Hopf/degree
argument).  On a large sphere F points radially outward (localization lemma
computation), so its degree is +1; around each positive charge F is locally
q_i*(x-x_i)/r^3, degree +1.  Hence if all equilibria are nondegenerate,

    sum over equilibria of sign det DF  =  1 - n.

A complete census with certified index signs must satisfy this; the tool
checks it and records the outcome.  A violation means a missed zero or a
wrong sign somewhere -- treat it as an error in the run, never as a finding.

Face-alignment note.  Certified zeros that lie exactly on a subdivision face
can never acquire an isolation certificate (they are interior to no leaf).
Symmetric configurations have zeros at symmetric points, and grids love to
put faces through those points, so the automatic search region is offset by
the non-dyadic rationals (1/97, 1/101, 1/103): every face coordinate is then
offset + dyadic, which no symmetric rational point of small height hits.

Output.  JSON certificate (see --help) with the configuration, lemma
instantiations, full leaf list including split paths (so an independent
checker can rebuild every box exactly and verify the tiling), isolation
boxes, det-signs, and the index-sum check.  Written under $MATHLAB_OUT
(default ./out).

Usage:
  equilibria.py --validate                      run the validation suite
  equilibria.py --config FILE.json [--out F]    count equilibria of a config
      config: {"charges": [{"q": "1", "pos": ["1", "0", "0"]}, ...],
               "region": optional [[lo,hi],[lo,hi],[lo,hi]] as strings}
Options: --min-width-bits K (default 30), --prec-bits P (default 128),
         --max-boxes N (default 400000).

Standard library only.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction
from math import isqrt

# ---------------------------------------------------------------------------
# outward-rounded dyadic interval arithmetic
# ---------------------------------------------------------------------------

PREC_BITS = 128
_SCALE = 1 << PREC_BITS


def set_precision(bits: int) -> None:
    global PREC_BITS, _SCALE
    PREC_BITS = bits
    _SCALE = 1 << bits


def _rdn(x: Fraction) -> Fraction:
    """Round down to the 2^-PREC_BITS grid."""
    return Fraction((x.numerator * _SCALE) // x.denominator, _SCALE)


def _rup(x: Fraction) -> Fraction:
    """Round up to the 2^-PREC_BITS grid."""
    return Fraction(-((-x.numerator * _SCALE) // x.denominator), _SCALE)


def _sqrt_dn(f: Fraction) -> Fraction:
    """Certified lower bound for sqrt(f), f >= 0."""
    if f < 0:
        raise ValueError("sqrt of negative")
    p, q = f.numerator, f.denominator
    s = isqrt(p * q * _SCALE * _SCALE)
    return Fraction(s, q * _SCALE)


def _sqrt_up(f: Fraction) -> Fraction:
    """Certified upper bound for sqrt(f), f >= 0."""
    if f < 0:
        raise ValueError("sqrt of negative")
    p, q = f.numerator, f.denominator
    t = p * q * _SCALE * _SCALE
    s = isqrt(t)
    if s * s == t:
        return Fraction(s, q * _SCALE)
    return Fraction(s + 1, q * _SCALE)


class IV:
    """Closed interval [lo, hi] with Fraction endpoints."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo: Fraction, hi: Fraction | None = None):
        if hi is None:
            hi = lo
        if lo > hi:
            raise ValueError("empty interval")
        self.lo = lo
        self.hi = hi

    # -- exact operations ---------------------------------------------------

    def __add__(self, o: "IV") -> "IV":
        return IV(self.lo + o.lo, self.hi + o.hi)

    def __sub__(self, o: "IV") -> "IV":
        return IV(self.lo - o.hi, self.hi - o.lo)

    def __neg__(self) -> "IV":
        return IV(-self.hi, -self.lo)

    def sub_point(self, p: Fraction) -> "IV":
        return IV(self.lo - p, self.hi - p)

    # -- rounded operations -------------------------------------------------

    def __mul__(self, o: "IV") -> "IV":
        cands = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return IV(_rdn(min(cands)), _rup(max(cands)))

    def scale(self, c: Fraction) -> "IV":
        if c >= 0:
            return IV(_rdn(self.lo * c), _rup(self.hi * c))
        return IV(_rdn(self.hi * c), _rup(self.lo * c))

    def sq(self) -> "IV":
        """Interval square, preserving nonnegativity (never a generic product)."""
        a, b = abs(self.lo), abs(self.hi)
        hi = _rup(max(a, b) ** 2)
        if self.lo <= 0 <= self.hi:
            return IV(Fraction(0), hi)
        return IV(_rdn(min(a, b) ** 2), hi)

    def div(self, o: "IV") -> "IV":
        """self / o, requiring o.lo > 0 (the only case this tool needs)."""
        if o.lo <= 0:
            raise ZeroDivisionError("divisor interval must be strictly positive")
        cands = (self.lo / o.lo, self.lo / o.hi, self.hi / o.lo, self.hi / o.hi)
        return IV(_rdn(min(cands)), _rup(max(cands)))

    def sqrt(self) -> "IV":
        return IV(_sqrt_dn(self.lo), _sqrt_up(self.hi))

    # -- predicates & helpers -----------------------------------------------

    def contains0(self) -> bool:
        return self.lo <= 0 <= self.hi

    def mid(self) -> Fraction:
        return _rdn((self.lo + self.hi) / 2)

    def width(self) -> Fraction:
        return self.hi - self.lo

    def inside_open(self, o: "IV") -> bool:
        return o.lo < self.lo and self.hi < o.hi

    def disjoint(self, o: "IV") -> bool:
        return self.hi < o.lo or o.hi < self.lo


IV0 = IV(Fraction(0))


# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------

class Config:
    """Charges with exact rational data."""

    def __init__(self, charges: list[tuple[Fraction, tuple[Fraction, Fraction, Fraction]]]):
        if len(charges) < 2:
            raise ValueError("need at least two charges")
        seen = set()
        for q, p in charges:
            if q == 0:
                raise ValueError("zero charge")
            if p in seen:
                raise ValueError("coincident charges")
            seen.add(p)
        self.charges = charges

    @property
    def n(self) -> int:
        return len(self.charges)

    def all_positive(self) -> bool:
        return all(q > 0 for q, _ in self.charges)

    @staticmethod
    def from_json(obj: dict) -> "Config":
        return Config([
            (Fraction(c["q"]), tuple(Fraction(v) for v in c["pos"]))
            for c in obj["charges"]
        ])

    def to_json(self) -> dict:
        return {"charges": [{"q": str(q), "pos": [str(v) for v in p]}
                            for q, p in self.charges]}


# ---------------------------------------------------------------------------
# interval field and Jacobian
# ---------------------------------------------------------------------------

def field_iv(cfg: Config, box: tuple[IV, IV, IV]) -> tuple[IV, IV, IV] | None:
    """Interval enclosure of F over box, or None if box may contain a charge."""
    out = [IV0, IV0, IV0]
    for q, p in cfg.charges:
        d = [box[k].sub_point(p[k]) for k in range(3)]
        r2 = d[0].sq() + d[1].sq() + d[2].sq()
        if r2.lo <= 0:
            return None
        r = r2.sqrt()
        r3 = r * r2
        for k in range(3):
            out[k] = out[k] + d[k].scale(q).div(r3)
    return tuple(out)


def jac_iv(cfg: Config, box: tuple[IV, IV, IV]) -> list[list[IV]] | None:
    """Interval enclosure of DF over box (DF_kl = dF_k/dx_l), or None."""
    J = [[IV0, IV0, IV0] for _ in range(3)]
    three = Fraction(3)
    for q, p in cfg.charges:
        d = [box[k].sub_point(p[k]) for k in range(3)]
        r2 = d[0].sq() + d[1].sq() + d[2].sq()
        if r2.lo <= 0:
            return None
        r = r2.sqrt()
        r3 = r * r2
        r5 = r3 * r2
        inv_r3 = IV(Fraction(1)).div(r3)
        for k in range(3):
            for l in range(3):
                if k == l:
                    term = inv_r3 - d[k].sq().div(r5).scale(three)
                else:
                    term = -((d[k] * d[l]).div(r5).scale(three))
                J[k][l] = J[k][l] + term.scale(q)
    return J


def det3_iv(J: list[list[IV]]) -> IV:
    a, b, c = J[0]
    d, e, f = J[1]
    g, h, i = J[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


# ---------------------------------------------------------------------------
# floating-point helpers (candidate machinery only -- never certifies)
# ---------------------------------------------------------------------------

def _field_float(cfg: Config, x: tuple[float, float, float]) -> list[float] | None:
    out = [0.0, 0.0, 0.0]
    for q, p in cfg.charges:
        d = [x[k] - float(p[k]) for k in range(3)]
        r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
        if r2 <= 0.0:
            return None
        r3 = r2 ** 1.5
        for k in range(3):
            out[k] += float(q) * d[k] / r3
    return out


def _jac_float(cfg: Config, x: tuple[float, float, float]) -> list[list[float]] | None:
    J = [[0.0] * 3 for _ in range(3)]
    for q, p in cfg.charges:
        d = [x[k] - float(p[k]) for k in range(3)]
        r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
        if r2 <= 0.0:
            return None
        r3 = r2 ** 1.5
        r5 = r3 * r2
        for k in range(3):
            for l in range(3):
                t = (1.0 / r3 if k == l else 0.0) - 3.0 * d[k] * d[l] / r5
                J[k][l] += float(q) * t
    return J


def _inv3_float(J: list[list[float]]) -> list[list[Fraction]] | None:
    a, b, c = J[0]
    d, e, f = J[1]
    g, h, i = J[2]
    A = e * i - f * h
    B = -(d * i - f * g)
    C = d * h - e * g
    det = a * A + b * B + c * C
    if det == 0.0 or abs(det) < 1e-300:
        return None
    adj = [[A, -(b * i - c * h), b * f - c * e],
           [B, a * i - c * g, -(a * f - c * d)],
           [C, -(a * h - b * g), a * e - b * d]]
    try:
        return [[Fraction(adj[r][s] / det) for s in range(3)] for r in range(3)]
    except (OverflowError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Krawczyk certificate
# ---------------------------------------------------------------------------

def krawczyk(cfg: Config, box: tuple[IV, IV, IV]):
    """Return ("unique", K) / ("empty", None) / ("unknown", None) for box."""
    m = tuple(box[k].mid() for k in range(3))
    Fm = field_iv(cfg, tuple(IV(mi) for mi in m))
    if Fm is None:
        return ("unknown", None)
    J = jac_iv(cfg, box)
    if J is None:
        return ("unknown", None)
    Jf = _jac_float(cfg, tuple(float(mi) for mi in m))
    if Jf is None:
        return ("unknown", None)
    Y = _inv3_float(Jf)
    if Y is None:
        return ("unknown", None)

    # K_k = m_k - (Y Fm)_k + sum_l (I - Y J)_kl (box_l - m_l)
    K = []
    for k in range(3):
        acc = IV(m[k]) - (Fm[0].scale(Y[k][0]) + Fm[1].scale(Y[k][1]) + Fm[2].scale(Y[k][2]))
        for l in range(3):
            # (I - Y J)_kl as an interval
            yj = J[0][l].scale(Y[k][0]) + J[1][l].scale(Y[k][1]) + J[2][l].scale(Y[k][2])
            ent = IV(Fraction(1 if k == l else 0)) - yj
            acc = acc + ent * box[l].sub_point(m[l])
        K.append(acc)
    K = tuple(K)

    if all(K[k].inside_open(box[k]) for k in range(3)):
        return ("unique", K)
    if any(K[k].disjoint(box[k]) for k in range(3)):
        return ("empty", None)
    return ("unknown", None)


def det_sign(cfg: Config, box: tuple[IV, IV, IV]) -> int:
    """Certified sign of det DF over box: +1, -1, or 0 (= undecided)."""
    J = jac_iv(cfg, box)
    if J is None:
        return 0
    det = det3_iv(J)
    if det.lo > 0:
        return 1
    if det.hi < 0:
        return -1
    return 0


# ---------------------------------------------------------------------------
# lemmas
# ---------------------------------------------------------------------------

def localization_radius_sq(cfg: Config) -> Fraction:
    """R_0^2 for the all-positive localization lemma (exact)."""
    if not cfg.all_positive():
        raise ValueError("localization lemma requires all-positive charges")
    return max(p[0] ** 2 + p[1] ** 2 + p[2] ** 2 for _, p in cfg.charges)


def domination_radii(cfg: Config) -> list[Fraction]:
    """Certified rho_i for the domination lemma, exact rational arithmetic."""
    rhos = []
    for i, (qi, pi) in enumerate(cfg.charges):
        d_lo = []  # certified lower bounds on distances to the other charges
        for j, (qj, pj) in enumerate(cfg.charges):
            if j == i:
                continue
            d2 = sum((pi[k] - pj[k]) ** 2 for k in range(3))
            d_lo.append((abs(qj), _sqrt_dn(d2)))
        rho = min(d for _, d in d_lo) / 2
        ok = False
        for _ in range(80):
            # |q_i|/rho^2 > sum |q_j|/(d_ij - rho)^2, all in Q
            if all(d - rho > 0 for _, d in d_lo):
                rhs = sum(aq / (d - rho) ** 2 for aq, d in d_lo)
                if abs(qi) / rho ** 2 > rhs:
                    ok = True
                    break
            rho = rho / 2
        if not ok:
            raise RuntimeError(f"no domination radius found for charge {i}")
        rhos.append(rho)
    return rhos


def box_inside_ball(box, center, rho_sq: Fraction) -> bool:
    """Exact test: is the box contained in the closed ball around center?"""
    s = Fraction(0)
    for k in range(3):
        s += max(abs(box[k].lo - center[k]), abs(box[k].hi - center[k])) ** 2
    return s <= rho_sq


def box_min_norm_sq(box) -> Fraction:
    s = Fraction(0)
    for k in range(3):
        if box[k].lo > 0:
            s += box[k].lo ** 2
        elif box[k].hi < 0:
            s += box[k].hi ** 2
    return s


# ---------------------------------------------------------------------------
# the census driver
# ---------------------------------------------------------------------------

# non-dyadic offsets so subdivision faces avoid symmetric rational points
_OFFSETS = (Fraction(1, 97), Fraction(1, 101), Fraction(1, 103))


def auto_region(cfg: Config) -> tuple[tuple[IV, IV, IV], Fraction]:
    """Search box covering the localization ball, offset off symmetric points."""
    r0_sq = localization_radius_sq(cfg)
    r_up = _sqrt_up(r0_sq)
    half = r_up + Fraction(1, 16)
    box = tuple(IV(_OFFSETS[k] - half, _OFFSETS[k] + half) for k in range(3))
    return box, r0_sq


def count_equilibria(cfg: Config, region=None, min_width_bits: int = 30,
                     kraw_start: Fraction = Fraction(1, 2),
                     max_boxes: int = 400_000, progress_every: int = 0) -> dict:
    """Run the certified census; return the certificate as a JSON-able dict."""
    t0 = time.time()
    complete = region is None
    r0_sq = None
    if region is None:
        region, r0_sq = auto_region(cfg)

    rhos = domination_radii(cfg)
    rho_sqs = [r ** 2 for r in rhos]
    min_width = Fraction(1, 1 << min_width_bits)

    leaves = []          # (path, kind, extra)
    isolated = []
    unresolved = 0
    excl_counts = {"excluded-sign": 0, "excluded-newton": 0,
                   "excluded-ball": 0, "excluded-outside": 0}
    stack = [(region, "")]
    processed = 0

    while stack:
        box, path = stack.pop()
        processed += 1
        if progress_every and processed % progress_every == 0:
            print(f"  ... {processed} boxes, {len(stack)} queued, "
                  f"{len(isolated)} isolated, {unresolved} unresolved, "
                  f"{round(time.time() - t0)}s", file=sys.stderr)
        if processed > max_boxes:
            raise RuntimeError(f"box budget exceeded ({max_boxes}); "
                               f"raise --max-boxes or shrink the region")

        # 1. outside the localization sphere (complete mode only)
        if r0_sq is not None and box_min_norm_sq(box) >= r0_sq:
            leaves.append((path, "excluded-outside", None))
            excl_counts["excluded-outside"] += 1
            continue

        # 2. inside a domination ball around a charge
        ball = next((i for i, (_, p) in enumerate(cfg.charges)
                     if box_inside_ball(box, p, rho_sqs[i])), None)
        if ball is not None:
            leaves.append((path, "excluded-ball", ball))
            excl_counts["excluded-ball"] += 1
            continue

        # 3. sign exclusion by interval evaluation
        F = field_iv(cfg, box)
        if F is not None and any(not F[k].contains0() for k in range(3)):
            leaves.append((path, "excluded-sign", None))
            excl_counts["excluded-sign"] += 1
            continue

        wmax = max(box[k].width() for k in range(3))

        # 4. Krawczyk isolation
        if F is not None and wmax <= kraw_start:
            verdict, K = krawczyk(cfg, box)
            if verdict == "empty":
                leaves.append((path, "excluded-newton", None))
                excl_counts["excluded-newton"] += 1
                continue
            if verdict == "unique":
                enc = tuple(IV(max(K[k].lo, box[k].lo), min(K[k].hi, box[k].hi))
                            for k in range(3))
                s = det_sign(cfg, enc)
                if s == 0:
                    s = det_sign(cfg, box)
                leaves.append((path, "isolated", None))
                isolated.append({"path": path,
                                 "box": _box_json(box),
                                 "enclosure": _box_json(enc),
                                 "det_sign": s})
                continue

        # 5. give up or split
        if wmax <= min_width:
            leaves.append((path, "unresolved", None))
            unresolved += 1
            continue
        axis = max(range(3), key=lambda k: box[k].width())
        mid = (box[axis].lo + box[axis].hi) / 2  # exact
        for side, part in enumerate((IV(box[axis].lo, mid), IV(mid, box[axis].hi))):
            child = tuple(part if k == axis else box[k] for k in range(3))
            stack.append((child, path + f"{axis}{side}"))

    # index-sum check (complete all-positive census, no unresolved leaves)
    index_sum = None
    index_sum_ok = None
    if complete and unresolved == 0:
        signs = [it["det_sign"] for it in isolated]
        if all(s != 0 for s in signs):
            index_sum = sum(signs)
            index_sum_ok = (index_sum == 1 - cfg.n)

    return {
        "tool": "equilibria.py",
        "config": cfg.to_json(),
        "prec_bits": PREC_BITS,
        "min_width_bits": min_width_bits,
        "region": _box_json(region),
        "complete": complete,
        "localization_r0_sq": str(r0_sq) if r0_sq is not None else None,
        "domination_rhos": [str(r) for r in rhos],
        "count": len(isolated),
        "unresolved": unresolved,
        "excluded": excl_counts,
        "isolated": isolated,
        "leaves": [{"path": p, "kind": k, "ball": b} for p, k, b in leaves],
        "index_sum": index_sum,
        "index_sum_ok": index_sum_ok,
        "boxes_processed": processed,
        "seconds": round(time.time() - t0, 3),
    }


def _box_json(box) -> list[list[str]]:
    return [[str(box[k].lo), str(box[k].hi)] for k in range(3)]


# ---------------------------------------------------------------------------
# validation suite
# ---------------------------------------------------------------------------

def _mk(charges) -> Config:
    return Config([(Fraction(q), tuple(Fraction(v) for v in p)) for q, p in charges])


def validate() -> int:
    failures = 0

    def check(name, cond):
        nonlocal failures
        print(("  ok   " if cond else "  FAIL ") + name)
        if not cond:
            failures += 1

    print("interval primitives")
    x = IV(Fraction(-1, 3), Fraction(1, 5))
    check("sq nonnegative on straddling interval", x.sq().lo == 0 and x.sq().hi >= Fraction(1, 9))
    two = IV(Fraction(2))
    r = two.sqrt()
    check("sqrt(2) enclosure", r.lo < r.hi and r.lo ** 2 <= 2 <= r.hi ** 2
          and r.hi - r.lo < Fraction(1, 1 << 100))
    q = IV(Fraction(1)).div(IV(Fraction(3)))
    check("1/3 enclosure", q.lo <= Fraction(1, 3) <= q.hi and q.hi - q.lo < Fraction(1, 1 << 100))

    print("two unit charges at (+-1,0,0): exactly 1 equilibrium, index -1")
    cert = count_equilibria(_mk([(1, (1, 0, 0)), (1, (-1, 0, 0))]))
    check("count == 1", cert["count"] == 1)
    check("no unresolved leaves", cert["unresolved"] == 0)
    check("index sum == 1-n", cert["index_sum_ok"] is True)
    if cert["count"] == 1:
        enc = cert["isolated"][0]["enclosure"]
        check("equilibrium at the midpoint (origin)",
              all(Fraction(enc[k][0]) <= 0 <= Fraction(enc[k][1]) for k in range(3)))

    print("three collinear unit charges at x = -1, 0, 1: exactly n-1 = 2")
    cert = count_equilibria(_mk([(1, (-1, 0, 0)), (1, (0, 0, 0)), (1, (1, 0, 0))]))
    check("count == 2", cert["count"] == 2)
    check("no unresolved leaves", cert["unresolved"] == 0)
    check("index sum == 1-n", cert["index_sum_ok"] is True)

    print("unit charges at an equilateral triangle (e1,e2,e3): exactly 4")
    cert = count_equilibria(_mk([(1, (1, 0, 0)), (1, (0, 1, 0)), (1, (0, 0, 1))]))
    check("count == 4", cert["count"] == 4)
    check("no unresolved leaves", cert["unresolved"] == 0)
    check("index sum == 1-n", cert["index_sum_ok"] is True)
    if cert["count"] == 4:
        c = Fraction(1, 3)
        centroid_hits = sum(
            1 for it in cert["isolated"]
            if all(Fraction(it["enclosure"][k][0]) <= c <= Fraction(it["enclosure"][k][1])
                   for k in range(3)))
        check("one equilibrium at the centroid", centroid_hits == 1)

    print("FAILURES:", failures)
    return 1 if failures else 0


# ---------------------------------------------------------------------------
# cli
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--config", help="JSON file with charges (and optional region)")
    ap.add_argument("--out", help="output path (default $MATHLAB_OUT/maxwell-cert-<t>.json)")
    ap.add_argument("--min-width-bits", type=int, default=30)
    ap.add_argument("--prec-bits", type=int, default=128)
    ap.add_argument("--max-boxes", type=int, default=400_000)
    ap.add_argument("--progress-every", type=int, default=0,
                    help="print a progress line to stderr every N boxes")
    args = ap.parse_args()

    set_precision(args.prec_bits)

    if args.validate:
        return validate()
    if not args.config:
        ap.error("need --validate or --config")

    with open(args.config, encoding="utf-8") as fh:
        obj = json.load(fh)
    cfg = Config.from_json(obj)
    region = None
    if obj.get("region"):
        region = tuple(IV(Fraction(lo), Fraction(hi)) for lo, hi in obj["region"])
    cert = count_equilibria(cfg, region=region,
                            min_width_bits=args.min_width_bits,
                            max_boxes=args.max_boxes,
                            progress_every=args.progress_every)

    out = args.out
    if not out:
        outdir = os.environ.get("MATHLAB_OUT", "./out")
        os.makedirs(outdir, exist_ok=True)
        out = os.path.join(outdir, f"maxwell-cert-{int(time.time())}.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(cert, fh, indent=1)
    print(f"count={cert['count']} unresolved={cert['unresolved']} "
          f"index_sum={cert['index_sum']} boxes={cert['boxes_processed']} "
          f"secs={cert['seconds']} -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
