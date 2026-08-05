#!/usr/bin/env python3
"""Independent checker for point-charge equilibrium certificates.

This tool re-verifies a certificate produced by the reference counting tool
in this directory.  It shares no code with it beyond fractions.Fraction, and
deliberately uses different machinery for every layer of the claim:

  arithmetic     directed-rounding decimal interval arithmetic (the decimal
                 module's ROUND_FLOOR / ROUND_CEILING contexts; Decimal
                 add/mul/div/sqrt are correctly rounded per IBM's decimal
                 specification), with a precision-doubling retry ladder --
                 instead of dyadic rational intervals;
  isolation      the preconditioned interval Newton operator with interval
                 Gaussian elimination:  N(B) = m - IGA(Y*J(B), Y*F(m)) for a
                 rational approximate midpoint-Jacobian inverse Y.  If the
                 interval Gauss algorithm succeeds (no pivot straddles 0)
                 and N(B) lies in the interior of B, then B contains exactly
                 one zero of F -- existence via the mean-value extension and
                 Brouwer, uniqueness because IGA success forces every matrix
                 in Y*J(B), hence in J(B), to be nonsingular; the
                 certificate is valid for ANY fixed Y, so the float-computed
                 Y is not part of the trust chain.  See A. Neumaier,
                 "Interval Methods for Systems of Equations", CUP 1990,
                 chapter 5 -- instead of the Krawczyk operator;
  tiling         an exact combinatorial check that the leaf boxes partition
                 the search region: leaf split paths must be prefix-free and
                 satisfy the Kraft equality  sum 2^-depth(leaf) = 1,  which
                 for a binary subdivision tree holds iff the leaves tile the
                 root; every leaf box is rebuilt exactly (Fraction midpoint
                 splitting) from its path;
  lemmas         re-derived, not trusted.  Localization (all-positive
                 charges): an equilibrium x satisfies |x| < R_0 = max|x_i|,
                 since x.F(x) = sum q_i x.(x-x_i)/|x-x_i|^3 is positive for
                 |x| >= R_0 (each term is q_i(|x|^2 - x.x_i)/... >= 0, not
                 all zero).  Domination (any signs): if 0 < rho < min_j d_ij
                 and |q_i|/rho^2 > sum_{j!=i} |q_j|/(d_ij - rho)^2 then no
                 zero of F has 0 < |x - x_i| <= rho, by the triangle
                 inequality on the field terms.  Both inequalities are
                 re-instantiated here in exact rational arithmetic with this
                 tool's own certified square-root bounds (integer isqrt).

Every leaf of the certificate is re-checked according to its kind, the
count is recomputed, and for complete all-positive censuses the index-sum
identity (sum of certified signs of det DF over the equilibria = 1 - n, a
classical Poincare-Hopf/degree computation) is re-verified with this tool's
own interval determinant.

Verdict: PASS means every leaf certificate was independently re-established.
FAIL lists every leaf that could not be re-established -- which does not by
itself mean the claim is false, but does mean the certificate cannot be
cited.  A certificate that PASSES here at a record-adjacent count is to be
reported as requiring escalation, not accepted (see the problem's
verification contract).

Usage:
  verify_equilibria.py --cert CERT.json          re-verify a certificate
  verify_equilibria.py --validate                run the validation suite

Standard library only.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
from math import isqrt

# ---------------------------------------------------------------------------
# directed-rounding decimal interval arithmetic
# ---------------------------------------------------------------------------

class DCtx:
    """A pair of decimal contexts rounding down / up at a given precision."""

    def __init__(self, digits: int):
        self.digits = digits
        self.dn = Context(prec=digits, rounding=ROUND_FLOOR)
        self.up = Context(prec=digits, rounding=ROUND_CEILING)

    def from_fraction(self, x: Fraction) -> "DIV":
        n = Decimal(x.numerator)
        d = Decimal(x.denominator)
        return DIV(self.dn.divide(n, d), self.up.divide(n, d), self)


class DIV:
    """Closed interval with Decimal endpoints and outward rounding."""

    __slots__ = ("lo", "hi", "ctx")

    def __init__(self, lo: Decimal, hi: Decimal, ctx: DCtx):
        if lo > hi:
            raise ValueError("empty interval")
        self.lo = lo
        self.hi = hi
        self.ctx = ctx

    def __add__(self, o: "DIV") -> "DIV":
        c = self.ctx
        return DIV(c.dn.add(self.lo, o.lo), c.up.add(self.hi, o.hi), c)

    def __sub__(self, o: "DIV") -> "DIV":
        c = self.ctx
        return DIV(c.dn.subtract(self.lo, o.hi), c.up.subtract(self.hi, o.lo), c)

    def __neg__(self) -> "DIV":
        return DIV(-self.hi, -self.lo, self.ctx)

    def __mul__(self, o: "DIV") -> "DIV":
        c = self.ctx
        lo = min(c.dn.multiply(a, b) for a in (self.lo, self.hi) for b in (o.lo, o.hi))
        hi = max(c.up.multiply(a, b) for a in (self.lo, self.hi) for b in (o.lo, o.hi))
        return DIV(lo, hi, c)

    def sq(self) -> "DIV":
        c = self.ctx
        a, b = abs(self.lo), abs(self.hi)
        hi = c.up.multiply(max(a, b), max(a, b))
        if self.lo <= 0 <= self.hi:
            return DIV(Decimal(0), hi, c)
        m = min(a, b)
        return DIV(c.dn.multiply(m, m), hi, c)

    def div(self, o: "DIV") -> "DIV":
        if o.lo <= 0:
            raise ZeroDivisionError("divisor interval must be strictly positive")
        c = self.ctx
        lo = min(c.dn.divide(a, b) for a in (self.lo, self.hi) for b in (o.lo, o.hi))
        hi = max(c.up.divide(a, b) for a in (self.lo, self.hi) for b in (o.lo, o.hi))
        return DIV(lo, hi, c)

    def sqrt(self) -> "DIV":
        # The decimal spec rounds squareRoot half-even REGARDLESS of the
        # context rounding, so directed bounds must be established by hand:
        # take the correctly-rounded root and step it by ulps until the
        # bound is proved by exact rational comparison of squares.
        c = self.ctx
        if self.lo < 0:
            raise ValueError("sqrt of negative")
        lo = c.dn.sqrt(self.lo)
        while _dec_to_frac(lo) ** 2 > _dec_to_frac(self.lo):
            lo = c.dn.next_minus(lo)
        hi = c.up.sqrt(self.hi)
        while _dec_to_frac(hi) ** 2 < _dec_to_frac(self.hi):
            hi = c.up.next_plus(hi)
        return DIV(lo, hi, c)

    def contains0(self) -> bool:
        return self.lo <= 0 <= self.hi

    def mig(self) -> Decimal:
        """Mignitude: min |y| over the interval (0 if it straddles 0)."""
        if self.contains0():
            return Decimal(0)
        return min(abs(self.lo), abs(self.hi))

    def inside_open_frac(self, lo: Fraction, hi: Fraction) -> bool:
        """Is this interval strictly inside the rational interval (lo, hi)?"""
        return (Fraction(self.lo.as_integer_ratio()[0], self.lo.as_integer_ratio()[1]) > lo
                and Fraction(self.hi.as_integer_ratio()[0], self.hi.as_integer_ratio()[1]) < hi)


def _zero(ctx: DCtx) -> DIV:
    return DIV(Decimal(0), Decimal(0), ctx)


def _dec_to_frac(d: Decimal) -> Fraction:
    return Fraction(*d.as_integer_ratio())


# ---------------------------------------------------------------------------
# certified rational square-root bounds (for the exact lemma re-derivations)
# ---------------------------------------------------------------------------

_LEMMA_BITS = 96


def sqrt_lower(f: Fraction) -> Fraction:
    """Certified rational lower bound for sqrt(f), f >= 0."""
    p, q = f.numerator, f.denominator
    s = 1 << _LEMMA_BITS
    r = isqrt(p * q * s * s)
    return Fraction(r, q * s)


def sqrt_upper(f: Fraction) -> Fraction:
    p, q = f.numerator, f.denominator
    s = 1 << _LEMMA_BITS
    t = p * q * s * s
    r = isqrt(t)
    if r * r != t:
        r += 1
    return Fraction(r, q * s)


# ---------------------------------------------------------------------------
# field and Jacobian in decimal intervals (independently written)
# ---------------------------------------------------------------------------

def d_field(charges, box, ctx: DCtx):
    """Enclosure of F over box; box is a triple of DIV.  None if r^2 may be 0."""
    acc = [_zero(ctx), _zero(ctx), _zero(ctx)]
    for q, pos in charges:
        qi = ctx.from_fraction(q)
        d = [box[k] - ctx.from_fraction(pos[k]) for k in range(3)]
        r2 = d[0].sq() + d[1].sq() + d[2].sq()
        if r2.lo <= 0:
            return None
        r3 = r2.sqrt() * r2
        for k in range(3):
            acc[k] = acc[k] + (qi * d[k]).div(r3)
    return acc


def d_jac(charges, box, ctx: DCtx):
    one = DIV(Decimal(1), Decimal(1), ctx)
    three = DIV(Decimal(3), Decimal(3), ctx)
    J = [[_zero(ctx) for _ in range(3)] for _ in range(3)]
    for q, pos in charges:
        qi = ctx.from_fraction(q)
        d = [box[k] - ctx.from_fraction(pos[k]) for k in range(3)]
        r2 = d[0].sq() + d[1].sq() + d[2].sq()
        if r2.lo <= 0:
            return None
        r3 = r2.sqrt() * r2
        r5 = r3 * r2
        for k in range(3):
            for l in range(3):
                if k == l:
                    t = one.div(r3) - (three * d[k].sq()).div(r5)
                else:
                    t = -((three * d[k] * d[l]).div(r5))
                J[k][l] = J[k][l] + qi * t
    return J


def d_det3(J):
    a, b, c = J[0]
    d, e, f = J[1]
    g, h, i = J[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


# ---------------------------------------------------------------------------
# interval Newton via interval Gaussian elimination
# ---------------------------------------------------------------------------

def interval_gauss_solve(A, b, ctx: DCtx):
    """Enclosure of {x : Mx = r, M in A, r in b} by interval Gaussian
    elimination with mignitude pivoting.  None if a pivot straddles 0
    (in which case nothing is certified)."""
    A = [row[:] for row in A]
    b = b[:]
    n = 3
    perm = list(range(n))
    for col in range(n):
        piv = max(range(col, n), key=lambda r: A[r][col].mig())
        if A[piv][col].mig() == 0:
            return None
        if piv != col:
            A[col], A[piv] = A[piv], A[col]
            b[col], b[piv] = b[piv], b[col]
            perm[col], perm[piv] = perm[piv], perm[col]
        for r in range(col + 1, n):
            m = _div_signed(A[r][col], A[col][col])
            for cc in range(col, n):
                A[r][cc] = A[r][cc] - m * A[col][cc]
            b[r] = b[r] - m * b[col]
    x = [None] * n
    for r in range(n - 1, -1, -1):
        s = b[r]
        for cc in range(r + 1, n):
            s = s - A[r][cc] * x[cc]
        x[r] = _div_signed(s, A[r][r])
    return x


def _div_signed(a: DIV, b: DIV) -> DIV:
    """a / b for an interval b not containing 0 (either sign)."""
    if b.lo > 0:
        return a.div(b)
    if b.hi < 0:
        return (-a).div(-b)
    raise ZeroDivisionError("pivot straddles zero")


def _float_jac(charges, x):
    """Plain floating-point midpoint Jacobian (candidate machinery only)."""
    J = [[0.0] * 3 for _ in range(3)]
    for q, pos in charges:
        d = [float(x[k]) - float(pos[k]) for k in range(3)]
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


def _float_inv3(J):
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
    return [[adj[r][s] / det for s in range(3)] for r in range(3)]


def newton_certify_unique(charges, box_frac, ctx: DCtx):
    """Preconditioned interval Newton certificate that box (rational
    endpoints) contains exactly one zero of F.  Returns True only on full
    success.

    Preconditioning is sound: for a fixed rational matrix Y, any zero x* of
    F in B satisfies Y*F(m) + (Y*M)(x* - m) = 0 for some M in J(B) (mean
    value theorem row by row), so x* lies in m - IGA(Y*J(B), Y*F(m)); and
    success of the interval Gauss algorithm on Y*J(B) makes every Y*M --
    hence every M -- nonsingular, which forbids two zeros in B.  The
    certificate does not depend on the quality of Y, only its success rate
    does."""
    N = _newton_image(charges, box_frac, ctx)
    if N is None:
        return False
    for k in range(3):
        if not N[k].inside_open_frac(box_frac[k][0], box_frac[k][1]):
            return False
    return True


def newton_certify_empty(charges, box_frac, ctx: DCtx):
    """Mean-value exclusion: every zero of F in B lies in N(B); if some
    component of N(B) is disjoint from B, then B contains no zero."""
    N = _newton_image(charges, box_frac, ctx)
    if N is None:
        return False
    for k in range(3):
        lo, hi = box_frac[k]
        if _dec_to_frac(N[k].hi) < lo or _dec_to_frac(N[k].lo) > hi:
            return True
    return False


def krawczyk_certify_empty(charges, box_frac, ctx: DCtx):
    """Krawczyk-image exclusion:  K(B) = m - Y*F(m) + (I - Y*J(B))(B - m)
    contains every zero of F in B (mean-value form, valid for any fixed Y,
    with no regularity requirement); if a component of K(B) is disjoint
    from B, then B has no zero.  This mirrors the reference tool's exclusion
    OPERATOR -- it is used here only for emptiness, where the interval
    Newton form is structurally weaker (a straddling pivot aborts the
    Gaussian elimination, while the Krawczyk image needs no linear solve).
    Isolation certificates never use it; independence for emptiness rests
    on the different arithmetic and implementation."""
    m = [(lo + hi) / 2 for lo, hi in box_frac]
    box = [DIV(ctx.from_fraction(lo).lo, ctx.from_fraction(hi).hi, ctx)
           for lo, hi in box_frac]
    Fm = d_field(charges, [ctx.from_fraction(mi) for mi in m], ctx)
    if Fm is None:
        return False
    J = d_jac(charges, box, ctx)
    if J is None:
        return False
    Jf = _float_jac(charges, m)
    if Jf is None:
        return False
    Yf = _float_inv3(Jf)
    if Yf is None:
        return False
    Y = [[ctx.from_fraction(Fraction(Yf[r][s])) for s in range(3)] for r in range(3)]
    one = DIV(Decimal(1), Decimal(1), ctx)
    for k in range(3):
        acc = ctx.from_fraction(m[k]) - (
            Y[k][0] * Fm[0] + Y[k][1] * Fm[1] + Y[k][2] * Fm[2])
        for l in range(3):
            yj = Y[k][0] * J[0][l] + Y[k][1] * J[1][l] + Y[k][2] * J[2][l]
            ent = (one - yj) if k == l else -yj
            acc = acc + ent * (box[l] - ctx.from_fraction(m[l]))
        lo, hi = box_frac[k]
        if _dec_to_frac(acc.hi) < lo or _dec_to_frac(acc.lo) > hi:
            return True
    return False


def _newton_image(charges, box_frac, ctx: DCtx):
    """N(B) = m - IGA(Y*J(B), Y*F(m)) as a list of DIV, or None."""
    m = [(lo + hi) / 2 for lo, hi in box_frac]
    box = [DIV(ctx.from_fraction(lo).lo, ctx.from_fraction(hi).hi, ctx)
           for lo, hi in box_frac]
    Fm = d_field(charges, [ctx.from_fraction(mi) for mi in m], ctx)
    if Fm is None:
        return None
    J = d_jac(charges, box, ctx)
    if J is None:
        return None
    Jf = _float_jac(charges, m)
    if Jf is None:
        return None
    Yf = _float_inv3(Jf)
    if Yf is None:
        return None
    Y = [[ctx.from_fraction(Fraction(Yf[r][s])) for s in range(3)] for r in range(3)]
    # precondition: YJ = Y * J(B), Yb = Y * F(m), as interval enclosures
    YJ = [[Y[r][0] * J[0][s] + Y[r][1] * J[1][s] + Y[r][2] * J[2][s]
           for s in range(3)] for r in range(3)]
    Yb = [Y[r][0] * Fm[0] + Y[r][1] * Fm[1] + Y[r][2] * Fm[2] for r in range(3)]
    try:
        delta = interval_gauss_solve(YJ, Yb, ctx)
    except ZeroDivisionError:
        return None
    if delta is None:
        return None
    return [ctx.from_fraction(m[k]) - delta[k] for k in range(3)]


# ---------------------------------------------------------------------------
# certificate re-verification
# ---------------------------------------------------------------------------

def _parse_box(bx) -> list[tuple[Fraction, Fraction]]:
    return [(Fraction(lo), Fraction(hi)) for lo, hi in bx]


def _rebuild_box(root, path: str) -> list[tuple[Fraction, Fraction]]:
    """Rebuild a leaf box from its split path, mirroring midpoint bisection."""
    box = list(root)
    if len(path) % 2:
        raise ValueError(f"odd-length path {path!r}")
    for i in range(0, len(path), 2):
        axis = int(path[i])
        side = int(path[i + 1])
        lo, hi = box[axis]
        mid = (lo + hi) / 2
        box[axis] = (lo, mid) if side == 0 else (mid, hi)
    return box


def _min_norm_sq(box) -> Fraction:
    s = Fraction(0)
    for lo, hi in box:
        if lo > 0:
            s += lo ** 2
        elif hi < 0:
            s += hi ** 2
    return s


def _max_dist_sq(box, center) -> Fraction:
    s = Fraction(0)
    for (lo, hi), c in zip(box, center):
        s += max(abs(lo - c), abs(hi - c)) ** 2
    return s


def verify(cert: dict, base_digits: int = 60, max_digits: int = 480,
           progress: bool = False) -> dict:
    t0 = time.time()
    charges = [(Fraction(c["q"]), tuple(Fraction(v) for v in c["pos"]))
               for c in cert["config"]["charges"]]
    n = len(charges)
    all_positive = all(q > 0 for q, _ in charges)
    root = _parse_box(cert["region"])
    complete = bool(cert.get("complete"))
    failures: list[str] = []

    # --- lemma re-derivations (exact) -------------------------------------
    r0_sq = None
    if complete:
        if not all_positive:
            failures.append("completeness claimed but charges are not all positive")
        else:
            r0_sq = max(p[0] ** 2 + p[1] ** 2 + p[2] ** 2 for _, p in charges)
            r_up = sqrt_upper(r0_sq)
            for k in range(3):
                if not (root[k][0] <= -r_up and r_up <= root[k][1]):
                    failures.append("search region does not contain the localization ball")
                    break

    rhos = [Fraction(r) for r in cert["domination_rhos"]]
    if len(rhos) != n:
        failures.append("wrong number of domination radii")
    else:
        for i, (qi, pi) in enumerate(charges):
            rho = rhos[i]
            ok = rho > 0
            rhs = Fraction(0)
            for j, (qj, pj) in enumerate(charges):
                if j == i:
                    continue
                d_lo = sqrt_lower(sum((pi[k] - pj[k]) ** 2 for k in range(3)))
                if d_lo - rho <= 0:
                    ok = False
                    break
                rhs += abs(qj) / (d_lo - rho) ** 2
            if not (ok and abs(qi) / rho ** 2 > rhs):
                failures.append(f"domination inequality fails for charge {i}")

    # --- structural: leaves tile the region -------------------------------
    leaves = cert["leaves"]
    paths = [lf["path"] for lf in leaves]
    if len(set(paths)) != len(paths):
        failures.append("duplicate leaf paths")
    spaths = sorted(paths)
    for a, b in zip(spaths, spaths[1:]):
        if b.startswith(a):
            failures.append(f"leaf {a!r} is a prefix of leaf {b!r}")
            break
    kraft = sum(Fraction(1, 1 << (len(p) // 2)) for p in paths)
    if kraft != 1:
        failures.append(f"Kraft sum is {kraft}, not 1: leaves do not tile the region")

    # --- per-leaf re-verification -----------------------------------------
    iso_by_path = {it["path"]: it for it in cert["isolated"]}
    unresolved = 0
    counted = 0
    index_signs = []
    kinds = {"excluded-sign": 0, "excluded-newton": 0, "excluded-ball": 0,
             "excluded-outside": 0, "isolated": 0, "unresolved": 0}

    for idx, lf in enumerate(leaves):
        if progress and idx % 2000 == 0 and idx:
            print(f"  ... {idx}/{len(leaves)} leaves", file=sys.stderr)
        kind = lf["kind"]
        kinds[kind] = kinds.get(kind, 0) + 1
        box = _rebuild_box(root, lf["path"])

        if kind == "excluded-outside":
            if r0_sq is None:
                failures.append(f"leaf {lf['path']!r}: outside-exclusion without "
                                f"a localization lemma in force")
            elif _min_norm_sq(box) < r0_sq:
                failures.append(f"leaf {lf['path']!r}: outside-exclusion is wrong")
            continue

        if kind == "excluded-ball":
            i = lf["ball"]
            if not (isinstance(i, int) and 0 <= i < n
                    and _max_dist_sq(box, charges[i][1]) <= rhos[i] ** 2):
                failures.append(f"leaf {lf['path']!r}: ball-exclusion is wrong")
            continue

        if kind == "excluded-sign":
            digits = base_digits
            ok = False
            while digits <= max_digits:
                ctx = DCtx(digits)
                dbox = [DIV(ctx.from_fraction(lo).lo, ctx.from_fraction(hi).hi, ctx)
                        for lo, hi in box]
                F = d_field(charges, dbox, ctx)
                if F is not None and any(not F[k].contains0() for k in range(3)):
                    ok = True
                    break
                digits *= 2
            if not ok:
                failures.append(f"leaf {lf['path']!r}: sign-exclusion not reproduced")
            continue

        if kind == "excluded-newton":
            digits = base_digits
            ok = False
            while digits <= max_digits:
                ctx = DCtx(digits)
                dbox = [DIV(ctx.from_fraction(lo).lo, ctx.from_fraction(hi).hi, ctx)
                        for lo, hi in box]
                F = d_field(charges, dbox, ctx)
                if F is not None and any(not F[k].contains0() for k in range(3)):
                    ok = True  # sign exclusion is an acceptable re-derivation too
                    break
                if newton_certify_empty(charges, box, DCtx(digits)):
                    ok = True
                    break
                if krawczyk_certify_empty(charges, box, DCtx(digits)):
                    ok = True
                    break
                digits *= 2
            if not ok:
                failures.append(f"leaf {lf['path']!r}: emptiness not reproduced")
            continue

        if kind == "unresolved":
            unresolved += 1
            continue

        if kind == "isolated":
            it = iso_by_path.get(lf["path"])
            if it is None:
                failures.append(f"leaf {lf['path']!r}: isolated leaf missing record")
                continue
            claimed = _parse_box(it["box"])
            if claimed != box:
                failures.append(f"leaf {lf['path']!r}: recorded box mismatches path")
                continue
            digits = base_digits
            ok = False
            while digits <= max_digits:
                if newton_certify_unique(charges, box, DCtx(digits)):
                    ok = True
                    break
                digits *= 2
            if not ok:
                failures.append(f"leaf {lf['path']!r}: interval Newton could not "
                                f"re-establish isolation")
                continue
            counted += 1
            # det sign over the claimed enclosure
            s_claim = it.get("det_sign", 0)
            enc = _parse_box(it["enclosure"])
            for k in range(3):
                if not (box[k][0] <= enc[k][0] and enc[k][1] <= box[k][1]):
                    failures.append(f"leaf {lf['path']!r}: enclosure not inside box")
            s_here = 0
            digits = base_digits
            while digits <= max_digits and s_here == 0:
                ctx = DCtx(digits)
                dbox = [DIV(ctx.from_fraction(lo).lo, ctx.from_fraction(hi).hi, ctx)
                        for lo, hi in enc]
                J = d_jac(charges, dbox, ctx)
                if J is not None:
                    det = d_det3(J)
                    if det.lo > 0:
                        s_here = 1
                    elif det.hi < 0:
                        s_here = -1
                digits *= 2
            if s_claim != 0 and s_here != 0 and s_claim != s_here:
                failures.append(f"leaf {lf['path']!r}: det sign disagrees")
            index_signs.append(s_here)
            continue

        failures.append(f"leaf {lf['path']!r}: unknown kind {kind!r}")

    if counted != cert["count"]:
        failures.append(f"recount: {counted} isolated leaves re-established, "
                        f"certificate claims {cert['count']}")
    if unresolved != cert["unresolved"]:
        failures.append(f"unresolved recount: {unresolved} vs claimed "
                        f"{cert['unresolved']}")

    index_sum = None
    index_sum_ok = None
    if complete and all_positive and unresolved == 0 and counted and \
            all(s != 0 for s in index_signs):
        index_sum = sum(index_signs)
        index_sum_ok = (index_sum == 1 - n)
        if not index_sum_ok:
            failures.append(f"index-sum identity fails: {index_sum} != {1 - n}")

    return {
        "verdict": "PASS" if not failures else "FAIL",
        "failures": failures,
        "count": counted,
        "unresolved": unresolved,
        "kinds": kinds,
        "index_sum": index_sum,
        "index_sum_ok": index_sum_ok,
        "seconds": round(time.time() - t0, 3),
    }


# ---------------------------------------------------------------------------
# validation suite
# ---------------------------------------------------------------------------

def validate() -> int:
    failures = 0

    def check(name, cond):
        nonlocal failures
        print(("  ok   " if cond else "  FAIL ") + name)
        if not cond:
            failures += 1

    print("decimal interval primitives")
    ctx = DCtx(40)
    third = ctx.from_fraction(Fraction(1, 3))
    check("1/3 outward", Fraction(str(third.lo)) <= Fraction(1, 3) <= Fraction(str(third.hi))
          and third.lo != third.hi)
    two = ctx.from_fraction(Fraction(2))
    r = two.sqrt()
    check("sqrt(2) outward", Fraction(str(r.lo)) ** 2 <= 2 <= Fraction(str(r.hi)) ** 2)
    x = DIV(Decimal("-1"), Decimal("0.5"), ctx)
    check("sq nonnegative on straddling interval", x.sq().lo == 0)

    print("certified rational sqrt bounds")
    check("sqrt_lower/upper bracket sqrt(2)",
          sqrt_lower(Fraction(2)) ** 2 <= 2 <= sqrt_upper(Fraction(2)) ** 2
          and sqrt_upper(Fraction(2)) - sqrt_lower(Fraction(2)) < Fraction(1, 1 << 90))

    print("interval Newton core: two unit charges at (+-1,0,0)")
    charges = [(Fraction(1), (Fraction(1), Fraction(0), Fraction(0))),
               (Fraction(1), (Fraction(-1), Fraction(0), Fraction(0)))]
    box = [(Fraction(-1, 40), Fraction(1, 48))] * 3  # off-center box at origin
    check("unique zero certified near origin",
          newton_certify_unique(charges, box, DCtx(60)))
    far = [(Fraction(1, 3), Fraction(1, 2))] * 3
    check("no false certificate away from the zero",
          not newton_certify_unique(charges, far, DCtx(60)))

    print("tamper detection on synthetic certificates")
    cfg2 = {"charges": [{"q": "1", "pos": ["1", "0", "0"]},
                        {"q": "1", "pos": ["-1", "0", "0"]}]}
    # (a) a region that truly contains the origin equilibrium, declared
    #     sign-excluded in one leaf: must be caught.
    bad = {
        "config": cfg2,
        "region": [["-9/64", "1/8"], ["-1/8", "1/8"], ["-1/8", "1/8"]],
        "complete": False,
        "domination_rhos": ["1/4", "1/4"],
        "count": 0, "unresolved": 0, "isolated": [],
        "leaves": [{"path": "", "kind": "excluded-sign", "ball": None}],
    }
    res = verify(bad, max_digits=120)
    check("wrong sign-exclusion detected", res["verdict"] == "FAIL"
          and any("sign-exclusion" in f for f in res["failures"]))
    # (b) a genuinely field-positive region away from zeros and charges:
    #     the same claim must be accepted.
    good = dict(bad)
    good["region"] = [["2", "9/4"], ["0", "1/4"], ["0", "1/4"]]
    res = verify(good)
    check("correct sign-exclusion accepted", res["verdict"] == "PASS")
    # (c) a leaf list with a hole: the Kraft check must catch it.
    gap = dict(good)
    gap["leaves"] = [{"path": "00", "kind": "excluded-sign", "ball": None}]
    res = verify(gap)
    check("coverage gap detected (Kraft sum)",
          any("Kraft" in f for f in res["failures"]))

    print("FAILURES:", failures)
    return 1 if failures else 0


# ---------------------------------------------------------------------------
# cli
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cert", help="certificate JSON to re-verify")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--digits", type=int, default=60,
                    help="starting decimal precision (doubles on retry)")
    ap.add_argument("--progress", action="store_true")
    args = ap.parse_args()

    if args.validate:
        return validate()
    if not args.cert:
        ap.error("need --cert or --validate")

    with open(args.cert, encoding="utf-8") as fh:
        cert = json.load(fh)
    res = verify(cert, base_digits=args.digits, progress=args.progress)
    print(json.dumps(res, indent=1))
    return 0 if res["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
