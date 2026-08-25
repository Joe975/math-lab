"""Exact arithmetic in a real multiquadratic field Q(sqrt d_1, ..., sqrt d_k).

The generators d_1..d_k must be pairwise coprime squarefree integers > 1.  Under
that hypothesis the 2^k numbers

    sqrt(prod_{i in A} d_i),   A a subset of {0..k-1},

are linearly independent over Q, so a field element is faithfully represented by
its coordinate vector and equality is a coordinate-wise integer comparison.

No floating point is used in any decision.  ``to_float`` and ``to_decimal``
exist only for reporting and for numeric cross-checks; nothing in this module
branches on them.

Representation: ``num`` is a tuple of 2^k ints indexed by the subset bitmask,
``den`` a positive int, kept in lowest terms.  Integer arithmetic throughout,
which is roughly an order of magnitude faster than Fraction for this workload.

An independently written reference implementation of the same arithmetic lives
in ``exact_field_poly.py``; ``selftest.py`` cross-checks the two.
"""

from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction
from math import gcd


def squarefree_part(m: int) -> tuple[int, int]:
    """Write m = f*f * s with s squarefree.  Returns (f, s).  m must be > 0."""
    if m <= 0:
        raise ValueError("squarefree_part expects a positive integer")
    f, s, d = 1, m, 2
    while d * d <= s:
        while s % (d * d) == 0:
            s //= d * d
            f *= d
        d += 1
    return f, s


def is_squarefree(m: int) -> bool:
    return m > 0 and squarefree_part(m)[1] == m


class MQField:
    """The field Q(sqrt d_1, ..., sqrt d_k) for pairwise coprime squarefree d_i."""

    __slots__ = ("gens", "dim", "_mask_prod", "_table")

    def __init__(self, gens):
        gens = tuple(int(g) for g in gens)
        for g in gens:
            if g <= 1 or not is_squarefree(g):
                raise ValueError(f"generator {g} must be squarefree and > 1")
        for i in range(len(gens)):
            for j in range(i + 1, len(gens)):
                if gcd(gens[i], gens[j]) != 1:
                    raise ValueError(
                        f"generators {gens[i]} and {gens[j]} are not coprime; "
                        "the subset basis is only valid for coprime generators"
                    )
        self.gens = gens
        self.dim = 1 << len(gens)
        prod = [1] * self.dim
        for mask in range(self.dim):
            p = 1
            for i, g in enumerate(gens):
                if mask >> i & 1:
                    p *= g
            prod[mask] = p
        self._mask_prod = tuple(prod)
        # table[i][j] = (result mask, rational multiplier)
        self._table = tuple(
            tuple((i ^ j, prod[i & j]) for j in range(self.dim))
            for i in range(self.dim)
        )

    def __repr__(self):
        return "Q(" + ", ".join(f"sqrt{g}" for g in self.gens) + ")"

    def __eq__(self, other):
        return isinstance(other, MQField) and self.gens == other.gens

    def __hash__(self):
        return hash(("MQField", self.gens))

    # -- constructors ----------------------------------------------------

    def zero(self):
        return MQ(self, (0,) * self.dim, 1)

    def one(self):
        return self.rational(1)

    def rational(self, q):
        q = Fraction(q)
        num = [0] * self.dim
        num[0] = q.numerator
        return MQ(self, tuple(num), q.denominator)

    def gen_sqrt(self, i):
        """sqrt(d_i)."""
        num = [0] * self.dim
        num[1 << i] = 1
        return MQ(self, tuple(num), 1)

    def sqrt_int(self, m):
        """sqrt(m) for a positive integer m, if it lies in this field.

        Raises ValueError when it does not, which is the point: it refuses to
        silently approximate a root the field cannot express.
        """
        f, s = squarefree_part(int(m))
        for mask in range(self.dim):
            if self._mask_prod[mask] == s:
                num = [0] * self.dim
                num[mask] = f
                return MQ(self, tuple(num), 1)
        raise ValueError(f"sqrt({m}) is not in {self}")

    def basis_element(self, mask):
        num = [0] * self.dim
        num[mask] = 1
        return MQ(self, tuple(num), 1)


class MQ:
    """An element of an :class:`MQField`.  Immutable and hashable."""

    __slots__ = ("field", "num", "den")

    def __init__(self, field, num, den, _normalized=False):
        if not _normalized:
            if den == 0:
                raise ZeroDivisionError("zero denominator")
            if den < 0:
                den = -den
                num = tuple(-x for x in num)
            g = den
            for x in num:
                g = gcd(g, x)
                if g == 1:
                    break
            if g > 1:
                num = tuple(x // g for x in num)
                den //= g
        self.field = field
        self.num = num
        self.den = den

    # -- structural ------------------------------------------------------

    def __eq__(self, other):
        if not isinstance(other, MQ):
            if other == 0:
                return self.is_zero()
            return self == self.field.rational(other)
        if self.field != other.field:
            raise ValueError("cannot compare elements of different fields")
        return self.num == other.num and self.den == other.den

    def __hash__(self):
        return hash((self.num, self.den))

    def is_zero(self):
        return not any(self.num)

    def is_rational(self):
        return not any(self.num[1:])

    def as_fraction(self):
        if not self.is_rational():
            raise ValueError("element is not rational")
        return Fraction(self.num[0], self.den)

    def coords(self):
        """Coordinates as Fractions, in bitmask order."""
        return tuple(Fraction(x, self.den) for x in self.num)

    def __repr__(self):
        parts = []
        for mask, x in enumerate(self.num):
            if x == 0:
                continue
            p = self.field._mask_prod[mask]
            parts.append(f"{x}" if p == 1 else f"{x}*sqrt{p}")
        body = " + ".join(parts) if parts else "0"
        return f"({body})/{self.den}" if self.den != 1 else f"({body})"

    # -- arithmetic ------------------------------------------------------

    def __add__(self, other):
        other = self._coerce(other)
        a, b = self.den, other.den
        g = gcd(a, b)
        la, lb = b // g, a // g
        num = tuple(x * la + y * lb for x, y in zip(self.num, other.num))
        return MQ(self.field, num, a * la)  # a*la == lcm(a, b)

    __radd__ = __add__

    def __neg__(self):
        return MQ(self.field, tuple(-x for x in self.num), self.den, _normalized=True)

    def __sub__(self, other):
        return self + (-self._coerce(other))

    def __rsub__(self, other):
        return self._coerce(other) + (-self)

    def __mul__(self, other):
        other = self._coerce(other)
        f = self.field
        dim = f.dim
        table = f._table
        acc = [0] * dim
        an, bn = self.num, other.num
        for i in range(dim):
            ai = an[i]
            if ai == 0:
                continue
            row = table[i]
            for j in range(dim):
                bj = bn[j]
                if bj == 0:
                    continue
                mask, mult = row[j]
                acc[mask] += ai * bj * mult
        return MQ(f, tuple(acc), self.den * other.den)

    __rmul__ = __mul__

    def conjugate(self, i):
        """Apply the automorphism sqrt(d_i) -> -sqrt(d_i)."""
        bit = 1 << i
        num = tuple(-x if (m & bit) else x for m, x in enumerate(self.num))
        return MQ(self.field, num, self.den, _normalized=True)

    def inverse(self):
        if self.is_zero():
            raise ZeroDivisionError("inverse of zero")
        # Multiply by every nontrivial conjugate; the result is rational.
        acc = self.field.one()
        for mask in range(1, self.field.dim):
            c = self
            for i in range(len(self.field.gens)):
                if mask >> i & 1:
                    c = c.conjugate(i)
            acc = acc * c
        norm = (self * acc)
        if not norm.is_rational():
            raise ArithmeticError("norm is not rational -- basis assumption violated")
        q = norm.as_fraction()
        if q == 0:
            raise ArithmeticError("zero norm for a nonzero element -- basis assumption violated")
        return acc * self.field.rational(Fraction(q.denominator, q.numerator))

    def __truediv__(self, other):
        other = self._coerce(other)
        return self * other.inverse()

    def _coerce(self, other):
        if isinstance(other, MQ):
            if other.field != self.field:
                raise ValueError("mixed fields")
            return other
        return self.field.rational(other)

    # -- reporting only (never used in a decision) -----------------------

    def to_float(self):
        f = self.field
        return sum(x * f._mask_prod[m] ** 0.5 for m, x in enumerate(self.num)) / self.den

    def to_decimal(self, prec=120):
        """High-precision Decimal value.  For cross-checks and reporting only."""
        ctx_prec = getcontext().prec
        getcontext().prec = prec + 20
        try:
            f = self.field
            total = Decimal(0)
            for m, x in enumerate(self.num):
                if x == 0:
                    continue
                p = f._mask_prod[m]
                total += Decimal(x) * (Decimal(1) if p == 1 else Decimal(p).sqrt())
            out = total / Decimal(self.den)
        finally:
            getcontext().prec = ctx_prec
        return +out


def field_for(*ints):
    """Smallest coprime-squarefree multiquadratic field containing sqrt of each arg.

    Splits the squarefree parts into pairwise coprime pieces, so e.g.
    ``field_for(3, 33)`` returns Q(sqrt3, sqrt11) and ``field_for(3, 27)``
    returns Q(sqrt3).
    """
    parts = []
    for m in ints:
        s = squarefree_part(int(m))[1]
        if s == 1:
            continue
        parts.append(s)
    gens = []
    for s in parts:
        cur = s
        new = []
        for g in gens:
            d = gcd(g, cur)
            if d == 1:
                new.append(g)
            else:
                if d != 1:
                    new.append(d)
                if g // d > 1:
                    new.append(g // d)
                cur //= d
        if cur > 1:
            new.append(cur)
        # dedupe
        gens = sorted(set(new))
    return MQField(tuple(gens))
