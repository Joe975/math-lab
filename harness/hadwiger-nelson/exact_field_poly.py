"""Reference implementation of multiquadratic arithmetic, written independently.

Where ``exact_field`` stores a coordinate vector of integers over a precomputed
subset multiplication table, this module stores a dict of monomials with
Fraction coefficients and multiplies polynomials generically, reducing modulo
x_i^2 - d_i afterwards.  Different data structure, different algorithm, no
shared code -- so agreement between the two is a real cross-check on the
arithmetic rather than a re-run.

It is deliberately slower and simpler.  Use it to check the fast path, not to
run searches.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product


class PolyQuad:
    """Element of Q[x_0..x_{k-1}] / (x_i^2 - d_i), i.e. Q(sqrt d_0, ...).

    Monomials are exponent tuples in {0,1}^k after reduction.
    """

    __slots__ = ("mods", "terms")

    def __init__(self, mods, terms):
        self.mods = tuple(mods)
        self.terms = {m: Fraction(c) for m, c in terms.items() if c != 0}

    # -- constructors ----------------------------------------------------

    @classmethod
    def constant(cls, mods, c):
        k = len(mods)
        return cls(mods, {(0,) * k: Fraction(c)})

    @classmethod
    def root(cls, mods, i):
        k = len(mods)
        e = [0] * k
        e[i] = 1
        return cls(mods, {tuple(e): Fraction(1)})

    # -- arithmetic ------------------------------------------------------

    def _co(self, other):
        if isinstance(other, PolyQuad):
            assert other.mods == self.mods
            return other
        return PolyQuad.constant(self.mods, other)

    def __add__(self, other):
        other = self._co(other)
        t = dict(self.terms)
        for m, c in other.terms.items():
            t[m] = t.get(m, Fraction(0)) + c
        return PolyQuad(self.mods, t)

    __radd__ = __add__

    def __neg__(self):
        return PolyQuad(self.mods, {m: -c for m, c in self.terms.items()})

    def __sub__(self, other):
        return self + (-self._co(other))

    def __rsub__(self, other):
        return self._co(other) + (-self)

    def __mul__(self, other):
        other = self._co(other)
        raw = {}
        for m1, c1 in self.terms.items():
            for m2, c2 in other.terms.items():
                m = tuple(a + b for a, b in zip(m1, m2))
                raw[m] = raw.get(m, Fraction(0)) + c1 * c2
        # reduce x_i^2 -> d_i
        out = {}
        for m, c in raw.items():
            coeff = c
            red = []
            for e, d in zip(m, self.mods):
                coeff *= Fraction(d) ** (e // 2)
                red.append(e % 2)
            red = tuple(red)
            out[red] = out.get(red, Fraction(0)) + coeff
        return PolyQuad(self.mods, out)

    __rmul__ = __mul__

    def __eq__(self, other):
        other = self._co(other)
        return self.terms == other.terms

    def is_zero(self):
        return not self.terms

    def __hash__(self):
        return hash(tuple(sorted(self.terms.items())))

    def __repr__(self):
        if not self.terms:
            return "0"
        bits = []
        for m, c in sorted(self.terms.items()):
            names = "".join(
                f"*sqrt{d}" for e, d in zip(m, self.mods) if e
            )
            bits.append(f"{c}{names}")
        return " + ".join(bits)

    def coord_dict(self):
        """Map from exponent tuple to Fraction, for comparison across modules."""
        k = len(self.mods)
        return {m: self.terms.get(m, Fraction(0)) for m in product((0, 1), repeat=k)}
