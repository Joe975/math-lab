---
title: Zaremba's Conjecture
short: Zaremba
order: 8
tagline: Every whole number should be the denominator of a fraction whose continued fraction never needs a digit bigger than 5.
posed: Stanisław Zaremba, 1971
---

## In plain terms

Any fraction m/n can be unpacked into a **continued fraction** — a tower of
nested divisions written [0; a₁, a₂, …], where the aᵢ are whole numbers called
partial quotients. Some fractions unpack gently, using only small digits;
others need a huge digit somewhere in the tower.

Zaremba's conjecture says: pick any denominator n whatsoever, and you can
always find a numerator m (sharing no factor with n) so that m/n unpacks using
digits **no larger than 5**. Every n, no exceptions. The 5 is sharp for at
least one case: with denominator 6, the best you can do is 5/6 = [0; 1, 5].

## What is known

Almost everything short of the conjecture itself. Bourgain and Kontorovich
proved in 2014 that the set of denominators that work has **density one** —
allowing digits up to 50. Huang then brought the digit bound down to the
conjectured 5: almost every n works exactly as Zaremba said. But "almost
every" leaves room for infinitely many exceptions, and closing that last gap
is the open problem.

Special families are settled completely: Niederreiter proved powers of 2 and 3
work with digits at most 3, and powers of 5 with digits at most 4.

A tempting stronger belief turned out to be **false**: Hensley conjectured
that any digit alphabet whose fractal of continued fractions is "big enough"
(dimension above ½) eventually captures every denominator. The alphabet
{2, 4, 6, 8, 10} is big enough by that measure, yet provably never produces a
denominator that is 3 more than a multiple of 4.

## Why it is hard

The fractions with bounded digits form a thin, fractal-like set, and asking
which denominators they hit is a **local–global problem**: the set obeys no
congruence obstruction, is large by every statistical measure, and yet nobody
can rule out sparse, structureless exceptions. The density-one methods count;
they cannot point at a specific n and certify it. Meanwhile the question is
completely concrete — for any single n, a finite computation settles it —
which makes the gap between "all but a vanishing fraction" and "all" feel
especially stark.

## What a breakthrough would mean

The conjecture began as applied mathematics: denominators with small partial
quotients are exactly the good moduli for lattice-based numerical integration,
so a proof would certify that an essentially optimal integration rule exists
at every resolution. Beyond that, it is the simplest clean test case for the
local–global philosophy on thin orbits — the same circle of ideas as the
Apollonian gasket — so techniques that close the last gap here would be
expected to travel.
