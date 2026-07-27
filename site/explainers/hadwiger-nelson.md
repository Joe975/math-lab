---
title: Hadwiger–Nelson Problem
short: Hadwiger–Nelson
order: 8
tagline: Colour every point of the plane so no two points exactly one inch apart match. How many colours do you need?
posed: Edward Nelson, 1950; Hugo Hadwiger, 1945
---

## In plain terms

Imagine painting every single point of an infinite sheet of paper. The only
rule: any two points exactly one inch apart must be different colours. How few
colours can you get away with?

Seven is enough — tile the plane with hexagons slightly under an inch across
and repeat seven colours in a pattern. Four is not enough, which you can see
from a small arrangement of seven points called the Moser spindle. So the
answer is 4, 5, 6 or 7, and for sixty-eight years nobody could narrow it.

Then in 2018 an amateur — Aubrey de Grey, better known for research on ageing —
found an arrangement of 1581 points that cannot be coloured with four. The
answer is now known to be **5, 6 or 7**. Which one is still anybody's guess.

## What is known

The problem is secretly finite. A compactness theorem says the whole infinite
plane needs *k* colours exactly when some finite set of points does. So
progress comes from building finite point sets — "unit-distance graphs" — that
resist colouring, and every one of them is an object you can hand to a computer
and check.

That is what happened after 2018. Polymath16, an open online collaboration,
spent three years shrinking de Grey's 1581-point graph. The record is now 509
points, found by Jaan Parts. Every step was a concrete object anyone could
verify independently.

Nobody has found a point set that forces six colours, and nobody has found a
colouring of the plane using fewer than seven.

## Why it is hard

The two ends of the problem resist for opposite reasons. Pushing the lower
bound to 6 means finding a finite graph that is 6-chromatic — but the search
space is the set of all finite point configurations in the plane, which is not
something you can enumerate. Pushing the upper bound below 7 means inventing a
genuinely new colouring of the whole plane, and the hexagonal one has stood
since the 1950s.

There is also a subtlety that makes the problem stranger than it looks: for
some ways of asking the question, the answer depends on which axioms of set
theory you adopt. Restricting to "reasonable" colourings — ones whose colour
classes are measurable — provably needs more colours than the general case
might.

## What a breakthrough would mean

This is the most *visual* problem in this library, and the one where an
amateur most recently moved a famous constant — which is itself the interesting
part. De Grey's graph was found by computer search guided by human insight
about which structures to spindle together, and then verified by SAT solvers.
The whole episode is a case study in what computer-assisted mathematics
currently looks like: a human picks the shape of the search, a machine settles
it, and independent machines re-check the answer.

Settling χ = 5 versus 6 versus 7 would not obviously cascade into other fields
the way a proof about primes might. Its value is different: it is one of the
cleanest test cases for whether machine search can push a genuinely open
geometric question, which is exactly the question this lab exists to gather
evidence about.
