---
title: Union-Closed Sets Conjecture
short: Union-closed sets
order: 1
tagline: Combine any two sets in your collection and you stay inside it. Must some single item then appear in half of them?
posed: Péter Frankl, 1979
---

## In plain terms

Take a collection of sets that is *union-closed*: whenever two of its sets are
combined, the result is also in the collection. The conjecture says some single
element must appear in at least half of the sets.

A small example. Take `{}`, `{1}`, `{1,2}`. Combining any two gives you
something already there, so it is union-closed. Element 1 appears in two of the
three sets — comfortably at least half.

It sounds like it should fall to a clever counting argument. It has resisted
for over forty years.

## What is known

Nothing at all was known about a constant fraction until 2022, when Justin
Gilmer proved that some element appears in at least ~1% of the sets, using an
argument from **information theory** rather than combinatorics. That was a
genuine surprise: entropy is not the tool anyone expected to crack this.

Within months, a series of papers sharpened the constant to
(3−√5)/2 ≈ 0.38197 — and then showed that the entropy method *cannot go past
that number* without a new idea. The current published record, 0.38271, comes
from a refinement by Liu.

## Why it is hard

The entropy argument uses union-closure in only one way: it assumes the
combination of two *typical* sets stays inside the collection. That is an
average-case fact, and Chase and Lovett built collections that are
*approximately* union-closed, satisfy that average-case fact, and still fail
the conjecture at 0.38197. So any advance has to use exact closure — the fact
that even weird, atypical, heavily-overlapping pairs stay inside — and nobody
knows how to charge for that.

## What a breakthrough would mean

Honestly: nothing outside mathematics. There is no application waiting on this.

Inside mathematics the interest is real, and it is about **method**. Gilmer's
proof was a demonstration that information-theoretic arguments can attack
extremal set theory, and the barrier at 0.38197 is unusually well understood —
we know precisely which step fails and why. A proof reaching 1/2 would almost
certainly carry a new technique for exploiting closure conditions, and that
technique would be the actual prize.

This problem is also unusual in being *quantitative*. Most famous conjectures
are open or closed. This one has a number attached that has moved from 0.01 to
0.38 in three years, which makes partial progress measurable — and makes it a
good target for a machine loop.
