---
title: Collatz Conjecture
short: Collatz
order: 7
tagline: Halve it if even, triple-plus-one if odd. Does every number fall to 1?
posed: Lothar Collatz, 1937
---

## In plain terms

Pick a whole number. If it is even, halve it. If it is odd, triple it and add 1.
Repeat.

Starting at 6: 6, 3, 10, 5, 16, 8, 4, 2, 1. Starting at 27, it wanders up past
9,000 and takes 111 steps to come down. The conjecture says every starting
number eventually reaches 1.

It is the most famous problem in this collection, and the one this lab
deliberately spent the least effort on.

## What is known

Verified by computer to around 2⁶⁸. The strongest theoretical result is Terence
Tao's, from 2019: *almost all* starting values eventually reach a value that is
almost bounded — in a precise density sense. That is remarkably close to the
conjecture and still, frustratingly, not it.

## Why it is hard

There is a hard warning sign. If you generalise the rule slightly — allowing
other multipliers and divisors — the resulting family of problems is
**undecidable**: no algorithm can settle them all. Conway proved this. So any
successful method must be specific to 3n+1 and must break for its near
neighbours, which rules out every general technique in advance.

The heuristic picture makes it worse. On average an odd step multiplies by 3/2
and the following even steps divide by more, so orbits drift downward like a
random walk with negative bias — which strongly suggests the conjecture is true
and gives no purchase on proving it, because a single exceptional orbit is all
it takes.

Erdős's assessment is still quoted because it is still accurate: mathematics is
not yet ready for such problems.

## What a breakthrough would mean

Nothing practical. No application depends on this.

What it would mean is a method for controlling the interaction between
multiplication and addition on the integers — the thing that makes the
arithmetic of 3n+1 so opaque. That is a genuine gap in mathematics, and it is
why the problem retains serious attention despite being a magnet for cranks.

**In this lab it is a deliberate long shot with a minority budget**, kept for a
different reason: the space of failed Collatz approaches is large and unusually
well documented, which makes it the best available test of whether a library of
recorded dead ends is useful at all. It was queued and never worked, so the
honest status is that nothing has been tried here.
