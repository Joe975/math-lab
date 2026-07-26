---
title: Erdős–Straus Conjecture
short: Erdős–Straus
order: 4
tagline: Can 4/n always be written as a sum of exactly three unit fractions?
posed: Paul Erdős and Ernst Straus, 1948
---

## In plain terms

The ancient Egyptians wrote fractions as sums of distinct unit fractions —
1/2, 1/3, 1/7 and so on. The conjecture says that for every whole number
n ≥ 2, the fraction 4/n can be written as a sum of exactly **three** of them.

For n = 5: 4/5 = 1/2 + 1/4 + 1/20. For most n it is easy, and there are simple
formulas that handle whole families of n at once.

## What is known

It has been verified by computer far past 10¹⁷ — nobody seriously doubts it is
true. It is enough to prove it for prime n, since composite cases follow.

The classical approach is **identities**: algebraic formulas that solve 4/n for
every n in some residue class. A dozen such families together handle 834 of the
840 classes modulo 840. The six that resist are exactly the perfect squares
among the units mod 840 — a fact traced to Mordell, and not a coincidence.

## Why it is hard

Mordell's obstruction is the crux, and it is a genuine barrier rather than a
gap in effort. Any solution built from finitely many polynomial identities
provably **cannot** cover those six classes. The identity method does not
merely fail there; it cannot succeed. So a proof needs a new mechanism for
those residues, and the classes that resist are precisely the ones defined by
quadratic-residue conditions — the same structure that makes so many questions
about primes hard.

## What a breakthrough would mean

Very little outside number theory, and it is worth being blunt about that.
There is no cryptographic angle here and no application waiting.

The value is in the mechanism. This is a problem where we know exactly which
tool fails and why, so a proof must introduce something genuinely new for
handling quadratic-residue obstructions in covering arguments. That machinery
would likely transfer to the wider family of Erdős–Graham questions about
Egyptian fractions.

There is also a quieter question the search data speaks to directly: not
*whether* 4/n has a representation, but **how many** it has, and why some
primes have so few. That statistical question is tractable now, and it is where
this lab's work on the problem actually landed.
