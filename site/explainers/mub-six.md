---
title: Mutually Unbiased Bases in Dimension Six
short: MUBs in d = 6
order: 17
tagline: A qubit has three perfectly "incompatible" ways to be measured. Does a qubit-plus-qutrit have four?
posed: Zauner, 1999 (in this form); the d = 6 question dates to the 1980s
---

## In plain terms

Measuring a quantum system means picking a *basis* — a set of perfectly
distinguishable outcomes. Two bases are **mutually unbiased** if knowing the
outcome in one tells you nothing at all about the other: every outcome is
equally likely. For a single qubit the three Pauli measurements X, Y and Z are
mutually unbiased; on the Bloch sphere they are three perpendicular axes.

In dimension d there can be at most d + 1 such bases. When d is a prime power
(2, 3, 4, 5, 7, 8, 9, …) that maximum is always reached. The first dimension
that is not a prime power is 6 — a qubit and a qutrit side by side — and there
nobody can find even **four** mutually unbiased bases. The conjecture is that
three is the most.

## What is known

Three bases are easy to build in dimension 6. Large numerical searches since
2007 have never found a fourth, and the best four bases anyone has found fall
just short: their average "distance" is about 0.9983 where perfect unbias
would be 1. Several families of candidate bases have been ruled out by proof,
but the space of candidates (6 × 6 complex Hadamard matrices) has never been
fully classified.

## Why it is hard

The objects are small — three 6 × 6 matrices of complex numbers — but the
space they live in is continuous and not classified, so there is no finite
list to check. Numerical searches can show that four bases are hard to find,
never that they cannot exist.

## What a breakthrough would mean

MU bases are the optimal measurement schemes for **quantum state tomography**
and the backbone of several **quantum key distribution** protocols, so the
answer decides how efficiently a qubit–qutrit system can be characterised.
Mathematically, a proof would be the first case showing that prime-power
dimensions are genuinely special for this problem, with knock-on effects for
related open questions: complete sets of equiangular lines (SIC-POVMs) and
the classification of complex Hadamard matrices. There is also a suggestive
parallel with finite geometry — complete MU sets in prime-power dimensions
behave like affine planes, and it is a classical theorem that no affine plane
of order 6 exists — but that parallel is only an analogy, not a proof
strategy anyone has made work. Finding four bases would be a genuine surprise and
would overturn nearly two decades of numerical consensus.
