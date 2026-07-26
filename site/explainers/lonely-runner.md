---
title: Lonely Runner Conjecture
short: Lonely runner
order: 5
tagline: Runners circle a track at different speeds. Does everyone eventually get some space to themselves?
posed: Jörg Wills, 1967; named by Luis Goddyn, 1998
---

## In plain terms

Several runners start together on a circular track of length 1 and run forever
at different constant speeds. Fix your attention on one of them. The conjecture
says that at some moment, that runner is **lonely**: at distance at least 1/k
from every other runner, where k is the total number of runners.

With two runners it is obvious — at some point the other one is on the far side
of the track. With seven or eight, it stops being obvious, and the bound 1/k is
exactly tight: with speeds 1, 2, …, k−1 the lonely moment happens, but with no
room to spare.

## What is known

Proved up to k = 7. A 2025 preprint by Rosenfeld claims k = 8, so that case is
likely settled pending review, leaving k ≥ 9 as the frontier.

The problem reduces neatly: speeds can be assumed to be whole numbers, one
runner can be parked at zero, and Tao showed it is enough to check speeds up to
a computable bound. So each k is, in principle, a finite calculation — just an
astronomically large one.

One belief about the problem turned out to be **false**. It was thought the
only tight configurations were the obvious speeds 1, 2, …, k−1; Goddyn and Wong
found others. This lab recovered those independently at k = 8, which is how it
checked its own search was working.

## Why it is hard

The finite bound is enormous, so brute force is hopeless beyond small k, and
the tight cases sit right on the boundary — you cannot get away with an
approximate argument, because the margin is exactly zero.

## What a breakthrough would mean

This one has the richest connections of any problem here. It is equivalent to a
question in **view obstruction** — whether a certain family of convex bodies
covers space — and it has been linked to **nowhere-zero flows** in graphs and
matroids, a genuinely applied corner of combinatorics. It also sits squarely in
Diophantine approximation, since it is really a statement about how badly
several numbers can be simultaneously approximated.

So a proof would not be an isolated curiosity. The techniques would be expected
to say something about the covering and flow questions it is tied to. That
breadth is unusual, and it is the honest reason to rate this problem above the
others here on implications.
