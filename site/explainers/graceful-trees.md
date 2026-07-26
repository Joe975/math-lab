---
title: Graceful Tree Conjecture
short: Graceful trees
order: 6
tagline: Can every tree be numbered so that its edge differences hit every value exactly once?
posed: Gerhard Ringel and Anton Kotzig, 1960s
---

## In plain terms

Take a tree with n points. Label the points 0 through n−1, using each label
once. Each edge then gets a value: the difference between the two labels at its
ends. The labelling is **graceful** if those edge values are exactly
1, 2, …, n−1, each appearing once.

A path of four points labelled 0, 3, 1, 2 gives edge differences 3, 2, 1 — every
value once. Graceful. The conjecture says *every* tree admits such a labelling,
no matter how it branches.

## What is known

Verified by computer for all trees up to at least 35 points, and proved for many
families: paths, stars, caterpillars, everything of diameter at most 5.
Lobsters — trees that sit within two steps of a central path — are a well-known
open subcase.

Nobody expects a counterexample. The difficulty is that the known proofs are all
*constructive* and family-specific: they build a labelling by exploiting the
shape of that particular family, and there is no general construction.

## Why it is hard

The constraint is unusually rigid. Every one of the n−1 differences must be hit
exactly once, so the labelling is a near-perfect packing with no slack anywhere.
Local greedy choices fail because a decision at one end of the tree silently
forbids a value at the other. Counting arguments run into the opposite problem:
most trees have *enormous* numbers of graceful labellings, so averaging says
almost nothing about whether the count can ever hit zero.

## What a breakthrough would mean

This is the one problem here with a clean, important consequence. By a theorem
of Rosa, a graceful labelling of a tree with n edges yields a decomposition of
the complete graph on 2n+1 points into 2n+1 copies of that tree. So the graceful
tree conjecture **implies Ringel's conjecture** on graph decompositions.

Ringel's conjecture was itself settled asymptotically by Montgomery, Pokrovskiy
and Sudakov — a celebrated result — but by a method that works for large graphs
and does not produce the labellings. A proof of gracefulness would give the
explicit, all-sizes, constructive version.

Applications to coding theory, radar pulse design and crystallography are often
cited in the literature on graph labelling. They are real in the sense that
someone has written them down, but the practical pull is thin, and it would be
overselling to present them as motivation.
