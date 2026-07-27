---
title: Conway's 99-Graph Problem
short: Conway 99-graph
order: 10
tagline: Ninety-nine people, every friendship in exactly one triangle, every pair of strangers with exactly two friends in common. Possible?
posed: John Horton Conway, 2014 ($1000 prize)
---

## In plain terms

Picture 99 people. Any two who know each other have exactly one friend in
common. Any two who do not know each other have exactly two friends in common.
Can such a group exist?

That is the whole question. Conway put $1000 on it, and nobody has collected.
Either you build the network, or you prove it cannot be built.

The conditions force a lot immediately. Everyone must have exactly 14 friends,
and each person's 14 friends must pair off into 7 couples who know each other
and nobody else in that circle. Everything about the structure is pinned down
locally — and yet it is unknown whether the pieces can be assembled globally.

## What is known

The frustrating part: **all the usual reasons for impossibility do not apply.**
Networks like this one have well-developed accounting rules — eigenvalue
counts, integrality conditions, bounds relating the numbers involved — and this
parameter set passes every one of them. There is no cheap contradiction and no
known construction. That is precisely why it is still open.

Two similar networks do exist: a 9-person one built from arithmetic mod 9, and
a 243-person one found by Berlekamp, van Lint and Seidel. Two other sizes in
the same family, with 6273 and 494019 people, are open just like 99.

What has been narrowed is symmetry. If the 99-network exists it must be quite
lopsided: it cannot look the same from every person's point of view, and a
chain of results — Makhnev and Minakova, then Behbahani and Lam, then Cesarz
and Woldar in 2023 — has whittled its possible symmetry group down to a handful
of tiny options. A 2026 paper reports that throwing modern SAT solvers at the
problem directly does not settle it, and analyses why the encoding defeats
them.

## Why it is hard

The search space is finite but beyond astronomical: there are vastly more ways
to wire 99 people together than there are atoms in the observable universe, and
the local rules do not prune it fast enough. Exhaustive search is out.

The other route — proving impossibility — needs an obstruction nobody has
found. The standard toolkit was already applied and came back clean, so a proof
would need a genuinely new invariant.

## What a breakthrough would mean

Modest, honestly, and this problem is in the library for a different reason.
Resolving it would settle one square in the table of possible networks of this
kind, and the *method* would be the real prize: a new obstruction would likely
close several open entries at once, and an explicit construction would be a new
combinatorial object worth studying on its own.

What makes it a good subject for this lab is the shape rather than the stakes.
It is a single, sharply-posed, finite yes-or-no question where the honest
expected outcome of any attempt is a precisely-described dead end — which is
exactly the kind of record this library is built to keep.
