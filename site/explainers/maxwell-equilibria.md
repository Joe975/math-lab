---
title: Maxwell's problem on points of equilibrium
short: Maxwell equilibria
order: 11
tagline: Fix a handful of electric charges in space. How many places can a test charge sit with no force on it?
posed: James Clerk Maxwell, 1873
---

## In plain terms

Pin a few positive electric charges at fixed points in space. Between them
the electric field pushes and pulls, and at certain special points the pushes
cancel exactly: a test charge placed there feels nothing. These are the
points of equilibrium. A famous theorem of Earnshaw says none of them is
stable — every one is a saddle, a mountain pass rather than a valley — which
is why you cannot levitate a charge with static fields alone.

Maxwell, in a footnote of his 1873 *Treatise on Electricity and Magnetism*,
claimed that n charges can never produce more than (n−1)² such points. Two
charges give at most one equilibrium; three should give at most four; five at
most sixteen. He sketched an argument and moved on. A century and a half
later, nobody had proved it — and for three charges, nobody has even proved
that the answer is four rather than five or six.

## What is known

Remarkably little, for how old and concrete the question is. It is not even
known that the number of equilibrium points is always finite. The modern
benchmark is a 2007 paper of Gabrielov, Novikov and Shapiro, which used deep
results about systems with few monomials to prove the first general bounds —
for three charges, at most twelve equilibria, against Maxwell's predicted
four. For three charges of *equal* strength the answer really is four
(Tsai, 2015), with the equilateral triangle attaining it: one equilibrium at
the center, three more near the edges.

Then, days before this problem was onboarded, the story broke open: a July
2026 preprint of Arathoon, Ball and Kvalheim constructs five charges with at
least twenty-four equilibria — eight more than Maxwell's sixteen. The
construction is delicate: an equilateral triangle of unit charges plus two
minuscule charges hovering just above and below its center. As the small
charges switch on, the central equilibrium shatters into twenty-one. If the
preprint holds up, the conjecture as Maxwell stated it is false, and the
true growth of the maximum is wide open. A companion preprint sharpens the
three-charge bound from twelve to six — leaving the original small case,
four versus five versus six, still unresolved.

## Why it is hard

The equilibrium equations look tame but are not algebraic — distances enter
through square roots — so the polynomial machinery that counts solutions of
polynomial systems does not apply directly, and what replaces it
(fewnomial theory) gives bounds that are almost certainly far from sharp.
The saddle points also interact: Morse theory constrains their indices
globally, so you cannot add equilibria one at a time; new ones must appear
in coordinated families, which is exactly what makes both counting them and
constructing many of them delicate.

And the new counterexample illustrates a second difficulty: it exists only
"for sufficiently small ε". Such proofs are rigorous as asymptotics yet name
no actual configuration, and pinning one down is genuinely hard — the
twenty-one new equilibria are born nearly merged, so any explicit instance
sits close to degeneracy, where both numerics and certification are at their
worst.

## What a breakthrough would mean

The equilibrium structure of Coulomb fields is textbook physics — ion traps,
charged-particle optics and levitation arguments all live in it — and it is
genuinely startling that the count of force-free points for five charges was
mispredicted for 150 years. Settling the three-charge case would close the
first open case of a question Maxwell thought obvious. An independently
certified explicit counterexample would turn a fresh asymptotic claim into a
concrete, checkable object. And any new record configuration would redraw
the map of what static charge arrangements can do.

**In this lab it carries a high budget.** Everything the contract needs is
decidable in exact rational arithmetic once the right variables are chosen,
the first targets are concrete and self-policing — a single certified
configuration with five equilibria from three charges, or twenty-four from
five, decides something — and the fresh refutation means the terrain has
just shifted, which is when independent certified searches earn the most.
Nothing has been tried here yet.
