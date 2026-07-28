---
title: Crouzeix's Conjecture
short: Crouzeix
order: 10
tagline: Feed a matrix into a polynomial. How much bigger can the answer get than the polynomial itself?
posed: Michel Crouzeix, 2004
---

## In plain terms

Take a square matrix A and a polynomial p, and form p(A) by substituting the
matrix for the variable. How large can p(A) be? There is an obvious region of
the complex plane attached to A — its *numerical range*, the set of all
quantities x*Ax over unit vectors x, which is a convex blob containing the
eigenvalues. If p is small everywhere on that blob, must p(A) be small too?

Not exactly: there is a factor to pay. Crouzeix conjectured the factor is 2 —
that p(A) can never exceed twice the largest size of p on the blob. Two is
certainly the right guess for the smallest possible factor, because a single
2×2 matrix already achieves it exactly: the matrix with a lone 1 above the
diagonal has size 1, while its numerical range is a disc of radius one half.

The conjecture says that little example is as bad as it ever gets, in any
dimension, for any polynomial.

## What is known

The factor is somewhere between 2 and 1 + √2 ≈ 2.414. Crouzeix proved 11.08 in
2007; Crouzeix and Palencia brought it down to 1 + √2 in 2017, and that has
stood since. Nobody has closed the remaining gap of about four tenths.

Special cases fall regularly. It is settled for 2×2 matrices, for normal
matrices, for matrices whose numerical range happens to be a disc, for matrices
close to a single Jordan block, for certain 3×3 tridiagonal matrices, and for
weighted shifts. Recent work has produced constants better than 1 + √2 for each
fixed dimension separately, and has reformulated the general conjecture as a
statement about a single differentiation operator — so the whole problem is
equivalent to one concrete question about analytic functions.

Numerically it has been hammered on. Greenbaum and Overton searched for a
counterexample by nonsmooth optimization over matrices and polynomials at once,
and found none: every search path they followed climbed to a ratio of 2 and
stopped at structures already understood.

## Why it is hard

The numerical range is a crude object. It knows the eigenvalues and something
about how far the matrix is from being normal, but it throws away most of the
matrix, and the conjecture asks it to control something — the size of p(A) —
that depends on the discarded detail. Proving that a small amount of
information suffices is a different kind of task from computing more
accurately.

There is also a shortage of hard examples. The known cases where the ratio
reaches 2 are all small and all closely related, so the conjecture has no rich
supply of near-misses to learn from, and the numerical searches keep returning
the same family. That is comforting if you believe the conjecture and useless if
you want to prove it.

And unlike most problems in this collection, the quantities in question are not
rational. Both the size of p(A) and the maximum of p over the numerical range
come from eigenvalue computations, so a computation that is merely accurate
proves nothing at all — it has to come with a proof of its own error bound.

## What a breakthrough would mean

The practical consequence is in numerical analysis. Bounds of exactly this shape
control how errors grow when iterative methods are applied to non-normal
matrices — the situation where the usual eigenvalue intuition fails and
computations misbehave for reasons that are hard to predict. A clean constant of
2 would make a family of convergence estimates sharp instead of approximate.

There is a physical reading too. The numerical range is the set of expectation
values an observable can take, so the conjecture is a statement about how far a
function of a non-normal operator can stray from the values the operator is
seen to take. Non-normality is exactly what distinguishes dissipative and open
quantum systems from the textbook self-adjoint ones.

**In this lab it carries a medium budget and a stated risk.** Every claim needs
a certified error bound built from scratch, which is more machinery than the
other problems here require before any mathematics happens. If that machinery
turns out to cost too much, the plan of record is to narrow to families where
the answer has a closed form and say so, rather than quietly switching to
uncertified numbers. Nothing has been tried here yet.
