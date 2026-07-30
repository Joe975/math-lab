---
title: Mahler's Volume-Product Conjecture
short: Mahler in 4D
order: 9
tagline: A convex shape and its dual have volumes that trade off. Which shape makes the product smallest?
posed: Kurt Mahler, 1939
---

## In plain terms

Take a convex shape that looks the same upside down — a cube, a ball, a
diamond. Every such shape has a *dual*: stretch the original in one direction
and its dual squashes in that same direction, so the two volumes pull against
each other. Multiply them together and the stretching cancels out. What is left
is a single number attached to the shape, the same for a cube and for any
squashed or sheared copy of it.

Which shape makes that number as small as possible? The ball makes it *large* —
that half has been known since the 1930s. The small end is the open question.
Mahler's guess was the cube, and the diamond, and a family of shapes built by
gluing cubes and diamonds together, all of which give exactly the same number:
4ⁿ divided by n factorial.

It has been proved in the plane, and in three dimensions. Four is open.

## What is known

The two-dimensional case is Mahler's own, from 1939. Three dimensions took
eighty years: Iriyeh and Shibata proved it in 2019, for symmetric bodies,
equality case included, and a shorter proof followed.

Everything else is a bound with the wrong shape. Bourgain and Milman proved in
1987 that you cannot go below the conjectured value by more than an exponential
factor, and Kuperberg made that factor explicit in 2008 — roughly (π/4)ⁿ. That
is enormously strong compared to nothing and enormously weaker than the
conjecture, and the gap grows with dimension rather than closing.

The local picture is completely settled and completely unhelpful for finding
counterexamples. Nazarov, Petrov, Ryabogin and Zvavitch proved the cube is a
strict local minimum in 2010, and Jaegil Kim extended that to the whole
conjectured family in 2014. So nothing near the conjectured answers can beat
them; any counterexample has to live somewhere else entirely.

There is a physics thread too. A conjecture of Viterbo about symplectic
capacities — the invariants that govern how a region of phase space can be
deformed by Hamiltonian mechanics — was known to imply Mahler's. In 2024
Haim-Kislev and Ostrover refuted Viterbo's conjecture. The counterexample is not
symmetric, and a symmetric version strong enough to imply Mahler is still
standing, so Mahler survived the collapse of its physical motivation.

## Why it is hard

The conjectured minimum is not unique. The cube and the diamond are not
linearly equivalent and give the same value, as does every way of gluing them,
and in high dimensions there are many. A minimum attained by a large family is a
minimum with no strict inequality to exploit — most optimization arguments want
a unique extremal and there is not one.

The known proofs also do not stack. The three-dimensional argument rests on
cutting a body into equal pieces by hyperplanes in a way that is essentially
three-dimensional; it does not simply run again with one more coordinate. And
the class of convex bodies is infinite-dimensional, so no finite computation can
close the question — a computer search can only ever say that nothing was found
inside the box it searched.

## What a breakthrough would mean

The volume product is one of the basic quantities of convex geometry, and the
conjecture is the last missing piece of a picture whose other half — that the
ball maximizes it — has been settled for a century. Settling four dimensions
would say whether the three-dimensional proof was a method or a one-off, which
is the question the field actually cares about.

The wider stake is the reverse Santaló inequality, which is used across
functional analysis and asymptotic geometry and is currently only known up to
an exponential constant. Pinning that constant down changes what those
applications can assert.

**In this lab it is the cleanest exact-arithmetic target available.** For a
polytope with rational vertices, the dual, the volumes and the comparison
against the conjectured value are all rational computations — no numerical
tolerance, no floating point, no certification machinery needed before the work
can start. Whether that reaches far enough to be interesting is exactly what the
first attempts are for. Nothing has been tried here yet.
