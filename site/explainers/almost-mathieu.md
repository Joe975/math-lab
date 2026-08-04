---
title: The Critical Almost Mathieu Operator
short: Almost Mathieu
order: 11
tagline: An electron on a grid in a magnetic field. Is its energy spectrum missing a piece at every possible place?
posed: Mark Azbel / Douglas Hofstadter (model, 1964/1976); "Dry Ten Martini" named after Mark Kac's bet
---

## In plain terms

Put an electron on an infinite grid of atoms and switch on a magnetic field.
The fraction of a flux quantum threading each little square of the grid — call
it α — controls the electron's allowed energies, and the picture of those
energies against α is the Hofstadter butterfly, one of the most famous fractal
images in physics. Each rational α = p/q gives q separate energy bands; each
irrational α gives an infinitely fragmented set.

The conjecture is about the fragmentation. A deep bookkeeping theorem ("gap
labelling") lists every place where the spectrum is *allowed* to have a gap.
The Ten Martini problem asked whether the spectrum really is fragmented for
every irrational field — that was proved in 2009. The **Dry** Ten Martini
problem asks for more: is every gap on the allowed list *actually open*, with
none collapsed to zero width? In 2023–24 this was settled for every strength
of the model except one — the physically central, self-dual, critical
coupling, which is exactly the Hofstadter case. There it remains open for
every single irrational α.

## What is known

The spectrum is a Cantor set for all irrational fields (Avila–Jitomirskaya,
2009). All labelled gaps are open away from critical coupling (Avila–You–
Zhou). At critical coupling the total size of the spectrum is zero, its
fractal dimension is at most one half (Jitomirskaya–Krasovsky), and for
rational flux p/q everything is exactly computable: the spectrum is where a
single degree-q polynomial stays between −4 and 4, all q − 1 gaps are open
except that for even q the two middle bands kiss at zero energy, and Thouless
conjectured — with striking numerics, still unproven — that q times the total
bandwidth tends to a specific constant, 32 times Catalan's constant over π.

## Why it is hard

Critical coupling is a phase transition point. Below it the model behaves
like free motion, above it like localized motion, and the two regimes come
with two different proof technologies; at the critical point both machines
break at once. The known gap-opening arguments run along one side of the
transition and lose all their room exactly there. Worse, the quantities that
survive at criticality (zero measure, dimension bounds) are the ones that
cannot see individual gaps.

## What a breakthrough would mean

The butterfly is now measured in real systems — moiré graphene, cold atoms,
photonic lattices, superconducting qubit arrays — and the same spectrum
governs the critical temperature of superconducting wire networks through the
de Gennes–Alexander equations, where the field-dependent T_c literally traces
the butterfly's edge. Whether the idealized spectrum has zero-width gaps is
the mathematical backbone behind which of those measured gaps are robust.
Inside mathematics it would close the last case of a program that has driven
quasi-periodic spectral theory for forty years.

**In this lab it carries a low budget and long-shot framing, deliberately.**
Nobody here expects to prove a gap open for irrational α. What the lab can
produce is exact: certified band edges, gap widths and bandwidth sums at
rational flux, over the real cyclotomic field with no floating point in any
certificate, pushed along Fibonacci approximants of the golden mean where the
conjectured constants have falsifiable shadows. Nothing has been tried here
yet.
