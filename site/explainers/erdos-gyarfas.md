---
title: Erdős–Gyárfás Conjecture
short: Erdős–Gyárfás
order: 2
tagline: If every junction has at least three roads, must there be a closed loop whose length is a power of two?
posed: Paul Erdős and András Gyárfás, 1995
---

## In plain terms

Draw a network where every point has at least three lines coming out of it. The
conjecture says you can always find a cycle — a closed loop back to where you
started — whose length is a power of 2: 4, 8, 16, 32, and so on.

Take a cube. Every corner has three edges, and each face is a loop of length 4.
Conjecture satisfied. The claim is that you can never draw such a network,
however large or strange, without one of these loops appearing somewhere.

## What is known

It holds for various restricted families of graphs, and computational searches
have confirmed it for all small cases. The natural hard case is **cubic**
graphs — where every point has exactly three lines, with no slack to spare.

The powers of 2 thin out fast, which is what makes the conjecture delicate. As
graphs get larger the gaps between 16, 32, 64 grow, so there is more room for a
graph to have plenty of cycles and still miss every target length.

## Why it is hard

To be a counterexample, a graph must dodge *every* power of 2 at once. Avoiding
short cycles is easy — force high girth. But a graph with no small cycles has to
be large, and large graphs tend to contain long cycles of many lengths,
including the powers of 2 further up. Nobody has found a way to make those two
pressures cancel, and nobody has proved they cannot.

## What a breakthrough would mean

This is the one problem in this set where a **counterexample is genuinely
findable**. It would be a single finite graph. You could publish the thing
itself, and anyone could verify it in seconds. That is a different kind of
target from the others here, and it is why this problem gets a large share of
the compute budget in this lab.

A proof, by contrast, would say something general about which cycle lengths are
unavoidable once a graph is even slightly dense — a question interesting in its
own right, and one where the surrounding theory (cycle spectra: which sets of
lengths can occur at all?) is thinner than you would expect.

Applications outside mathematics: none that anyone has identified.
