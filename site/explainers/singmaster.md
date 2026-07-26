---
title: Singmaster's Conjecture
short: Singmaster
order: 3
tagline: Does any number appear more than a handful of times in Pascal's triangle?
posed: David Singmaster, 1971
---

## In plain terms

Write out Pascal's triangle and look for repeats. Most numbers appear once or
twice. Some appear more: 120, 210, 1540, 7140 each appear six times.

The record holder is **3003**, which appears eight times. Singmaster conjectured
there is some fixed ceiling — that no number appears a thousand times, or a
million. Nobody knows whether even 3003's record has ever been beaten.

## What is known

The best general bound is due to Kane, and it grows with the size of the number
rather than staying fixed — so it does not settle the conjecture. In 2022,
Matomäki, Radziwiłł, Shao, Tao and Teräväinen proved it holds in the
**interior** of the triangle, using analytic number theory. The near-edge
region remains open.

There are infinitely many numbers appearing at least six times, coming from a
Fibonacci-parameterised family. So the ceiling, if it exists, is at least six.

## Why it is hard

Each way a number can repeat is a separate Diophantine equation — asking when
C(n, j) = C(m, k) — and each such equation defines an algebraic curve. Deep
results bound how many rational points a single curve can have. What
Singmaster's conjecture needs is a bound that is **uniform across the whole
infinite family of curves at once**, and that kind of uniformity is exactly
what current technology struggles to deliver.

## What a breakthrough would mean

A proof would very likely arrive as a by-product of progress on uniform bounds
for rational points on families of curves — which is a central problem in
Diophantine geometry with reach far beyond Pascal's triangle. In that sense
this conjecture is a readable test case for machinery that matters elsewhere,
which is a better reason to care about it than the triangle itself.

There is also a real, if unlikely, **discovery payoff** at the computational
end: finding a number with multiplicity 9 or more would refute nothing but
would instantly rewrite what everyone assumes. It is a well-defined search with
a clear target, which is why this lab ran a complete census rather than a
sampled one — a census can say "there is nothing here", where a search can only
say "we did not find one".
