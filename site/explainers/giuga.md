---
title: Giuga's Conjecture
short: Giuga
order: 9
tagline: One sum of powers should tell you whether a number is prime. Does it ever lie?
posed: Giuseppe Giuga, 1950
---

## In plain terms

Take a number n. Raise every smaller number to the power n−1, add them all up,
and divide by n. Giuga noticed that when n is prime the remainder is always
n−1 — the sum comes out at "−1" in clock arithmetic. He conjectured the reverse
too: **only** primes do this.

That would make one arithmetic test a perfect primality test. The easy half is
proved. The open half is that no composite number sneaks through.

## What is known

Nobody has found a cheat, and the constraints on what one would have to look
like are extraordinary. A composite counterexample must be simultaneously two
rare kinds of number — a Carmichael number and a Giuga number — and that forces
its prime factors into a tightly interlocking pattern: each prime p dividing n
must divide n/p − 1 exactly.

Those constraints are strong enough to search against. In 1996 Borwein,
Borwein, Borwein and Girgensohn showed any counterexample must have at least
3459 distinct prime factors, and therefore more than 13,000 digits. In 2013
Borwein, Maitland and Skerritt pushed the bound past **19,908 digits**. Nobody
tested numbers one at a time to get there — the searches enumerate possible
sets of prime factors and rule out whole branches at once.

There is a second face to the problem. Agoh's reformulation says the same thing
using Bernoulli numbers, a completely different piece of machinery. The two
statements are known to be equivalent, which is unusually convenient: a claim
can be checked two ways that share no arithmetic.

## Why it is hard

The bounds grow, but they never close. Every search establishes "no
counterexample smaller than this", and there is no known reason the answer
should be findable at any particular size — a counterexample, if one exists, is
under no obligation to be near the current frontier.

Meanwhile the proof side has no purchase. The condition is a statement about
every prime factor at once, and there is no known mechanism forcing that
pattern to be impossible rather than merely astronomically rare.

## What a breakthrough would mean

A proof would hand number theory a clean characterisation of primality by a
single congruence — the sum-shaped companion to Wilson's theorem, which does
the same job with a product. That is aesthetically satisfying rather than
practically useful: computing the sum is far slower than the primality tests
everyone actually uses.

A counterexample would be more interesting than a proof. It would be a specific
integer with thousands of prime factors, exhibiting a coincidence that all
current heuristics say should not happen — and heuristics that are wrong about
this would likely be wrong about neighbouring questions on Carmichael numbers
too.

For this lab the appeal is different and more practical: it is the problem here
where a modest computation, done carefully, produces a real result with a
stated range — and where the published frontier can be recovered from scratch
first, which is how any search here earns the right to be believed.
