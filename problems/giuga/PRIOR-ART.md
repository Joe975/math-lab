# Giuga — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`.

## Attempts

**None yet.** Added 2026-07-27; the first attempts are in flight. Until one
lands there is no prior art to be informed by, so `blind` and `informed` mode
are equivalent here.

## Editorial view of the attack surface

Why this problem was added: it is the cheapest problem in the library to get a
real result from. The counterexample conditions are so restrictive that the
search is a prunable enumeration over sets of primes, in exact integer
arithmetic, with no external dependencies — and there is a published numeric
frontier to recover from scratch first, which is how this lab validates a
pipeline before trusting it.

The mechanism family (constrained prime-chain enumeration with divisibility
pruning) is not represented anywhere else in the library. The nearest neighbour
is the Erdős–Straus covering-congruence work, and the two share no tooling.

Concrete lines, if you want them:

- Recover the small Giuga numbers and the published factor-count bound with an
  independent implementation. Do this before anything else; a search whose
  calibration is untested is not evidence.
- Attack the pruning soundness explicitly: write down what each pruning rule
  discards and why it cannot contain a counterexample. This is where a wrong
  result would come from, and a clean statement of the rules is worth more than
  a bigger number.
- The Agoh/Bernoulli formulation is an independent second implementation for
  free on small n. Use it as the skeptic's re-implementation rather than
  re-running the same enumerator.
- Odd Giuga numbers, and the transfer of technique to Lehmer's totient problem,
  are natural adjacent targets if the enumerator works.
