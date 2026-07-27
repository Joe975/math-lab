# Giuga's Conjecture

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** Giuga (1950): for n > 1,

    n is prime  ⟺  1^(n−1) + 2^(n−1) + ... + (n−1)^(n−1) ≡ −1  (mod n)

The forward direction is elementary. The open direction is that **no composite
n satisfies the congruence**.

Agoh's equivalent formulation, in Bernoulli numbers: n is prime iff
n·B_{n−1} ≡ −1 (mod n). The two statements are known to be equivalent, which
gives the problem two computational faces that fail independently — useful,
because a claim checkable two ways is checkable by two implementations that
share no arithmetic.

## Published status

Open. The structure of a counterexample is heavily constrained:

- A composite n satisfies the congruence **iff it is both a Carmichael number
  and a Giuga number**. A *Giuga number* is a composite n with p | (n/p − 1)
  for every prime p | n; a Carmichael number is squarefree composite with
  (p−1) | (n−1) for every prime p | n. So a counterexample is squarefree and
  satisfies, for every prime p | n, both p | (n/p − 1) and (p−1) | (n/p − 1).
- Borwein, Borwein, Borwein and Girgensohn (1996) showed any counterexample
  has at least 3459 distinct prime factors and hence at least ~13,800 decimal
  digits.
- Borwein, Maitland and Skerritt (2013) improved the bound to at least 19,908
  decimal digits. Neither bound came from testing integers one at a time: both
  exploit the factor structure above to prune an enumeration over prime sets.
- Giuga numbers themselves are known and small ones are listed in OEIS
  A007850 (30, 858, 1722, 66198, …). Every known Giuga number is even; whether
  an odd one exists is open. The counterexample question is strictly harder,
  since it needs the Carmichael condition simultaneously.
- Lehmer's totient problem (does φ(n) | n−1 have a composite solution?) has the
  same shape — a constrained prime-chain search with published factor-count
  bounds — and results transfer between them in the literature.

**Sources.** Figures above come from encyclopaedia summaries and abstracts
rather than the primary papers; treat them as machine transcribed `[T]`. The
1996 and 2013 bounds are reported inconsistently across secondary sources
(13,800 vs 13,887 digits); pin the number from the paper before claiming to
have matched or beaten it.

## Verification contract

The published results here are statements about a **search**, not about a range
of integers, and that changes what has to be shown.

- **Exact integer arithmetic throughout.** No floating point anywhere in a
  divisibility or factor-count argument.
- **State the enumeration and prove the pruning sound.** A claim of the form
  "no counterexample has fewer than m prime factors" is worth nothing without
  an argument that every pruned branch was incapable of extending to a
  counterexample. A pruning bug produces a *stronger*-looking result and leaves
  no trace in the output — this is the single most likely way to be wrong here.
- **Recover the frontier before extending it.** The small Giuga numbers and the
  published factor-count bound are a calibration set: rediscover them from
  scratch with your own code, and say so, before claiming anything past them.
- **A verified congruence for a specific n is EVIDENCE about that n.** Say
  which n. A search over a family is EVIDENCE bounded by the family; state the
  family exactly, including how many prime factors and what size limit.
- **Cross-check across the two formulations.** A count computed from the power
  sum and from the Bernoulli/Agoh side should agree; when they can both be run,
  disagreement is a bug and agreement is real (if partial) evidence.

## Harness (tier 0)

None yet. A contributor adding one should put the Giuga/Carmichael condition
checkers and the prime-chain enumerator here — they verify the objects
themselves and are not specific to any route.
