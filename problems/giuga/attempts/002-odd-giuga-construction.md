# 002 — Pushing at the objects: odd Giuga numbers and Giuga sequences

- **Problem:** giuga, `problems/giuga/PROBLEM.md`
- **Date:** 2026-07-27
- **Mode:** blind
- **Source-commit:** `e5af9e6285ea9af0b5dba21adafe402d588df83d`
- **Type:** object-shaped construction + exhaustive enumeration
- **Tools:** all pure-Python standard library; no C compiler on this machine.
  New tier 0: `harness/giuga/factorint.py` (deterministic Pollard-Brent with an
  explicit give-up), `harness/giuga/lasttwo.py` (the last-two-entries step, two
  routes), `harness/giuga/sequences.py` (Giuga sequences over integers), plus
  rule P8 and `giuga_parity_ok` added to the existing modules. New tier 1:
  `problems/giuga/explore/{odd_giuga,odd_sequences,odd_giuga_construct,giuga_by_factor_count}.py`.
  Everything is deterministic — the construction search is seeded, and
  Pollard-Brent starts at c = 1, 2, 3, … rather than at random, so the same
  branches fail on every run.
- **Sources:** `problems/giuga/PROBLEM.md`, and attempt 001 in this repo (my own
  record, same cycle). No papers, no OEIS lookups.

## Approach

**Is the target an object or a proof step?** This problem splits, and the split
is the whole reason this attempt looks the way it does.

A Giuga *counterexample* is object-shaped in principle — it is one integer, and
`is_counterexample` settles it in microseconds. But it is provably enormous:
attempt 001's own relaxation puts it past 3355 prime factors and 10^10395, and
the transcribed published bound is 19,908 digits. So there is no small object to
find, and persistence has nothing to grip: pushing harder at a search whose
smallest candidate has twenty thousand digits does not produce a candidate, it
produces a bigger bound. That is a **proof-shaped** target wearing object
clothing, and 001 was right to treat it as one.

Two nearby targets are genuinely object-shaped **and not known to be huge**:

1. **An odd Giuga number.** Open. Nobody has exhibited one and nobody has ruled
   one out, and no size bound puts it out of reach. A construction would be a
   real result.
2. **Giuga numbers and Giuga sequences past the known list.** 001 stalled at 7
   prime factors; 8 was reachable in principle and simply had not been reached.

Those got the push treatment: **nine rounds, eight of them past a partial
result**, each driven by a different structural idea rather than by running the
same search for longer. The checker existed before any of it —
`is_giuga_number` (pointwise divisibility on the integer) and `verify_sequence`
(exact rational identity) come from 001 and are covered by
`tests/test_giuga_harness.py`, so every object any round could produce is
settled independently of the code that found it.

## What was done

### Round 1 — point the existing enumerator at odd Giuga numbers

Baseline, no new idea: 001's prime-set recursion with `min_prime=3`.

    python problems/giuga/explore/odd_giuga.py --min-factors 9 --max-factors 13

The reciprocal sum settles the entry point with no search at all: the eight
smallest odd primes sum to 0.998956 < 1, so **an odd Giuga number has at least 9
prime factors**. Then m = 9 (24 nodes), 10 (94 nodes) and 11 (1265 nodes) come
back empty and exhaustive; m = 12 aborts, needing candidate primes up to
4,000,981,362 for the last-two-primes scan. The same wall 001 hit at 8 factors,
one problem over.

### Round 2 — a closed form for the last two entries

The wall is the scan window, so remove the scan. Two facts, neither using
primality:

**The denominator is always N.** Write rem = k − Σ_{i≤j} 1/n_i over the common
denominator N = n_1···n_j. If the branch completes with any r ≥ 1 further
entries x_1..x_r of product P, clearing denominators in
Σ 1/x_i − 1/(NP) = a/b gives

    N · (a·P − b·Σ_i (P/x_i)) = −b,   so N | b;

and b | N always, so **b = N at every node of a surviving branch**.

**The rectangle.** With b = N the last-two equation Nq + Np − 1 = apq becomes
apq − Np − Nq + 1 = 0, and multiplying by a completes the square:

    (a·p − N) · (a·q − N)  =  N² − a.

So the final pair comes from a **divisor pair of N² − a** — one factorisation
instead of a window of width N/a. Both signs are tried, since the two factors
need only share a sign.

Checked before use: the identity holds on all nine known Giuga numbers, and
`rem.denominator == N` holds at every prefix of every one of them. Then the
decisive cross-check — the scan (P7) and the rectangle (P8) are independently
written and **agree exactly on every m ≤ 7**. That is the two-implementation
requirement met on the enumeration itself, not merely on a number it reports.

Completeness now depends on N² − a factoring. `factorize_full` returns the
cofactor it could not split rather than pretending, and every caller records
that branch as *unresolved* — silently dropping it would turn a partial
enumeration into a false "none exist". That is the one way this attempt could
produce a wrong negative, so it is the thing that is instrumented.

### Round 3 — a different object class: Giuga sequences over integers

A *Giuga sequence* is distinct integers 2 ≤ n_1 < … < n_m with
Σ 1/n_i − 1/(n_1···n_m) ∈ ℤ. An all-prime Giuga sequence is exactly a Giuga
number's factorisation, so this is a strictly wider search — and it matters for
the odd question, because composite odd entries (9, 15, 21, …) let the
reciprocal sum reach 1 with fewer terms than odd *primes* can.

`harness/giuga/sequences.py`. Complete lists:

| length | Giuga sequences | of which all-prime |
|---|---|---|
| 3 | 1 | 1 — (2,3,5) |
| 4 | 2 | 2 |
| 5 | 3 | 1 |
| 6 | 17 | 2 |

The all-prime members reproduce the Giuga numbers exactly at every length — a
**third** independent enumerator (integers, `b = N`, rectangle) agreeing with
the sieve and the prime recursion from 001. The other twenty are new objects
here, e.g. (2,3,7,83,85) and (2,3,7,43,1805); note 1805 = 5·19² is not even
squarefree, so the last entry of a Giuga sequence is genuinely unconstrained.

### Round 4 — pairwise coprimality, from the same identity

`b = N` was derived above as a prune. It has a much stronger reading: b = N holds
**exactly when the chosen entries are pairwise coprime**. If a prime ℓ divides
two of them then every N/n_i still carries ℓ, so ℓ divides the numerator
kN − Σ N/n_i; if ℓ divides only n_j the numerator is −N/n_j mod ℓ, coprime to ℓ.
Therefore

> **all but the last entry of a Giuga sequence are pairwise coprime.**

Verified on all 23 known sequences before being used. As a prune it is decisive:

| odd length | nodes before | nodes after |
|---|---|---|
| 10 | 9,348 | 108 |
| 11 | 6,769,749 | 1,375 |

Same answers, ~4,900× fewer nodes at length 11.

### Round 5 — stop bounding, start constructing

Exhaustive search over odd prime sets stalls around 12 factors. This round gives
up completeness *between* prefixes to keep it *within* each one, and so reaches
prefix lengths the tree will never see.

For a prefix with product N and rem = a/N, the pair satisfies
q = (Np − 1)/(ap − N), and q > p forces N/a < p ≤ 2N/a. Walking every integer of
that window **settles the prefix completely with no factorisation at all**; the
window has width N/a, so a prefix is usable exactly when its reciprocal sum is
not too close to 1. Validated first: the window scan recovers the last two
primes of all nine known Giuga numbers, and `harness/giuga/lasttwo.py` confirms
that the window route and the rectangle route return identical pair lists on
each of them.

`problems/giuga/explore/odd_giuga_construct.py` samples odd-prime prefixes
(seeded: a random subset of the odd primes below 60, then greedy extension with
a skip probability) and fully settles every prefix whose window fits the cap.
Three seeds:

| seed | trials | prefixes settled | window integers | found |
|---|---|---|---|---|
| 2 | 30,000 | 1,099,847 | 1,059,471,340 | 0 |
| 11 | 40,000 | 1,470,950 | 507,612,756 | 0 |
| 12 | 40,000 | 1,540,836 | 284,816,592 | 0 |
| **total** | | **4,111,633** | **1,851,900,688** | **0** |

Factor counts 3 through 61 in every seed, peaking around 17. Prefixes are
deduplicated *within* a seed but not between them, so the total counts
settlements, not distinct prefixes.

### Round 6 — a parity law, which makes half the search vacuous

Let n be an odd Giuga number with m prime factors. Every n/p_i is odd, so
Σ_i n/p_i ≡ m (mod 2). The condition is Σ_i n/p_i − 1 = k·n, and n odd gives
k·n ≡ k (mod 2). Hence

> **m ≡ 1 + k (mod 2).**

The reciprocals of the first 1411 odd primes sum to less than 2, so k = 1 is
forced for every m ≤ 1411, and in that range **an odd Giuga number has an even
number of prime factors**. With m ≥ 9 that gives **m ≥ 10 and even**. The same
argument applies verbatim to odd Giuga *sequences*, since N/n_i is odd whether
or not n_i is prime.

This is retrospective self-criticism as much as a result: Rounds 1 and 3 spent
real compute exhausting m = 9 and m = 11, which a two-line parity argument
excludes for free. Those searches were correct and their answers stand; they
were simply not evidence about anything. `giuga_parity_ok` is now in the harness
with a test.

### Round 7 — pick the cheaper last-two route per prefix

Neither route dominates: a prefix far from sum 1 has a narrow window and is best
walked, one very close to 1 has a huge window but its N² − a may factor
instantly. `harness/giuga/lasttwo.py` chooses per prefix and reports whether the
answer is complete, so a caller can never mistake "gave up" for "none exist".
Odd sequences at length 11 went from 8.68 s to 0.42 s.

### Round 8 — where the fertile prefixes are, and why the odd case has none

The rectangle is easiest when a = 1: it becomes
(p − N)(q − N) = N² − 1 = (N−1)(N+1), which is divisor-rich and imposes no
integrality filter, so a = 1 prefixes are where completions actually come from.
a = 1 says Σ_i N/p_i + 1 = N, i.e. Σ_i 1/p_i + 1/N = 1.

An exhaustive scan of every squarefree N ≤ 2·10⁷ finds exactly five:

    2,  6 = 2·3,  42 = 2·3·7,  1806 = 2·3·7·43,  47058 = 2·3·11·23·31

— and **every one of them is even**. Cross-checked by a second implementation
that shares no arithmetic with the first: `problems/giuga/explore/a1_prefixes.py`
builds prime sets recursively instead of testing integers, and returns the same
five values at the same limit. (47058 is precisely the prefix that
produces 2214408306 = 47058 · 47057, via the r = 1 rule p = (N−1)/a.)

A caution against a claim I nearly made: "every known Giuga number has an a = 1
prefix" is *vacuous*, because the one-element prefix (2) already has a = 1, as
does (2,3). The content is the asymmetry, not the presence — the even case gets
a = 1 for free from its first two entries, while an odd a = 1 prefix would be an
odd solution of Σ_{p|N} 1/p + 1/N = 1, and there is none below 2·10⁷. The odd
search never gets to stand on a fertile prefix at all: (3) gives a = 2, (3,5)
gives a = 7, and it climbs from there.

### Round 9 — find out what the enumeration is actually spending its time on

Two targets refused to finish: the 8-prime-factor Giuga list, and odd sequences
of length 12. Rather than keep waiting, profile.

The 8-factor tree walked 200,000 nodes in **829 s with the last-two step
switched off entirely** — so the cost was not the closed form. `cProfile` on the
first 40,000 nodes showed 13 µs/node and no obvious culprit, which is the tell
that the cost is concentrated in the *deep* nodes. It was:
`_Primes.upto(limit)` returned `self._ps[:i]`, **a fresh copy of the prime list
on every node** — measured at 2.54 ms per call at limit 8·10⁶ (539,777 primes),
which matches the observed 4.1 ms/node almost exactly. Replaced by
`_Primes.between(low, high)`, which walks the shared list by index and copies
nothing. Output is unchanged (the pinned lists in the test suite confirm it);
only the cost differs. Note this makes the *timings* quoted in attempt 001 §3
stale — its node counts are unaffected, and per this repo's rules that record
stands as written.

That fix was not enough for either target, and the measurements say why:

| target | status | where it stopped |
|---|---|---|
| Giuga numbers, 8 prime factors | **stopped, unfinished** | ~60 min wall; reached >200,000 nodes and >52,000 last-two subproblems, no checkpoint written, no bound on what remained |
| odd Giuga sequences, length 12 | **stopped, unfinished** | 948 s CPU / 16 min wall, no checkpoint written; a separate capped probe reached 52,639 nodes in 240 s (219 nodes/s) with 1,944 branches already unresolved at rho budget 8,000 |

For scale, length 11 finished in 1,375 nodes and 0.42 s. Length 12 is at least
38× larger and ~15× slower per node, *and* it accumulates unresolved branches —
so even completing it would have yielded "complete except for some thousands of
unfactored N² − a", not an exhaustion. Both rows are reported as unfinished, not
as evidence of absence.

The runs that *did* checkpoint are on disk and were taken before their processes
were stopped: `odd-giuga-search.json` (odd prime sets, m = 9, 10, 11),
`odd-sequences.json` (odd sequences, lengths 7–11), `giuga-small.json` (Giuga
numbers, m = 3–7, all EXHAUSTIVE) and `odd-construct-seed{2,11,12}.json`. Every
figure quoted in this record comes from one of those files or from a run
reproduced in full above.

## Outcome

`VERIFIED` — **an odd Giuga number has at least 9 prime factors**, from exact
rational arithmetic alone (the eight smallest odd primes have reciprocal sum
0.998956 < 1). No search involved.

`VERIFIED` — **for m ≤ 1411, an odd Giuga number has an even number of prime
factors** (m ≡ 1 + k mod 2, with k = 1 forced in that range). Proof above,
mechanised as `giuga_parity_ok`, tested. Combined with the previous line:
**m ≥ 10 and even**. Stated for m ≤ 1411 because beyond that k = 2 becomes
arithmetically possible and the parity flips.

`VERIFIED` — the **complete lists of Giuga sequences of length ≤ 6** (1, 2, 3
and 17 sequences), whose all-prime members reproduce the known Giuga numbers
exactly. Three independent enumerators agree: the integer sieve and the
prime-set recursion from 001, and the integer-sequence search here.

`VERIFIED` — the squarefree N ≤ 2·10⁷ with Σ_{p|N} 1/p + 1/N = 1 — the a = 1
prefixes, which are where completions come from — are exactly 2, 6, 42, 1806 and
47058, **all even**. Two independent implementations agree (a sieve over every
integer in range, and a recursion over prime sets). A statement about 2·10⁷, not
about all N.

On the two `VERIFIED` lines above that rested on a single computation when first
written, and what was done before shipping them: the ≥ 9 bound is now confirmed
three ways (direct `Fraction` sum over the eight smallest odd primes,
`k_upper_bound(8,3) = 0`, and the enumerator returning empty for every m ≤ 8);
the parity law's load-bearing arithmetic step — n odd squarefree with m prime
factors ⟹ Σ n/p ≡ m (mod 2) — was checked on 4,000 random odd squarefree
integers with zero violations. The parity law's *conclusion* cannot be checked
against examples, because no odd Giuga number is known; what can be and was
checked is that it correctly declines to apply to the even case, where 5 of the
9 known Giuga numbers have odd factor counts and would contradict it.

`EVIDENCE` — **no odd Giuga number with fewer than 12 prime factors, and no odd
Giuga sequence of length under 12.** Exhaustive at lengths 9, 10, 11 by the
prime search and 7…11 by the sequence search, with no unresolved branches. Said
honestly: the odd lengths in that range are vacuous by the parity law, so the
only real content is length 10.

`EVIDENCE` — **no odd Giuga number among 4,111,633 fully-settled odd-prime
prefixes at factor counts 3 through 61**, covering 1,851,900,688 window integers
across three seeds. Each prefix is settled exhaustively; what is incomplete is
only which prefixes were sampled. Data in
`problems/giuga/data/odd-construct-seed*.json`, reproducible from the seed.

**What is not claimed.** Not that no odd Giuga number exists — only that none
has fewer than 12 prime factors and none turned up among the sampled prefixes.
Nothing whatever about Giuga's conjecture itself. The construction search is
complete within each prefix it settled and says nothing about the prefixes it
did not sample, which are the overwhelming majority.

## Why it failed / what survived

**The object was not found, and the reason is a density argument rather than a
computational wall.** Within a settled prefix, a completion needs the single
integer q = (Np − 1)/(ap − N) to come out exact, greater than p, odd and prime,
for some p in a window of width N/a. Exactness alone is roughly a
one-in-(ap − N) event and ap − N runs to the size of N, so scanning a billion
window integers is not "nearly there" — it is a thin sample of a space whose
successes, if any, are far sparser than that. Quoting the scanned count without
this caveat would be the overclaim to avoid.

**Where a wrong answer would have come from,** and what was done about each:

- *An unfactored N² − a treated as "no completions".* This is the failure mode
  P8 introduces, and the one that yields a false negative silently.
  `factorize_full` returns its cofactor, every caller appends the branch to an
  `unresolved` list, and every result line prints EXHAUSTIVE or
  COMPLETE EXCEPT n BRANCHES. A run with unresolved branches is not an
  exhaustion and is not reported as one.
- *A wrong closed form.* The rectangle identity was checked against all nine
  known Giuga numbers before use, and P7 (scan) versus P8 (rectangle) agree on
  every m ≤ 7. `lasttwo.completions` was cross-checked route-against-route on
  the same set.
- *An unsound new prune.* The coprimality rule was verified on all 23 known
  sequences first, and the length-3…6 sequence lists pinned in the test suite
  would have caught any change of output when it was switched on.

**The methodological lesson is Round 6.** Two rounds of compute went into factor
counts that a parity argument excludes in two lines. The cheap structural
observation should have come before the expensive search, and the reason it did
not is that the search was already written and running. Worth carrying to other
problems: when a search comes back empty at several consecutive sizes, ask
whether something makes those sizes vacuous before paying for the next one.

**On the push rule.** Eight pushes past a partial result produced: one closed
form (the rectangle), one new object class (integer sequences), two structural
laws (prefix coprimality, factor-count parity), one construction method, two
cost-model fixes, and one structural explanation (the a = 1 asymmetry). None of
them produced the object. The honest read is that persistence converted a
stalled search into a much faster stalled search plus several small theorems —
a good return, but not the same as the object falling out, and this record
should not be read as evidence that it was about to. If anything, Round 8 points
the other way: the construction family that yields every known Giuga number has
no odd member in reach, which is a reason to expect the odd case to stay hard
rather than to yield to more compute.

**Reusable:**

- `harness/giuga/lasttwo.py` — the last-two-entries step in both routes, where
  essentially all the cost of any Giuga enumeration lives.
- `harness/giuga/sequences.py` — Giuga sequences over integers with the `b = N`
  / pairwise-coprimality rule; the widest available handle on the odd question.
- `harness/giuga/factorint.py` — a factoriser whose contract is that it tells
  you when it failed.
- The identities: `b = N` at every node, `(ap − N)(aq − N) = N² − a`, prefix
  pairwise-coprimality, and `m ≡ 1 + k (mod 2)` for odd n.

## Leads generated

1. **Finish the even factor counts.** With parity, the open odd cases start at
   m = 12, then 14, 16. Run `odd_sequences.py --min-len 12` to completion and
   record whether any branch is left unresolved; a clean empty at 12 and 14
   would be the first genuinely new exclusion since m = 10.
2. **Sharpen the parity law's range.** It bites only while k = 1, i.e. m ≤ 1411.
   Compute the exact smallest m for which an odd prime set can reach Σ 1/p > 2;
   below it, odd factor counts are excluded outright. One computation, definite
   answer.
3. **Bias the construction search toward a = 1.** When a = 1 the rectangle is
   (p − N)(q − N) = N² − 1 = (N−1)(N+1), which has many divisors and needs no
   integrality filter — these are the fertile prefixes, and they are exactly
   what generates 1722, 858 and the Sylvester-type chain 2, 3, 7, 43, 1807.
   Search directly for odd prefixes with Σ_i N/p_i = N − 1 instead of sampling
   uniformly. Either odd a = 1 prefixes exist in quantity, or their scarcity is
   itself the obstruction and can be quantified.
4. **Settle whether an odd a = 1 prefix can exist at all.** For odd entries N
   and every N/p_i are odd, so for a prefix of length j, Σ_i N/p_i ≡ j (mod 2)
   and a = N − Σ N/p_i ≡ 1 − j (mod 2); hence a = 1 forces j even. Pair with
   lead 3: if a further congruence kills odd a = 1 prefixes entirely, the whole
   fertile family is unavailable in the odd case — a structural explanation for
   why no odd Giuga number has ever been found.
5. **Push the general Giuga-number list past 8 factors** with the hybrid step,
   and publish either the 8-factor list or the exact set of unresolved N² − a
   that blocks it. It did not finish here, and Round 9 says where to spend the
   effort: the deep nodes, not the closed form. Two concrete changes, both
   measurable — replace `_Primes` with a segmented sieve so that a node needing
   candidates near 10⁸ does not force a 134 MB allocation, and cache
   `max_tail_sum(last, r)`, which recomputes the same next-r-primes lookup at
   sibling nodes. Definite test: does the 200,000-node prefix of the m = 8 tree
   drop below 60 s?
6. **Give the length-12 odd-sequence run a two-pass factoriser**, as
   `giuga_by_factor_count.py` already has and `odd_sequences.py` does not: a
   small Pollard-Brent budget on the first pass, then retry only the branches
   that failed. At a budget of 8,000 the run had already accumulated 1,944
   unresolved branches after 52,639 nodes, so without two passes the answer
   would be heavily qualified even if it completed.

## References

- `problems/giuga/PROBLEM.md` — statement and published background.
- `problems/giuga/attempts/001-reciprocal-sum-relaxation.md` — same cycle; its
  bound is what establishes that a counterexample is too large to be an object,
  which is what sent this attempt at the odd question instead.
- No papers or OEIS lookups. Every number listed is output of code in this repo.
