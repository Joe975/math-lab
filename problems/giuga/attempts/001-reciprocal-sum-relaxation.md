# 001 — Reciprocal-sum relaxation of the counterexample conditions

- **Problem:** giuga, `problems/giuga/PROBLEM.md`
- **Date:** 2026-07-27
- **Mode:** blind
- **Source-commit:** `e5af9e6285ea9af0b5dba21adafe402d588df83d`
- **Type:** elementary reduction + exhaustive computational search
- **Tools:** written for this attempt, all pure-Python standard library (no C
  compiler on this machine: `where gcc` / `where cc` / `where tcc` all miss, so
  the C-kernel route in the README was unavailable and nothing here uses it).
  Tier 0: `harness/giuga/conditions.py` (three formulations of the congruence),
  `harness/giuga/enumerate_giuga.py` (two Giuga-number enumerators),
  `harness/giuga/admissible.py` (the factor-count / product bound search).
  Tier 1: `problems/giuga/explore/*.py`. Everything is deterministic; longest
  completed run 336 s (one further run was stopped unfinished at 1250 s — see
  the obstruction below).
- **Sources:** `problems/giuga/PROBLEM.md` only. The published bounds quoted
  there (3459 prime factors, ≥13,800 and ≥19,908 decimal digits) are used purely
  as a target to compare against; they are marked `[T]` there as machine
  transcriptions and I did not consult the papers. No OEIS lookup was made — the
  Giuga numbers below are output of the code in this repo, and the four named in
  `PROBLEM.md` are the only externally-supplied values I compared against.

## Approach

The obvious line is to enumerate Giuga numbers by factor count and filter them
for the Carmichael condition. I built that first (`enumerate_by_prime_sets`) and
it does work — but it dies at **8 prime factors**, while a counterexample needs
thousands. The two scales never meet, and no amount of tuning closes a gap of
that size. So that enumeration is worth building only as calibration.

What I did instead is throw information away on purpose. The counterexample
conditions imply three much weaker facts that are *purely additive and
pairwise*:

- every prime factor is odd;
- no prime factor divides (another prime factor − 1);
- the reciprocals of the prime factors sum to more than 1.

Call a finite set of primes with those three properties **admissible**. Then

    (number of prime factors of a counterexample)  ≥  min { |S| : S admissible }

and the right-hand side is a self-contained combinatorial optimisation with no
divisibility conditions in it at all. It is finite, it is searchable, and — this
is the reason to prefer it — its pruning rules are simple enough that their
soundness can be *stated and tested*, which the verification contract in
`PROBLEM.md` singles out as the way to be wrong here.

The cost of the relaxation is that it can never see the exact Giuga divisibility
conditions again. I measured that cost too (§6), rather than leaving it unknown.

## What was done

### 1. Deriving the counterexample conditions rather than assuming them

`PROBLEM.md` states the characterisation; I re-derived it so the search rests on
something I checked. Write S(n) = Σ_{k=1}^{n−1} k^(n−1).

**L1.** For any prime p | n, the residues k mod p run over 0…p−1 exactly n/p
times as k runs over 0…n−1, and Σ_{r=1}^{p−1} r^t ≡ −1 (mod p) when (p−1) | t
and ≡ 0 otherwise. Hence

    S(n) ≡ −(n/p) (mod p)  if (p−1) | (n−1),      S(n) ≡ 0 (mod p) otherwise.

So S(n) ≡ −1 (mod n) iff for **every** p | n both (p−1) | (n−1) and
n/p ≡ 1 (mod p).

**L2 (squarefree).** If p² | n then p | (n/p), so the displayed congruence gives
S(n) ≡ 0 (mod p); but −1 ≢ 0 (mod p). No non-squarefree n satisfies the
congruence. (A proof, not a search.)

**L3 (odd).** If 2 | n, pick an odd p | n; then (p−1) | (n−1) forces an even
number to divide the odd number n−1.

**L4 (pairwise).** If p, q | n are distinct and q | (p−1), then q | (p−1) |
(n−1), so q divides both n and n−1, so q | 1.

**L5 (reciprocal sum).** Multiplying Σ_i 1/p_i − 1/n by n gives Σ_i n/p_i − 1,
which is ≡ n/p_j − 1 ≡ 0 (mod p_j) for each j because p_j divides every other
term. So Σ_i 1/p_i − 1/n = k ∈ ℤ, and k ≥ 1 because the quantity is positive for
m ≥ 2. Hence **Σ_{p|n} 1/p = k + 1/n > 1**.

L3, L4 and L5 are exactly the admissibility conditions; L2 makes the prime set
determine n.

`problems/giuga/explore/check_lemmas.py --limit 3000000` re-checks L3–L5
empirically over every Carmichael and Giuga number below 3·10⁶: 63 Carmichael
numbers, none even, none with a factor dividing another factor minus one; all
four Giuga numbers below 10⁶ have reciprocal sum > 1. It also reports the
quantity the whole approach turns on — the largest reciprocal sum any actual
Carmichael number below 3·10⁶ achieves is **0.565846**, at 62745 = 3·5·47·89.
Real Carmichael numbers are nowhere near the 1 a counterexample needs.

### 2. Calibration: three formulations, checked against each other

`harness/giuga/conditions.py` implements the congruence three ways that share no
arithmetic: the literal power sum (`power_sum_mod`), the same residue rebuilt
from the factorisation by CRT using L1 (`power_sum_mod_structural`), and Agoh's
Bernoulli form n·B_{n−1} ≡ −1 (mod n) with exact `Fraction` Bernoulli numbers
(`agoh_holds`).

    python problems/giuga/explore/calibrate_faces.py --direct-limit 4000 --agoh-limit 260

- direct vs structural: **2432 exact residues** compared (the squarefree n in
  2…4000), all equal; the booleans agree for all n in 2…4000;
- both agree with primality for every n in 2…4000 — no composite below 4000
  satisfies the congruence;
- Agoh agrees with primality for n in 2…261;
- no even n in 4…260 passes Agoh, which is the Bernoulli side independently
  reproducing L3 (B vanishes at odd index ≥ 3).

### 3. Calibration: recovering the Giuga numbers from scratch

Two enumerators sharing nothing (`harness/giuga/enumerate_giuga.py`): a
smallest-prime-factor sieve over the integers, and a recursion over *prime sets*
using Σ 1/p_i − 1/n = k ∈ ℤ.

    python problems/giuga/explore/calibrate_giuga_numbers.py --scan-limit 1000000 --max-factors 7

- sieve to 10⁶ → 30, 858, 1722, 66198 — the four values named in `PROBLEM.md`,
  recovered independently;
- prime-set recursion, complete by factor count:
  - 3 factors: 30 = 2·3·5
  - 4 factors: 858 = 2·3·11·13, 1722 = 2·3·7·41
  - 5 factors: 66198 = 2·3·11·17·59
  - 6 factors: 2214408306 = 2·3·11·23·31·47057,
    24423128562 = 2·3·7·43·3041·4447
  - 7 factors: 432749205173838 = 2·3·7·59·163·1381·775807,
    14737133470010574 = 2·3·7·71·103·67213·713863,
    550843391309130318 = 2·3·7·71·103·61559·29133437
- the two agree exactly on the overlap (everything ≤ 10⁶ with ≤ 7 factors);
- every set returned is re-checked by the pointwise divisibility form
  `p | (n/p − 1)` on the integer, which is different arithmetic from the
  rational recursion that produced it.

**8 factors does not finish.** A branch reaches a state where the closed form for
the last two primes needs candidates up to **54,164,630,208**, which is not
sievable. That is the precise failure point, and it is why this route cannot be
pushed toward a counterexample.

Pruning rules P1–P7 are listed with justifications in the module docstring; each
is a *necessary* condition, so nothing that could extend to a solution is
discarded. Tested rather than asserted:

    python problems/giuga/explore/validate_giuga_prunes.py --max-factors 7

Disabling P4 or P7 leaves the output identical at every m ≤ 7 and changes only
the cost (m = 7: 1775 nodes with P7, 674182 without).

### 4. The extension: a certified lower bound on the factor count

`harness/giuga/admissible.py`. Fix a branch limit X. Any admissible S splits as
T = S ∩ [3, X] and B = S ∩ (X, ∞). Every member of B is > X and T-admissible, so
its reciprocal sum is at most that of the |B| smallest T-admissible primes above
X. That gives the least |B| that could close the gap to 1, hence a bound
|T| + t_min(T) on |S|; the minimum over all admissible T ⊆ [3, X] is a valid
lower bound L(X) for **every** admissible set, including those whose large primes
are nothing like the ones the bound used.

No floating point anywhere: reciprocals are integers scaled by 10²⁴ with
*directed* rounding (`hi[i] ≥ SCALE/p`), so every comparison the search makes
errs toward reporting a weaker bound.

    python problems/giuga/explore/run_bound.py --pool 3000000 --cap 50000 --xs 3,13,23,53,73,101,151,199,251,307,401,503,601,701

| X | L(X) ≥ | nodes | time |
|---|---|---|---|
| 3 | 27 | 2 | 0.00 s |
| 13 | 202 | 10 | 0.00 s |
| 23 | 554 | 40 | 0.01 s |
| 53 | 704 | 72 | 0.03 s |
| 73 | 965 | 117 | 0.06 s |
| 101 | 1510 | 325 | 0.30 s |
| 151 | 1670 | 449 | 0.67 s |
| 199 | 2109 | 921 | 1.57 s |
| 251 | 2308 | 1281 | 2.10 s |
| 307 | 2431 | 1651 | 3.31 s |
| 401 | 2724 | 3180 | 7.41 s |
| 503 | 2927 | 5935 | 18.5 s |
| 601 | 3204 | 17746 | 69.9 s |
| 701 | **3355** | 43707 | 336 s |
| 809 | — | — | stopped, unfinished after 1250 s |

Full table including every X between 3 and 701 in
`problems/giuga/data/bound-by-X.json`.

The same search with the product as objective (`run_product_bound.py`, exact
bigint products, no logarithms) gives a size bound:

| X | n ≥ | time |
|---|---|---|
| 13 | 10^577 | 0.00 s |
| 53 | 10^2522 | 0.70 s |
| 101 | 10^6059 | 4.25 s |
| 199 | 10^8853 | 11.6 s |
| 251 | 10^9802 | 16.6 s |
| 307 | **10^10395** | 21.9 s |

The product search has a `cap` on how many tail primes a leaf bound may use.
Running X = 101 and X = 199 at cap 4000, 8000 and 20000 gives identical answers
(6060 and 8854 digits), so the cap is not binding on the reported minimum.

The minimising configuration is the same at every X: **T never contains 3.**
Including 3 buys 1/3 of the required sum but forbids every prime ≡ 1 (mod 3) —
half of all primes — and the trade is never worth it. The witness T is exactly
the greedy admissible set started at 5.

### 5. Verification of the bound

Three handles, in decreasing strength:

1. **A separate naive implementation.** `explore/independent_bound.py` shares
   nothing with the search: `itertools.combinations` over every subset of the odd
   primes ≤ X, no pruning, no bitmasks, no scaled integers, exact `Fraction`s
   throughout. It reproduces L(X) exactly at X = 13, 23, 31, 43, 53 (202, 554,
   554, 704, 704). At X = 73 it had not finished after 14 minutes and was
   stopped — 2²⁰ raw subsets, each surviving one costing a `Fraction` sum of
   ~10³ terms — so X ≤ 53 is the honest reach of the fully independent check.
2. **Pruning on vs off, plus an exact-Fraction recheck of the winner.**
   `explore/validate_pruning.py` runs the search with `prune=False` (a complete
   enumeration of admissible T) and recomputes the winning leaf with
   `leaf_bound_exact`, which uses plain division and `Fraction`. All three agree
   at every X ≤ 73 — at X = 73 the pruned search visits 117 nodes and the
   unpruned one 45588, and both return 965. This is the direct test of the
   failure mode `PROBLEM.md` warns about: an unsound prune would delete the
   branch carrying the *smaller* bound and report a larger, better-looking
   number.
3. **The lemmas the reduction rests on**, re-checked over every Carmichael number
   below 3·10⁶ (§1).

**What is not cross-checked:** L(X) for X > 73 rests on a single search program.
For the headline X = 701 the winning leaf is re-verified exactly, but the claim
"no other T does better" comes from that one program, with its pruning validated
only at smaller X. The 3355 figure should be read with that caveat.

### 6. The ceiling: how far this relaxation can *ever* go

Because the bound is a minimum over admissible sets, **any explicit admissible
set is a ceiling on it**. Greedily taking the smallest still-allowed prime
starting from 5 produces one:

    python problems/giuga/explore/method_ceiling.py --pool 5000000 --skips "3" --dump problems/giuga/data/ceiling-set.json
    python problems/giuga/explore/verify_ceiling_set.py problems/giuga/data/ceiling-set.json

**8135 primes**, largest 321053, reciprocal sum 1.000000919, product 40733
digits. Re-checked by a separate script working from the definition rather than
the sieve that built it: all 8135 elements odd, distinct, increasing, prime by
Miller-Rabin; the pairwise rule verified over all 66,170,090 ordered pairs; the
reciprocal sum > 1 as an exact `Fraction`. The set is stored at
`problems/giuga/data/ceiling-set.json` (61 KB) so it can be re-checked in one
command.

So min{|S| : S admissible} ≤ 8135, and **no amount of pushing X can take this
method past 8135 prime factors**.

For contrast, greedy *including* 3 fails: using all 75651 admissible primes below
5·10⁶ it reaches only 0.989052.

## Outcome

`VERIFIED` — the three formulations of the congruence agree, and agree with
primality, for **n = 2…4000** (2432 exact residues compared on the squarefree n;
Agoh's form checked for n = 2…261). A statement about that range.

`VERIFIED` — the Giuga numbers with **at most 5 prime factors** are exactly 30,
858, 1722, 66198, recovered by two implementations sharing no arithmetic
(integer sieve to 10⁶, prime-set recursion).

`EVIDENCE` — the complete lists at **6 and 7 prime factors** are the five numbers
in §3. Membership is double-checked (rational form and pointwise divisibility);
*exhaustiveness* rests on one enumerator whose pruning rules were validated by
toggling them off, not by a second program.

`EVIDENCE` — **a composite solution of Giuga's congruence has more than 3355
distinct prime factors, and exceeds 10^10395.** Exact family searched: sets of
odd primes, pairwise satisfying "no member divides another member minus one",
with reciprocal sum > 1; branch limit X = 701 for the factor count and X = 307
for the size; tail primes drawn from a pool of the odd primes below 3·10⁶ (10⁶
for the size bound). Exact integer arithmetic throughout, rounding directed so
that error can only weaken the result.

`VERIFIED` — **the method's ceiling is 8135.** An admissible set of exactly that
size is exhibited and independently re-checked, so this relaxation cannot prove a
factor-count bound above 8135 at any branch limit.

**What is not claimed.** Nothing about Giuga's conjecture itself. The bound is
about a *relaxation*: every counterexample gives an admissible set, but almost no
admissible set is a counterexample, so a large admissible set is not a near miss.
3355 does not reach the transcribed published figure of 3459 prime factors, and
10^10395 is well short of the transcribed 19,908 digits — I did not match either
and did not try to. No claim that the enumerations above 5 prime factors are
exhaustive beyond what §5 supports.

## Why it failed / what survived

**The obstruction is the branching factor over small primes, and it is
quantified.** L(X) is a minimum over admissible subsets T of the primes ≤ X, and
the number of those grows exponentially in π(X). The pruning knocks the visited
count down by three orders of magnitude (117 vs 45588 nodes at X = 73) but does
not change the growth rate: nodes went 5935 → 17746 → 43707 for X = 503 → 601 →
701, and time 18 s → 70 s → 336 s. Each further step of X costs roughly 4–5×.
**X = 809 was still running after 1250 s of CPU and was stopped**, which is
where this attempt actually ends. `SPECULATION`: on that curve, reaching the
transcribed 3459 needs X somewhere around 800–900 and an hour or two in this
implementation, so the published factor count is close but not free; reaching
the 8135 ceiling would need X in the thousands and is out of reach in pure
Python.

**The second obstruction is the ceiling, and it is the more interesting one.**
Even with unlimited compute this relaxation stops at 8135 prime factors, because
an admissible set of that size demonstrably exists. That single number explains
the shape of the published results. `SPECULATION`, but a testable one: the 1996
factor-count bound of 3459 sits comfortably inside the reciprocal-sum
relaxation's reach, so it is plausibly this argument or a close relative, while
the 2013 improvement to 19,908 digits cannot come from here and must re-use the
exact divisibility conditions that L1 discards.

**Where a wrong answer would have come from,** as `PROBLEM.md` predicts: the
pruning. Rules R2 and R3 lower-bound a subtree; if either over-estimated, the
branch carrying the smaller bound would be deleted silently and the search would
report a bigger, better-looking number with nothing in the output to show for it.
That is why §5 exists, and why the unpruned run and the separate `itertools`
implementation matter more than the headline. The same hazard is present in the
Giuga enumerator (P4, P7) and is tested the same way.

**A near-miss worth recording.** The scaled-integer arithmetic rounds every
reciprocal *upward*. That is the direction that makes the search believe a branch
is more capable than it is, so it prunes less and reports a smaller bound. Had
the rounding been directed the other way the bound would have come out
marginally too large, and nothing in §5 would necessarily have caught it, because
both scaled-integer implementations would have shared the convention. The
`Fraction`-based checks are the only guard against that class of error.

**Reusable:**

- `harness/giuga/conditions.py` — three mutually-checking formulations of the
  congruence, the Giuga/Carmichael/counterexample predicates, and exact Bernoulli
  numbers. Any future computational claim here can be cross-checked with the face
  it did not use.
- `harness/giuga/enumerate_giuga.py` — both Giuga-number enumerators, with
  per-rule toggles, so a successor can re-run the calibration and the soundness
  test rather than trusting this record.
- `harness/giuga/admissible.py` — the bound search in both objectives, with the
  pruning switch that makes it auditable.
- L1–L5 above as a self-contained derivation of the counterexample conditions.
- `problems/giuga/data/ceiling-set.json` — the 8135-prime witness. Anyone
  claiming a factor-count bound above 8135 from an argument that uses only
  oddness, the pairwise rule and Σ1/p > 1 is wrong, and this file shows it.

## Leads generated

1. **Close the gap to the published 3459 and locate the crossover.** Finish
   X = 809 (stopped here after 1250 s) and run 907, 1009. If L(X) passes 3459 at
   X ≈ 900, that is concrete support for the speculation that the 1996 bound is
   this relaxation; if L(X) flattens well below 3459, the 1996 argument uses
   something this relaxation does not, and the difference is worth finding.
   Definite outcome either way. The obvious enabling optimisation is to carry
   the filtered candidate list plus its prefix sums down the recursion, turning
   the ~110 `t_min` scans in each node's child loop into binary searches; that
   is where the 4 ms per node goes.
2. **Compute min{|S| : S admissible} exactly rather than bounding it.** The
   greedy ceiling is 8135; a better construction lowers it. Run local search
   (swap one prime for two smaller admissible ones, and so on) from the witness
   set. If the true minimum lands near 3459, the published bound sits *at* this
   method's ceiling and the whole family is closed — a strong negative result
   about the approach.
3. **Add back one divisibility condition and re-measure the ceiling.** The
   cheapest is the Giuga condition at the smallest prime: if 5 | n then the
   product of the other factors ≡ 1 (mod 5), while every one of them is ≢ 1
   (mod 5) by L4. That constrains the multiset of residues mod 5 and is checkable
   at a leaf. Rebuild the ceiling set subject to it: if the ceiling rises above
   8135 the extra condition has real force and points at how the 19,908-digit
   bound is reached; if it barely moves, the family is capped near 8135.
4. **Push the Giuga enumeration to 8 factors with a two-sided last step.** The
   blocker is a closed form needing primes to 5.4·10¹⁰. Meet-in-the-middle on the
   final two primes — factorising b(Np−1) instead of scanning p — removes the
   scan. Definite outcome: either the 8-factor list is produced, or the
   factorisations needed are themselves too hard, which is itself a clean
   statement of the barrier.
5. **Rule out odd Giuga numbers at small factor counts by exhaustion.**
   `k_upper_bound(m, 3)` is 0 for m ≤ 8 — the reciprocals of the eight smallest
   odd primes sum to 0.998 < 1 — so no odd Giuga number has fewer than 9 prime
   factors, immediately. Running the prime-set enumeration with `min_prime=3` for
   m = 9, 10, 11 extends that by exhaustion and is cheap. It is a fact about
   Giuga numbers rather than counterexamples, but it uses tooling already here.

## References

- `problems/giuga/PROBLEM.md` — the statement, Agoh's formulation, and the
  published figures (marked `[T]` there as machine transcriptions; not
  independently checked here, and deliberately not used as inputs to anything).
- No papers were consulted. The Giuga numbers listed in §3 beyond the four in
  `PROBLEM.md` are output of `harness/giuga/enumerate_giuga.py`, not lookups.
- Tooling and data produced by this attempt: `harness/giuga/`,
  `problems/giuga/explore/`, `problems/giuga/data/`,
  `tests/test_giuga_harness.py`.
