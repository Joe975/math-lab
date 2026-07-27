# 001 — The forced pair model, and the partner-regular case split

- **Problem:** conway-99, `problems/conway-99/PROBLEM.md`
- **Date:** 2026-07-27
- **Mode:** blind
- **Source-commit:** `e5af9e6285ea9af0b5dba21adafe402d588df83d`
- **Type:** formalization + calibrated computational search + barrier analysis
- **Tools:** all standard-library Python 3.13, all deterministic.
  - `harness/conway-99/srg.py` — exact feasibility conditions and a
    definition-based strongly-regular checker (tier 0)
  - `harness/conway-99/constructions.py` — graphs built from their definitions,
    including the Berlekamp–van Lint–Seidel graph from the ternary Golay code
  - `harness/conway-99/local_model.py` — the local decomposition and its checks
  - `problems/conway-99/explore/pair_search.py` — the constraint search
  - `problems/conway-99/explore/soundness.py` — propagator soundness control
  - `problems/conway-99/explore/control_ab.py`, `partner_regular.py` — the A/B control
  - `problems/conway-99/explore/autos.py` — automorphism restrictions
  - `problems/conway-99/explore/crosscheck.py` — second implementations of every number
  - `problems/conway-99/explore/make_evidence.py` — regenerates the data file
  - `tests/test_conway99_srg.py` — 30 calibration tests, 0.5 s
  - No C compiler, `nauty` or SAT solver was available on this machine
    (`where gcc`, `where cc` both empty); everything is pure Python. No SAT
    encoding was used, so no UNSAT claim appears anywhere below.
- **Sources:** only `problems/conway-99/PROBLEM.md`, which is itself marked
  machine-transcribed `[T]`. No papers were consulted — this was a blind,
  offline attempt, and the automorphism results below were re-derived from
  scratch rather than read.

## Approach

The obvious attack is to throw the 99-vertex adjacency matrix at a solver:
C(99,2) = 4851 binary variables with λ/μ counting constraints. `PROBLEM.md`
records that this has been done and does not settle the problem. Repeating it
with a hand-written DPLL would have been strictly worse, and its UNSAT answers
would have been worthless without a DRAT checker I also did not have.

So I did the opposite: **collapse the search space to its true size first, and
calibrate every step on a parameter set whose answer is published.** The
decisive feature of this problem is that (99,14,1,2) sits in the λ=1, μ=2
family, and that family has *two realised members* — (9,4,1,2) and
(243,22,1,2). An open problem with a positive control in the same family is
rare, and it means every structural claim and every piece of tooling can be
tested against a graph that actually exists before being pointed at 99. That
control turned out to be the single most informative thing in the attempt: it
is what stopped me reporting a fake obstruction (see *Why it failed*).

## What was done

### 1. Exact feasibility, and a correction to the tier-0 statement

`srg.py` decides the standard conditions in integers only (`exact_isqrt`
returns `None` on non-squares; every comparison is between integers). The
multiplicities are obtained by solving

    f + g = n − 1        and        k + f·r + g·s = 0   (tr A = 0)

rather than by quoting a closed form, so the result is a re-derivation.

For (99,14,1,2) this gives eigenvalues r = 3, s = −4 with multiplicities
**f = 54, g = 44**.

> **`PROBLEM.md` is wrong here.** It states the multiplicities are "56 and 42".
> Those fail both trace identities: 14 + 3·56 − 4·42 = 14 ≠ 0, and
> 14² + 9·56 + 16·42 = 1372 ≠ 99·14 = 1386. The correct pair 54, 44 satisfies
> both. Cross-checked three independent ways in `crosscheck.py`: the linear
> system, the closed form ½[(n−1) ∓ (2k+(n−1)(λ−μ))/√D], and a brute-force
> integer scan over all (f, g) with f+g = 98. **I did not edit `PROBLEM.md`** —
> tier-0 edits were out of scope for this cycle. It should be corrected, since
> every future blind agent reads it.

The pipeline was then calibrated *in both directions* on eleven published
parameter sets — it must admit the realisable ones and reject the impossible
ones **for the right reason**. `(28,9,0,4)` must die by the absolute bound,
`(33,8,1,2)` by multiplicity integrality, `(21,10,4,5)` by the conference-graph
sum-of-two-squares condition. All eleven behave correctly
(`tests/test_conway99_srg.py`).

Re-deriving the family: SRG(n,k,1,2) forces n = (k²+2)/2 and 4k−7 = t², and
multiplicity integrality then leaves exactly five sets with n ≤ 2·10⁶:

    (9,4,1,2)  (99,14,1,2)  (243,22,1,2)  (6273,112,1,2)  (494019,994,1,2)

This confirms the tier-0 five-set claim, `VERIFIED` for n ≤ 2·10⁶.

### 2. Graphs built from scratch, including the positive control

`constructions.py` builds ten graphs from their definitions. The important one
is **BvLS(243,22,1,2)**, constructed as the coset graph of the perfect ternary
Golay code: the degree-5 divisor of x¹¹−1 over GF(3) is found by brute force
over all 3⁵ monic quintics (no generator polynomial is quoted), the code's
perfectness is *checked* (the 729 cosets of weight ≤ 2 representatives cover
all 3¹¹ vectors exactly), and the graph is verified to be SRG(243,22,1,2) by
the definition-based checker. This is the positive control.

### 3. The forced pair model

Fix v₀. With λ=1 the neighbourhood N(v₀) is 1-regular, i.e. a perfect matching
(7K₂ for k=14); write x′ for the partner of x. Then:

- a distance-2 vertex d has exactly μ=2 neighbours in N(v₀), and they are never
  partners (a partner pair plus v₀ would give an adjacent pair two common
  neighbours);
- conversely a non-partner pair {x,y} ⊂ N(v₀) is non-adjacent, so has exactly
  two common neighbours — v₀ and one other, at distance 2.

So **d ↦ {its two neighbours in N(v₀)} is a bijection from the distance-2 layer
onto the non-partner pairs of N(v₀)**, and the counts agree identically:
C(k,2) − k/2 = n − k − 1. For k = 14: 91 − 7 = 84 = 99 − 15.

Consequently, for each x ∈ N(v₀) the 12 distance-2 vertices whose pair contains
x induce a perfect matching M_x, readable as a perfect matching on
N(v₀) \ {x,x′}.

Everything outside D×D is therefore forced up to isomorphism, and the entire
search space is the C(84,2) = **3486** binary variables on D×D — down from 4851
on the raw adjacency matrix, and, more importantly, with all of the easy
structure already consumed rather than left for a solver to rediscover.

This reduction is not assumed: `local_model.py` checks every step and it was
run at **every vertex** of Paley(9), the 3×3 rook graph and BvLS(243) — 261
independent decompositions, all passing. It also correctly *refuses* Petersen
(λ=0), so the check is not vacuous.

### 4. The partner-regular case split

Running the reduction on the realised graphs turned up a sharp regularity:
**in both Paley(9) and BvLS(243), M_x is the partner matching, at every vertex
and for every x.** (In the cycle-type language of `local_model.cycle_type`,
the type is all-1s.) This is not forced by the SRG axioms — it follows from
those graphs being Cayley graphs on elementary abelian groups with connection
set closed under negation, where d = s+t and M_s(t) = −t. But it means the
configuration is exactly the one both existing graphs realise.

Call v₀ **partner-regular** if M_x is the partner matching for every x ∈ N(v₀).
This gives an auditable split:

- **Case A:** some vertex is partner-regular.
- **Case B:** no vertex is partner-regular.

The split is exhaustive by construction (B is the negation of A). Case A is
the case containing every candidate that resembles the two graphs that exist.

Under Case A the structure tightens a long way. Index N(v₀) by 7 *positions*
(the matching edges) × 2 signs; D then falls into 21 *blocks* of 4, indexed by
pairs of positions, each block inducing a 4-cycle. Two derivations, done
independently, agree exactly:

- **by hand:** the vertex {x,y} has exactly one neighbour in A_{x′} and it lies
  in its own block, so *no* out-of-block neighbour may lie in a block that
  meets its own; hence all 10 remaining neighbours lie in the 10 blocks
  disjoint from it, and each label (k,c) is used exactly twice;
- **by the propagator:** imposing partner-regularity and closing under the
  constraints forces exactly the same 1680 non-adjacencies.

The resulting closure, for k = 14:

| quantity | value |
|---|---|
| total D-pair variables | 3486 |
| decided by partner-regularity + propagation | 1806 (51.8%) |
| edges forced | 84 (two per D-vertex) |
| **left open** | **1680** |
| per D-vertex: neighbours still needed | 10 |
| per D-vertex: candidates | 40 |

The multiplicity profile n_B (how many of a vertex's 10 out-neighbours land in
each of the 10 disjoint blocks) satisfies Σ n_B = 10, n_B ≤ 2, and per-position
degree 4 — a 4-regular multigraph on K₅. Brute force over all 3¹⁰ multiplicity
functions gives **73 labelled solutions in exactly 4 orbits** under relabelling
the 5 positions: all-ones on K₅; multiplicity 2 on a 5-cycle; and two mixed
shapes. Re-enumerated by a second algorithm (choose the multiplicity-2 support,
then solve for the rest) with identical output. This case split is exhaustive.

In the all-ones orbit the 10 out-neighbours must pair into 5 edges joining
*disjoint* blocks, i.e. a perfect matching of the Petersen graph K(5,2) — of
which there are 6.

### 5. The search, and the control that matters

`pair_search.py` is a backtracking search over the 3486 variables with an
exact-count propagator (not a clausal encoding): for each pair it maintains
certain and possible common-neighbour counts against the forced target, and it
enforces the layer cardinality constraints |N(d) ∩ A_z| = 2 − [z′ ∈ pair(d)].

Two controls were run before any measurement was taken from it:

- **Soundness** (`soundness.py`): the real graphs are planted into the search
  state and must survive. Paley(9), rook(3) and BvLS(243) are all accepted at
  several base vertices and reconstruct to genuine SRGs. So the propagator does
  not prune away real solutions.
- **A/B** (`control_ab.py`): the *identical* search under the *identical*
  hypothesis is run on k=22 (satisfiable — BvLS is partner-regular) alongside
  k=14 (open).

The A/B control is what makes the rest of this section honest, and it came out
badly for the search — see below.

### 6. Automorphism restrictions, re-derived blind

Four statements, derived from the pair bijection alone. Their arithmetic is
checked exactly in `autos.py` and `crosscheck.py`; **the proofs themselves have
not had a skeptic pass, so they are `LIVE`, not `VERIFIED`.**

- **L1 (rigidity).** If σ fixes v and fixes N(v) pointwise then σ = id, because
  every other vertex is *named* by its pair in N(v) under the bijection.
- **L2.** If σ has prime order p > k/2 and fixes a vertex, then σ permutes the
  k/2 matching edges of N(v) with orbits of size 1 or p, so p > k/2 fixes every
  edge setwise; σ|N(v) then has order dividing both 2 and the odd prime p, so
  is trivial, and L1 gives σ = id. Hence such a σ is fixed-point-free and p | n.
  For n = 99 = 3²·11 and k/2 = 7: **no prime ≥ 13 divides |Aut(Γ)|**, leaving
  only {2,3,5,7,11}.
- **L3 (order 7).** 7 = k/2, so a non-identity σ of order 7 permutes the 7
  edges in a single cycle and fixes no vertex of N(v) for any fixed v; Fix(σ)
  is therefore an independent set, and it is closed under common neighbours
  (odd order acting on a set of size ≤ 2 fixes it pointwise), which forces
  |Fix| ≤ 1. Since |Fix| ≡ 99 ≡ 1 (mod 7), **|Fix(σ)| = 1 exactly**.
- **L4 (order 5).** |Fix| ≡ 99 ≡ 4 (mod 5) so Fix ≠ ∅. For v ∈ Fix the 7 edges
  split as one 5-orbit plus 2 fixed edges, giving |Fix ∩ N(v)| = 4 inducing
  2K₂; Fix is closed under common neighbours, so Fix induces a graph with
  k=4, λ=1, μ=2, and the counting identity forces exactly 9 vertices. Hence
  **an order-5 automorphism forces Fix(σ) to be the 3×3 rook graph**, and a
  short count then forces **every one of the 90 moved vertices to have exactly
  one neighbour in Fix**. I could not push this to a contradiction: the edge
  count closes consistently (18 + 90 + 45 + 180 + 360 = 693), so no
  contradiction is available at this level.

L1 is checked computationally on all three realised graphs; L2's consequence is
checked against |Aut(Paley 9)| = |Aut(rook 3)| = 72, computed by exhaustive
backtracking. Petersen and Clebsch are included as **negative controls**: they
have λ=0, the bound duly fails for them, which is what shows the passing checks
are testing something.

### Reproducing every number

```bash
python -m pytest tests/test_conway99_srg.py -q          # 30 tests, 0.5 s
python problems/conway-99/explore/crosscheck.py         # second implementations
python problems/conway-99/explore/soundness.py          # propagator soundness
python problems/conway-99/explore/partner_regular.py    # the closure table
python problems/conway-99/explore/autos.py              # automorphism results
python problems/conway-99/explore/control_ab.py 20000   # the A/B control
python problems/conway-99/explore/make_evidence.py      # regenerates the data file
```

## Outcome

**`MAP`**, with a small number of `VERIFIED` sub-results whose scope is stated,
and one correction to the problem statement.

`VERIFIED`:

- (99,14,1,2) satisfies every standard feasibility condition, in exact integer
  arithmetic, with spectrum {14¹, 3⁵⁴, (−4)⁴⁴}. **The multiplicities in
  `PROBLEM.md` (56 and 42) are wrong**; 54 and 44 are correct, by three
  independent derivations.
- The λ=1, μ=2 family has exactly five feasible sets — scope: **n ≤ 2·10⁶**.
- The pair reduction holds — scope: **proved in general, and machine-checked at
  every vertex of all three realised graphs in the family (261 decompositions)**.
- Under partner-regularity at v₀ with k=14, exactly **1806 of 3486** pair
  variables are forced and **1680 remain open**, each D-vertex needing 10
  neighbours from 40 candidates. Derived by hand and by propagation
  independently, agreeing exactly.
- The n_B multiplicity profile has exactly **73 labelled solutions in 4
  orbits** — scope: exhaustive over all 3¹⁰ multiplicity functions, by two
  algorithms.

`LIVE` (proved here, not yet adversarially reviewed): the four automorphism
statements L1–L4, in particular *no prime ≥ 13 divides |Aut(Γ)|* and *an
order-7 automorphism has exactly one fixed point*.

**What is not claimed.** No existence or nonexistence result, and no
elimination of any case. Case A is *not* eliminated — it is only reduced. The
search proved nothing: it ran to a node budget, and, as the next section
explains, its stopping point is not evidence about 99. Nothing here is a
nonexistence proof, no SAT encoding or UNSAT answer is involved, and the four
automorphism claims are proof steps that still need a skeptic pass.

## Why it failed / what survived

### The obstruction, stated precisely

The route dies at a specific, quantifiable place: **the forced local structure
consumes only about half the search space, and the λ/μ conditions transmit
almost no information into the other half.**

Concretely, under the strongest structural hypothesis available
(partner-regularity, the configuration both realised graphs satisfy), closing
under every constraint the propagator has leaves each of the 84 D-vertices
needing to choose **10 neighbours out of 40 candidates** — C(40,10) ≈ 8.5·10⁸
per vertex — with 1680 variables still free and *zero* vertices fully settled.
The same measurement on the satisfiable control is worse in proportion: at
k=22, knowing **all 22 matchings** determines only 8250 of 24090 pairs (34%),
and pure propagation adds essentially nothing beyond what the matchings state
outright.

That is the quantity going the wrong way. The λ=1, μ=2 conditions are *local*:
they constrain pairs and their common neighbourhoods, and the pair model has
already extracted everything local. What remains — the adjacency between
disjoint pairs — is where the global obstruction must live, and none of the
standard conditions reach it. This is a structural statement about the problem,
not about my implementation, and it is the same reason a SAT encoding strains:
after unit propagation on the forced structure there is nothing left to
propagate, and the residual instance is a large, nearly-unconstrained design
problem.

### The honest negative: my search failed its own positive control

I measured where the extension search explodes, and then the A/B control told
me the measurement was worthless as evidence about 99:

| k | n | satisfiable? | settled prefix, 3000 nodes | settled prefix, 20000 nodes |
|---|---|---|---|---|
| 4 | 9 | yes (Paley 9) | 4/4 — solved in 1 node | 4/4 |
| 22 | 243 | **yes (BvLS 243)** | 14/220 | **14/220 — unchanged** |
| 14 | 99 | open | 6/84 | 7/84 |

**The search cannot reconstruct a graph that provably exists.** At k=22 it
settles 14 of 220 vertices and then thrashes: the extra 17000 nodes between the
two columns bought *nothing*, even though BvLS is sitting in the space it is
searching and the propagator has been proved sound against it.
Therefore the point at which it stalls on k=14 is a fact about this search, not
about (99,14,1,2), and I am not reporting it as an obstruction. Per the
verification contract — "a pipeline that cannot recover known
strongly-regular existence and nonexistence results is not evidence about 99" —
this route fails that bar for the search component specifically, and passes it
for the feasibility and reduction components.

That failure is itself the most transferable finding here. A plain
exact-count propagator plus DFS is not merely slow on 99; it is **too weak to
be diagnostic at all** in this family, because the residual instance has almost
no propagation available. Any future attack needs either isomorph rejection
strong enough to quotient the block structure, or a genuinely different
global invariant.

I also caught and fixed a real gap mid-attempt: the first propagator omitted
the μ-condition between a D-vertex and a P-vertex, which made it miss a lemma I
had already proved by hand and confirmed on BvLS. Adding those layer
constraints roughly doubled the forced structure (966 → 1806 pairs). Worth
recording as a methodological lesson: **the propagator and the hand derivation
disagreeing was the signal**, and it would have been invisible without the
hand derivation to compare against.

### What survived, and is reusable

- The **pair model** itself: an exact, validated reduction of any SRG(n,k,1,2)
  to a D×D adjacency problem on the non-partner pairs, with the forced seed
  fully determined up to the hyperoctahedral group B_{k/2} of order
  2^{k/2}·(k/2)!. For 99 this is 3486 variables, and it is the right object for
  any future attack.
- The **positive control methodology**: BvLS(243) built from scratch, and the
  soundness harness that plants a real graph into a search state. Any future
  search on this problem can and should be run through `soundness.py` and
  `control_ab.py` before its output is believed.
- The **exact feasibility calculator**, calibrated in both directions, and the
  correction to the tier-0 spectrum.
- The **partner-regular case split** with its exhaustive 4-orbit sub-split, and
  the block/Petersen-matching structure inside it.
- The four **automorphism lemmas**, re-derived blind.

## Leads generated

1. **Eliminate Case A, all-ones orbit.** Under partner-regularity with n_B ≡ 1,
   each D-vertex's out-neighbourhood is one vertex per disjoint block, pairing
   into a perfect matching of the Petersen graph K(5,2) — 6 choices per vertex,
   84 vertices, plus consistency. This is a far smaller object than the raw
   1680 variables and may be finitely checkable. Outcome either way: Case A
   all-ones is eliminated, or a candidate emerges. Note the analogous
   configuration at k=22 has n_B ≡ 1 on only 18 of 36 disjoint blocks, so 99 is
   in a *different* regime from BvLS here — 10 out-neighbours across exactly 10
   disjoint blocks — which is a specific, checkable rigidity.

2. **Push L4 to a contradiction or a construction.** An order-5 automorphism
   forces Fix(σ) = the 3×3 rook graph, with all 90 moved vertices having
   exactly one neighbour in Fix, class sizes 10, and prescribed inter-class
   degrees (1 between adjacent classes, 2 between non-adjacent). Edge counts
   close consistently, so any contradiction must come from λ=1 inside the
   classes. Definite outcome: either 5 ∤ |Aut(Γ)| falls out, or an explicit
   σ-invariant candidate does.

3. **Attack the order-7 case using |Fix| = 1.** L3 pins the fixed-point count
   exactly, so the quotient is 14 orbits of size 7 plus one fixed vertex, and
   the fixed vertex's 7 matching edges are permuted in a single 7-cycle. Build
   the quotient explicitly and search *it* rather than the 99-vertex graph;
   the object is roughly 14 variables' worth of orbit-adjacency, not 3486.
   This is by far the smallest search this problem admits.

4. **Test whether partner-regularity is forced at some vertex.** If some vertex
   of any SRG(99,14,1,2) must be partner-regular, Case B vanishes and lead 1
   becomes the whole problem. Concretely: count the M_x cycle types compatible
   with the layer constraints at a single vertex and check whether any non-all-
   ones type survives. Definite outcome either way, and cheap.

5. **Correct `PROBLEM.md`.** The multiplicities 56 and 42 should read 54 and
   44. Mechanical, and it currently misleads every blind agent.

## References

- `problems/conway-99/PROBLEM.md` — the tier-0 statement, itself marked machine
  transcribed `[T]`. The only source read during this attempt. Its
  multiplicities are incorrect; see above.
- No papers were consulted. The automorphism restrictions of Makhnev–Minakova,
  Behbahani–Lam and Cesarz–Woldar, and the SAT analysis of Keramatipour, are
  named in `PROBLEM.md` but were **not** read; L1–L4 above are independent
  re-derivations and I make no claim about how they relate to those results
  beyond what `PROBLEM.md` states.
- Berlekamp–van Lint–Seidel graph: constructed here from the perfect ternary
  Golay code [11,6,5]₃ directly, not taken from a table.
- Data: `problems/conway-99/data/pair-model-evidence.json`, regenerated by
  `problems/conway-99/explore/make_evidence.py`.
