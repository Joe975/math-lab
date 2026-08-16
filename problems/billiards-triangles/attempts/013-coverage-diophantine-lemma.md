# 013 — The coverage Diophantine lemma: windows-cover-(0,1) is unconditional

- **Problem:** billiards-triangles, `problems/billiards-triangles/PROBLEM.md`
- **Date:** 2026-08-15
- **Mode:** informed
- **Type:** elementary proof + exact computational verification
- **Tools:** `problems/billiards-triangles/explore/coverage_diophantine.py`
  (stdlib only, deterministic — seeded RNG; `selftest` runs 73,542 exact
  checks in ~1 s)
- **Sources:** none from the literature (two targeted searches found no
  citation for the lemma; see Novelty below). Starting point: an external
  ChatGPT session's handoff supplied a proof sketch of 006's lead 2; every
  step was independently re-derived here before being trusted (provenance
  note in Approach).

## Approach

Attempt 006 (lead 2) and queue item 12 reduce "every obtuse angle has an
alive W member" — *assuming 006's SPECULATION birth law and exact-window
claim* — to an elementary Diophantine statement:

> **LEMMA.** For every real t ∈ (0, 1) there exist integers a ≥ 1, b ≥ 1
> with a + b < t·a(b+1) < a + b + 1.

(For an obtuse apex angle gamma, t = (180 − gamma)/90, and (a, b) names
the family member W(a, b) whose window [gamma_birth, gamma_d] strictly
contains gamma under the window law.)

006 expected this to need continued fractions ("the needed member length
at gamma is governed by the continued-fraction structure of t"). The
observation here is that *existence* needs nothing of the sort: fixing
q = b + 1 and letting a run turns the condition into an overlapping chain
of open rational intervals whose union telescopes to (1/q, 2/q), and
those cover (0, 1). Continued fractions remain relevant only to the
*optimal* (minimal a+b, hence minimal word length) witness.

Provenance: the proof below was proposed in a user-supplied handoff from
an external ChatGPT session. Per lab rules the handoff was treated as
unverified input: each algebraic step was re-derived from scratch here
(and the equivalences additionally machine-checked exactly, section
"What was done"), the constructive witness formulas and their bound
proofs were added, and the sharpness observation at a = q − 1 is new to
this record.

## What was done

### The proof, re-derived

**Step 1 (change of variable).** Put q = b + 1 ≥ 2. For integers
a ≥ 1, q ≥ 2, dividing the chain a + q − 1 < t·aq < a + q by aq > 0
gives the equivalent statement that t lies in the open interval

    I(a, q) = ( 1/q + (q−1)/(aq) ,  1/q + 1/a ),

using (a + q − 1)/(aq) = 1/q + (q−1)/(aq) and (a + q)/(aq) = 1/q + 1/a.

**Step 2 (chain overlap).** Both endpoints of I(a, q) strictly decrease
in a. Consecutive intervals overlap iff the upper endpoint of I(a+1, q)
exceeds the lower endpoint of I(a, q):

    1/(a+1) > (q−1)/(aq)   ⟺   aq > (q−1)(a+1)   ⟺   a > q − 1.

So for a ≥ q the chain I(q, q), I(q+1, q), … is connected. Sharpness: at
a = q − 1 the two endpoints are *equal* (1/q = (q−1)/((q−1)q)), so the
chain touches without overlapping — a ≥ q is exactly the right cutoff,
and starting the chain at a = q loses nothing (I(q−1, q) ⊂ (1/q, 2/q)
contributes no new points… its upper endpoint 1/q + 1/(q−1) exceeds 2/q,
so for completeness: the union below is over a ≥ q and is already all of
(1/q, 2/q); smaller a only adds points ≥ 2/q, which the next q handles).

**Step 3 (union over a).** The upper endpoint of I(q, q) is 2/q; the
lower endpoints 1/q + (q−1)/(aq) decrease to 1/q as a → ∞. With Step 2,

    ⋃_{a ≥ q} I(a, q) = (1/q, 2/q).

**Step 4 (union over q).** Given t ∈ (0, 1), the open interval
(1/t, 2/t) has length 1/t > 1, so it contains an integer q, and q > 1/t
> 1 forces q ≥ 2. Then 1/q < t < 2/q, so by Step 3 some a ≥ q has
t ∈ I(a, q); with b = q − 1 the lemma follows. ∎

**Constructive witness** (proved, and used by the checker):

    q = ⌊1/t⌋ + 1,     a = ⌊(q−1)/(qt−1)⌋ + 1,     b = q − 1.

For q: q > 1/t by construction, and q ≤ 1/t + 1 < 2/t since 1/t > 1
(strict even when 1/t is an integer). For a: t ∈ I(a, q) is equivalent to
(q−1)/(qt−1) < a < q/(qt−1) (divide the endpoint inequalities through by
t − 1/q = (qt−1)/q > 0); that open window has length 1/(qt−1) > 1
because t < 2/q, so the floor-plus-one lands inside it (strict at both
ends by the same 1/(qt−1) > 1); and a > (q−1)/(qt−1) > q − 1 because
qt − 1 < 1, so a ≥ q holds automatically.

### Exact verification, two independent arithmetic layers

`coverage_diophantine.py selftest` (deterministic, ~1 s, 73,542 exact
checks, zero failures):

1. **Form equivalence** — the interval form (Fraction arithmetic) agrees
   with the original inequality in pure integer arithmetic
   (r(a+b) < p·a(b+1) < r(a+b+1) for t = p/r) at 20,000 random
   (a, b, p) triples.
2. **Overlap criterion** — "consecutive intervals overlap iff a > q−1"
   checked exactly for q ≤ 120, a < 2q + 60 (21,539 pairs), including
   the touching equality at a = q − 1.
3. **Farey sweep** — for *every* reduced rational t = p/r with
   r ≤ 300 (27,397 values), the constructive witness satisfies the
   strict inequalities in both layers; for r ≤ 60 an independent
   brute-force scan (minimal a+b, integer layer only, sharing no
   formulas with the witness) also succeeds and its witness is never
   longer.
4. **Endpoint adversaries** — 2,494 exact boundary rationals: the
   endpoints of I(a, q) themselves (where I(a, q) fails by strictness —
   verified to fail) and the points 1/q, 2/q; every one is recovered by
   a different member, as the proof requires.
5. **Edge stress** — denominators to 2·10^50 near both ends of (0, 1),
   plus 2,000 random rationals with denominators to 10^40.
6. **Informational** — on 006's 149 half-degree arcs, the constructive
   witness's a+b is within ×15.33 of the brute-force minimum (mean
   ×1.88); the minimum governs the true minimal word length 4(a+b)+2.

Reproduce:

    python problems/billiards-triangles/explore/coverage_diophantine.py selftest
    python problems/billiards-triangles/explore/coverage_diophantine.py gamma 135

Consistency spot-check with prior art: gamma = 135° gives t = 1/2, and
the constructive witness is (a, b) = (5, 2) — precisely the W(5,2),
length 30, that 006 certified alive at exactly 135° (and 007
re-verified). The witness for gamma = 144° (t = 2/5) is W(4, 3), 006's
pinch-gap member family; for t → 0 (gamma → 180°) witnesses grow like
a ≈ q²/(numerator scale), matching 001's observation that angle reach
costs length.

### Novelty check

Two targeted literature searches (2026-08-15; general web + arXiv-heavy
results) for the lemma and for its unit-fraction reformulation
(0 < 1/a + 1/q − t < 1/(aq), i.e. a two-unit-fraction over-approximation
with error below the product of the denominators) found no statement of
this lemma. Nearby standard material: Dirichlet's approximation theorem,
greedy/Sylvester unit-fraction expansions, Erdős–Stein on sums of
distinct unit fractions — none is this statement, though the proof
technique (overlapping mediant-style interval chains) is entirely
classical. **The lemma is recorded as elementary and likely folklore,
not as new mathematics**; the contribution is closing 006's labelled
sublemma with a proof and an exact constructive checker.

## Outcome

- **VERIFIED (proof + exact computation, scope as stated):** the lemma —
  for every t ∈ (0,1) there exist integers a, b ≥ 1 with
  a+b < t·a(b+1) < a+b+1 — has an elementary proof, re-derived here
  independently of the handoff that proposed it, every algebraic
  equivalence in it additionally machine-checked in exact arithmetic,
  and the constructive witness verified exactly for all 27,397 reduced
  rationals with denominator ≤ 300, 2,494 exact boundary adversaries,
  and denominators to 2·10^50 (two independent arithmetic layers; the
  brute-force layer shares no formulas with the construction).
  Per-lab-convention caveat: this record's own proof has not yet had an
  adversarial skeptic pass (lead 1); treat the proof's VERIFIED as
  carrying that standing obligation.
- **CONDITIONAL, and only conditional:** *if* 006's birth law
  gamma_birth(a,b) = 180 − 90(a+b+1)/(a(b+1)) and exact-window claim
  hold (both still SPECULATION), *then* every obtuse gamma ∈ (90°, 180°)
  lies strictly inside the window of the explicit member
  W(a, b) above. The Diophantine side of queue item 12 is closed;
  the load-bearing open problem is now entirely the birth/sufficiency
  theorem (queue item 11).
- **NOT claimed:** the birth law; the exact-window (interior aliveness)
  claim; the coverage conjecture itself (do not upgrade it — it
  inherits SPECULATION from the window law); anything about arcs as
  intervals (the certificates in 006/007 are pointwise; per-triangle
  arc coverage is open — 007's C2 stands); optimality of the
  constructive witness (it is provably suboptimal, ×15 on one grid
  arc); any claim about the a = 1 column being needed or not (the
  witness always has a ≥ q = b+1 ≥ 2, so the lemma never needs a = 1).

## Why it failed / what survived

Nothing failed. What the result changes: 006 lead 2 guessed the
existence question was governed by continued fractions; it is not —
existence is a two-line interval-chain argument, and continued
fractions matter only for the *optimal* witness. The interesting
residue, made precise by check 6: the constructive witness overshoots
the minimal a+b by up to ×15 on 006's own grid, exactly at arcs just
above a window corner 90/j (t just above 1/q), where the constructive
chain enters at huge a while a much smaller member from a *different* q
column covers the same t. So the minimal-length staircase (006 lead 6)
is genuinely a different, still-open computation — this lemma bounds it
above but does not compute it.

Reusable: the checker (exact witness for any rational t or rational
gamma, integer-layer verifier usable as a component in any future
window-arithmetic tool); the sharpness fact that the a ≥ q cutoff is
exact (chain touches at a = q − 1) — any future tightening of the
window law that shifts an endpoint by even one lattice step will break
the chain, so the checker's endpoint-adversary suite is the regression
test to keep.

## Leads generated

1. **Skeptic pass on this record** (default stance: refute). Attack
   surface, in order: (a) the strictness bookkeeping in Step 4 and in
   the witness-bound proofs (the 1/t-integer and (q−1)/(qt−1)-integer
   edge cases); (b) re-implement the checker's witness from the record's
   formulas alone and diff against the committed one on the Farey sweep;
   (c) check the Step 2/Step 3 union argument covers interval endpoints
   interior to the union (a point equal to some lower endpoint must lie
   in the *next* interval — verify the inequality used is the right
   one); (d) the novelty claim (find a citation; if found, re-file as
   rediscovery).
2. **Birth law theorem (= queue item 11(ii)), now the sole blocker for
   coverage.** With this lemma, a proof of the birth law + interior
   aliveness upgrades the coverage conjecture for W windows immediately.
   Falsifiable as in 006 lead 1.
3. **Minimal-witness staircase.** Compute min a+b over ALL valid (a, b)
   per arc (the checker's brute-force column does this for rational t)
   and derive the continued-fraction law 006 lead 2 guessed — now a
   clean standalone question about the interval chains, decoupled from
   existence. Concrete start: prove or refute that the minimal witness
   always has q ∈ {⌊1/t⌋+1, ⌈2/t⌉−1} or one of the two neighboring
   columns.
4. **a = 1 column irrelevance** (a one-line fact for the reviewer to
   confirm, not an open lead): a = 1 requires 1 + b < t(b+1), and since
   1 + b = b + 1 this forces t > 1 — so no t ∈ (0,1) ever has an a = 1
   witness, and the lemma's a ≥ 1 hypothesis is effectively a ≥ 2
   (indeed a ≥ q ≥ 2 in the construction).

## References

- `problems/billiards-triangles/attempts/006-design-family-past-135.md`
  (lead 2, the target sublemma; birth law SPECULATION).
- `problems/billiards-triangles/attempts/007-skeptic-review-of-006.md`
  (C2: pointwise-vs-arc caveat inherited here).
- `problems/billiards-triangles/attempts/005-complete-death-law-theorem.md`
  (death law, the proven half of the window).
- External handoff: user-supplied ChatGPT session output proposing the
  proof (unpublished; treated as unverified input and re-derived).
- Literature consulted in the novelty check (none contains the lemma):
  Dirichlet approximation theorem (standard); Erdős & Stein, *Sums of
  distinct unit fractions*, Proc. AMS 14 (1963).
