# 014 — Skeptic review of 013 (coverage Diophantine lemma): adversarial verification

- **Problem:** billiards-triangles, `problems/billiards-triangles/PROBLEM.md`
- **Date:** 2026-08-16
- **Mode:** informed
- **Type:** adversarial verification of `013-coverage-diophantine-lemma.md`
  (default stance: refute), taking exactly the attack surface 013's lead 1
  pre-registered: (a) strictness bookkeeping in Step 4 and the witness-bound
  proofs, (b) checker re-implementation from the record's formulas alone,
  (c) the Step 2/Step 3 union argument at interval endpoints, (d) the
  novelty claim. Plus 013 lead 4 (the a = 1 one-liner) and every
  consistency spot-check in the record.
- **Outcome in one line:** the lemma, its proof, its constructive witness,
  and the committed checker all SURVIVE — every step re-derived by hand,
  every count reproduced by a from-scratch implementation, witness-for-
  witness identical on all 27,397 Farey rationals — but the record's prose
  contains two errors: the γ = 144° consistency spot-check names a member
  (W(4,3)) that is not a witness there at all, and the Step 2 parenthetical
  about I(q−1, q) asserts a false inclusion. Neither is load-bearing.
- **Tools:** `explore/cdsk_review.py` (stdlib only, deterministic, seeded
  RNG; `selftest` ~50 s, `diff` ~5 s). Independence protocol: everything
  except the diff harness was written from 013's prose alone, *before*
  reading `explore/coverage_diophantine.py`; the diff harness (which
  imports the committed module) was added only after the independent
  selftest passed. Committed selftest also re-run as published (PASS,
  73,542 checks).
- **Sources:** two fresh literature searches (2026-08-16, queries disjoint
  from 013's: the two-unit-fraction over-approximation form, and the raw
  inequality form); no citation found.

## Claims attacked

1. **The proof, Steps 1–4** — in particular the strictness bookkeeping:
   Step 4's "an open interval of length > 1 contains an integer" and the
   two integer edge cases in the witness-bound proofs (1/t an integer;
   (q−1)/(qt−1) an integer).
2. **The union argument (Step 2/Step 3)** — does ⋃_{a≥q} I(a,q) really
   equal (1/q, 2/q) *including* the points that sit on a lower endpoint of
   some I(a, q) interior to the union (013 lead 1c)?
3. **The constructive witness formulas** and the claimed automatic a ≥ q.
4. **The checker** — re-implement witness + both arithmetic layers + the
   brute-force column from the record's formulas alone, re-run every
   numbered check with different code (and different seeds where the check
   is randomized), diff the witnesses on the full Farey sweep r ≤ 300.
5. **The record's headline numbers** — 27,397 Farey values, 21,539 overlap
   pairs, the ×15.33 / ×1.88 overshoot statistics of check 6.
6. **The consistency spot-checks** — γ = 135° ↦ W(5,2), γ = 144° ↦ W(4,3).
7. **The novelty claim** ("elementary and likely folklore, no citation
   found").
8. **013 lead 4** — no t ∈ (0,1) ever has an a = 1 witness.

## Refutations found

### R1. The γ = 144° consistency spot-check is false: W(4,3) is not a witness there

013 ("Consistency spot-check with prior art"): *"The witness for
gamma = 144° (t = 2/5) is W(4, 3), 006's pinch-gap member family."*

Three independent computations refute the sentence:

- **W(4,3) fails the lemma's inequality at t = 2/5 outright:**
  a + b = 7 but t·a(b+1) = (2/5)·16 = 32/5 = 6.4 < 7. Its interval is
  I(4,4) = (7/16, 1/2), i.e. its (conditional) window is
  γ ∈ (135°, 140.625°) — it cannot contain 144°.
- **The actual constructive witness at t = 2/5 is (a,b) = (11,2)** (word
  length 54): q = ⌊5/2⌋+1 = 3, a = ⌊2/(1/5)⌋+1 = 11; check
  13 < 66/5 = 13.2 < 14. The *committed checker itself agrees*:
  `python explore/coverage_diophantine.py gamma 144` prints W(11,2) — so
  the error is purely in 013's prose, not in the code or the proof.
- **The minimal witnesses at t = 2/5 have a+b = 9** ((6,3) and (4,5), by
  the independent brute-force scan), so (4,3) is not the minimal witness
  either.

Sharper version, from the same computation: **W(4,3) is never the
constructive witness for any t.** The construction needs q = 4, i.e.
t ∈ (1/4, 1/3], while a = 4 = ⌊3/(4t−1)⌋+1 needs t ∈ (7/16, 1/2] — the
ranges are disjoint. Likely cause of the slip: W(4,3) is 006's
pinch-gap member, whose window starts at exactly 135° — it belongs to the
γ = 135° discussion one sentence earlier, not to 144°. Not load-bearing:
nothing downstream cites the sentence.

### R2. The Step 2 parenthetical about I(q−1, q) asserts a false inclusion

013 Step 2: *"(I(q−1, q) ⊂ (1/q, 2/q) contributes no new points… its
upper endpoint 1/q + 1/(q−1) exceeds 2/q, …)"*. The two clauses
contradict each other, and the first is false: the lower endpoint of
I(q−1, q) is 1/q + (q−1)/((q−1)q) = 2/q exactly, so

    I(q−1, q) = ( 2/q , 1/q + 1/(q−1) )

is **disjoint from** (1/q, 2/q), not a subset of it (machine-checked
exactly for q ≤ 120; e.g. q = 5: I(4,5) = (2/5, 9/20)). The correct
statement — smaller a only adds points ≥ 2/q, which the choice of q in
Step 4 makes irrelevant — is what the rest of the same parenthesis
already says. The committed checker's docstring carries a shard of the
same garble ("the chain TOUCHES at a = q-1 (upper endpoint of
I(q, q-1)...)" — mis-ordered arguments). Not load-bearing: the proof
only ever uses a ≥ q, and the touching fact it *does* use (equality of
the relevant endpoints at a = q−1) is true and confirmed.

No other refutations. In particular the proof itself, all four steps,
survives unmodified.

## Claims that survive

1. **The lemma and the proof (attacks 1–3): CONFIRMED.** All four steps
   re-derived by hand from the record's prose. Step 1's algebra checks
   (interval width is exactly 1/(aq) > 0). Step 2's overlap criterion
   re-derived: 1/(a+1) > (q−1)/(aq) ⟺ a > q−1, with exact touching at
   a = q−1 (re-checked exactly at 21,539 (a,q) pairs — same count as
   013). Step 3 was attack (c): the union argument is sound, and a
   cleaner form of it closes the endpoint worry — for t ∈ (1/q, 2/q)
   take the *minimal* a₀ ≥ q with lower(a₀) < t (exists since
   lower(a) ↓ 1/q); if a₀ = q then t < 2/q = upper(q,q), else
   t ≤ lower(a₀−1) < upper(a₀) by the overlap inequality, so
   t ∈ I(a₀, q) either way. Machine side: 27,810 exact points including
   every lower endpoint of I(a,q), a ∈ [q, q+80], q ≤ 60 — each is
   caught by a strictly later interval, zero failures. Step 4's
   bookkeeping survives both integer edge cases: any open interval of
   length > 1 contains ⌊x⌋+1 (strict at both ends even when x ∈ ℤ);
   q ≤ 1/t + 1 < 2/t is strict even when 1/t ∈ ℤ because t < 1; and
   a = ⌊L⌋+1 with L = (q−1)/(qt−1) satisfies L < a ≤ L+1 < L + 1/(qt−1)
   strictly even when L ∈ ℤ, because qt−1 < 1. The automatic a ≥ q
   (L > q−1 since qt−1 < 1) also re-derived.
2. **The checker (attacks 4–5): CONFIRMED, witness-for-witness.**
   `cdsk_review.py selftest` — written from the record's formulas before
   reading the committed code, different structure (Fraction-floor via
   integer division rather than `__floor__`, brute force scanning b
   ascending rather than a, different seeds) — passes with **zero
   failures** and reproduces 013's counts exactly: 27,397 Farey
   witnesses (r ≤ 300) valid in both layers with a ≥ q; brute-force
   cross-check r ≤ 60 never beaten; overlap criterion at 21,539 pairs;
   my own endpoint-adversary suite (3,196 boundary rationals, a larger
   net than 013's 2,494) all recovered by a different member with
   strictness failing on the generating member; edge stress to
   denominator 2·10⁵⁰ and 2,000 fresh random rationals to 10⁴⁰. The
   diff harness then compared witnesses on **all 27,397 reduced
   rationals with r ≤ 300: zero mismatches** against the committed
   `witness()`. The committed selftest also re-runs green as published
   (73,542 checks). Reproduce:
   `python problems/billiards-triangles/explore/cdsk_review.py selftest`
   then `... diff`.
3. **Check 6's overshoot statistics: CONFIRMED.** On the 149 half-degree
   arcs, worst constructive/minimal a+b ratio ×15.33 at γ = 134.5°
   (constructive (91,1) vs minimal (4,2)), mean ×1.88 — matching 013 to
   the printed digits, from the independent implementation. The
   mechanism 013 names (t just above 1/q forces the constructive chain
   in at huge a) is visible in the worst case: t = 91/180, barely above
   1/2.
4. **The γ = 135° spot-check: CONFIRMED.** t = 1/2 ↦ (5,2), exactly
   006/007's certified-alive W(5,2).
5. **Lead 4 (a = 1 irrelevance): CONFIRMED.** a = 1 needs
   1 + b < t(b+1), i.e. t > 1 — impossible on (0,1); 5,000 randomized
   integer-layer trials found no spurious witness. The lemma's a ≥ 1
   hypothesis is indeed effectively a ≥ q ≥ 2.
6. **The novelty label: SURVIVES AS WRITTEN.** Two fresh searches (the
   0 < 1/a + 1/q − t < 1/(aq) unit-fraction form; the raw a+b <
   t·a(b+1) < a+b+1 form with interval-chain/mediant keywords) found no
   statement of the lemma — nearest hits are standard Farey/mediant and
   Egyptian-fraction material. 013's calibrated claim ("elementary and
   likely folklore, not new mathematics; the contribution is the closed
   sublemma + exact checker") is the right label, and this pass found
   nothing requiring a re-file as rediscovery.
7. **Scope discipline: CONFIRMED.** 013's conditional framing (the
   billiards meaning rests entirely on 006's SPECULATION birth/window
   law; pointwise-vs-arc caveat inherited from 007 C2) is stated
   correctly everywhere it matters, including in the committed checker's
   user-facing output.

Status consequence: with this pass, 013's `VERIFIED` no longer carries
the "skeptic pass pending" obligation. The Diophantine side of queue 12
is closed for real; the sole blocker for conditional coverage remains
the birth/sufficiency theorem (queue 11), and the minimal-witness
staircase (013 lead 3) stays open — R1's data point (constructive (11,2)
vs minimal (6,3)/(4,5) at t = 2/5, where the minimal witness sits in a
*different* q column than either constructive candidate) is a concrete
seed for it.

## References

- `problems/billiards-triangles/attempts/013-coverage-diophantine-lemma.md`
  (the record under review).
- `problems/billiards-triangles/attempts/006-design-family-past-135.md`
  (lead 2, the sublemma's origin; W(4,3) pinch-gap context for R1).
- `problems/billiards-triangles/attempts/007-skeptic-review-of-006.md`
  (C2 pointwise-vs-arc caveat, confirmed correctly inherited).
- `problems/billiards-triangles/explore/coverage_diophantine.py` (committed
  checker, re-run and diffed against).
- `problems/billiards-triangles/explore/cdsk_review.py` (this review's
  independent implementation).
- Literature consulted in the novelty re-check (none contains the lemma):
  standard Farey/mediant expositions; Erdős & Stein, *Sums of distinct
  unit fractions*, Proc. AMS 14 (1963) — same nearest-neighbor 013 found.
