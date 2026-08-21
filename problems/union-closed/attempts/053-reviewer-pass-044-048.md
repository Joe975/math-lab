# 053 — Reviewer pass on 044–048: the certificate and every checkpoint hold; eight corrections, two of them mathematical (N2-CONC's tightness clause, L1's boundary identity)

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** reviewer pass on 044, 045, 046, 047 and 048 (reviewer batch 5;
  default stance REFUTE).
- **Numbering note:** this review was commissioned as "049" and its
  re-implementation is therefore named `uc_reviewer049_reimpl.py`, but
  050–052 were written and committed while it ran, so
  `scripts/new_attempt.py` assigned 053. The file name is left as
  committed rather than renamed; 049 is an unused id.
- **Independence disclosure:** fresh-context subagent with zero shared
  conversation state, but spawned from the author's session with a
  reviewing brief, and the same model family — the caveat recorded in
  029/034/039/043. Whether that meets 024's fresh-session bar is for the
  human owner to judge. Everything below was re-derived or re-computed
  here; nothing was accepted on the records' say-so.
- **Tools:** `explore/uc_reviewer049_reimpl.py` (new; 78 checks, exit 0
  iff all pass; imports nothing from `uc_hu_ownconst`, `uc_hu_bestorder`,
  `uc_hu_n2`, `uc_hu_n2_dichotomy`, `uc_hu_L1`, their skeptics, or the
  `uc_hu_order2` / `uc_hu_certify` / `uc_reviewer036` stack — the HU
  coupling, CR, the canonical and rollout rules, c\*, the n = 2 closed
  form and branch L1 are all rebuilt from the records' prose, with an
  exact-rational history-pair evaluator whose h() is 80-digit
  `Decimal.ln()`, i.e. neither of the repo's two kits). Runtime ≈ 1 min.
  The ten pipeline re-runs (below) were done outside that file.
- **Sources:** none beyond the repo.

## Approach

Batch 5 covers one certified refutation (044), one survival claim (045)
and three proof-probe records (046/047/048) whose content is mostly
hand argument. So the effort went, in order: re-derive the five hand
proofs from scratch before reading how the records argue them;
recompute the 044 certificate at 50+ digits with an independent
evaluator; re-run every committed pipeline and byte-compare; and audit
every quoted number against the checkpoints.

The reason to re-derive first rather than read first: four of the five
statements are short enough that an independent derivation is cheaper
than checking someone else's, and it is the only way to catch a
*clause* that is wrong inside a lemma that is otherwise right — which
is exactly what turned up twice.

## What was done

**P1. The five hand statements, re-derived independently.**

- **046 N2-ONE-BAD — correct.** If u₀, u₁ < 1/2 then each of 1−u₀,
  1−u₁ > 1/2, so f₁ = x(1−u₀)+(1−x)(1−u₁) > 1/2 for every x ∈ [0,1]
  (convex combination of two numbers each > 1/2; degenerate x ∈ {0,1}
  included). Confirmed on 20,000 in-regime samples, 0 exceptions.
- **047 N2-ONE-ABOVE — correct.** x·p₀ + (1−x)·p₁ = f₁ ≤ q exactly, so
  p₀, p₁ > q is impossible for x ∈ (0,1); at x ∈ {0,1} one conditional
  is unconstrained but carries zero weight, so the statement is
  degenerate rather than false. 0 exceptions in 20,000 samples.
- **046 N2-CONC — the inequality and the averaging step are correct;
  the tightness clause is FALSE.** ψ(t) = h(min(1/2,t)) is concave as
  h|[0,1/2] is concave nondecreasing and min(1/2,·) is concave, so
  ψ(p_a+p_b) ≥ (ψ(2p_a)+ψ(2p_b))/2, which is exactly
  s(p_a,p_b) ≥ (G(p_a)+G(p_b))/2. The averaging step is right for the
  reason the record gives: the cell weights (z₀, x−z₀, x−z₀, 1−2x+z₀)
  have row marginals (x, 1−x) *and* column marginals (x, 1−x)
  (re-checked numerically), so Σ w_ab (G(p_a)+G(p_b))/2 collapses to
  x G(p₀) + (1−x) G(p₁). **But equality does not require p_a = p_b:**
  ψ is *constant* on [1/2, ∞), so whenever min(p_a,p_b) ≥ 1/4 all three
  ψ-values are 1 and the bound is tight with p_a ≠ p_b. Worked
  instance: p_a = 0.30, p_b = 0.40 gives s = 0.073879… =
  (G(0.30)+G(0.40))/2 exactly. The correct clause is *equality iff
  p_a = p_b **or** min(p_a,p_b) ≥ 1/4*, verified on 200,000 random
  pairs with zero mismatches.
- **046's exact decomposition — correct, in both cases, and its
  σ(p) = 0 above 1/2 convention is right.** Independent derivation:
  with t := min(x−1/2, 1−x) the off-diagonal cell weight and
  w₀₀ = x−t, w₁₁ = (1−x)−t, the Case-A sum
  Σ w_ab ψ(p_a+p_b) regroups as xψ(2p₀)+(1−x)ψ(2p₁)+t·Δ, giving
  F = (c\*(f₀)−c\*(q))h(f₀) + [xσ(p₀)+(1−x)σ(p₁)] + tΔ −
  c\*(q)[x h(p₀)+(1−x)h(p₁)]. In Case B (say p₁ > 1/2) the clip's
  `min(·, u_a, u_b)` arm fires, the (1,1) cell contributes s = 0 and
  the off-diagonals (h(u₁)−h(u₀))/2, and the same identity holds with
  σ(p₁) = 0 and Δ = h(p₁) − ψ(2p₀) — i.e. **Δ must be read as the
  ledger form 2s(u₀,u₁) − s(u₀,u₀) − s(u₁,u₁), not as the ψ closed
  form, which is Case-A only.** Verified to 6.5e-16 over 20,000
  in-regime samples spanning both cases (9,952 of them Case B).
- **048's boundary identity — exact, but only for q ≤ 1/4.** At
  (p₀,p₁) = (0,1/2): Δ = 2h(1/2) − h(0) − h(1/2) = 1 so C = t, and
  −B = q·c\*(q) since h(0) = 0, h(1/2) = 1, c\*(1/2) = 0. Hence
  ratio = t/(q c\*(q)), and t = min(1/2−q, q) equals q **iff
  q ≤ 1/4**. Above 1/4 the ratio is (1/2−q)/(q c\*(q)): at q = 0.3 it
  is 4.94930, not 1/c\*(0.3) = 7.42395. Checked in exact
  rational/60-digit arithmetic at eight values of q.

**P2. The c\* correction of 046 F, both directions.**

- (a) The stated form (h(min(2p̄,1))−h(p̄))/h(p̄) does differ from the
  computed one above p̄ = 1/4 and is strictly smaller there, and all
  seven numbers 046 F quotes reproduce. **But it is negative only above
  p̄ = 1/3**, not above 1/4: h(2p) = h(p) ⟺ 2p = 1−p ⟺ p = 1/3.
  On (1/4, 1/3) the stated form is positive and merely wrong —
  stated(0.30) = +0.101737 against the correct +0.134707 — so the
  record's "above it … it is **negative**" and "would make (HU-TAX)
  vacuous exactly on the range 034 re-posed it for" overshoot on that
  subinterval.
- (b) The reassuring half is clean. Every `cstar`/`cst`/`zstar`/`zst`/
  `sig`/`sigma` definition in `explore/` (17 scanned) computes the
  max(1/2, 1−2p) form. A repo-wide scan for the mis-stated form finds
  exactly one live occurrence — a diagnostic f-string in
  `uc_reviewer034_reimpl.py:733` — which feeds no checkpoint and no
  assertion (the assertion beside it uses
  `z = min(max(0.5, 1−2p), 1−p)`). **No computed number in the line
  depends on the mis-stated formula**, as 046 says.

**P3. Ten committed pipelines re-run; six checkpoints byte-compared.**
`uc_hu_ownconst.py` (4 s), `uc_hu_ownconst_certify.py` (0 s),
`uc_hu_bestorder.py` (55 s), `uc_hu_bestorder_skeptic.py` (2 s, exit 0),
`uc_hu_n2.py` (819 s), `uc_hu_n2_skeptic.py` (exit 0),
`uc_hu_n2_dichotomy.py` (7 s), `uc_hu_n2_dichotomy_skeptic.py` (13 s,
exit 0), `uc_hu_L1.py` (11 s), `uc_hu_L1_skeptic.py` (1 s, exit 0).
`hu_ownconst.json`, `hu_ownconst_certify.json`, `hu_bestorder.json`,
`hu_n2.json`, `hu_n2_dichotomy.json`, `hu_L1.json` all came back
**byte-identical** (md5 against a pre-run backup, and `git status`
clean). No nondeterminism anywhere in the batch.

**P4. The 044 certificate, recomputed at 50+ digits.** Own exact
evaluator on the same `limit_denominator(1e7)` rationalization:

    CR_roll = -0.0009596936656816684080411028274699116639225834836287452
    CR_best = +0.01610189393146184153930553062753534528530869411207410
    H       =  3.137336715625876600621222310634776141735427921025340

Both values lie inside **both** kits' enclosures in
`hu_ownconst_certify.json`. The rollout order (1,0,2,3) was re-derived
here from the 035/037 rules (canonical = greedy min H(A_i|A_S), ties to
lowest index; rollout = maximise full CR of the canonical completion)
without importing the engine's; the best order (0,3,2,1) reproduces by
enumeration; max marginal 0.49499665 < 497/1000; 18 of the 24 orders
are negative. The witness **is** already present in 037's committed
checkpoint (`hu_order2.json` → `D_bestorder_cap_0.497` → start
`floor:windowkill`, 9 atoms, n = 4), exactly as 044 claims. Beyond the
record: all 120 `hu_order2.json` endpoints were re-audited here with an
own rollout derivation and an own evaluator — max deviation from
`hu_ownconst.json` 6.4e-16, and the re-audit finds **the same single
violation and no other**.

**P5. The numbers, re-implemented.** 045: all 67 endpoint rows
re-scored by full 24/120-order enumeration through this file's
evaluator (max |Δratio| 7.2e-16), zero kills, zero sharp violations,
the three global floors exact (+2.190730e-4 / +2.190730e-4 /
+1.295407e-4), twelve rows saturating their own constant at
+2.6e-9…+2.7e-9, the (d,4)-family perturbation at own margin exactly 0.
046: identity-order margin ≥ 0 on an own 91³ grid (worst −2.8e-16, at a
product), B&B volume bookkeeping self-consistent and matching the
record's 88.4% / 11.6% / zero residue, (**) failing on 21.1% of an own
sampler's Case-A points with the original margin positive there. 047:
the dichotomy census re-run on an own sampler and evaluator —
**needing both = 0, failing = 0** — with the q = f₀ branch's C-failures
0 on 7,548 cases. 048: ratio floor 1.1764 on an own 400k branch sample,
P4's 0/N domination, and P5's two-regime minimiser table reproduced
row-for-row on an own 140² mesh.

**P6. Reporting audit.** Everything checked against the committed
checkpoints; the discrepancies are in "Outcome" below. 046's own
self-correction (G′) is accurate: the committed engine really does use
a largest-box-first heap and really does produce the quoted volumes,
and an own float-enclosure B&B of the same geometry at a 200,000-box
budget reproduces the lesson qualitatively — depth-first certifies
0.001% of the root box and leaves 99.9% unprocessed, largest-box-first
certifies 40.6% and leaves 12.8%.

## Outcome

- **VERIFIED: 044's refutation, in full and at 50+ digits.** The
  certificate, the witness provenance in 037's checkpoint, the order
  derivations, the marginal bound, the positive control, and the
  212-endpoint audit (independently re-run over the 120 rows whose
  measures are recoverable) all hold. Rollout-order HU positivity above
  cap ≈ 0.495 stays REFUTED.
- **VERIFIED: 045 in full.** Every endpoint, floor, kill list and
  structural observation reproduces through an independent evaluator.
- **VERIFIED: 046's certified and computational content** — the closed
  form, N2-ONE-BAD, the decomposition identity (in both cases), the
  equality-set check, the B&B accounting, part E's measured failure of
  (**), part D's first-order transversal growth, part G's boundary
  family (+7.0998e-2 / +2.0970e-3 / +7.7494e-11, and CR = 0 exactly at
  the two disjoint singletons), and the G′ method lesson.
- **VERIFIED: 047's dichotomy and N2-ONE-ABOVE**, and **048's** P1/P4/P5
  and the logarithmic approach of c\*(q) → 1 (checked in exact
  arithmetic down to q = 10⁻¹⁰⁰).
- **Correction 1 (046, Lemma N2-CONC, propagated to the index's
  `one_line`): "tight iff p_a = p_b" is false.** Equality also holds
  whenever min(p_a,p_b) ≥ 1/4, because ψ is constant on [1/2,∞) —
  e.g. (0.30, 0.40). The lemma's inequality, its averaging step and
  every use made of it survive; the clause does not, and 048's "the
  concavity route gives a bound tight only at p₀ = p₁" inherits the
  same imprecision (its conclusion — that the bound is *not* tight at
  L1's extremal (0,1/2) — is correct: Δ = 1 there).
- **Correction 2 (048, part P2 and the index's `one_line`): the
  boundary identity holds only for q ≤ 1/4.** ratio(q,0,1/2) =
  t/(q c\*(q)) with t = min(1/2−q, q); above q = 1/4 that is
  (1/2−q)/(q c\*(q)) ≠ 1/c\*(q) (4.94930 vs 7.42395 at q = 0.3). Both
  the engine (six q ≤ 0.1) and its skeptic (q ≤ 0.2) only probed
  q ≤ 1/4, so the gap survived the skeptic pass. Everything the record
  concludes from the identity survives: the corner still satisfies L1
  strictly at those q, and the infimum-1 story lives at q → 0.
- **Correction 3 (046, part F and the index's `corrections`): the
  mis-stated constant is negative only above p̄ = 1/3**, not above 1/4;
  on (1/4, 1/3) it is positive but too small (+0.101737 vs +0.134707 at
  0.30). The correction's substance — wrong above 1/4, right below,
  and used by no computation — stands.
- **Correction 4 (047, scope): the record's "A ≥ 0 … C ≥ 0 … only B can
  be negative" and its whole census are Case-A only.** The engine draws
  p₀,p₁ ∈ (0,1/2); outside Case A — 50.2% of in-regime measures under
  the natural uniform (x,u₀,u₁) parametrization — the ledger
  interaction C = tΔ is *negative* on 63% of sampled Case-B points, and
  047's closed-form A+B+C is not the margin there at all (deviation up
  to 0.5). The index's `range` field does say "Case A"; the record's
  §A and its Outcome bullet ("the deficit structure at n = 2 is now
  completely determined") do not. Not a threat to the n = 2 conjecture:
  the margin stays ≥ 0 over 200,000 Case-B samples (worst +1.6e-5), and
  046's B&B covers the whole box.
- **Correction 5 (048, part P3 and the index's `one_line`): "375
  versus 0" is a `--fast` run.** The committed checkpoint says
  **5,698 versus 0**; an own full re-run of the same generator gives
  5,698. The qualitative claim (no tight configuration on the
  t = x−1/2 arm) is unaffected.
- **Correction 6 (046, parts A and B, and the index's `one_line` and
  `range`): the "1.38M-point grid" is a superseded run.** The committed
  skeptic's S4 sweeps N = 60 — 219,539 points, worst −6.939e-16, which
  is what 046 part C quotes as −6.9e-16. 1,379,840 = 110·112·112 is
  exactly the N = 111 grid, i.e. an earlier skeptic. Part A's
  "−6.7e-16" matches neither committed run.
- **Correction 7 (044, part B and the index's `one_line`/`leak_terms`):
  "rollout ranks 13/24" is not reproducible** under a best-to-worst
  ranking. All six orders that reveal coordinate 1 first give *exactly*
  the same CR (−9.596936657e-4), so rollout sits in a six-way tie at
  descending ranks 7–12; 13 is that block's first position counting
  from the *worst* end, under which rollout is above the median — the
  opposite of the record's gloss "worse than the middle".
- **Correction 8 (044, part A): the endpoints sitting at exactly 0 are
  the `hu_blocks.json` block-tensor rows** (`B_attack:[(d,2),(d,2)]`,
  `B_attack:[(d,3),(d,3)]` to 5e-20, `C_sat497`), **not "the
  family-saturated anneal endpoints"** — the closest anneal endpoint is
  at +1.53e-11, near the family but not on it. The audit's substance
  (211 of 212 clear) is unaffected.
- **Reporting nits (no number moves).** 047's parts B and C are
  different samples, so their deficit totals do not agree (42,575 vs
  7,435 + 34,650 = 42,085) and the prose reads as if they were one
  census. 046 part E's "identity verified to 2.6e-15 over 99,970
  in-regime samples" corresponds to no committed check (the identity
  is true — verified here to 6.5e-16 — but the committed cross-check is
  047's skeptic S1, 2.89e-15 over 30,000 samples). 045's "+3.07e-3" for
  the 044 witness is the cap-0.499 row; at caps 0.495/0.497 it is
  +5.13e-3. 048's P5 calls the q = 0.3 minimiser (0.42857, 0.0)
  "interior" though p₁ = 0 is on the branch boundary.
- **Not refuted:** nothing. No claim in the five records is
  contradicted by a number; the two mathematical corrections are a
  false tightness clause inside a true lemma and a missing hypothesis
  on a true identity.
- **Not claimed here:** any statement at n ≥ 3 (050/051 landed while
  this ran and were not reviewed); that the descent adversaries of 045
  are strong; that 046's B&B covers more than its measured 88.4%; that
  the six-way tie at the 044 witness is generic.

## Why it failed / what survived

The batch survives its adversarial pass on every number, which is the
part that matters for 044's REFUTED verdict — that one is now checked
twice by disjoint code at 50+ digits, plus an independent re-derivation
of the order rule that names it.

What did *not* survive is a pattern worth naming, because it repeated
three times in one batch: **a correct statement carrying an incorrect
sharpness clause or a missing hypothesis.** N2-CONC's "tight iff
p_a = p_b" is wrong on a positive-measure set; L1's boundary identity is
right only below q = 1/4; 046 F's own correction of 034 over-reaches by
a factor of the interval (1/4, 1/3). Each is invisible to the engines,
because each engine samples exactly the region where the clause happens
to be true — Case A for the concavity clause, q ≤ 0.2 for the identity,
p̄ ≥ 0.38 for the constant — and each skeptic inherited the sampling
region from the record it was checking. The lesson for this repo's
skeptic protocol is specific: **a skeptic that reuses the engine's
sampling region cannot see a missing hypothesis**, so the first thing a
skeptic should vary is the domain, not the arithmetic.

The two ancestry problems (a `--fast` number and a superseded grid size
quoted as if from the committed run) are the mirror image: both are
cases where the prose was written from a transcript that the final
committed run then superseded. Since prose numbers are what the index
and the website carry, quoting from the checkpoint rather than the
scrollback is the cheap fix.

## Leads generated

1. **Redo 047's census over Case B.** The decomposition there is
   F = G(f₀) + xσ(p₀) + t·[h(p₁) − ψ(2p₀)] − c\*(q)H, the interaction
   term is frequently negative, and the dichotomy as stated does not
   apply. Either a Case-B analogue of the dichotomy exists (then n = 2
   is one theorem) or the two cases need different arguments — worth
   knowing before L1 is invested in further.
2. **Re-probe L1's corner at q > 1/4.** P5 says the minimiser leaves
   the corner above q ≈ 0.03, and the corner's value there is
   (1/2−q)/(q c\*(q)), not 1/c\*(q); the "two regimes" picture should
   be redrawn with the correct corner curve before lead 1 of 048
   (a monotonicity argument in p₁) is attempted.
3. **Sharpen N2-CONC and use it.** The corrected equality set of the
   cell bound — p_a = p_b *or* min(p_a,p_b) ≥ 1/4 — is itself
   structure: the interaction Δ vanishes identically on the whole
   region min(p₀,p₁) ≥ 1/4, so on that region the theorem must be
   carried by A and B alone. That is a clean sub-case someone can close
   by hand.
4. **Standing skeptic rule (proposed):** every skeptic must sample at
   least one region the engine does not, and say which. Three of this
   batch's corrections would have been caught by that one line.

## References

- This repo: 044/045/046/047/048 (the records reviewed), 043/039/036/034/029
  (the reviewer-pass pattern and the independence caveat this repeats),
  042/041/040/038/037/035 (the endpoints, orders and families re-scored
  here), 031 (the HU coupling definition used to rebuild the evaluator),
  009 (CR). Checkpoints: `data/hu_ownconst.json`,
  `data/hu_ownconst_certify.json`, `data/hu_bestorder.json`,
  `data/hu_n2.json`, `data/hu_n2_dichotomy.json`, `data/hu_L1.json`,
  `data/hu_order2.json`.
- No external sources.
