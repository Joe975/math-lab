# 043 — Reviewer pass on 041/042: both survive on the numbers; 042's "endpoints on the family" story is wrong for exactly the two campaigns where the polish moved

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-20
- **Mode:** informed
- **Type:** reviewer pass on 041 and 042 (queue: "joins the next reviewer
  batch"), stance refute. Covers the two hand re-derivations both records
  queued for review — 041's count-pair-DP collapse and 042's block
  additivity (with 040's DIAG identity re-derived as its ingredient) —
  plus an independent re-implementation, byte-compare pipeline re-runs,
  and a reporting audit.
- **Independence disclosure:** fresh-context subagent with zero shared
  conversation state, but spawned from the author session with the review
  brief, same model family — the caveat recorded in 029/034/039, one
  level below 036's genuinely-fresh bar. Whether that meets the 024
  fresh-session bar is for the human owner to judge.
- **Tools:** `explore/uc_reviewer043_reimpl.py` (new; imports nothing
  from the uc_hu_interp / uc_hu_blocks / uc_hu_order2 / uc_hu_canon /
  uc_hu_roll_anneal stack — the HU coupling, CR, c*, the canonical
  completion and the rollout rule rebuilt from 031/009/034/037 prose,
  with two own evaluator structures: a posterior-pair atom evaluator and
  a count-pair DP written from this review's own re-derivation, no
  sub-threshold pruning; deterministic, ~10 s; exit 0). Pipeline re-runs
  in a shell, unmodified, against backed-up checkpoints.
- **Sources:** none beyond the repo.

## Approach

Default stance refute, per the verification standard. Priority order per
the brief: (P1) hand re-derivation of 041's claim that the HU coupling on
a mixture of K product measures collapses exactly to a
(ones-in-a-prefix, ones-in-b-prefix) count-pair DP; (P2) hand
re-derivation of 042's every-order block additivity and its corollary
(every tensor of diag blocks and Bern(p) factors at common p has
CR/H = c*(p), via 040's DIAG identity, also re-derived); (P3) an
independent re-implementation reproducing sampled 041 rows, 042 part A,
and every B/C endpoint; (P4) unmodified byte-compare re-runs of all four
committed pipelines; (P5) an audit of every number in the two records'
prose and prior-art entries against the checkpoints.

## What was done

**P1 — count-pair DP re-derived; holds.** For μ = Σ_k λ_k Bern(q_k)^⊗n,
the posterior over components given a revealed prefix of length k with
j ones is ∝ λ_k q_k^j (1−q_k)^(k−j), so the predictive
P(next = 0 | prefix) depends on the prefix only through (k, j). Both
sides of a pair-history cell therefore carry the same (x, y, z) whenever
their ones-counts agree, and the 4-way HU transition weights
(z, x−z, y−z, 1−x−y+z) depend only on the count pair — the cell process
merges exactly onto (ja, jb) with no approximation. H(μ) =
−Σ_m C(n,m) p_m log₂ p_m with p_m = Σ_k λ_k q_k^m (1−q_k)^(n−m), since
all C(n,m) atoms in a ones-count class share the weight p_m. Edge cases
attacked: components with q ∈ {0,1} (path A) rely on 0⁰ = 1 in the
posterior weight — the implementation's `q ** j * (1-q) ** (k-j)` does
this correctly, and any cell of positive coupling weight has positive
prefix probability on both margins (the coupling's cell margins are μ's
conditionals), so the `den > 0` guard never masks a live cell. The
committed `hu_cr_pmix` matches this derivation line by line; its one
liberty is dropping transitions with cw ≤ 1e-18, a mass leak bounded by
~4·(k+1)²·1e-18 per step (≲1e-13 cumulative at n = 64, consistent with
the −5.3e-14 endpoint noise; my re-implementation prunes nothing and
matches to <1e-9). `comps_A`/`comps_B` reproduce the record's path
definitions; both keep every coordinate marginal exactly p
(Σ λq = p algebraically; checked to 1e-15 at the component level and
1e-12 at the atom level), and path B hits Bern(p)^n at t = 0 and exactly
the two-point diagonal at t = 1 (λ = p identically, since b − a = t).

**P2 — block additivity re-derived; holds, interleavings included.**
For μ = μ₁ ⊗ μ₂ and any revelation order π: when a block-1 coordinate is
revealed, x = μ(A_i = 0 | a) = μ₁(A_i = 0 | a↾S₁) by independence, same
for y, so z and the four transition weights depend only on the block-1
sub-cell. By induction the joint cell distribution stays the product of
the two block cell distributions (a block-1 step multiplies the block-1
factor by the exact transition the block-1 coupling alone would apply
under its induced sub-order π↾S₁, and leaves the block-2 factor fixed);
the h(z) charge of a block-j step marginalizes to the block-j process
because the other factor sums to 1. Hence
CR_HU(μ₁⊗μ₂, π) = CR_HU(μ₁, π↾S₁) + CR_HU(μ₂, π↾S₂) and H adds — the
interleaved case is not special. DIAG identity (three lines,
re-derived): on (1−p)δ∅ + pδ_S the first revealed live coordinate has
x = y = 1−p > 1/2, so z = min(max(1/2, 1−2p), 1−p) = max(1/2, 1−2p);
conditioning collapses the posterior to δ∅ or δ_S, after which
x, y ∈ {0,1} force z ∈ {0,1} and h(z) = 0 at every later step; so
CR = h(max(1/2, 1−2p)) − h(p) and H = h(p) for every order and every
n ≥ 1, both c* branches. Bern(p) is the m = 1 case. Corollary: a tensor
of B blocks (diag or Bernoulli) at common p has CR = B·(h(max(1/2,
1−2p)) − h(p)) and H = B·h(p) under every order (each block's induced
sub-order is irrelevant because the per-block value is order-free), so
CR/H = c*(p) exactly, independent of the partition. 042's recorded
proof is correct as written.

**P3 — independent re-implementation** (`uc_reviewer043_reimpl.py`,
exit 0): 132 sampled 041 rows (both paths, n ∈ {2, 8, 64} ×
p ∈ {0.38271, 0.49}, every 4th grid point) match `data/hu_interp.json`
margins to 1e-9 through my own no-pruning DP; own DP vs own
posterior-pair atom evaluator ≤ 1e-10 at n ≤ 6; order invariance under
random permutations ≤ 1e-12; all 041 summaries (worst margin, argmin,
zero negatives, refined minima, 1968 rows) recompute from the stored
rows. 042 part A re-done with my own tensor builder: 5 cases including
the p = 0.20 clamp branch, under identity, my own rollout, a fixed
non-contiguous interleaving and a random order — |CR/H − c*(p)| ≤ 1e-12
everywhere; DIAG spot-checked at random orders, n ∈ {3,6},
p ∈ {0.05, 0.20, 0.35, 0.494}. Every endpoint row of
`data/hu_blocks.json` (B_attack, C_crash8, C_sat497) re-scored through
my own rollout + evaluator to 1e-8, all in-regime; C_sat497's floor
equals c*(0.494) to 1e-12 on the exact two-atom measure
{∅: 0.506, S: 0.494}.

**P4 — pipelines re-run unmodified, byte-compare:** `uc_hu_interp.py`
and `uc_hu_blocks.py` regenerate their checkpoints **byte-identically**
(md5 bfe37ca… and b705558… before and after); both skeptics exit 0.

**P5 — reporting audit:** 041's 1968-point count, the −5.3e-14 worst
margin at (n=64, p=0.45, t=0), all six interior minima at p = 0.49
(+2.20e-4/+2.1e-3/+9.1e-3 path A, +1.41e-7/+9.8e-7/+8.7e-6 path B, all
at t = 0.025, the first grid point) and their global-minimum locations
(n=2, p=0.49) check against the rows. 042's 4.4e-16 over 40 = 8×5
case-order pairs, the five B-table floors, crash8's +0.000952/+0.000310,
and the extremal constants c*(0.49) = +0.000289, c*(0.497) = +0.000026
all check. What does **not** check is 042's characterization of the
B-attack endpoints — the corrections below.

## Outcome

- **VERIFIED: 041 in full.** The count-pair-DP collapse is proved (P1),
  the committed implementation matches the derivation, sampled margins
  reproduce independently at all three n up to 64, the checkpoint
  summaries recompute, the pipeline re-runs byte-identically, and every
  number in the record's prose is right. The record's claims are stated
  within its evidence.
- **VERIFIED: 042's additivity proof and equality family** — the
  every-order block additivity and the c*(p) corollary are re-derived
  here in full (P2), and the machine check reproduces through an
  independent stack including non-contiguous interleavings and the
  clamp branch. 042's part A "VERIFIED" status is properly earned.
- **VERIFIED with corrections: 042's parts B/C numbers** — every floor
  reproduces to 1e-8, crash8 and the 0.497 saturation close-outs are
  right, zero violations confirmed (and strengthened: I re-checked every
  endpoint against its OWN constant c*(max marginal), all ≥ 0 to float
  noise — the committed scripts' `floor < 0` violation flag alone is
  weaker than (HU-TAX)). But:
- **Correction 1 (042, part B prose + outcome bullet):** "the polished
  endpoints are single diagonals at p closer to the cap — e.g.
  +0.000364 ≈ c*(0.4887)" and "every endpoint ON the family (single
  diagonals nearer the cap)" are **false for exactly the two campaigns
  where the polish improved on the anneal** (cap 0.49 n=6 d2⊗d2⊗d2,
  floor +0.000364; cap 0.497 n=6 d2⊗d2⊗d2, floor +0.0000610). Those
  endpoints are 8-atom measures on the d2⊗d2⊗d2 support with all six
  marginals equal (0.489001 and 0.495947) but the blocks pairwise
  **correlated** (P(two blocks both on) = 0.238143 vs 0.239122 for the
  product at cap 0.49): they are not single diagonals, not block
  tensors, and not on the equality family at all — their CR/H exceeds
  c*(own marginal) by +1.47e-5 and +1.36e-5. The +0.000364 ≈ c*(0.4887)
  numerology also fails on its own terms: c*(0.4887) = 3.686e-4, not
  3.639e-4 (the q solving c*(q) = floor is 0.48877, but the endpoint is
  not a diagonal at any q). The three unpolished endpoints ARE exact
  family members (d3(0.487)⊗d3(0.487), d2(0.494)⊗d2(0.494), and
  diag_2(0.487) padded by two marginal-0 coordinates — the last also
  not literally "at p closer to the cap", and the padding coordinates
  mean the equality set extends beyond the common-p family as stated).
- **Clarification (042, not an error in the numbers):** the corrected
  picture is *weaker but still positive*: the polish stalls strictly
  ABOVE the family (own-margin ~+1.4e-5) rather than converging onto
  it, and the cap-0.49 B floor (+0.000364) is therefore a descent stall
  point, not a family value — it lies below the projection-boundary
  family value c*(0.487) = +0.000488 only because the polish, unlike
  the anneal, is not projection-bounded and can raise the marginal
  toward the cap. "Every attacked cap's floor sits on the diagonal
  family" survives only in the anneal-with-projection sense used by
  040/C_sat497, not for the part-B polished floors.
- **Not refuted:** zero violations anywhere (confirmed independently,
  and against own constants); both records' EVIDENCE statuses; 041's
  status and every claim in it.
- **Not claimed:** anything the records did not claim (n > 64 / n > 8,
  non-exchangeable paths, whole-equality-set characterization); that
  the anneal/polish search itself was re-run seed-for-seed beyond the
  byte-identical pipeline re-execution; own-constant checks for
  intermediate (non-endpoint) measures visited by the anneal, which the
  committed instrumentation does not record.

## Why it failed / what survived

The refutation attempt failed on everything load-bearing: both hand
re-derivations came out clean, the DP and the additivity proof are
correct as recorded, and every committed number reproduces through an
independent implementation. What it caught is narrative drift at the
exact spot 042 leaned past its data: the adversary does NOT "find no
direction off the family" — in the two polished campaigns it ends OFF
the family, slightly above it, and the record dressed those stall
points as family members with a spurious c*(0.4887) match. The
corrected reading is friendlier to (HU-TAX) in one way (even off-family
stall points respect the own-constant bound) and less tidy in another
(the "adversary converges onto the family" story now has two counter-
instances; whether the polish would reach the family with more rounds,
or is stuck at a genuine local structure — correlated blocks at equal
marginals — is open and mildly interesting). Reusable: the reviewer
stack now contains a fifth independent HU/CR implementation
(posterior-pair state) and an own count-pair DP; the off-family stall
measures themselves are candidate seeds for a sharper attack.

## Leads generated

1. **The correlated-block stall measures** (8 atoms, equal marginals,
   blocks pairwise correlated, own-margin +1.4e-5) are the closest
   off-family points of record. Descending FROM them with a finer step
   floor, or characterizing the d2⊗d2⊗d2-support slice analytically
   (it is a 3-parameter exchangeable-block family — small enough for
   exact work), would either close the +1.4e-5 gap onto the family or
   exhibit a genuine off-family local minimum of CR/H.
2. **Instrument own-constant margins in future attacks:** the
   `floor < 0` violation flag should be
   `CR/H < c*(max marginal) − tol` per candidate, or violations
   strictly between 0 and c*(p̄) go unflagged (none occurred at any
   endpoint here, but intermediates are unchecked).
3. The equality family should be stated with padding: marginal-0
   coordinates (and by symmetry any deterministic coordinate) preserve
   equality, so the set is block-tensors at common p on a sub-cube,
   tensored with point masses — relevant to 042 lead 1's "iff".

## References

- This repo: 041, 042 (under review), 040 (DIAG identity, anneal
  machinery), 039/036/034/029/024 (reviewer-pass precedents and the
  independence-bar discussion), 037 (rollout rule), 031/030 (HU
  coupling), 009 (CR). `explore/uc_reviewer043_reimpl.py`,
  `data/hu_interp.json`, `data/hu_blocks.json`.
- No external sources.
