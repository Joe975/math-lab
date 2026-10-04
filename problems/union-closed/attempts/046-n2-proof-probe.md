# 046 — n = 2: the order quantifier is vacuous, the equality set is exact, and the pair-interaction term is load-bearing

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** proof probe + certified computation (045 lead 3: the first
  order-quantified positive statement on the HU line, attempted at the
  smallest nontrivial size).
- **Tools:** `explore/uc_hu_n2.py` (new; closed form, structural
  lemmas, exact-rational branch-and-bound certification, local
  structure, the reduction test; deterministic; checkpoint
  `data/hu_n2.json`); `explore/uc_hu_n2_skeptic.py` (new; own
  history-recursion evaluator in nats, no import of the engine or the
  037 stack; re-derives every structural claim and re-samples the
  residue boxes; exit 0). Reproduce: run the two in that order.
- **Sources:** none.

## Approach

At n = 2 the half-union coupling is a three-parameter object, so the
conjecture promoted in 045,

    CR_HU(μ) ≥ c*(q)·H(μ),   q = max marginal,
    c*(p) = (h(z*(p)) − h(p))/h(p),  z*(p) = max(1/2, 1−2p),

becomes a concrete inequality on a box. Parametrize by x = P(A₀=0) and
u_a = P(A₁=0 | A₀=a); the marginals are f₀ = 1−x and
f₁ = x(1−u₀) + (1−x)(1−u₁), and "in regime" means
q = max(f₀,f₁) ∈ (0,1/2).

Why n = 2 rather than the 2-block genre of 041: the block-tensor
genre is already covered by 042's additivity, so it cannot produce the
missing averaging step. n = 2 with *arbitrary* dependence is the
smallest case where 031's deficit cells actually appear, and 033/035
both pointed at it as "the minimal concrete instance of the averaging
problem".

## What was done

**A. The closed form, and two structural lemmas.**

With z₀ = clip(x,x) and cell weights (z₀, x−z₀, x−z₀, 1−2x+z₀),

    CR = h(z₀) + Σ_{a,b} w_ab h(clip(u_a,u_b)) − H(μ),
    H(μ) = h(x) + x·h(u₀) + (1−x)·h(u₁),

which reproduces the general evaluator to **2.4e-15** over 19,886
in-regime samples.

- **Lemma N2-ONE-BAD (proved).** At most one conditional can be out of
  regime: if u₀ < 1/2 and u₁ < 1/2 then
  f₁ = x(1−u₀)+(1−x)(1−u₁) > 1/2, so μ is out of regime. Hence the
  clip has only two cases, and every cell has an explicit closed form
  — in particular a deficit cell (both arguments below 1/2) can only be
  the (1,1) cell, where s(u₁,u₁) = 0 exactly, and the off-diagonal
  cells contribute (h(u₁) − h(u₀))/2.
- **The order quantifier is vacuous at n = 2 (EVIDENCE).** The
  identity order is strictly worse than the best order on 8,207 of
  19,886 samples, yet its margin never goes below 0 — worst
  **+3.4e-05** over the samples and **−6.7e-16** (i.e. 0) over a
  1.38M-point grid, attained at products. Since the two orders are
  related by relabelling and the regime is relabel-closed, positivity
  for the identity order on all μ is equivalent to positivity for
  *every* order. So at n = 2 the conjecture needs no max over orders
  — the quantifier that 033/035/044 showed is essential at n ≥ 4 is
  empty at n = 2.

**B. The equality set, exactly.** Three families, each an identity
(derived by hand in "Why it failed", re-checked here in exact rational
interval arithmetic under **both** kits alone, 24 points, every
enclosure of F′ containing 0):

    products Bern(p)^⊗2   (x = u₀ = u₁),
    diagonals {∅:1−p, {1,2}:p}   (u₀ = 1, u₁ = 0),
    the n = 1 degenerations (one coordinate constant).

Everything else in a 1.38M-point grid is strictly positive. This is
the n = 2 shadow of 042's block-tensor family, and it confirms that
family is the whole equality set here — no fourth family.

**C. Certified branch-and-bound.** F′ := CR·h(q) − (h(z*(q))−h(q))·H
(the division-free form of the claim) enclosed in exact rational
interval arithmetic on sub-boxes, subdividing the longest edge:

    budget 2,000,000 boxes, minimum edge 1/512, queue largest-box-first
    root box volume            0.500000
      certified F' > 0         0.235977   (88.4% of the in-regime volume)
      wholly out of regime     0.233095
      residue at minimum edge  0.000000   (no box resisted to the floor)
      unprocessed (budget)     0.030930   (11.6% of in-regime)

  **Zero residue** is the informative part: no box anywhere in the
  regime resisted refinement down to edge 1/512 — the certification
  never hit an obstruction, it only ran out of budget, and what is
  left is the neighbourhood of the equality set and the q = 1/2 face
  where boxes must be split finest. The claim this supports is
  therefore "F′ ≥ 0 is certified on 88.4% of the in-regime volume in
  exact rational arithmetic, with no counterexample region", not a
  theorem.

  A numeric certificate can never close a neighbourhood of a point
  where the inequality is an equality, so the equality set (B) and the
  q → 1/2 face (G) are exactly where the unprocessed volume sits. The
  skeptic's own 60³ grid over the whole regime box (independent
  evaluator, no shared code) finds no negative margin anywhere — worst
  −6.9e-16, i.e. zero, at a product — which is the sampling
  counterpart to the certified fraction. (Its residue-sampling check
  is vacuous here: there are no residue boxes.)

**D. Local structure at the equality families.** Transversal probes at
ε = 10⁻² and 10⁻³ give margin growth of order ε¹ in every transversal
direction tested (products: +2.4e-2 → +2.4e-3; diagonals: +2.5e-2 →
+3.7e-3, order ≈ 0.8–1.3). So the equality set is a *first-order*
ridge, not a degenerate valley — the inequality is not tight to second
order anywhere, which is why the residue volume shrinks so fast under
refinement, and which tells a future proof it only needs a first-order
argument transversally.

**E. The proof route that fails, and exactly why.** Two steps are
clean:

- since c* is decreasing and q ≥ f₀, the first coordinate's surplus
  can be discarded: F ≥ Σ_ab w_ab·s(u_a,u_b) − c*(q)·[x h(p₀) +
  (1−x) h(p₁)], with s(a,b) = h(clip(a,b)) − (h(a)+h(b))/2 the cell
  ledger of 037 and p_a = 1−u_a;
- **Lemma N2-CONC (proved).** ψ(t) := h(min(1/2,t)) is concave
  (h is concave increasing on [0,1/2], min(1/2,·) is concave
  nondecreasing), so in the both-good case
  s(p_a,p_b) ≥ (G(p_a)+G(p_b))/2 where G(p) = c*(p)h(p) is the
  diagonal-cell surplus, with equality iff p_a = p_b. Averaging over
  the cells (whose two marginals are both (x, 1−x)) turns this into
  Σ w_ab s ≥ x G(p₀) + (1−x) G(p₁).

Together they reduce the theorem to the scalar inequality

    (**)  x G(p₀) + (1−x) G(p₁) ≥ c*(q)·[x h(p₀) + (1−x) h(p₁)].

**(**) is FALSE.** It fails on **85,380 of 400,000** Case-A points
(21.3%), worst deficit **−9.72e-02** at x = 0.7587, p₀ = 0.3180,
p₁ = 0.0010 — while the original margin at those very points is
**positive** (+9.3e-2 there, +7.9e-2…+1.2e-1 across the worst five).
So the reduction is lossy, and the exact decomposition says where:

    F = (c*(f₀) − c*(q))·h(f₀)
        + [x σ(p₀) + (1−x) σ(p₁)]
        + t·Δ
        − c*(q)·[x h(p₀) + (1−x) h(p₁)],
    t = min(x − 1/2, 1 − x),  Δ = 2s(u₀,u₁) − s(u₀,u₀) − s(u₁,u₁),
    σ(p) = G(p) for p ≤ 1/2 and 0 otherwise

— an identity verified to **2.6e-15** over 99,970 in-regime samples.
Δ ≥ 0 in the both-good case is exactly N2-CONC (0 negative samples,
min −1.1e-16), and the scalar part **without** t·Δ is negative on
**5.3%** of Case-A points. **The pair-interaction term t·Δ is
load-bearing**: no argument that averages the cell interaction away
can prove even the n = 2 case.

**F. A correction found in passing: the stated form of c\* is wrong
above p = 1/4.** Every engine in the line computes

    c*(p) = (h(z*(p)) − h(p))/h(p),   z*(p) = max(1/2, 1−2p)

(the clip evaluated at a product), and every recorded constant matches
it: +0.041739 at 0.38271, +0.007278 at 0.45, +0.000289 at 0.49,
+0.0000260 at 0.497. But 034's prose states the constant as
`c*(p̄) = (h(min(2p̄,1)) − h(p̄))/h(p̄)`, and 037 §Leads and STATUS.md
repeat it. Those two agree only for p ≤ 1/4; above it the stated form
is not merely different, it is **negative** — −0.181287 at p = 0.38271,
−0.527591 at 0.45, −0.858519 at 0.49 — so as written it would make
(HU-TAX) vacuous exactly on the range 034 re-posed it for. The
mis-statement is prose-only: no computation in 030–045 ever evaluated
it (all of them call the `max(1/2, 1−2p)` form), so no number moves.
Root cause: 034 wrote the p < 1/4 branch (`h(2p)`, where the Fréchet
floor binds) as if it were the general formula, with `min(·,1)`
guarding the wrong endpoint; the correct single expression is
`h(max(1/2, 1−2p))`, which equals h(2p) below 1/4 and 1 above it.

**G′. A method error caught by measuring it — depth-first B&B is not
coverage.** The first two certification runs used a plain DFS stack and
reported "residue volume 1.0e-2 against a root volume of 0.5", which
reads as near-complete coverage. Adding volume accounting for the
*unprocessed* stack showed the truth: DFS dives into one corner and
subdivides it to the minimum edge, so at a 2,000,000-box budget it had
certified **0.045%** of the box by volume and left **99.8%
unprocessed**. Any residue figure quoted without the unprocessed
volume beside it is meaningless. Fixed by making the queue
largest-box-first (heap on −volume), which grows coverage before
depth; the numbers in C are from that version. Lesson for every future
certified sweep in this repo: report certified / out-of-regime /
residue / **unprocessed** volumes together, or do not report a
fraction at all.

**G. A fourth, limiting equality structure — found by probing this
record's own lead.** Minimising the margin subject to staying a fixed
distance from the three interior families does NOT settle at an
interior point: it drifts to the regime boundary, to
μ = {{1}:1/2, {2}:1/2} — two disjoint singletons at marginals exactly
1/2 — where the coupling gives CR = 0 exactly while c*(1/2) = 0, so
the inequality is tight trivially (0 ≥ 0). Along the in-regime family
{{1}:1/2−ε, {2}:1/2−ε, ∅:2ε} the **ratio** margin
CR/H − c*(q) is strictly positive at every ε > 0 (+7.10e-2 at ε = 0.2,
+2.10e-3 at 0.01, +7.7e-11 at 1e-6) and vanishes only in the limit.
So this is an *asymptotic* equality at the q → 1/2 boundary, not a new
interior one — B's characterisation stands as stated for the open
regime.

  The methodological point is the useful part: the absolute margin
  F = CR − c*(q)H is degenerate as q → 1/2 because **both** sides
  vanish (c*(q) → 0), so any adversary descending on F drifts to that
  boundary and reports vanishing margins that are not near-violations.
  This is the exact mirror of 023's H(μ) → 0 corner lesson — there the
  entropy vanished, here the constant does — and it means a q ≤ cap
  guard, or the scale-free ratio form, belongs in every future
  (HU-TAX) adversary. It also accounts for part of C's residue: boxes
  hugging the q = 1/2 face cannot be certified in the F form no matter
  how finely they are split.

## Outcome

- **PROVED (hand, skeptic-re-derived): Lemma N2-ONE-BAD** — at n = 2
  at most one conditional is out of regime, so deficit cells are
  confined to the (1,1) cell and contribute exactly 0.
- **PROVED (hand, skeptic-re-derived): Lemma N2-CONC** — concavity of
  ψ(t) = h(min(1/2,t)) gives s(p_a,p_b) ≥ (G(p_a)+G(p_b))/2, tight iff
  p_a = p_b. Reusable at every n: it is the sharp cell-wise bound.
- **PROVED (hand): the exact decomposition identity** above, and the
  three equality identities.
- **NOT PROVED: the n = 2 case of (HU-TAX).** The natural reduction
  dies at (**), by a measured margin, and the obstruction is now
  named: the pair-interaction term.
- **EVIDENCE (certified, partial): F′ ≥ 0 on the n = 2 regime box** —
  exact-rational B&B certifies the overwhelming majority of the box by
  volume, with residue confined to a neighbourhood of the equality set
  and no negative point found in any residue box.
- **EVIDENCE: the order quantifier is vacuous at n = 2** (identity
  order always suffices) — in contrast to n ≥ 4, where 033/035/044
  show it is essential.
- **CORRECTION to 034 (reporting-level, propagated to 037 and
  STATUS.md):** the stated constant `(h(min(2p̄,1))−h(p̄))/h(p̄)` is
  wrong for every p̄ > 1/4 and negative there; the constant every
  engine and every recorded number uses is
  `(h(max(1/2,1−2p̄))−h(p̄))/h(p̄)`. Prose only — no computation used
  the stated form, so nothing measured moves. Per the repo rule 034
  is left as written.
- **Not claimed:** a proof at n = 2; anything about n = 3; that the
  B&B exhausts the box (it does not — the budget and the equality set
  both bound it); that (**) failing means (HU-TAX) fails (it does not
  — the original margin is positive exactly there).

## Why it failed / what survived

The reduction failed for a reason worth carrying: **the cell ledger's
pair interaction is not slack, it is the mechanism.** Concavity gives
the sharp cell-wise bound and it is *exactly* tight at products —
which is why products are an equality family — but the moment the two
conditional marginals separate, the true surplus
h(min(1/2,p₀+p₁)) − (h(p₀)+h(p₁))/2 runs far above the average of the
diagonal surpluses, and the theorem lives on that excess. In the
decomposition, that excess is t·Δ with t = min(x−1/2, 1−x): it
vanishes at x = 1/2 and at x = 1, and peaks at x = 3/4. So the
extremal structure is a competition between the first coordinate's
own surplus (large when f₀ is small) and the interaction weight t
(largest at f₀ = 1/4) — which is exactly the p = 1/4 threshold that
034 found governs the clamp. That coincidence is worth a look before
anything larger is attempted.

For the route: the n = 2 case being *unproved but certified* and
*order-quantifier-free* is the honest state. It says the difficulty in
(HU-TAX) is not the order quantifier (empty here, essential at n ≥ 4)
and not the deficit cells (they contribute 0 here) but the averaging
of the interaction — the same thing 031 §"Why it failed" named, now
isolated to a two-variable function on a three-parameter box.

## Leads generated

1. **Prove n = 2 by keeping t·Δ.** The target inequality is now
   explicit and scalar: with p₀ ≥ p₁ and t = min(x−1/2, 1−x),
   show (c*(f₀)−c*(q))h(f₀) + xσ(p₀) + (1−x)σ(p₁) + tΔ ≥
   c*(q)(x h(p₀) + (1−x)h(p₁)). Δ has the closed form
   2h(min(1/2,p₀+p₁)) − h(min(1/2,2p₀)) − h(min(1/2,2p₁)) in Case A.
   Falsifiable and small.
2. **Where is the interior binding configuration?** The guess in
   "Why it failed" — f₀ = 1/4, where the interaction weight t peaks
   and the clamp switches branch — was NOT confirmed: a
   distance-constrained descent on the absolute margin leaves the
   interior entirely and runs to the q → 1/2 boundary (part G). Redo
   it on the **ratio** margin CR/H − c*(q), which is scale-free, and
   with a q ≤ 0.49 guard; only then does "which interior μ binds"
   have a well-posed answer.
3. **Complete the certification** with a tighter enclosure (centered
   forms, or the case split of N2-ONE-BAD applied before enclosing) so
   the B&B exhausts the box outside explicit ε-balls at the equality
   set — that would upgrade this record's EVIDENCE to a certified
   theorem-modulo-the-equality-neighbourhoods.
4. **n = 3 next, not the block genre**: 042 already owns block
   tensors; the value of n = 2 was the interaction term, and n = 3 is
   where a *second* interaction layer appears.

## References

- This repo: 045 (the promoted conjecture and its lead 3), 044 (the
  own-constant standard used for the margin), 042/041 (the equality
  families this recovers at n = 2), 037 (the cell-ledger surplus
  form s), 034 (the corrected constant c* and the p = 1/4 threshold),
  031 (the averaging problem this isolates), 029 (single-kit
  certification), 026/022/016 (kits). `data/hu_n2.json`.
- No external sources.
