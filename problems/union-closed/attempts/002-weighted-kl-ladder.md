# 002 — Weighted-KL ladder (Idea B): evaluated and closed

- **Problem:** union-closed sets conjecture (Frankl), `problems/union-closed.md`
- **Date:** 2026-07-25
- **Type:** rigorous evaluation of candidate idea B from
  `001-entropy-barrier-map.md` (formalization + adversarial test). Outcome:
  **dead end, definitively** — the refutation is at the level of *truth* of
  the needed statement, not of proof technique.
- **Tools:** `tools/uc_weighted_kl.py` (written and run here, stdlib-only,
  deterministic, < 1 s); primary sources re-read in full for the two decisive
  facts: Sawin arXiv:2211.11504 (Proposition 6 and its proof, transcribed
  from ar5iv), Ellis arXiv:2211.12401.

Notation as in 001: `h` = binary entropy (bits), `d(a‖b)` = Bernoulli KL,
`ψ = (3−√5)/2 ≈ 0.381966`, `φ = 1−ψ`. `A, B` i.i.d. from a distribution `μ`
on `2^[n]`, `U = A∪B`. `D(U‖μ)` means KL of the law of `U` against `μ`
(in the family application `μ = Unif(F)`, so this is 001's `D(U‖A)`).
Record to beat: Liu's ≈ 0.38271 (arXiv:2306.08824).

## Approach

Idea B proposed the family-level inequality `H(U) + c·D(U‖Unif(F)) ≤ log|F|`
as a "strengthened use of exact closure", with a computed product-obstruction
ladder `p(c)` reaching 0.3926 already at c = 0.1. Plan: (1) formalize and
locate where the actual mathematical content sits; (2) re-derive and verify
the `p(c)` curve independently; (3) adversarially test the needed statement
against Sawin's counterexample to Gilmer's Conjecture 1 (and relatives),
computing its full c-profile; (4) record the outcome precisely enough that
the idea is never naively re-attempted.

## What was done

### 1. Formalization: the identity empties the family level

For **any** random variable `X` supported inside `F`:

    H(X) + D(X‖Unif(F)) = log|F|        (one line: D(X‖Unif) = Σ P(x)·log(P(x)|F|)).

Hence for exact union-closed `F` (which forces `supp(U) ⊆ F`):

    H(U) + c·D(U‖Unif(F)) = log|F| − (1−c)·D(U‖Unif(F)).

- **c ≤ 1:** the inequality `≤ log|F|` is *identically true* — it is
  equivalent to `(1−c)D ≥ 0`. It is not a new constraint on `F`: it re-uses
  the single fact `supp(U) ⊆ F` that Gilmer's `H(U) ≤ log|F|` already used,
  and nothing else. No additional exactness is injected by the KL term.
- **c > 1:** the inequality is *false for every union-closed `F` with
  `|F| ≥ 2`*. (Lemma, proved here: `U` is never uniform on such `F` — for a
  minimal-under-inclusion `M ∈ F`, `U = M` forces `A = B = M`, so
  `Pr[U=M] = 1/|F|² < 1/|F|`; hence `D > 0` and `H(U)+cD = log|F| + (c−1)D >
  log|F|`.) So the admissible weights are exactly `c ∈ [0,1]`.

**Where the difficulty moves.** Since the family side is a tautology, all
content must come from the matching *distributional* lower-bound statement
(the c-analogue of what AHS/CL/Sawin actually prove for c = 0):

    (S_c) at threshold p:  for every μ on 2^[n] with H(A) > 0 and all
    marginals Pr[i∈A] < p:   H(U) + c·D(U‖μ) > H(A).

(S_c) at p, applied to `μ = Unif(F)`, contradicts the identity above and
yields Frankl with constant p. The ladder is thus an interpolation between
two known endpoints: **c = 0 is the AHS/CL/Sawin distributional theorem
(true, sharp at ψ)** and **c = 1 is precisely Gilmer's Conjecture 1
(refuted by Sawin and Ellis)**. Note the equality structure at c = 1:
uniform on any union-closed `F` satisfies `H(U)+D = log|F| = H(A)` *exactly*,
so Conjecture 1 (strict) would have implied the full conjecture. Strictness
is never the loophole: we also checked (small lemma) that no non-degenerate
"union-stationary" μ exists (`A∪B ~ A` for i.i.d. A,B forces a point mass,
via `|A∪B| ≥ |A|` with equality iff `B ⊆ A`, plus symmetry).

**Bookkeeping that makes the evaluation mechanical.** For fixed μ,

    f_μ(c) := H(U) + c·D(U‖μ) − H(A)   is LINEAR in c with slope D ≥ 0.

So each violating μ is summarized by one number, `c*(μ) = (H(A)−H(U))/D`
(defined when `H(U) < H(A)`, `0 < D < ∞`): μ refutes (S_c) exactly for
`c ≤ c*(μ)`, at every threshold above its max marginal. Consequences:
counterexamples propagate *downward* in c; the true ceiling
`P(c) := sup{p : (S_c) holds at p}` is nondecreasing in c; `ψ ≤ P(c) ≤ p(c)`
(lower bound because `f_μ(c) ≥ f_μ(0)` and the c = 0 theorem is true at ψ;
upper bound = product obstruction below).

**Why the chain-rule machinery was never going to tensorize** (recorded for
completeness; moot given §3): the KL chain rule gives
`D(U‖μ) = Σᵢ E_{u_{<i}~U} d(P(Uᵢ|u_{<i}) ‖ μ(Aᵢ | A_{<i}=u_{<i}))` — the
reference conditionals are `A`-conditionals evaluated at *U-histories*,
events the `A`-process reaches with different (possibly zero) probability;
they are a third, F-dependent quantity unrelated to the `x, y` of the key
one-bit inequality, and the "drop conditioning" step (concavity of entropy)
has no KL analogue. The c-weighted one-bit inequality
`h(xy) + c·d(xy‖z) ≥ λ_c(x·h(y) + y·h(x))` is true on the diagonal for
`x = y = z = 1−p`, `p ≤ p(c)` — but it cannot be assembled into (S_c),
because (S_c) is *false* (§3), not merely hard.

### 2. The threshold curve p(c), verified

`p(c)` = first root of `h(2p−p²) + c·d(2p−p²‖p) = h(p)` — equivalently the
level where the Bernoulli-product family `Bern(p)^⊗n` stops violating (S_c);
computed by bisection on the closed form `c*_prod(p) = (h(p)−h(q))/d(q‖p)`,
`q = 2p−p²` (`tools/uc_weighted_kl.py`, part A). All values in 001's table
reproduce to the printed precision:

    c    : 0        0.05     0.1      0.2      0.3      0.5      0.7      0.9      1.0
    p(c) : 0.381966 0.387275 0.392644 0.403565 0.414733 0.437818 0.461919 0.487048 0.500000

001's claims verified: c = 0 → 0.38197 (= ψ), c = 0.1 → 0.39264. The weight
that would beat Liu's record is tiny: `c_needed = c*_prod(0.38271) ≈ 0.00704`.
But this curve is a ceiling computed *from products only*; §3 shows products
are not the extremal distributions for any c > 0.

### 3. The decisive test: known counterexamples and their c-profiles

**Chase–Lovett's approximate family refutes no c > 0.** For `μ = Unif(F_CL)`
(slice + top), `U` escapes `F_CL` with positive probability, so
`D(U‖μ) = +∞` and (S_c) holds vacuously for every c > 0 on it. The KL term
genuinely detects inexact closure — this is the (correct) intuition that
motivated Idea B. It fails anyway, because the detection is only
*logarithmic in planted mass*, hence removable by smoothing:

**Sawin's Proposition 6 kills every c at every p > ψ.** Construction
(arXiv:2211.11504, transcribed and independently verified): fix `ū`, sample
`k ~ Geometric(θ)`, then `A | k ~ Bernoulli(1−(1−ū)^{k+1})^{⊗n}`. The
component class is *closed under the union-convolution*:
`(1−p_a)(1−p_b) = (1−ū)^{a+b+2}`, so `U` is a mixture of the *same* products
with shifted index `k' = k_A + k_B + 1`. Sawin proves: marginals ≤ ū + O(θ);
`H(A) = Θ(n)`; `H(U) ≤ (h(2ū−ū²)/h(ū) + o_θ(1))·H(A)`; and — the killer —
`D(U‖A) = O(1)` *uniformly in n*, because each component `S_{k'}` of `U` is
`A` conditioned on an event of probability `(1−θ)θ^{k'}`, so
`D(S_{k'}‖A) ≤ k'·log(1/θ) + O(1)`, and convexity of KL finishes.

Consequence (the whole point): for any `p > ψ` and ANY fixed `c ≥ 0`, take
`ū ∈ (ψ, p)`; then `h(2ū−ū²) < h(ū)` (since `(1−ū)² < ū` iff `ū > ψ`), so

    f_μ(c) ≤ (h(2ū−ū²)/h(ū) − 1)·Θ(n) + c·O(1) → −∞.

`c*(μ_n) → ∞`. Every rung of the ladder — indeed every c, not just c ≤ 1 —
is refuted at every threshold above ψ. Combined with the lower bound from
§1: **P(c) = ψ exactly, for every c ≥ 0. The ladder is flat.**

**Finite-n certificates** (exact evaluations, not asymptotics — both laws are
exchangeable mixtures of products, so all quantities reduce to weight
profiles; `tools/uc_weighted_kl.py` part C):

    n      ū       θ      max marg   H(A)      H(U)      D       c*(μ)
    2000   0.390   0.05   0.402269   1927.70   1867.07    4.350   13.9
    20000  0.386   0.02   0.390799  19239.59  18974.63    5.655   46.9
    60000  0.3823  0.001  0.382536  57578.99  57527.34    9.966    5.2

The third row has **every marginal ≤ 0.382536 < 0.38271** and violates (S_c)
for all `c ≤ 5.2` (f(0.007) = −51.6, f(0.1) = −50.7, f(1) = −41.7 bits). So
the ladder cannot even recover the current record, let alone beat it; its
exact ceiling is ψ = 0.381966.

**Gap analysis (the question 001 posed).** 001 asked for "the largest c₀
surviving known counterexamples". Answer: the apparent gap was an artifact
of looking at small counterexamples. Ellis's n = 2 example
(`μ = {∅:0.3, {1}:0.2, {2}:0.2, {12}:0.3}`; part B) has `c* = 1.148` at
marginals exactly 1/2, and its perturbations lose their bite quickly
(marginals 0.48 → c* = 0.965; 0.45 → 0.678; 0.40 → 0.169) — this family
alone would have left every `c ≲ 0.17` alive at p ≤ 0.40, seemingly leaving
room for the `c ≈ 0.007` needed. Sawin's large-n family closes the gap
*completely*: no c survives at any p > ψ.

**Independent cross-check** (to rule out an artifact of Sawin's specific
mixture; part D): a "regularized Chase–Lovett" distribution built here —
`μ = (1−t)·Unif(slice at w₀ ≈ p·n) + t·Bern(q)^{⊗n}`, `q = 1−(1−p)²`,
`t = 2^{−εn}` — has full support (D finite) and satisfies the *certified*
closed-form bounds `H(U) ≤ n·h(q)`, `CrossEnt(U,μ) ≤ εn + n·h(q)` (μ
dominates `t·Bern(q)^n` pointwise; cross-entropy against a product depends
only on U's marginals), `H(μ) ≥ (1−t)log₂C(n,w₀) + t·n·h(q)`. This gives
`f(c) ≤ F(c) = n·h(q) + cεn − H_lb` for all c ∈ [0,1]; computed:
F(1) = −2543 at (n = 2·10⁵, p = 0.3927), F(1) = −6253 at
(n = 2·10⁷, p = 0.3822 < record). Same phenomenon, different construction.

## Outcome

**Dead end — closed with a proof-level obstruction, not a technique-level
one.** The precise final statement:

> For every `c ≥ 0`, the distributional statement (S_c) holds at threshold
> ψ (implied by the c = 0 theorem, since `f_μ` is increasing in c) and fails
> at every threshold `p > ψ` (Sawin's Proposition 6 distribution with
> `ū ∈ (ψ, p)`: entropy deficit Θ(n), KL cost O(1)). Hence the weighted-KL
> ladder's true value is `P(c) ≡ ψ`: no admissible weight `c ∈ [0,1]` yields
> any constant beyond ψ = 0.381966. The `p(c)` table of 001 is the ceiling
> for Bernoulli products only; products are not extremal for any c > 0.

## Why it failed / what survived

- **The two-part obstruction.** (i) *Family level:* the identity
  `H(X) + D(X‖Unif F) = log|F|` makes `H(U)+cD ≤ log|F|` (c ≤ 1) a
  tautology — the KL term adds zero new information about F; the only
  closure fact used is still `supp(U) ⊆ F`. (ii) *Distributional level:*
  KL against μ charges an "escape" region only by its log-likelihood under
  μ, not by its probability under the pair distribution. A union-target
  distribution planted as a mixture component of mass δ costs at most
  `log(1/δ)` bits of divergence while enabling a Θ(n) entropy drop; taking
  δ constant-in-n but small (Sawin: geometric weights θ^k) caps D at O(1).
  Any functional bounded by a fixed multiple of `D(U‖μ)` is defeated the
  same way.
- **Generalized no-go (for future cycles — do not re-attempt variants
  naively):** any strengthening of Gilmer's interface of the form
  "`Φ(law(U), μ) ≤ log|F|` for exact F, plus a marginal-threshold lower
  bound on Φ for all μ", where Φ depends only on the two laws and increases
  by at most `O(log(1/δ))` when μ is smoothed by a δ-mass full-support
  component, has ceiling exactly ψ. This covers: c·D with any c, D against
  smoothed references, Rényi-divergence weightings of bounded order, and
  any `Φ ≤ H + O(D)`. To pass ψ, a correction must charge closure
  violations by their *probability* (worst-case/pairwise closure), or use
  information that is not a functional of `(law(U), μ)` at all (couplings,
  conditioning on F-dependent events — Ideas A/C/D).
- **Survived:**
  - The `p(c)` curve, verified and correctly re-labeled as a product-only
    ceiling; `c_needed ≈ 0.00704` documents how tantalizingly cheap the
    (false) statement would have been.
  - The `c*(μ) = (H(A)−H(U))/D` linearity calculus: a one-number test that
    instantly profiles any proposed strengthened functional against any
    candidate counterexample. Implemented and reusable.
  - Sawin's component family (geometric mixtures of Bernoulli products,
    closed under union-convolution) as a reusable adversarial gadget, now
    with an exact finite-n evaluator (exchangeability → weight profiles).
  - The two small lemmas: no nontrivial union-stationary distribution;
    `U` never uniform on `F` when `|F| ≥ 2` (so c > 1 fails family-side
    everywhere).
  - `tools/uc_weighted_kl.py` (all tables above re-runnable in < 1 s;
    `--big` for larger n).

## Leads generated

1. **Adversarial test-first protocol** for any future strengthened
   distributional inequality on this problem: before any proof effort, run
   its c*-style profile against (i) Bernoulli products, (ii) Sawin geometric
   mixtures, (iii) the smoothed-slice family of part D. All three are in
   `uc_weighted_kl.py`. A proposal that survives all three is already
   interesting.
2. **Conditioning beats reweighting** (sharpened redirect toward Idea C/D):
   for exact F and uniform A, closure gives `H(A∪B | (A,B) ∈ E) ≤ log|F|`
   for *every* event E in the pair space — constraints that charge escapes
   by probability. They detect the CL family's inexactness directly: E = the
   high-overlap pairs makes U sweep the deleted middle band, whose log-size
   exceeds log|F_CL|. And they are not defeated by the smoothing gadget,
   because they are genuinely non-distributional: the right-hand side
   log|F| equals H(A) only because μ is *uniform on its support* — exactly
   the quantifier every (S_c)-type statement relaxed away, and the reason
   the no-go above does not apply.
   SPECULATION as a route, but it is exactly the licensing structure behind
   Sawin's and Liu's +0.0007, and it is untouched by this cycle's no-go.
3. **The no-go's boundary is the coupling class.** Liu's conditionally-iid
   scheme changes the sampler, not the functional; the analogue of c* for
   coupling families (which coupling classes have product/mixture
   obstructions above 0.383?) remains the open quantitative frontier — 001's
   Idea C first target stands, now with better test tooling.
4. Bookkeeping correction to 001: its lead (1) ("determine the c-profile of
   Sawin's counterexample; if some c₀ > 0 survives...") is resolved
   negatively and should not be re-queued: the profile is total
   (c*(μ_n) → ∞), and the same paper's Proposition 6 — not just its
   Conjecture-1 refutation — was already, in hindsight, a proof that Idea B
   as stated cannot work.

## References

- Sawin, arXiv:2211.11504 — Proposition 6 (geometric mixture of Bernoulli
  products; `D(A∪B‖A) = O(1)` with linear entropy deficit); also Examples
  4–5 (sharpness of the c = 0 scheme).
- Ellis, arXiv:2211.12401 — n = 2 counterexample to Gilmer's Conjecture 1.
- Gilmer, arXiv:2211.09055 — Conjecture 1 (= the c = 1 endpoint).
- Alweiss–Huang–Sellke arXiv:2211.11731; Chase–Lovett arXiv:2211.11689 —
  the true c = 0 theorem at ψ; CL Example 1.4 (the family whose smoothing
  underlies part D).
- Liu, arXiv:2306.08824 — the 0.38271 record the ladder needed to beat.
