# 025 — Gap 2's non-Sinkhorn branches attacked: both fail as stated, and one lemma repairs both

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-19
- **Mode:** informed
- **Type:** adversarial attack (023 lead 2 / queue 1b pre-step) on the 009
  recipe's half-mixing and block-adaptive branches — the two branches no
  CR adversary had exercised — plus the CR/H floor measurement that queue
  1(b) asks for before proof effort. Run inside a human-opened persistence
  window; verification bar unchanged.
- **Tools:** `explore/uc_branch_attack.py` (new; stdlib only; **no RNG
  anywhere** — grid scans and analytic constructions, so fully
  deterministic with no seed to record; ~30 s; checkpoint
  `data/branch_attack.json`); `explore/uc_branch_skeptic.py` (new;
  independent re-implementation, natural-log arithmetic, own coupling
  builders, shares no code with the engine or `uc_mitax.py`; output
  `data/branch_skeptic_out.txt`, exit 0 = no refutation). Reproduce:
  `python3 problems/union-closed/explore/uc_branch_attack.py` then
  `python3 problems/union-closed/explore/uc_branch_skeptic.py`.
  Python 3.11; no swarm, no external workers.
- **Sources:** none.

Notation as in 009/011: `CR(μ,π) = Σᵢ E_π[h(z_i)] − H(μ)`,
`ST = Gain − CR`, in-regime = all element marginals < 0.38271,
`ψ = (3−√5)/2`. The 009 recipe (Gap 3, as currently stated): λ-swept
Sinkhorn tilt on generic/slice μ; block-adaptive shared-latent coupling on
product mixtures; half-mixing on the `{0}∪[1/2,1)` genre.

## Approach

023 attacked only the Sinkhorn branch and left lead 2 explicit: exercise
the other two branches before trusting Gap 2's survival. Both branches
have closed-form structure, so instead of an anneal (023's tool) this
attack is analytic-first: locate the failure mechanism in the closed form,
construct explicit in-regime instances, and only then confirm at
full-history precision. That choice found failures an anneal seeded at
generic geometries would likely have missed (the sliver has ph ≈ 1−10⁻³
and grid-scan density ~0.1% at n = 8) — and it produces the repair, not
just the kill. Cost of the avoided alternative: an anneal over the branch
parameter spaces would have cost more compute and given only float
endpoints with no mechanism.

## What was done

### 1. The extended mixing lemma (proved here; skeptic-passed at
implementation level, reviewer pass pending)

**Lemma EM.** Let ν be *any* law on 2^[n], (ε_A, ε_B) *any* joint law on
{0,1}² with equal marginals P(ε=0) = P₁, independent of s ~ ν, and
q = P(ε_A = ε_B = 0). For the coupling A = ε_A·s, B = ε_B·s (both
marginals μ = P₁δ∅ ⊕ (1−P₁)ν):

    ST = 0 identically, hence CR = Gain = F(q) − F(P₁),
    where F(w) = H(w·δ∅ + (1−w)·ν).

*Proof.* If U_{<i} = u ≠ 0 then ε_A∨ε_B = 1 and s_{<i} = u are forced;
each (a,b) cell refining u fixes (ε_A, ε_B) but, since s ⊥ (ε_A, ε_B),
the conditional law of U_i = s_i in every such cell is law(s_i | s_{<i}=u)
— the same as conditioning on u alone. If U_{<i} = 0, the event
{A_{<i} = B_{<i} = 0} *is* the event {U_{<i} = 0}: the same σ-field atom.
So σ(A_{<i},B_{<i}) refines σ(U_{<i}) without changing the conditional
law of U_i, and every ST summand vanishes. The U-law is qδ∅ ⊕ (1−q)ν. ∎

009's HM lemma (011-verified) is the special case ν = Bern(ph)^⊗n,
q = P₁/2. Numerics: |ST| ≤ 3.2e-14 over a 12-point (n, ph, P₁, q) battery
(engine, A1) and ≤ 5.7e-15 on independently built couplings including
mixture s-laws (skeptic, S2).

**Corollary (exact optimality frame).** F is strictly concave along the
δ∅–ν segment, q ranges over [max(0, 2P₁−1), P₁], and so

    sup_q CR > 0  ⟺  F′(P₁) < 0  ⟺  −log₂ μ(∅) < E_ν[−log₂ μ] :

the branch helps iff ∅'s surprisal under μ is below the ν-average
surprisal. Verified as an iff against 401-point q-sweeps at 1096 in-regime
grid points, 1096/1096 agreement (A4).

### 2. Half-mixing branch (q = P₁/2) REFUTED as stated, at every n

The stated branch's CR is F(P₁/2) − F(P₁): a fixed cost
h(P₁) − h(P₁/2) > 0 against a shared-entropy gain ≈ (P₁/2)·n·h(ph). As
ph → 1 the gain vanishes while the cost stays — and the marginal cap
(1−P₁)ph < 0.38271 is *satisfiable* there (P₁ ∈ (0.61729, 2/3]). Explicit
in-regime failure family at P₁ = 0.64 (engine A2; every line has max
marginal ≤ 0.360, H(μ) ≈ 0.97):

    n=1 : ph=0.991334  CR_hm = −0.02939
    n=2 : ph=0.996214  CR_hm = −0.01539
    n=6 : ph=0.998943  CR_hm = −0.01532
    n=16: ph=0.999653  CR_hm = −0.01532   (…same to 5 dp for n ≥ 3:
    n=100: ph=0.999955 CR_hm = −0.01532    the family converges in the
                                           n·h(ph) parametrization)

So **for every n there are in-regime genre measures where the stated
branch has CR < 0** — (TAX at p) as instantiated by the 009 recipe is
false. This is a branch-assignment failure, not a Gap-2-core failure: see
§4. Failure-region size (A3 grid, 98 521 in-regime points per n): 3.91%
of the genre grid at n = 1, 0.15% at n = 4, 0.03% at n = 16; worst
CR_hm ≈ −0.067 at (ph ≈ 0.9985, P₁ ≈ 0.617). The skeptic reproduced the
three spot instances from 004's *verbal* coupling description to 4e-15
agreement (S1).

### 3. Block branch: iid-given-k REFUTED (Θ(n), the ψ threshold);
adaptive variant fails strictness only

Per coordinate, the shared-latent iid-given-k coupling turns a
p-component into union marginal 1−(1−p)², so its gain is
h((1−p)²) − h(p), which **changes sign exactly at p = ψ** (skeptic S5:
the bisection root equals (3−√5)/2 to 1e-12 — the AHS constant surfacing
as a *component-wise* threshold inside mixtures). Consequences, all
in-regime (B1/B2; profile arithmetic cross-checked against direct 2^n
enumeration to 1e-8 at n = 8, 12 by the skeptic, S3):

    2block  {0:.62, .95:.38}         (marg .361): Gain = −0.0993n + O(1)
    3block  {0:.60, .90:.30, .99:.10}(marg .369): Gain = −0.1244n + O(1)
    hi+lo   {.10:.50, .60:.50}       (marg .350): Gain = −0.0521n + O(1)

CR ≤ Gain always (009 fact (i)), so the iid-block branch is Θ(n)-negative
on any in-regime mixture carrying weight above ψ — even `hi+lo`, where
the high component is only 0.60. The **adaptive** variant (009 part D's
m_k rule: p ≤ 1/4 → union 2p; p ∈ (1/4, 1/2] → union 1/2; p > 1/2 →
comonotone) does better: on `hi+lo` its full-history CR is positive at
every n tested (+0.071 at n = 2 rising to +0.638 at n = 10, B2b), but on
mixtures whose every nondegenerate component exceeds 1/2 it degenerates
to the comonotone coupling, and comonotone CR = 0 *exactly* (identity:
B_{<i} adds nothing to A_{<i}, so CR = −T_A = 0; numerically −2.4e-15,
skeptic −4.8e-15). Zero is not > 0: the adaptive branch fails the strict
(TAX at p) inequality on the all-above-half corner rather than going
negative.

### 4. The repair: one branch covers everything the two broken ones owned

**Lemma EM-coverage (proved modulo two float-level 1-D sweeps).** For
every μ = P₁δ∅ ⊕ (1−P₁)Bern(ph)^⊗n with ph ∈ [1/2, 1), P₁ ∈ (0, 1) and
(1−P₁)ph < 0.38271: sup_q CR > 0, i.e. the ∅-mixing coupling with
optimized q always certifies. Proof structure (constants recomputed
independently by the skeptic, S4):

- With b = (1−ph)^n, m₀ = μ(∅), the criterion F′(P₁) < 0 reads
  G := (1−b)·log₂(m₀/(1−P₁)) + n·h(ph) + b·log₂ b > 0, and G is
  strictly increasing in P₁ (both terms are).
- **P₁ ≥ 1/2, any n:** m₀ ≥ P₁ ≥ 1−P₁ makes the first term ≥ 0, and
  n·h(ph) + b·log₂ b > 0 for ph ∈ (0,1). Rigorous.
- **P₁ < 1/2, n ≥ 3:** in-regime forces P₁ > 1−2(0.38271) = 0.23458 and
  ph < 0.76542, so h(ph) > h(0.76542) = 0.78591, the first term is
  ≥ −log₂(0.76542/0.23458) = −1.7062, and −b·log₂ b ≤ 0.53074; margin
  3(0.78591) − 0.53074 − 1.7062 = +0.1208 > 0. Rigorous.
- **P₁ < 1/2, n ∈ {1, 2}:** by monotonicity in P₁ it suffices on the
  boundary P₁ = 1 − 0.38271/ph; min G = +0.3448 (n = 1) and +0.8658
  (n = 2) over 200 001-point sweeps (float-level; skeptic re-checked the
  sign by central-differencing brute-force entropies, 1001/1001 negative
  F′). The margins are ~10¹² times float error, but this leg is
  computer-checked, not hand-proved.

On the failure instances: CR(q = 1/2) ≈ +0.067…+0.073 everywhere the
stated branch is negative (A2), and the A3 grids found **zero** in-regime
points with sup_q CR ≤ 0. The 3-block all-above-half mixture — where the
adaptive block branch flatlines at 0 — is also covered, because Lemma EM
never needed ν to be a product: ∅-mixing with s ~ the conditional mixture
gives CR = +0.693 at q* = 0.20 (closed form; cr_eval confirms to 1e-9;
skeptic reproduces from its own triple-enumeration builder). The Sinkhorn
sweep also rescues every constructed failure instance (sup_λ CR ∈
[+0.056, +0.224], λ* including 0), so the recipe correction is
overdetermined: route the `{0}∪[1/2,1)` genre and all-above-half mixtures
to ∅-mixing-with-optimized-q (provable, closed form), or fall back to the
λ-sweep (measured).

### 5. Floor ratios (queue 1(b) pre-work)

From `data/cr_attack.json` (023's committed run): over the in-regime,
H-bounded endpoints, CR/H(μ) = 0.0834, 0.0387, 0.0434, 0.0522, 0.0695 —
min **0.0387 at the mmabskill n = 5 floor**. Endpoint H vs seed H: 2.68
vs 2.54, 1.92 vs 2.65, 2.33 vs 2.65 — the anneal cut CR ~5× while moving
H by ≤ 1.4×, so 023's floors are genuine CR floors, not entropy
artifacts; the one H-destroyed endpoint (windowkill_n6, H → 0) is the
known degenerate corner. Any H-scaled target CR ≥ c·H(μ) on the generic
genre must have c ≤ 0.0387 on current evidence (per-trajectory caveat
from 024 applies: these are floors of *some* runs).

## Outcome

- **REFUTED (the two stated branch assignments of the 009 recipe).**
  Half-mixing (q = P₁/2): in-regime CR < 0 instances at every n
  (constructed for n ≤ 100, mechanism analytic, margins ~1.5e-2 ≫ float
  error), confirmed by two independent implementations agreeing to
  ≤ 4e-15. Block iid-given-k: Θ(n)-negative in-regime via the ψ component
  threshold, slopes matching the per-coordinate prediction to 3 decimal
  places at n = 40. Block adaptive: CR = 0 exactly (not > 0) on
  all-above-half mixtures. Float-level, not exact-rational — exact
  certification of one instance per family is queued (the 022 log₂ kit
  applies; lead 3).
- **Lemma EM and its coverage corollary: proved here** (hand derivations
  in §1/§4, the two n ∈ {1,2} legs computer-checked at float level);
  implementation-level skeptic pass done same-session by
  `uc_branch_skeptic.py` (natural-log, own builders, S1–S6, zero
  refutations after two reporting-level catches — a truncated-input
  comparison and a misquoted constant 0.78603 → 0.78591, both fixed).
  **Reviewer-level independence pending** per 024's standard; until a
  fresh-session pass re-derives EM/EM-coverage, treat them as
  candidate-VERIFIED.
- **EVIDENCE (Gap 2 survival, strengthened in scope):** with the recipe
  corrected to ∅-mixing on the affected genres, no in-regime μ with
  CR ≤ 0 is known: A3's 5 × 98 521-point genre grids, every constructed
  block instance, and 023/024's Sinkhorn-branch endpoints all end
  positive under the corrected assignment.
- **Not claimed:** any lower bound on CR; any statement about
  non-exchangeable mixtures or components varying per coordinate (the
  block instances are exchangeable); exact-rational certification of any
  number here; that q-optimized ∅-mixing covers *product mixtures with a
  below-half component* (B2b's positive adaptive values cover the one
  family tested, `hi+lo`, and nothing more); reviewer-level independence
  of the new lemmas.

## Why it failed / what survived

Gap 2's candidate survived again, but the *recipe* — Gap 3's current
draft — did not: both non-Sinkhorn branches fail exactly where their
coupling stops buying entropy. The mechanism is the same in both kills:
a fixed O(1) bookkeeping cost (mixing-weight entropy h(P₁) − h(P₁/2);
mixture-boundary terms) must be paid from a per-coordinate entropy gain
that the branch's *fixed* coupling choice lets collapse to zero
(h(ph) → 0 as ph → 1; h((1−p)²) − h(p) < 0 past ψ). The repair works
because Lemma EM exposes the whole q-family at zero tax, and concavity
hands over the optimum. Two structural facts worth keeping:

- **ψ is a component-wise threshold.** Inside any mixture, iid-coupling a
  component above ψ is Θ(n)-toxic no matter how light the component or
  how small the overall marginals. Any total recipe must never iid-couple
  above-ψ mass — a sharp, checkable design rule for Gap 3.
- **The surprisal criterion.** sup_q CR > 0 for ∅-mixing iff ∅'s
  surprisal is below the ν-average — a one-line test any future recipe
  can apply before choosing this branch, and the natural boundary object
  where ∅-mixing hands off to the tilt branch.

Reusable: `uc_branch_attack.py` (closed forms F/G, deterministic grid
engines, adaptive/iid block builders, profile arithmetic with empty-atom
mixing); `uc_branch_skeptic.py` (independent CR accounting, B-major
grouping, natural-log — a second full-history evaluator for future
skeptic passes); the failure families as hard instances for any
quantitative TAX statement (the sliver instances have H(μ) ≈ 0.97 and
CR under the *best known* coupling ≈ +0.068 — tight test points for
CR ≥ c·H).

## Leads generated

1. **Fresh-session reviewer pass on EM/EM-coverage** (the 024 standard):
   re-derive both lemmas, re-run the two engines, independently rebuild
   one violating instance per family. Until then the lemmas are not
   VERIFIED.
2. **Restate the recipe (Gap 3) with the corrected assignments** and
   re-pose (TAX at p) against it; then resume queue 1(b) proof effort
   with the c ≤ 0.0387 envelope and the sliver instances added to the
   test battery.
3. **Exact certification of one kill per family** (022 log₂ kit +
   dyadic accumulation): the half-mixing closed form needs only certified
   log₂ at rational (ph, P₁), which the kit already does.
4. **The uncovered cell: mixtures with both below-half and above-ψ
   components** (e.g. `hi+lo`). Adaptive-block is positive on the one
   family tested; either prove it (per-coordinate gain ≥ 0 with the
   O(1) cross-terms bounded by mixture entropy?) or find its failure
   corner. This is the last branch region with neither a proof nor a
   kill.
5. **ψ-threshold ripple:** the component-wise ψ rule may bite other
   routes that couple mixtures (008's assembly budgets; 021's
   transfer-operator lead) — a cheap `/ripple` scan candidate.

## References

- This repo: 009/011 (candidate, HM lemma, recipe statement), 023/024
  (Sinkhorn-branch attack, floor instances, determinism standard),
  004 (half-mixing coupling R1), 020/022 (in-regime convention, log₂
  kit), `docs/SWARM.md` §independence (skeptic-pass rules).
- No external sources.
