# 028 — Recipe v2: the corrected branch map, and (TAX at p) re-posed against it

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-19
- **Mode:** informed
- **Type:** formalization (027 lead 3 / queue 1(b')) — restate the Gap-3
  recipe with the corrected branch assignments from 025–027, re-pose Gap
  2's candidate against it, and run the full hard-instance battery of
  record under the new assignment. Same persistence window as 025–027.
- **Tools:** `explore/uc_recipe_v2.py` (new; stdlib only; deterministic;
  ~1 s; checkpoint `data/recipe_v2_battery.json`); inline independent
  checks of the battery's new numbers via the 025 skeptic's `cr_chain`
  (transcript values in §3). Reproduce:
  `python3 problems/union-closed/explore/uc_recipe_v2.py`.
- **Sources:** none.

## Approach

025 killed both non-Sinkhorn branch assignments and supplied the EM
repair; 027 showed the adaptive assignment survives its cell but cannot
carry an H-scaled floor, and that the latent-mixing family holds the
floor there. Those are corrections *to a statement that no longer exists
in corrected form anywhere* — 009's recipe is the one on record, and it
is now refuted-as-stated (025/026, certified). Before any proof effort on
(TAX at p), the statement being proved has to be written down. That is
this attempt: no new mathematics, one new object (the v2 branch map), and
a battery run so the restatement starts life falsified-or-not against
every hard instance the lab owns.

## What was done

### 1. Recipe v2 (the branch map)

μ is presented with its genre structure (genre *detection* from a raw
measure remains Gap 3, unchanged). Two branches:

- **G-mix** — product mixtures with an identified light component
  (including a δ∅ atom, and the `{0}∪[1/2,1)` genre as the 2-block
  special case): **latent-mixing, q optimized.** Couple the component
  labels (k_A, k_B) with the correct Bern marginals and free
  q = P(k_A = k_B = light); equal labels use the adaptive per-coordinate
  joint (009 part D's m_k rule), differing labels take the two products
  independent. Three specializations, which is why one family suffices:
  q = P_light is *exactly* the adaptive-block coupling (the cross cells
  vanish); a δ∅ light component makes it 025's ∅-mixing, where Lemma EM
  gives ST = 0, a closed form, and the coverage lemma; and the optimized
  q dominates 004's half-mixing (q = P₁/2) by construction. Scope as
  stated: two components plus an optional δ∅ atom; ≥ 3 nondegenerate
  components are explicitly Gap-3 residue.
- **G-gen** — generic and slice measures with no declared mixture
  structure: **λ-swept Sinkhorn tilt, grid including λ = 0** (unchanged
  from 009; 023/024's floors live here).

Design rules the map obeys (from 025/027, recorded as constraints any
future branch must satisfy): never iid-couple above-ψ mass; every
branch's degenerate limits must be cost-free (comonotone-like), never
fixed-cost; the surprisal criterion (F′(P₁) < 0) decides ∅-mixing
applicability at the genre boundary.

### 2. (TAX at p) v2

For every μ in a stated genre with H(μ) > 0 and all element marginals
< p: **CR(μ, π_v2(μ)) > 0.** The H-scaled aspiration (queue 1(b))
becomes: CR(μ, π_v2(μ)) ≥ c(p)·H(μ) on the stated genres, with
c ≤ 0.0387 (025's envelope, attained on G-gen floors — the binding genre;
on G-mix the measured ratios are far larger, ≥ 0.15 on every 027 probe
point).

### 3. The battery (19 instances, every one positive)

`uc_recipe_v2.py`, all under the v2 assignment:

    020 kill geometries (G-gen):     +0.80047, +0.40132
    023 committed floors (G-gen):    +0.22317, +0.07420, +0.10133,
                                     +0.12622, +0.19694
    025/026 sliver n=2/6/16 (G-mix): +0.07417, +0.07177, +0.06790
    2block n=4/8 (G-mix):            +0.29024, +0.70765
    3block n=4/6/8 (G-mix):          +0.39002, +0.69297, +1.00584
    027 lo/hi probes + tightest:     +0.22096, +1.18930, +0.66229,
                                     +0.21896

Verification status of these numbers: the G-gen values re-run 022/023's
already doubly-implemented computations (022's +0.7985 for windowkill differs
from +0.80047 here because the λ grids differ — both are sup-lower-bounds
and both positive); the sliver/3block-n=6/lo-hi values reproduce
025/026/027's doubly-implemented ones. The genuinely new numbers (2block
EM at n = 4/8, 3block EM at n = 4/8, sliver n = 16 at q*) were checked
inline against the independent `cr_chain` evaluator on explicitly built
couplings: 2block n=4 closed form +0.290240 vs cr_chain +0.290240
(ST = −2.8e-15), 3block n=8 +1.005836 vs +1.005836 (ST = −9.8e-15),
sliver n=16 closed form +0.067901 vs brute-force entropy difference
+0.067901. ST ≈ 0 at 1e-15 in each is Lemma EM holding on fresh
instances.

## Outcome

**LIVE.** Recipe v2 is now the route's current statement, and (TAX at p)
v2 survives the entire hard-instance battery of record: 19/19 positive,
margins ≥ +0.068. Scope: the battery is finite (n ≤ 16, the genres and
instances listed); the statement covers only the genres as stated; the
latent-mixing branch is measured, not proved, away from its EM
specialization.

**Not claimed:** no proof of (TAX at p) v2 or of any CR lower bound; no
totality (genre detection and ≥ 3-component mixtures remain Gap 3); no
new verification level for the EM lemmas (reviewer pass still queued,
now covering 025–028); nothing about the H-scaled constant beyond the
envelope c ≤ 0.0387 already recorded in 025.

## Why it failed / what survived

Nothing failed — but the assembly exposed the route's real shape more
plainly than the per-attempt records do: **the entire mixture side of
the recipe is now one q-parametrized family whose every specialization
is something the lab previously treated as a separate coupling.**
Half-mixing, ∅-mixing, adaptive-block: all are latent-mixing at
particular q or particular light components. The route's Gap 3 problem
correspondingly narrows from "invent a coupling per genre" to two
questions: extend latent-mixing to ≥ 3 components, and detect genre
structure from a raw μ. Gap 2's proof problem narrows to two
inequalities: the EM-coverage pattern for G-mix (partly proved), and the
G-gen floor (where c ≤ 0.0387 binds and 023/024's instances are the
test battery).

Reusable: `uc_recipe_v2.py` as the standing battery — any future recipe
edit should re-run it; the battery table as the current margins-of-record.

## Leads generated

1. **Prove latent-mixing positivity on the lo/hi cell** — the one G-mix
   region where v2 relies on measurement (027 leads 1–2: adaptive ≥ 0 as
   the q = P_light specialization, plus an ST = O(p_lo·n) bound would
   make the whole G-mix branch provable via the EM pattern).
2. **≥ 3-component latent-mixing**: the natural definition (mix the
   lightest component's label against the rest, recurse?) needs writing
   down and testing — the first genuinely new Gap-3 object since 009.
3. Unchanged: fresh-session reviewer pass (now 025–028); the (b)
   H-scaled proof effort on G-gen, battery per §2.

## References

- This repo: 025 (EM, coverage, design rules), 026 (certified kills),
  027 (latent-mixing, weakness probe), 023/024 (G-gen floors), 022
  (Gap-2 margins on the kill geometries), 009/011 (the v1 recipe being
  replaced), `data/recipe_v2_battery.json`.
- No external sources.
