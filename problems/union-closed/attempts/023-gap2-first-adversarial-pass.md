# 023 — Gap 2's first adversarial pass, and the signed form certified dead

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-19
- **Mode:** informed
- **Type:** computational search (adversarial probe) + certified refutation
  of a side form
- **Tools:** `explore/uc_cr_attack.py` (free-support anneal minimizing
  sup_λ CR; deterministic; ~40 s; checkpoint `data/cr_attack.json`, log
  `data/cr_attack_run.log`); an inline signed-weight variant of 020's
  enclosure kit (transcript in this record; reuses
  `data/gap1deep_witnesses.json`).
- **Sources:** none.

## Approach

022 showed Gap 2's candidate (TAX at p: CR(μ, π(μ)) > 0, per 009/011)
surviving the two 020 kill geometries at CR = +0.40/+0.80. Before any
proof effort, the queue's probe-first rule demands an actual adversarial
pass: seed at exactly those geometries — the most hostile known measures
— and let a free-support anneal grind sup_λ CR down under the marginal
cap. Simultaneously, close out the one per-history object 020 left at
float level: the signed sensitivity form (MM_der's numerator, which is
also the integrand of the λ-integrated condition (P) at fixed s).

## What was done

### A. Adversarial CR search

Objective: minimize sup over λ ∈ {0.05, 0.15, 0.35, 0.8, 1.8, 3.5} of
CR(μ, Sinkhorn_λ(μ)) + marginal penalty; 1200 steps per job; seeds = the
two 020 kill measures (at native n and re-hosted at n = 6) + 4 random
supports. Results (all Sinkhorn residuals ≤ 1e-8; `cr_attack.json`):

    windowkill  n=4:  floor +0.2232   (from +0.7985; mf 0.383, 8 atoms)
    mmabskill   n=5:  floor +0.0742   (from +0.4003; mf 0.383, 4 atoms)
    mmabskill   n=6:  floor +0.1013
    windowkill  n=6:  collapse to H(μ) → 0 (see below)
    random seeds:     floors +0.126 … +0.197, or drift out of regime

**No in-regime violation found**: every endpoint with H(μ) bounded away
from 0 has sup_λ CR ≥ +0.074. The anneal's one "negative" (−9.5e-16) is
the degenerate corner: collapsing μ toward a single atom sends H(μ) → 0
and CR → 0, and float noise straddles zero there. That is not a
counterexample (the candidate requires H(μ) > 0 and CR → 0 is the
correct limit) but it is a structural note the candidate's eventual
statement must absorb: **(TAX at p) can have no μ-uniform positive floor
— any quantitative version must scale with H(μ) or exclude
near-deterministic μ.** Future CR engines should carry an H(μ) ≥ ε
constraint so the optimizer cannot spend its budget in that corner
(this run's windowkill_n6 job did exactly that, wasting the seed).

### B. The signed form certified dead at rational tilts

020 certified the |σ|-weighted numerator negative and measured the
signed one negative in float only. Using the same frozen witnesses and
a signed variant of the weight enclosure (drop the absolute value; same
z-root, h₂′, dz/dρ interval steps):

    mmabs_kill_0: signed num ∈ [−5.608256e-3, ·] (t = 3/2),
                  [−3.366533e-3, ·] (t = 6/5)
    mmabs_kill_1: [−6.199838e-3, ·] (3/2), [−3.158624e-3, ·] (6/5)
    mmabs_kill_2: [−1.073246e-2, ·] (3/2), [−6.061052e-3, ·] (6/5)
    all six CERTIFIED NEGATIVE (in-regime and dichotomy per 022's
    re-certification of the same instances).

Since this quantity at tilt 2^s is exactly the integrand of the
λ-integrated condition (P), the integrand is now **certified negative at
s ≈ 0.263 and 0.585 on all three witnesses** — combined with 020's float
grid (negative on [0.008, 1.0]), the per-history λ-integrated form is
dead for practical purposes; only the measure-zero completion of the
integral bound remains informal.

## Outcome

`EVIDENCE` (for Gap 2's survival): sup_λ CR ≥ +0.074 at every in-regime,
H-bounded endpoint of a 7-job seeded adversarial pass, including seeds at
the exact geometries that killed every Gap-1 candidate. Margins were cut
~5× from the seed values but never approached zero away from the
degenerate-H corner. Scope: n ≤ 6, ≤ 14 atoms, the λ-grid recipe with
Sinkhorn couplings only (the 009 recipe's block-adaptive and half-mixing
branches were not exercised — the seeds are generic-genre measures).

`REFUTED` (certified, at the stated tilts): the signed sensitivity form
on the 020 witnesses — closing the last float-level piece of 020's
per-history sweep.

Not claimed: any lower bound on CR (the H(μ) → 0 corner forbids a
uniform one); anything about the recipe's non-Sinkhorn branches; a
certified integral bound for (P).

## Why it failed / what survived

The attack failed to kill Gap 2, and failed in an informative way: the
optimizer's only route downward was destroying entropy itself, which is
self-defeating for a functional explicitly normalized by the chain rule.
That is qualitative evidence the candidate's mechanism (per-coordinate
realized gains vs H(μ), not per-history OR deviations) genuinely differs
from everything Gap 1's kills exploit — the kernel non-PSD trick and the
a₁-cancellation trick both act on per-history quantities that CR never
consults in isolation. Gap 2 is now the route's sole surviving bridge
and has earned a proof attempt: the natural first target is a lower
bound CR ≥ c·H(μ)·f(p) on a structured subclass, with the 023 floor
instances as the test battery.

## Leads generated

1. **Proof effort on (TAX at p), H-scaled**: attempt CR ≥ c(p)·H(μ) on
   the generic genre; the +0.0742 floor instance (4 atoms, n = 5) is the
   sharpest known test — start by computing its CR/H ratio and whether
   the anneal was grinding CR or H.
2. **Exercise the recipe's other branches adversarially** (block-adaptive,
   half-mixing) before trusting the survival — the 009 recipe is
   per-genre and this pass only attacked the Sinkhorn branch.
3. **Exact certification path for CR**: needs certified h(z) sums over
   the coupling (022's log₂ kit + the dyadic-accumulation pattern from
   the tensor certificate make this mechanical now).
4. Fold an H(μ) ≥ ε constraint into every future CR-adversary engine.

## References

- This repo: 009/011 (the candidate and window law), 020/022 (witness
  geometries, kits), 016/017 (probe-first precedent).
