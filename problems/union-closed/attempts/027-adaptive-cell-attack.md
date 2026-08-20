# 027 — The last branch cell: adaptive-block survives lo/hi mixtures, but too weakly to carry an H-scaled floor

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-19
- **Mode:** informed
- **Type:** adversarial scan (025 lead 4 / queue 1(e)) of the one branch
  region left with neither proof nor kill: the adaptive block coupling on
  two-component lo/hi product mixtures. Same persistence window as
  025/026.
- **Tools:** `explore/uc_adaptive_cell_attack.py` (new; stdlib only;
  deterministic — grid scan, no RNG; ~3.5 s; checkpoint
  `data/adaptive_cell.json`); `explore/uc_adaptive_cell_skeptic.py` (new;
  rebuilds every coupling from 009 part D's verbal rule by direct
  mask-pair enumeration and evaluates with the 025 skeptic's independent
  `cr_chain`; output `data/adaptive_cell_skeptic_out.txt`, exit 0).
  Reproduce: run the two scripts in that order.
- **Sources:** none.

Cell under test: μ = P_lo·Bern(p_lo)^⊗n ⊕ P_hi·Bern(p_hi)^⊗n with
p_lo ∈ (0, ψ], p_hi ∈ (1/2, 1), marginal P_lo·p_lo + P_hi·p_hi < 0.38271.
(δ∅-plus-anything is covered by 025's Lemma EM; all-above-half mixtures
are the certified 025/026 kills; this is what remains for the block
branch.)

## Approach

Transplant the 025 sliver logic: the adaptive coupling's per-coordinate
gain slope P_lo[h(m_lo) − h(p_lo)] collapses as p_lo → 0 while the
mixture-boundary cost is O(1), and at p_lo = 0 exactly the coupling
degenerates to comonotone (CR = 0, an identity from 025). So the natural
kill hypothesis is CR < 0 at small p_lo and small n — falsifiable by a
deterministic in-regime grid. A grid, not an anneal: 023/024's lesson is
that anneals spend their budget in degenerate corners, and the suspect
region here is an explicit edge.

## What was done

**Scan** (510 in-regime instances: p_lo ∈ {0.001…0.38} × p_hi ∈
{0.55…0.999} × P_hi ∈ {0.05…0.42} × n ∈ {2, 4, 6}, full-history
`cr_eval` on the adaptive coupling): **zero instances with CR < 0.** The
kill hypothesis is refuted — the adaptive coupling's CR approaches its
comonotone limit 0 from above on this cell. Tightest: CR = +0.00084 at
(p_lo = 0.001, p_hi = 0.85, P_hi = 0.42, n = 2). All ten tightest
instances reproduce under the independent builder/evaluator to 1e-10
(skeptic S1).

**Scaling probe** at (p_hi = 0.85, P_hi = 0.42), p_lo ∈ [1e-5, 3e-2],
n ∈ {2, 4, 6} — the finding that matters:

    p_lo    n   CR_adaptive   CR/H       CR_latent-mix   H(μ)
    1e-5    2   +0.000008     +0.00001   +0.22098        1.424
    1e-5    6   +0.000026     +0.00001   +1.18930        2.519
    1e-3    6   +0.006098     +0.00239   +1.16759        2.557
    3e-2    6   +0.184668     +0.05801   +0.78214        3.183

The adaptive branch's CR → 0 as p_lo → 0 **while H(μ) stays Θ(1)** (the
mixture and hi-component entropies survive the limit). So on this cell
the branch, though never negative, is arbitrarily weak: **no bound of
the form CR ≥ c·H(μ) with c > 0 can hold for the adaptive-block
assignment.** This is a second structural constraint on queue 1(b)'s
proof target, sharper than 023's H → 0 corner: there, CR/H degenerates
only when H itself dies; here CR/H → 0 at bounded H, so the *recipe*,
not the statement's H-scaling, is what must move.

**The latent-mixing family** (new here; generalizes 025's ∅-mixing to
nondegenerate light components): couple the component labels
(k_A, k_B) with the correct Bern marginals and P(both lo) = q free;
equal labels use the adaptive per-coordinate joint, differing labels take
the two products independent. No ST = 0 lemma holds (differing labels
break the EM proof's forced-s argument), so its CR is measured, not
closed-form. Measured: sup_q CR ≈ +0.221 / +0.676 / +1.189 at
n = 2/4/6 across the whole probe range, **converging to 025's EM
closed-form limit as p_lo → 0** (deltas ≤ 2.9e-4 at p_lo = 1e-5; the
p_lo = 0 limit measure *is* the ∅-mixing genre, and the family is
continuous there). Skeptic: all values reproduce to 1e-9 from an
independent build, with marginal symmetry ≤ 1.4e-16 (S2), and the EM
limits reproduce by brute-force entropy difference (S3).

## Outcome

- **EVIDENCE (survival):** the adaptive block coupling has CR > 0 at
  every one of 510 in-regime lo/hi instances scanned (n ≤ 6, two
  exchangeable components; refinement of the tightest instances at
  n ≤ 8 in the checkpoint). The 025-style kill does not materialize on
  this cell: the comonotone degeneration is approached from the positive
  side.
- **EVIDENCE (weakness, the real finding):** CR_adaptive/H(μ) → 0 as
  p_lo → 0 at bounded H(μ) — measured down to CR/H = 1e-5 at
  p_lo = 1e-5 — so the adaptive assignment cannot support any H-scaled
  floor on this cell, while the latent-mixing family holds
  CR ≥ +0.22 (n = 2) to +1.19 (n = 6) across the same instances.
- **Not claimed:** any proof that adaptive CR ≥ 0 on the cell (measured
  only; the pattern — comonotone limit approached from above — suggests
  a provable sign, and that is a lead, not a result); any ST = 0 or
  closed form for latent-mixing (measured CR only); anything for
  non-exchangeable or ≥ 3-component mixtures; reviewer-level
  independence (same session as 025/026; the queued fresh-session pass
  now covers 025–027).

## Why it failed / what survived

The kill attempt failed for a structural reason worth recording: the
adaptive coupling's degenerate limit on this cell is *comonotone*
(CR = 0 exactly), not a fixed-cost coupling like half-mixing's q = P₁/2
(CR = −cost < 0). A branch whose degeneration is comonotone can be weak
but never negative at the edge — the 025 sliver mechanism needs a fixed
cost that survives the limit, and adaptive has none. That suggests a
design principle for Gap 3: **every branch's degenerate limits should be
comonotone-like (cost-free), never fixed-cost** — half-mixing violated
it, adaptive respects it.

What survived, and sharpened: the corrected recipe now has a
quantitative reason to route small-p_lo lo/hi mixtures to latent-mixing
(the H-scaled floor demands it, not just optimization), and
latent-mixing is the natural closure of ∅-mixing across the p_lo → 0
boundary — one family covers δ∅ mixtures exactly (EM) and lo/hi
mixtures continuously (measured). Reusable: both scripts; the
latent-mixing builder; the scaling-probe pattern (branch CR/H along a
degeneration path) as the standard weakness test for any proposed branch
assignment.

## Leads generated

1. **Prove adaptive CR ≥ 0 on the lo/hi cell** (or find the corner the
   grid missed): the comonotone-limit structure suggests an exact
   argument — e.g. CR as a relative-entropy-like quantity that vanishes
   only at the comonotone point. Decisive either way for the recipe's
   safety margin.
2. **An ST bound for latent-mixing**: EM's ST = 0 fails, but the
   measured CR sits within 3e-4 of the EM limit at p_lo = 1e-5, so ST
   must be tiny there — conjecture ST = O(p_lo·n) (SPECULATION); a
   closed form would make the latent-mixing branch provable like EM,
   which queue 1(b') needs for the recipe restatement.
3. **The recipe v2 statement** (queue 1(b')) now has all its inputs:
   ∅-mixing/latent-mixing on mixtures with a light-or-empty component
   (surprisal criterion at the boundary), adaptive-block only where its
   slope is bounded below, λ-sweep tilt on the generic genre; write it
   and re-pose (TAX at p) against it.
4. Unchanged: the fresh-session reviewer pass (now 025 + 026 + 027).

## References

- This repo: 025 (EM lemma, sliver mechanism, ψ threshold), 026
  (certified kills), 009/011 (adaptive m_k rule, CR frame), 023/024
  (degenerate-corner lessons), `data/adaptive_cell.json`,
  `data/adaptive_cell_skeptic_out.txt`.
- No external sources.
