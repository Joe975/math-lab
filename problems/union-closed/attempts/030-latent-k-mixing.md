# 030 — K-component latent-mixing: the mixture branch generalizes, and near-ψ mixtures pin the floor

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-19
- **Mode:** informed
- **Type:** construction + adversarial battery (028 lead 2 / queue 1(b')
  residue): define latent-mixing for K ≥ 3 components and test it on the
  genre v2 left uncovered. Persistence window (round 2), same author as
  025–028; the fresh-context reviewer pass (029,
  VERIFIED_WITH_CORRECTIONS — every load-bearing 025–028 claim
  survived) ran earlier this window and does not cover this record.
- **Tools:** `explore/uc_latentK_probe.py`, `explore/uc_latentK_skeptic.py`
  (both new; stdlib only; deterministic — grid/refinement, no RNG;
  checkpoint `data/latentK_probe.json`). Reproduce: run the two scripts
  in that order.
- **Sources:** none.

## Approach

028's recipe v2 stated the mixture branch only for two components plus an
optional δ∅ atom; ≥ 3 nondegenerate components had no assigned coupling
at all. The natural generalization: couple the component labels
(k_A, k_B) by a K×K matrix Q with both marginals P — diagonal cells use
the adaptive per-coordinate joint (009 part D), off-diagonal cells take
the two products independent. Q = diag(P) is exactly the adaptive-block
coupling, so the family can only improve on it; parametrize by symmetric
transfers t_jk ≥ 0 (q_jk = t_jk, diagonal reduced accordingly) and
optimize deterministically (coarse grid + refinement — no anneal, no
RNG). A K = 2 anchor run must reproduce 027's latent-mixing values.

## What was done

**Anchor (K = 2).** At 027's probe points the optimizer reproduces the
recorded values exactly: +0.22098 / +1.18930 / +1.16759 (Δ < 1e-5 from
027's q-sweep values), with the diagonal matching the recorded adaptive
values to 1e-6. The transfer parametrization and 027's q-parametrization
agree where they overlap, as they must.

**Battery (K = 3, 15 rows, all in-regime, n ∈ {2, 4, 6}).** Five
instance shapes; every row positive after optimization:

    lo+mid+hi  {(.5,.01),(.3,.55),(.2,.9)}:   CR* +0.26→+1.49, diag +0.007→+0.038
    spread     {(.6,.05),(.25,.7),(.15,.99)}: CR* +0.22→+0.95, diag +0.041→+0.284
    near0+2hi  {(.62,.001),(.22,.9),(.16,.99)}: CR* +0.20→+0.62, diag ~+0.01
    all-low+1mid {(.4,.2),(.35,.3),(.25,.55)}: CR* +0.19→+0.61 ≈ diag
    near-psi   {(.4,.36),(.3,.39),(.3,.40)}:  CR* = diag = +0.083→+0.248

Optimal transfers always mix the lightest label against the heavy ones
(t01 dominant, saturating its cap on near0+2hi); on near-ψ mixtures the
optimum is exactly the diagonal — transfers never help there.

**Degeneration probe.** On near0+2hi with the light component's p → 0
(1e-2 … 1e-4, n = 4): the K = 3 optimum holds at +0.4037/+0.3970/+0.3950
while the adaptive diagonal collapses (+0.038 → +0.0003) — and the K = 3
optimum **strictly beats the 2-block EM limit (+0.3523)** in which the
two heavy components are fused into one s-law. Pairing the light label
fully against one specific heavy component keeps structure that fusing
destroys: K-component latent-mixing is strictly stronger than
∅-mixing-with-mixture-s-law, not just equal in the limit.

**The near-ψ observation that became a theorem.** On the near-ψ battery
rows CR/H = 0.04303 exactly, constant in n, because CR = n − H(μ) to
machine precision (≤ 6e-14 at n = 6). That identity is provable:

  **Mini-theorem HU-mix (proved here; skeptic-passed numerically,
  reviewer pass pending).** Let μ = Σ_k P_k Bern(p⃗_k)^⊗n with every
  p_{k,i} ∈ [1/4, 1/2] and all element marginals ≤ p̄ < 1/2. The
  shared-label coupling whose per-coordinate joints have both-zero
  probability exactly 1/2 satisfies

      CR(μ, π) = n − H(μ) ≥ n(1 − h(p̄)) ≥ ((1 − h(p̄))/h(p̄)) · H(μ).

  *Proof.* In any history cell (a, b), the conditional both-zero
  probability is a posterior mixture over k of the per-component value
  1/2, hence exactly 1/2; so E[h(z_i)] = h(1/2) = 1 and Σᵢ E[h(z_i)] = n
  exactly. (Also U ~ uniform on 2^[n], so Gain = n − H(μ) = CR and
  ST = 0.) And H(μ) ≤ Σᵢ h(margᵢ) ≤ n·h(p̄) since h is increasing on
  [0, 1/2]. ∎

  At p̄ = 0.38271 the constant is (1−h)/h = 0.041739 — the first
  **proved** H-scaled floor on any sub-genre, and it sits just below the
  measured G-gen floor ratios (0.0387–0.083, 025) and the near-ψ battery
  value (0.04303 at marginal 0.381, matching (1−h(0.381))/h(0.381) up to
  the sub-cap slack). Near-ψ measures pin the constant in every genre
  tested so far. This theorem seeded the 031 line (the per-history
  totalization of the m = 1/2 rule).

## Outcome

**EVIDENCE (coverage):** K-component latent-mixing is defined, anchors
to 027 exactly, and is positive on all 15 battery rows (n ≤ 6, K = 3,
exchangeable components) with CR*/H ∈ [0.043, 0.395]; the ≥ 3-component
hole in recipe v2 is now covered at the tested instances by the same
one-family mechanism.

**Proved (pending reviewer pass): mini-theorem HU-mix** — CR = n − H(μ)
with the H-scaled floor (1−h(p̄))/h(p̄)·H(μ) on the all-components-in-
[1/4,1/2] sub-genre.

**Not claimed:** optimality of the transfer parametrization (only
symmetric Q explored); anything for non-exchangeable components or
K > 3; any lower bound for the optimized family outside the sub-genre
theorem; reviewer-level independence (029's scope is 025–028; this
record joins the next reviewer batch).

## Why it failed / what survived

Nothing failed. Structure gained: (i) transfers help exactly where
components leave [1/4, 1/2] — the sub-genre theorem explains why (inside
it, the diagonal already achieves the per-coordinate entropy maximum
E h(z) = 1, so no label-shuffling can add anything); (ii) the extremal
role of near-ψ mixtures is now a theorem-backed anchor for the H-scaled
constant c(p̄) = (1−h(p̄))/h(p̄); (iii) fusing components loses real CR —
recipe totality for mixtures should keep components separate and couple
labels, not collapse them into one s-law.

## Leads generated

1. Prove positivity of the optimized K-family outside the sub-genre
   (the 027 leads 1–2 now generalize: adaptive ≥ 0 as the diagonal
   specialization plus a transfer-monotonicity argument).
2. Non-exchangeable components (p⃗_k varying per coordinate): the
   HU-mix proof survives verbatim if every p_{k,i} ∈ [1/4, 1/2] —
   battery it.
3. The 031 line (opened same-window): totalize the m = 1/2 rule
   per-history to arbitrary μ.

## References

- This repo: 027 (K = 2 family, anchor values), 028 (recipe v2, the
  uncovered cell), 025 (EM lemma; the fused-s-law construction beaten
  here), 009 (adaptive rule, part D), `data/latentK_probe.json`.
- No external sources.
