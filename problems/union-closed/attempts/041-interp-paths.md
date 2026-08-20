# 041 — Interpolating between (HU-TAX)'s two equality families: no dip, margins grow with n

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-20
- **Mode:** informed
- **Type:** falsification sweep (040 lead 1).
- **Tools:** `explore/uc_hu_interp.py` (new; both interpolation paths
  as exchangeable product mixtures, evaluated by an exact count-pair
  DP — (k+1)² cells at step k — cross-checked against the general
  evaluator at n ≤ 6 and against the rollout order for
  exchangeability; deterministic, no RNG; checkpoint
  `data/hu_interp.json`); `explore/uc_hu_interp_skeptic.py` (new;
  sampled small-n margins through the independent 037-skeptic stack,
  endpoint identities at every n, summary consistency; exit 0).
  Reproduce: run the two in that order.
- **Sources:** none.

## Approach

040 left (HU-TAX) with two proved equality families at opposite
correlation extremes — products (independent) and two-point diagonals
(perfectly correlated), both attaining c*(p̄) = (h(max(1/2, 1−2p̄)) −
h(p̄))/h(p̄). If the conjecture fails in the mixture world, the natural
place is between them. Two paths, marginals exactly p at every point:

- **A (convex):** μ_ρ = (1−ρ)·Bern(p)^n ⊕ ρ·diag(p) — a mixture of the
  Bern(p), Bern(0), Bern(1) products;
- **B (product-mixture):** μ_t = (1−λ)Bern(a)^n ⊕ λBern(b)^n with
  a = p(1−t), b = p+(1−p)t, λ = (p−a)/(b−a) — t = 0 the product, t = 1
  exactly the diagonal, the whole path inside the mixture-of-products
  genre the recipe's latent-mixing branch owns.

Both are exchangeable, so the coupling is order-free (checked), and
both are mixtures of ≤ 3 products, so the HU coupling collapses to an
exact DP over (ones-in-a-prefix, ones-in-b-prefix) count pairs — the
sweep runs to n = 64 where the generic evaluator is hopeless.

Margin compared against the instance's OWN constant c*(p) (034's
lesson); any negative margin refutes the re-posed (HU-TAX).

## What was done

Grids: n ∈ {2, 4, 8, 16, 32, 64} × p ∈ {0.30, 0.38271, 0.45, 0.49} ×
41 path points, both paths (1968 evaluations), plus 21-point
refinements around each path's global minimum.

    path        negatives   worst margin   interior min (t ∈ (0,1))
    A convex        0        −5.3e-14       +2.20e-04  (n=2, p=0.49, t=0.025)
    B prodmix       0        −5.3e-14       +1.41e-07  (n=2, p=0.49, t=0.025)

- The worst margins are the **equality endpoints at float noise**
  (t = 0 products; the −5.3e-14 is the n = 64 entropy sum's rounding).
  Refined minima sit at the endpoints, not inside.
- **Interior margins are strictly positive and grow with n**: at
  p = 0.49 the interior minimum (always at the first grid point off
  the product) rises +2.2e-4 → +2.1e-3 → +9.1e-3 (path A, n = 2, 8,
  64) and +1.4e-7 → +9.8e-7 → +8.7e-6 (path B). Leaving the product
  family in either direction increases CR/H − c*; there is no valley
  between the two equality families on these paths, at any n tested.
- Cross-checks: DP vs the general evaluator ≤ 1e-10 at n ≤ 6 with
  exact marginals; identity vs rollout order ≤ 1e-12 (exchangeability);
  the skeptic reproduces 50 sampled margins through the independent
  stack and confirms the endpoint identities at every n — the latter
  is a closed-form validation of the DP at n = 16/32/64, where no
  independent evaluator exists.

## Outcome

- **EVIDENCE: (HU-TAX, re-posed) survives the interpolation test** —
  zero violations on 1968 grid points spanning both paths between its
  two equality families, n up to 64, caps up to 0.49, with margins
  strictly positive in the interior and minimized exactly at the known
  equality cases.
- **EVIDENCE (structural): both equality families are locally rigid on
  these paths** — the margin leaves zero with positive slope in t on
  both sides at every (n, p) tested, and the slope grows with n. On
  this evidence the products and diagonals are isolated equality
  cases, not edges of a flat region.
- **Not claimed:** anything off these two paths (non-exchangeable
  interpolations are untested); n > 64; that the paths' positivity is
  monotone in t (only that no grid/refined point dips); any statement
  about K ≥ 3 mixture components beyond path A's three.

## Why it failed / what survived

The falsification attempt failed cleanly: the conjectured bound is not
merely unrefuted between its equality families — it is locally
reinforced there, with the deficit direction nowhere in sight. For the
proof effort this is the useful shape: 040 lead 2's per-cell argument
should expect strict convexity-like behavior along mixture paths, with
equality characterized by cells sitting at the clamp boundary
(product: z = 1/2 everywhere at p ≥ 1/4; diagonal: one live cell, then
determinism).

## Leads generated

1. **Non-exchangeable interpolation:** both paths preserve
   exchangeability, which the equality families share. The cheapest
   symmetry-breaking test: tensor products diag(p)^⊗(n/m) of small
   diagonals (blocks of perfectly-correlated coordinates, independent
   across blocks) — CR is additive over independent blocks, so these
   are equality cases too (each block contributes c*(p)·h(p)·...);
   verify, then perturb BETWEEN blocks. If block-tensors are also
   equality cases, the equality set is much bigger than two families —
   that changes the proof target's shape and deserves its own record.
2. **The count-pair DP is reusable** for any exchangeable
   product-mixture question at large n (the K = 3 battery of 030, the
   latent-mixing family of 027/028) — it is exact and O(n³) per
   evaluation.
3. Unchanged from 040: the per-cell proof shape (lead 2), the 0.497
   anneal saturation (lead 3), and — from 039 — the crash8
   roll-descent close-out.

## References

- This repo: 040 (the DIAG identity and the lead), 034 (the corrected
  constant and the own-marginal lesson), 037/038 (the rollout rule and
  ROLL-DOM), 031/030 (HU, HU-mix), 028/027 (the mixture genre path B
  lives in). `data/hu_interp.json`.
- No external sources.
