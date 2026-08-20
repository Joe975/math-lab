# 040 — Anneal vs rollout-HU at n = 6/7: zero kills, and the adversary converges onto a second equality family that pins every order rule's ceiling

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-20
- **Mode:** informed
- **Type:** adversarial attack + one elementary identity (038 lead 3:
  the n ≥ 6 story needed a stronger adversary than 035's pattern
  descent, which stalled 40× above the n ≤ 5 floors).
- **Tools:** `explore/uc_hu_roll_anneal.py` (new; simulated annealing
  over measures — weight kicks, atom add by random mask / union /
  complement, atom drop, pairwise transfer — with δ∅ regime
  projection, then a 035-style pattern-descent polish from the
  annealed best; deterministic, seeds 947000+; checkpoint
  `data/hu_roll_anneal.json`); `explore/uc_hu_roll_anneal_skeptic.py`
  (new; every endpoint re-scored through the 037-skeptic independent
  stack; exit 0). Reproduce: run the two in that order.
- **Sources:** none.

## Approach

038's honest caveat was that its n = 6/7 "zero violations" came from a
visibly stalling monotone descent. This attempt upgrades the adversary
in the two ways 038 lead 3 asked for: annealing (accepts worsening
moves, escapes the descent's basins) and structural moves (the support
itself is searched, not just the weights). Hostile seeding embeds the
known-hard geometry — the two 035 cap-0.49 kill measures tensored with
independent near-cap Bernoulli coordinates (top-14 atoms kept) — plus
fresh random supports.

## What was done

**A. Anneal + polish campaigns** (4000 steps at n = 6, 1200 at n = 7;
6 runs per block; polish = 035's pattern descent from the annealed
best):

    n=6 cap 0.49 : global floor +0.000488   0 violations   [extremal +0.000289]
    n=6 cap 0.497: global floor +0.000180   0 violations   [extremal +0.000026]
    n=7 cap 0.49 : global floor +0.000488   0 violations   [extremal +0.000289]

  The anneal is a genuinely stronger adversary: it lands ~60× below
  038's stalled descent floors (0.0305 → 0.000488 at n = 6, 0.49).
  All endpoints reproduce under the independent skeptic stack, all
  in-regime.

**B. The attractor is identified, and it is not the product.** The
sharpest endpoints at both n = 6 and n = 7 (identical floors to
3e-12) are, up to vanishing atoms, the **two-point diagonal measure**
μ_p = (1−p)·δ_∅ ⊕ p·δ_S at p = cap − 0.003 (the projection
boundary), S the full coordinate set: every S-coordinate has marginal
exactly p, all perfectly correlated. The measured floor equals the
family value to 2e-13.

**C. The DIAG identity (elementary, order-free).** For μ_p and ANY
revelation order: the first revealed S-coordinate has x = y = 1−p, so
z = min(max(1/2, 1−2p), 1−p) = max(1/2, 1−2p) (as 1−p > 1/2), and
conditioning on it makes everything deterministic — every later cell
has z ∈ {0,1}. Hence, for every n and every order,

    CR_HU(μ_p) = h(max(1/2, 1−2p)) − h(p),    H(μ_p) = h(p),

so CR/H = (1 − h(p))/h(p) for p ∈ [1/4, 1/2) and
(h(2p) − h(p))/h(p) for p < 1/4 — **exactly the corrected (HU-TAX)
product constant c\*(p̄) of 034, on both branches**. Verified
numerically at p ∈ {0.05 … 0.4999} × n ∈ {2,4,6} × sampled orders, all
to 1e-12; the identity is three lines and is recorded here as proved
(same-session; reviewer re-derivation joins the 036–038 batch).

**Consequences.** (i) (HU-TAX)'s constant has a SECOND equality
family besides products — two-point diagonals attain c\*(p̄)·H at
every fixed n (in the p → p̄ limit), so any proof must handle both,
and they are structurally opposite (independent vs perfectly
correlated). (ii) Since the identity is order-free, **no order rule —
including the best-order oracle — can guarantee more than
c\*(p̄)·H(μ) at marginal cap p̄**: 037 lead 3's ceiling question is
now pinned analytically, closing it. The order-engineering program has
exactly c\*(p̄) as its target, and rollout empirically saturates it
(the anneal converged onto this family and could not go below).
(iii) 037 part D's best-order floors (+0.000694 at 0.49) were above
the ceiling only because its descents never found the diagonal.

## Outcome

- **EVIDENCE: rollout-HU survives a genuinely stronger adversary at
  n = 6 (caps 0.49, 0.497) and n = 7 (0.49)** — zero violations, and
  this time the floors are not stall points: the adversary saturated
  an identified extremal family (value matched to 2e-13).
- **PROVED (elementary, same-session): the DIAG identity** — on
  two-point diagonal measures CR_HU = h(max(1/2, 1−2p)) − h(p) for
  every order and every n, so CR/H = c\*(p̄) there, pinning the
  guaranteed ratio of EVERY order rule at cap p̄ to at most c\*(p̄).
- **CLOSED: 037 lead 3** (best-order ceiling): the ceiling is
  c\*(p̄), by construction, not by search.
- **Not claimed:** sharpness of the 0.497 floor (+0.000180 sits above
  the diagonal value 0.000104 at p = 0.494 — that block's anneal did
  not fully reach the family); anything at n > 7; that roll *attains*
  c\*(p̄)·H for all μ (that is exactly the (HU-TAX, roll form)
  conjecture, SPECULATION as always).

## Why it failed / what survived

The attack failed to kill, and failed in the most informative way
available: it converged, from random supports at two different n, onto
a closed-form family that achieves the conjectured extremal ratio
exactly. Together with 034's product analysis, the extremal landscape
of (HU-TAX) now has two named equality families at opposite
correlation extremes, and the constant is provably unimprovable by any
order rule. What survives unproven is the whole middle: measures
between independence and full correlation, where every attack so far
lands strictly above c\*(p̄).

## Leads generated

1. **Interpolation family:** ρ-correlated mixtures between Bern(p̄)^n
   and the diagonal (e.g. (1−ρ)·product ⊕ ρ·diagonal, or Gaussian-
   copula-style). If CR_HU/H dips below c\*(p̄) anywhere along the
   path, (HU-TAX) dies; if it is monotone-above on the whole path,
   that is the first structured evidence the two equality families are
   the only ones. Cheap, falsifiable, closed-form endpoints.
2. **Proof shape suggestion:** both equality cases make every history
   cell hit its clamp boundary (product: z = 1/2 idle everywhere at
   p ≥ 1/4; diagonal: one live cell then determinism). A per-cell
   convexity/majorization argument between those boundary behaviors is
   the natural attack on (HU-TAX, roll form) — 031's ledger, but now
   with both extremal patterns known.
3. **The 0.497 gap:** close the anneal onto the diagonal at cap 0.497
   (seed it there directly) so every cap's floor is family-saturated;
   trivial run, tidies the evidence table.
4. Unchanged: 038's sandwich framing (best ≥ roll ≥ canon), now with
   the top of the sandwich capped by c\*(p̄).

## References

- This repo: 038 (the lead and the stalled-descent caveat), 037 (the
  rollout rule under attack), 035 (kill measures as embedded seeds,
  descent move set as polish), 034 (the corrected two-branch constant
  c\*(p̄) that DIAG realizes), 031/030 (HU, the ledger).
  `data/hu_roll_anneal.json`.
- No external sources.
