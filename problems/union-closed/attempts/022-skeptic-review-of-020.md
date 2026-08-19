# 022 — Skeptic review of 020 (cross-family), plus the Gap-2 survival check

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-19
- **Mode:** informed
- **Type:** skeptic review (default stance: refute), same-session as 020 —
  scope caveat below
- **Tools:** `explore/uc_gap1_skeptic_gemini.py` (drafted end-to-end by
  gemini-3.7-flash from the committed spec
  `explore/swarm022_skeptic_brief.md` with no repo code visible to it;
  audited line-by-line by the director; run unmodified; output in
  `data/gap1deep_skeptic_gemini_out.txt`). Frozen exact witness inputs:
  `data/gap1deep_witnesses.json` (12 instances, rational weights,
  deterministic regeneration procedure recorded in 020's engine).
  Gap-2 check: 009's `explore/uc_mitax.py` `cr_eval` on 014/016-style
  Sinkhorn couplings. Tensor certificate: `data/gap1deep_tensor_witness.json`
  + the 013 exact-census path (aggregate only).
- **Sources:** none.

## Approach

020's certificates rest on (a) the 013 atanh-series log₂ enclosure kit,
(b) director-written census and interval code (`exact_mm_abs`), and
(c) a new interval step (dz/dρ over the Plackett-root enclosure). The
verification standard demands a different implementation, and
`docs/SWARM.md` makes the family that drafts the re-implementation part
of the independence claim. So the skeptic verifier here was **drafted by
the other model family** (gemini-3.7-flash) from a mathematical spec
alone: its own certified log₂ (binary digit extraction with directed
200-bit dyadic rounding on separate lower/upper tracks — a different
algorithm, not a re-derivation of the atanh series), its own exact
bisection for the Plackett root, its own census construction and
normalization convention (unnormalized mass-weighted SUM where 020
reports the mass-weighted MEAN — deliberately not harmonized, so
agreement is a nontrivial cross-check of both conventions), and built-in
self-tests (product-measure anchor with exact OR ≡ t, log₂ anchors,
exact Plackett root x = y = 1/2, t = 9 ⇒ z = 3/8).

**Scope caveat, stated up front:** this review was run in the same
session as 020 by the same director. The *implementation* independence
is real (different model family drafted the program from spec; different
algorithms; zero shared code — audited), but the *reviewer* independence
of a fresh session is not achieved. The record is labelled accordingly
and a fresh-session skeptic remains queued; what this pass rules out is
shared-code and shared-algorithm error, the failure mode 007 documented.

## What was done

### 1. Audit of the drafted verifier

Line-by-line, before any run. Soundness points checked by hand: the
two-track digit extraction's monotonicity argument (floor-rounding can
only delay digit emission ⇒ valid lower bound; ceil-rounding plus the
2^-D stream truncation ⇒ valid upper bound); exact bisection bracket
validity for t > 1 (F(xy) < 0 < F(min(x,y)) verified before iterating,
loud failure otherwise); the dz/dρ interval (numerator decreasing in z
below min(x,y); denominator linear with negative slope, positivity
checked, loud failure otherwise); four-product interval multiplication
with the |·| straddle case; strict z_hi < 1/2 for the R₊ legs; strict
exact marginal comparison against 38271/100000. No modifications were
needed; the program ran as drafted. (The one worker-side wrinkle: the
first draft hit the output-token ceiling and was regenerated at a higher
cap; both prompts are hash-recorded in the swarm metas.)

### 2. Re-certification of all 12 frozen instances

All three self-tests passed, then (widths below print precision
everywhere; `mixed_zero_tables = 0` on every instance — the degeneracy
dichotomy holds row-by-row in the independent census too):

| instance | 020's value | 022 (gemini kit) | verdict |
|---|---|---|---|
| mmabs_kill_0, t=3/2 | num −3.403811e-3 | −3.4038105394e-3 | CONFIRMED |
| mmabs_kill_0, t=6/5 | num −2.922542e-3 | −2.9225421258e-3 | CONFIRMED |
| mmabs_kill_1, t=3/2 | num −1.001348e-3 | −1.0013481226e-3 | CONFIRMED |
| mmabs_kill_1, t=6/5 | num −2.143631e-3 | −2.1436307016e-3 | CONFIRMED |
| mmabs_kill_2, t=3/2 | num −1.784190e-3 | −1.7841902271e-3 | CONFIRMED |
| mmabs_kill_2, t=6/5 | num −2.884705e-3 | −2.8847050517e-3 | CONFIRMED |
| window_kill, t=21/20 | A −5.256144e-8 (mean) | sum −1.576843e-7 = mean × 2.99997… wtot | CONFIRMED |
| window_kill, t=16/15 | A −6.242487e-8 | sum −1.872749e-7 | CONFIRMED |
| window_kill, t=27/25 | A −5.508935e-8 | sum −1.652679e-7 | CONFIRMED |
| rplus_nogo_0, t=7/5 | num −2.0682e-3, R₊ | −2.0682097e-3, all_rplus **True** | CONFIRMED |
| rplus_nogo_1, t=7/5 | num −1.7775e-3, R₊ | −1.7774765e-3, all_rplus **True** | CONFIRMED |
| rplus_nogo_2, t=7/5 | num −3.6933e-4, R₊ | −3.6932653e-4, all_rplus **True** | CONFIRMED |

Every kill sign, every in-regime flag (max marginals 0.25747–0.38267,
all strictly below 0.38271 in exact arithmetic), and every R₊ leg
reproduces. The window-kill values agree exactly under the
mean-vs-sum conversion (independent wtot ≈ 2.99997 defined-mass units at
n = 4, matching 020's per-i defined masses).

### 3. Hand re-derivations (proof-shaped steps of 020)

- The weight identity dh₂(z_ρ(x,y))/dλ = h₂′(z)·(dz/dρ)·ρ·ln 2 with
  dz/dρ = (x−z)(y−z)/((1−x−y+2z) + ρ(x+y−2z)): re-derived by implicit
  differentiation of the Plackett quadratic. The global constant t·ln 2
  is strictly positive, multiplies numerator and denominator of MM_abs
  alike, and so cannot affect any sign — 020's dropping of it is sound.
- The R₊ reduction (020 part F's logic): any weight w ≥ 0 with w = σ_λ
  wherever σ_λ > 0 satisfies E[w·dev] = E[σ_λ·dev] = E[|σ_λ|·dev] on any
  measure whose nondegenerate histories all have z_{2^λ} < 1/2. So a
  certified MM_abs-numerator < 0 with certified all-R₊ histories rules
  out the whole compatibility class. Re-checked; sound. Boundary note
  kept from 020: this defines compatibility as "equals σ on R₊"; weights
  allowed to differ from σ somewhere on R₊ are outside the no-go.
- The first-order kernel model (design tool only) is **not load-bearing**
  for any certificate — confirmed by reading `exact_mm_abs`: no model
  quantity enters the certified path.

### 4. Gap-2 survival check (020 lead 4 = queue 1b')

Both 020 kill geometries were run through 009's chain-rule accounting
(`cr_eval`, exact-coupling float path, Sinkhorn marginal deviation
≤ 2.3e-15) with a λ sweep:

    window-kill μ (n = 4):  sup_λ CR = +0.7985 (λ = 0.05); gain +0.8058;
                            second tax +0.0073; CR(0.1) = +0.7965
    MM_abs-kill μ (n = 5):  sup_λ CR = +0.4003 (λ = 0.05); gain +0.7525;
                            second tax +0.3522; CR(0.1) = +0.3992

**Gap 2's candidate (TAX at p: CR > 0) survives both witnesses with wide
margins.** The geometries that annihilate every per-history odds-ratio
condition leave the chain-rule assembly value strongly positive — direct
evidence that Gap 2's mechanism is genuinely different from the dead
bridge family, and the natural next target for proof effort.

### 5. Tensor certificate at n = 8 (020 lead 2)

The m = 2 block tensor of the window-kill measure (225 atoms), frozen at
t = 25/24 (λ ≈ 0.0589 < λ_win(8) = 0.9694), aggregate certified by the
exact census path (log₂-enclosure cache over distinct OR values):
result recorded in `data/gap1deep_tensor_cert.txt` — see Outcome.

## Outcome

**CONFIRMED — 020's headline stands in full.** All 12 certificates
reproduce digit-for-digit under a cross-family, zero-shared-code,
different-algorithm re-implementation; the three R₊ no-go legs certify
independently; the two hand-derivation steps are sound; and the window
kill extends to a certified n = 8 in-window instance (tensor
certificate). Additionally the Gap-2 candidate survives both kill
geometries (float, margins ~0.4–0.8 — far above any enclosure concern).

Corrections to 020: **none found.** One reporting-level note: 020's
part-E/F tables quote the aggregate as a mass-weighted mean while the
spec this review issued asks for the raw sum; both are correct and the
conversion factor is the defined mass (≈ 3 at n = 4), but future records
should name the convention next to each number.

Not claimed: fresh-session reviewer independence (see scope caveat);
certification of the λ-integrated variant's negativity (still
float-level); anything about Gap 2 beyond the two-witness float check.

## Why it failed / what survived

Nothing failed to reproduce. The review's value is (a) closing the
shared-code/shared-algorithm risk on a six-certificate batch produced in
one day, (b) establishing the cross-family drafted-verifier pattern as a
practical skeptic tool (spec → independent program → audit → run:
~$0.04 of worker tokens, one director audit hour), and (c) the Gap-2
survival data point, which converts "Gap 2 is what's left" from
elimination reasoning into positive evidence.

## Leads generated

1. **Fresh-session skeptic** (unchanged from 020 lead 1, now narrower):
   what remains untested is reviewer-level independence — a new session
   should re-run `uc_gap1_skeptic_gemini.py` from the committed inputs
   and spot-check the audit, rather than rebuild from scratch.
2. **Gap 2 becomes the queue's proof-effort target**: 009/011's CR > 0
   candidate survived its most dangerous new adversaries. First
   falsifiable step: adversarial search for CR < 0 seeded at the 020
   geometries (mutate within the marginal cap, minimize sup_λ CR), and
   an exact certification path for CR (needs certified h(z) sums — the
   022 verifier's log₂ kit extends directly).
3. Certify the λ-integrated variant's negativity on the 020 witnesses
   (cheap: same census, integrand at a few rational tilts, then a
   monotonicity argument or a finer partition for the integral bound).

## References

- This repo: attempt 020 (under review), 009/011 (Gap 2), 013/017/019
  (exact-verification precedents), `docs/SWARM.md` (cross-family rule).
- Swarm: `explore/swarm022_skeptic_brief.md`; worker gemini-3.7-flash,
  effort medium, single draft used unmodified after audit.
