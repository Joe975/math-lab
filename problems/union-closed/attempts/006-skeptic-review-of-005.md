# 006 — Skeptic review of 005 (odds-ratio control refuted): adversarial verification

- **Problem:** union-closed, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-07-30
- **Mode:** informed
- **Type:** adversarial verification of `005-odds-ratio-control-refuted.md`
  (default stance: refute). Every proposition re-derived by hand from the
  003/004 definitions; every number re-computed with an independent
  implementation (`explore/uc_or_skeptic.py`, no imports from
  `uc_odds_ratio.py`; column-first Sinkhorn with both-margin residuals vs
  005's row-first with row residual; potential symmetry *checked* rather than
  assumed; agreement is meaningful because the scaling coupling of a fixed
  positive kernel is unique). The memory-cited Karlin–Rinott step was checked
  against the literature this cycle, and independently re-derived.
- **Outcome in one line:** 005 holds. All six propositions and the crash
  family survive re-derivation and independent recomputation (including at
  parameter values 005 never ran); the corrections found are small — one
  equality-characterization edge case (λ = 0), two reporting slips in the
  part-B prose, two wrong sentences in the tool's docstring — none
  load-bearing. The headline (Gap 1 refuted as stated, both directions; route
  LIVE; gap restated as averaged control) stands, and the restated gap
  survives its own first falsifiable test, run here.
- **Tools:** `explore/uc_or_skeptic.py` (parts S1–S8; stdlib only;
  deterministic, fixed seeds; runtime ~11 s; checkpoints
  `data/006_partS[1-8].json`, log `data/006_full_run.log`). Command
  reproducing every 006 number below:
  `python problems/union-closed/explore/uc_or_skeptic.py`
- **Sources:** 003/004/005 (this repo). Fallat–Lauritzen–Sadeghi–Uhler–
  Wermuth–Zwiernik, *Total positivity in Markov structures*,
  arXiv:1510.01290 (ar5iv HTML fetched 2026-07-30 [T]) — its Prop 3.4(ii)
  states MTP₂ closure under marginalization and attributes it to
  Karlin–Rinott 1980, Prop 3.2; used to check 005's memory citation. The
  K–R original was not fetched (paywalled); the binary-cube case 005 needs is
  re-proved from the Ahlswede–Daykin four-functions theorem below, so nothing
  rests on the transcription.

Notation as in 005: `π_λ(A,B) ∝ u(A)u(B)2^{λ|A∩B|}`, histories `a = A_{<i}`,
`b = B_{<i}`, `m = n − i`, `OR_i(a,b)` the conditional odds ratio,
`W = ⊗_{j>i} [[1,1],[1,2^λ]]`, `F_α(x) = u(a,α,x)`, `G_β(y) = u(b,β,y)`.

## Claims attacked

1. The reduction `OR_i(a,b) = 2^λ·R`, R the cross-ratio of `⟨F,WG⟩`
   inner products — plus the two lemmas it silently uses (equal potentials on
   both sides; "any positive u is its own marginal's potential").
2. Prop 2 (diagonal Cauchy–Schwarz: `OR_i(a,a) ≥ 2^λ`, equality iff
   proportional slices), including W ⪰ 0 and the identification `G_β = F_β`.
3. Prop 3 (4-atom crash family: `OR = 2^{λ(3−n)}` at a positive-mass history,
   marginals ≤ 0.37), the potential-cancellation argument, the exponent
   arithmetic, and the "where the crash bites" inequality chain.
4. Prop 4 (sharp universal range `[2^{λ(1−m)}, 2^{λ(1+m)}]`, vertex
   extremality, attainment at both ends).
5. Prop 5 (potential-level MTP₂ ⟹ `OR ≥ 2^λ` everywhere) — with the
   Karlin–Rinott marginalization step cited from memory in 005.
6. The Sinkhorn-breaks-MTP₂ instance (`005_partF.json`,
   `worst_case_instance`, claimed violating ratio 0.9922).
7. Prop 6 (`log₂OR = λ + O(δ²)` near product μ) and the δ²-scaling numerics;
   part-E census figures; part-B prose vs its own checkpoint.
8. The headline framing itself: whether the refuted statement is really
   003's Gap 1 verbatim, whether both directions are refuted *by admissible
   μ* in the record-relevant marginal regime, and whether the restated gap
   (averaged control) is already dead on 005's own hard instances.

## Refutations found

Nothing load-bearing. Four genuine but small defects:

### R1. Prop 2's equality characterization fails at λ = 0
Prop 2 is stated for λ ≥ 0 with "equality iff the two future slices are
proportional". At λ = 0, `W` is the all-ones matrix (each factor
`[[1,1],[1,1]]` has det 0), `⟨F,WG⟩ = (ΣF)(ΣG)`, so `R = 1` and `OR = 1 = 2⁰`
**for every μ** — equality with generically non-proportional slices.
Verified: a random full-support potential at λ = 0 has max |OR − 1| =
1.1e-15 over all histories with slices nowhere near proportional (S2). The
inequality half of Prop 2 is untouched (and at λ = 0 is trivially tight);
the "iff" needs the hypothesis λ > 0, where each W-factor has det
2^λ − 1 > 0 and W is strictly PD, making Cauchy–Schwarz equality ⟺
proportionality genuine. Cosmetic, nothing downstream uses the λ = 0 case —
but the proposition as printed is false at its boundary point.

### R2. Part-B prose misreports its own checkpoint (two slips)
- "minimum of OR/2^λ over **all** diagonal histories = 1.0000000011": per
  `005_partB.json` that number is `global_min_diag_or_over_2lam_i_lt_n`; the
  all-histories minimum is 0.9999999999999992 (the i = n rows, exactly 2^λ up
  to float, several of them a few ulps *below* 1 — so "≥ 1 as proved" is
  float-false for the quantity the sentence names). The sentence's own
  parenthetical shows the i<n reading was intended; the number is real,
  the label on it is wrong.
- "median per-trial minimum excess 5.4e-4": recomputed from
  `005_partB.json`, the median of `min_diag_or_over_2lam_i_lt_n − 1` over
  the 60 trials is **5.2e-4** (max diag ratio 3.84 ✓, 0/60 off-diagonal
  ORs below 1 ✓). Reporting slip only.

### R3. Two wrong sentences in the tool's own docstring
`uc_odds_ratio.py` header, part F line: "if the potential u is
log-supermodular then OR >= 1 at every history" — the code, checkpoint, and
record all establish the stronger `OR ≥ 2^λ`. Part G line: "overall
deviation like delta" — the code's own output (and `005_partG.json`, where
`max_overall_dev ≡ max_diag_excess`) shows δ² for both. The *record* states
both points correctly; only the docstring is wrong. Flagged because the
docstring is what the next agent reads first.

### R4. "About 83% of history mass has OR > 2^λ" is the λ = 0.5 figure
Part-E checkpoint: 83.1% at λ = 0.5 but 81.1% at λ = 1.5. Rounding-level
overstatement in a sentence covering both runs; the qualitative claim
(essentially all deviation mass sits above 2^λ, none below 1) is exactly
right (S5 reproduces every field to 1.1e-14).

## Claims that survive (and what was done to break them)

### The reduction (claim 1) — re-derived by hand, then stress-tested off 005's grid
Hand derivation from 003's definitions: split `A = (a, α, x)`,
`B = (b, β, y)`; `|A∩B| = |a∧b| + αβ + ⟨x,y⟩` (coordinates < i, = i, > i);
so `P(α,β|a,b) ∝ 2^{λ|a∧b|}·2^{λαβ}·⟨F_α, WG_β⟩`, the prefix factor and the
normalization cancel in the odds ratio, and the 2^{λαβ} factor contributes
2^λ from the (1,1) cell alone: `OR = 2^λ·R`. Two silent lemmas checked:

- *Equal potentials.* Sinkhorn gives `π = r(A)c(B)K`; by symmetry of K and
  equality of the marginals, `π^T` is also a scaling coupling with the same
  marginals, so uniqueness forces `π = π^T`, hence `r = t·c` and
  `u = √(rc)` works. Verified computationally: max relative deviation of
  r/c from constant ≤ 5.7e-12 across every fit in S3/S5/S6 (`sym_dev`).
- *Any positive u is its own marginal's potential* — immediate from
  uniqueness of the scaling coupling; this is what makes the four slices
  free and is used by S1/S2/S4/S6b to build test instances without running
  Sinkhorn.

Numerical kill attempt (S1): 15,336 nondegenerate histories over random
potentials at n = 4, 5, 6, λ ∈ {0.6, 1.7}, dense AND sparse random supports
(005's sweeps were full-support only): census OR vs `2^λ·R` computed
directly from the slices — max relative error 4.7e-15. CONFIRMED.

### Prop 2 (claim 2) — inequality airtight; see R1 for the λ = 0 edge
W ⪰ 0 re-derived: each factor `[[1,1],[1,2^λ]]` has positive trace and det
`2^λ − 1 ≥ 0`; tensor products of PSD matrices are PSD; strict PD for λ > 0.
On the diagonal `a = b` the equal-potentials lemma gives `G_β = F_β`
legitimately, and `R = ‖F₁‖²_W ‖F₀‖²_W / ⟨F₁,F₀⟩²_W ≥ 1` is Cauchy–Schwarz
(the cross term is > 0 at a nondegenerate history, so the division is
legal). Kill attempts: 60 random potentials at n = 6 (005 used n = 5),
λ ∈ {0.25, 1, 3}, including sparse supports: min diag OR/2^λ = 1 − 7e-16
(float-exact ≥ 1; sparse supports produce genuine equality rows —
point-mass slices are proportional, consistent with the characterization,
not violations of it). A constructed proportional-slice instance
(F₁ = 3F₀ at one prefix, λ = 1.3) hits OR = 2^λ to 0.0e+00 (S2).
CONFIRMED (with R1's λ > 0 amendment to the "iff").

### Prop 3, the crash family (claim 3) — the record's central claim; it is exact
Hand re-check of everything:

- Marginals: element 1 ∈ S₃,S₄ only: 0.04 + 0.30 = 0.34; element 2 ∈ S₁,S₃:
  0.33 + 0.04 = 0.37; elements ≥ 3 ∈ S₂,S₃: 0.37. All ≤ 0.37 < 0.38271 ✓.
  H(μ) = 1.7625 bits ✓. Masses sum to 1.00 ✓.
- Cell classification at i = 2, (a,b) = (A₁=0, B₁=1): A₁=0 selects {S₁,S₂},
  B₁=1 selects {S₃,S₄}; coordinate 2 then splits S₁ (has 2) from S₂ and S₃
  (has 2) from S₄ — one atom pair per cell: (1,1)=(S₁,S₃), (0,0)=(S₂,S₄),
  (1,0)=(S₁,S₄), (0,1)=(S₂,S₃) ✓.
- Cancellation: each uᵢ appears once in numerator and once in denominator
  (u₁u₂u₃u₄ both) ✓ — the OR is kernel-only for ANY potentials, so Sinkhorn
  convergence is irrelevant to the value. Exponent:
  |S₁∩S₃| = |{2}| = 1, |S₂∩S₄| = 0, |S₁∩S₄| = 0, |S₂∩S₃| = n−2, so
  `OR = 2^{λ(1+0−0−(n−2))} = 2^{λ(3−n)}` ✓ (n ≥ 4 keeps the atoms distinct).
- "Where the crash bites": for fixed zero-margins x = y = 1−p, z(ρ) is
  increasing in ρ with z(1) = xy, so ρ < 1 ⟹ z < (1−p)², m = 1−z > q;
  q ≥ 1−p ⟺ p² − 3p + 1 ≤ 0 ⟺ p ≥ ψ (equality at ψ: (1−ψ)² = ψ), and
  1−p ≥ 1/2; h decreasing on [1/2,1] gives h(m) < h(q) ≤ h(p) ✓.

Kill attempts (S3): own Sinkhorn at n ∈ {4,5,8,11} and λ ∈ {0.7, 1.0, 1.9}
— 005 only ever ran λ = 1; the identity is linear in λ and the engine
confirms `log₂OR = λ(3−n)` to 1.8e-15 at all twelve points, marginal-fit
residuals ≤ 1e-13 on both margins. Crash-history mass at λ = 1: 0.1661,
0.1588, 0.1441, 0.1376 (n = 4, 5, 8, 11) — 005's "0.13–0.16" ✓ (mass is
λ-dependent: ~0.09 at λ = 1.9, still Θ(1)). Full-support ε-mix re-run at
λ = 1.0 (matches 005: −4.6204/−4.9587/−4.9958) and additionally at λ = 1.5
(−6.5553/−7.3859/−7.4883 → −7.5), so the robustness statement is not
λ = 1-specific either. CONFIRMED — and this alone refutes Gap 1 as stated,
since 2^{λ(3−n)} < 1 for every λ > 0, n ≥ 4, with marginals in-class.

### Prop 4 (claim 4) — vertex argument re-derived; both ends re-attained
Re-derivation: at a nondegenerate history the vectors WG_β are entrywise
positive, so R is, in F₁ (others fixed), a ratio of two positive linear
forms; `(Σaᵢtᵢ)/(Σbᵢtᵢ) ≤ maxᵢ aᵢ/bᵢ` (mediant inequality) puts the max at
a point mass, and four sequential replacements — each not decreasing R —
land on all-point-mass configurations, which are realizable as 4-atom
supports. There `log₂R = λ⟨x₁−x₀, y₁−y₀⟩ ∈ λ[−m, m]` since entries of the
difference vectors lie in {−1,0,1} ✓. Kill attempts (S4): 30 random
potentials at n = 6, λ = 0.8, dense and sparse: no violation (min slack
−1.2e-15, i.e. the i = n boundary where the range collapses to {2^λ}).
Attainment re-checked at n = 5 AND n = 7 (005: n = 5 only): crash family
hits `λ(1−m)` and the mirror family `λ(1+m)` to 1e-15. CONFIRMED.

### Prop 5, MTP₂ (claim 5) — memory citation checked, step re-proved, bound holds
- *The citation.* arXiv:1510.01290 (Fallat et al.), Prop 3.4(ii): "If X has
  an MTP₂ distribution, then for every A⊆V, the marginal distribution X_A of
  X is MTP₂", attributed there to Karlin–Rinott 1980, Prop 3.2 [T]. 005's
  memory-cited statement is the standard one.
- *Independent proof for the case 005 needs* (binary cube), so nothing rests
  on a transcription: it suffices to sum out one coordinate. Let f be MTP₂
  on {0,1}^{k+1}, g(x) = Σ_t f(x,t). Fix x, x′ and set a(s) = f(x∨x′, s),
  b(s) = f(x∧x′, s), c(s) = f(x, s), d(s) = f(x′, s); MTP₂ of f applied to
  the pairs (x,s), (x′,t) gives a(s∨t)b(s∧t) ≥ c(s)d(t) for all s, t ∈
  {0,1}, and the Ahlswede–Daykin four-functions theorem (on the two-element
  lattice) yields (Σa)(Σb) ≥ (Σc)(Σd), i.e. g(x∨x′)g(x∧x′) ≥ g(x)g(x′).
  Iterate. Numerically corroborated (S6a): 120 random ferromagnetic f at
  n = 5–7, random coordinate subsets summed out — worst log-supermodularity
  ratio of the marginal 1.0065 ≥ 1, zero violations.
- *The assembly of Prop 5.* Slices of MTP₂ are MTP₂ (fix coordinates);
  products over disjoint variable sets add log-supermodularity; each
  W-factor 2^{λx_jy_j} is a pairwise interaction with coefficient λ ≥ 0,
  hence MTP₂ — so Ψ = F_α(x)G_β(y)W(x,y) is MTP₂ on {0,1}^{2m+2} and
  marginalizing (x,y) leaves Φ(α,β) TP₂: R ≥ 1, OR ≥ 2^λ ✓. Note the
  hypothesis is genuinely potential-level and full-support, as 005 says.
- *Kill attempt* (S6b): 20 random log-supermodular potentials at n = 5
  (005: n = 4), λ ∈ {0.5, 2.0}: min OR/2^λ over all histories = 1 − 8e-16.

CONFIRMED.

### The Sinkhorn-breaks-MTP₂ instance (claim 6) — real, and correctly interpreted
Refit of `worst_case_instance` from `005_partF.json` with this review's
Sinkhorn (S6c): the stored μ **is** MTP₂ (worst log-supermodularity ratio
1.00211 ≥ 1 — this mattered: a non-MTP₂ μ would have made the instance
vacuous), fit residual 3.2e-14 on both margins, and the fitted potential's
worst log-supermodularity ratio is **0.992152021** — matching 005's claimed
0.9921520206 to 1.8e-15 from independent code. All conditional ORs of the
instance still ≥ 2^λ (min ratio 1 − 3e-16), as 005 reports. CONFIRMED.

### Prop 6 and the δ² numerics (claim 7) — survives, with a caveat worth recording
First-order cancellation re-derived: with u₀ product, the future slices at
any prefix are proportional, `F⁰_α = c_α h_A`, `G⁰_β = d_β h_B`; the O(δ)
term of `log⟨F_α, WG_β⟩` is `T^A_α + T^B_β` with `T^A` independent of β and
`T^B` independent of α, and the alternating sum (1,1)+(0,0)−(1,0)−(0,1)
annihilates any such additive form — so `log R = O(δ²)` at every history,
uniformly at fixed n and full support. The one unproved step is smoothness
of the potential in μ, exactly as 005 labels it. Kill attempts (S7): 005's
family reproduced (max dev 1.53e-5 … 3.02e-3, log-log slopes 1.97 → 1.81);
a SECOND family with a >1/2 component, `(1−δ)Bern(0.35)⁵ + δBern(0.70)⁵`,
initially looked like a refutation — slopes 1.75 → 1.06 on 005's δ-range —
but extending down to δ = 1.56e-4 gives slopes 1.992, 1.983, 1.966, …: the
δ² law holds with a much larger δ³ coefficient and a correspondingly narrow
asymptotic window. CONFIRMED, with the caveat: **the O(δ²) constant and the
onset of the asymptotic regime are family-dependent**; 005's lead 3
(perturbative assembly with budget c·δ²) must treat c as worst-case over
directions, and 005's own part-G family is among the gentlest.

### Part-E census (claim 7 continued) — reproduced to 1.1e-14
Independent engine, same mixture: every field of `005_partE.json` (min/max
log₂OR, mass-weighted means, mass above 2^λ / below 1, per coordinate, both
λ) agrees to 1.1e-14 (S5). The record's qualitative reading — deviation
budget sits above 2^λ, zero mass below 1, means decreasing to exactly λ at
i = n — is faithful to the data (R4's 83%-vs-81% aside).

### The headline framing (claim 8) — checked against 003's text and the data
- The refuted statement is 003 §3 item 1 / Gap 1 verbatim ("λ ≥ 0 forces
  all conditional odds ratios into [1, 2^λ]", general μ, SPECULATION label).
  No goalpost was moved: 005 refutes exactly the recorded sentence.
- Both directions are witnessed by admissible μ inside the record-relevant
  marginal regime: lower by the crash family (marginals ≤ 0.37 < 0.38271);
  upper by real measures, not just abstract potentials — the part-E mixture
  (marginals 0.225) has 80%+ of history mass strictly above 2^λ. So even a
  weakened Gap 1 restricted to marginals < 0.38271 is dead in both
  directions.
- "Route stays LIVE" touches nothing this review found: the licensing lemma,
  separations, no-go evasion, and 004's 0.431496 ceiling are untouched by
  005's content, and 005 claims no bound.
- **The restated gap survives its own first falsifiable test, run here**
  (005 lead 1 names it; S8 executes it): on the crash family the
  mass-weighted mean `M_i = E_π[log₂OR_i]` over defined histories is +1.85,
  +2.09 (n = 5) and +2.28, +3.54 (n = 8) at λ = 1 — the 2^{λ(n−1)}-diagonal
  histories outweigh the crash history; on the mirror family M₁ = 1.000 =
  λ exactly (the conjectured bound is attained, so `M_i ≥ λ` is sharp if
  true) and M₂ = λ(n−1). No `M_i < λ` found. The averaged-control
  conjecture remains SPECULATION — this is two families, not a proof — but
  it is not already dead on the instances that killed the pointwise claim.
  (Bookkeeping note: on 4-atom supports, histories at i ≥ 3 are all
  degenerate — after coordinates 1–2 the atom is determined — so `M_i` there
  is defined on zero mass; any averaged-control statement needs to say what
  it averages when the conditional table degenerates.)

## Verdict

| # | Claim | Verdict |
|---|-------|---------|
| 1 | Reduction `OR = 2^λ·R` + equal-potential and own-potential lemmas | **CONFIRMED** (hand re-derivation; 15,336 histories incl. sparse supports, err ≤ 4.7e-15) |
| 2 | Prop 2 diagonal Cauchy–Schwarz `OR ≥ 2^λ` | **CONFIRMED**; equality "iff" **CORRECTED** to λ > 0 (false at λ = 0, where OR ≡ 1 for every μ) |
| 3 | Prop 3 crash family `OR = 2^{λ(3−n)}`, marginals ≤ 0.37, mass Θ(1) | **CONFIRMED** (exact by hand; engine-exact also at λ = 0.7, 1.9 which 005 never ran; ε-mix robustness also at λ = 1.5) |
| 4 | Prop 4 sharp range `[2^{λ(1−m)}, 2^{λ(1+m)}]` + attainment | **CONFIRMED** (mediant/vertex argument re-derived; both ends re-attained at n = 5 and 7) |
| 5 | Prop 5 MTP₂ ⟹ `OR ≥ 2^λ`; Karlin–Rinott step | **CONFIRMED** (citation checked via arXiv:1510.01290 Prop 3.4(ii) [T]; binary case re-proved from Ahlswede–Daykin; 120-marginal direct test; n = 5 OR sweep) |
| 6 | Sinkhorn-breaks-MTP₂ instance, ratio 0.9922 | **CONFIRMED** (μ verified MTP₂; independent refit reproduces 0.992152021 to 1.8e-15) |
| 7 | Prop 6 `λ + O(δ²)` + numerics; part-E census | **CONFIRMED** (second, harsher family also δ² once δ small enough; caveat: family-dependent asymptotic onset). Part-B/E prose: **CORRECTED** (R2, R4 reporting slips; docstring bugs R3) |
| 8 | Headline: Gap 1 refuted as stated, both directions, in-regime; route LIVE; gap restated as averaged control | **CONFIRMED** (matches 003's text; both directions witnessed by admissible in-regime μ; restated gap passes its own first falsifiable test, run here — S8) |

**Net assessment:** 005 is the cleanest record in this problem's chain so
far: the central refutation is an exact, potential-free identity that
survives every attack tried here, the partial positive results are real
theorems (with one boundary-case wording fix), and its self-labeling
(SPECULATION/EVIDENCE markers, "not claimed" list) accurately matches what
the data supports. Gap 1 (`plackett-odds-ratio-control`) should be
considered REFUTED as stated, and the route's gap list correctly reads
`averaged-odds-ratio-control` in its place.

## Residual risk

- **Prop 6's smoothness step** (differentiability of the Sinkhorn potential
  in μ) remains taken-as-standard; 005 labels it and this review did not
  close it. All numerics are consistent with it.
- **The averaged-control conjecture (`M_i ≥ λ`)** was tested here on exactly
  two structured families plus 005's realistic instance; S8 is EVIDENCE for
  the restated gap being well-posed, not support for its truth. The
  degenerate-history bookkeeping issue noted above must be resolved before
  anyone attempts a proof.
- **Karlin–Rinott** was verified against a secondary source and re-proved
  for the binary cube; the original 1980 paper itself was still not read.
  For the binary-cube use in Prop 5 this risk is now closed by the
  Ahlswede–Daykin derivation; anyone extending Prop 5 beyond product binary
  lattices should fetch the original.
- **Shared-reduction risk.** Both engines represent atoms as bitmasks and
  couplings as atom-pair dictionaries; a conceptual error in the *census
  definition* itself (what "conditional odds ratio given histories" means)
  would evade both. Mitigated by S1: the census is checked against the
  independently-derived slice/inner-product formula at every history, which
  is a different mathematical route to the same object — but not by a third
  formalization.

## References

- 005-odds-ratio-control-refuted.md; 004-skeptic-review-of-003.md;
  003-dependent-couplings.md (this repo).
- `explore/uc_or_skeptic.py`; checkpoints `data/006_partS[1-8].json`,
  `data/006_full_run.log`; reviewed artifacts `explore/uc_odds_ratio.py`,
  `data/005_part[A-G].json`.
- Fallat, Lauritzen, Sadeghi, Uhler, Wermuth, Zwiernik, arXiv:1510.01290,
  Prop 3.4 [T: ar5iv HTML, fetched 2026-07-30] — MTP₂ closure under
  marginalization, attributing it to S. Karlin, Y. Rinott, J. Multivariate
  Anal. 10 (1980) 467–498, Prop 3.2 (original not fetched).
- Ahlswede–Daykin four-functions theorem (1978) — used here in its
  elementary two-element-lattice form to re-prove the marginalization step.
