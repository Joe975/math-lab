# 019 — Skeptic review of 018 (surviving controls vs the ladder): adversarial verification

- **Problem:** union-closed, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-06
- **Mode:** informed
- **Type:** adversarial verification of `018-surviving-controls-vs-ladder.md`
  (default stance: refute). The decisive check is a **from-scratch exact
  re-certification of the headline kill** with a kit sharing no code and no
  algorithms with 018's path: own scaled-integer census written from 007 §1's
  prose, Plackett-root enclosure by exact-rational Newton with an exact
  sign-checked bracket (018 used blind dyadic bisection; my bracket endpoint
  signs are proved symbolically and re-checked exactly per call), and a
  dyadic interval-squaring log₂ enclosure (the idea of 017's kit,
  re-implemented — NOT 013's atanh-series kit, which 018 imported). The
  secant identity, the product-measure zero, and the degenerate-history
  convention one-liner were re-derived by hand (below) and the last two also
  engine-checked as exact identities.
- **Outcome in one line:** 018's headline is REAL — my independent
  certificate `MM_sec ∈ [−0.022254658205655785176146,
  −0.022254658205655785176145]` (width 5.9e-30) sits strictly inside 018's
  enclosure and is negative, in-regime, dichotomy-exact — and no surviving
  signed variant of the margin-modulated control was found in the fork 018
  left untested; but 018 ships one prose sentence its own checkpoint
  contradicts (raw-ladder MM_sec/MM_der signs at n ≥ 96), one false
  golden-ratio parenthetical repeated twice, a literal `PENDING-F`
  placeholder for certification jobs that had not landed, and two smaller
  imprecisions. VERIFIED_WITH_CORRECTIONS.
- **Tools:** `explore/uc_gap1c_skeptic.py` (parts S0–S4; stdlib only;
  deterministic, fixed seeds; runtime ~15 s; checkpoints
  `data/gap1csk_s[0-4].json`, log `data/gap1csk_run.log`). Command
  reproducing every number below:
  `python problems/union-closed/explore/uc_gap1c_skeptic.py 2>&1 | tee problems/union-closed/data/gap1csk_run.log`
  Independence contract: no imports from `uc_gap1_candidates.py`,
  `uc_or_avg.py`, `uc_or_avg_skeptic.py`, `uc_agg_ctrl_probe*.py`; witness
  atoms, θ\* weights and search endpoints copied as DATA; the MU(n, r)
  ladder and the dilution tuner re-implemented from 016 §1's prose and
  checked atom-for-atom (`MU(7,1) ==` the 10 witness atoms, exact dict
  equality). Two WebSearch queries for the novelty check (K5).
- **Sources:** 018 (record under attack, its engine read for claims only),
  007 §1 + lead 2 (the definitions 018 formalizes), 013 (the
  exact-arithmetic standard), 016/017 (ladder spec, θ\*, window law),
  gap1c checkpoints/logs (data).

Notation as in 007/018: tilt `π ∝ u(A)u(B) t^{|A∩B|}`, `t = 2^λ`; at a
nondegenerate history the conditional 2×2 table has zero-margins
`x = P(A_i = 0)`, `y = P(B_i = 0)`, realized both-zero `z̃`, odds ratio OR;
`z_ρ(x, y)` the Plackett both-zero probability at odds ratio ρ; MM_sec /
MM_der / MM_abs as in 018's part B; record threshold 0.38271.

## Claims attacked

1. **(K1, headline)** The certified kill: MM_sec on 007's 10-atom witness at
   exactly-rational t = 181/16 is negative
   (`∈ [−0.02225465820566…, −0.02225465820566…]`), in-regime (all
   elementwise marginals < 38271/100000), dichotomy exact — plus the three
   supporting derivations: the secant identity `h₂(z̃) = h₂(z_OR(x, y))`,
   `MM_sec ≡ 0` on product measures, and the degenerate-histories-
   contribute-0 convention one-liner.
2. **(K2)** The definitional fork: that 007 lead 2 admits exactly the three
   readings 018 tests, and that no natural *signed* variant survives —
   attacked by constructing and testing the readings 018 ignored
   (sensitivity evaluated at the **realized** OR instead of the target,
   signed and unsigned; `|dh/dρ|` weighting).
3. **(K3)** The float survival numbers (EVIDENCE level, spot-checks): witness
   MM_abs +2.354328 and MM_der −0.895542 at λ = 3.5; the θ\* ladder at
   n = 96, λ = 2 (plain −0.000760, MM_sec −3.819e-4, MM_abs +0.969); the
   part-D window value A = +3.747e-3 at n = 96, λ_win = 4.847/93, and the
   a₁/a₂ fit arithmetic and `a₁·n` growth claims; the part-R perturbation
   battery (re-run with different RNG seeds).
4. **(K4)** Scope honesty: Outcome claims vs what the logs/checkpoints show;
   the free-climb +0.0101 endpoint's regime/nondegeneracy status; the
   prior-art entry's leak_terms/gaps.
5. **(K5)** Novelty: is the margin-modulated statement (and its kill)
   plausibly a rediscovery?

## Refutations found

None headline-level: the certified kill, both survival claims, and the fork
analysis all stand. Four corrections and one cosmetic:

### C1. The raw-ladder MM_sec/MM_der sign sentence is contradicted by 018's own checkpoint
Location: 018 part B, "On the ladders at λ = 2 (raw witness weights and
016's θ\*), MM_sec is negative at every n ∈ {16, …, 160} (…), MM_der
negative throughout".

False for the **raw** ladder at n ∈ {96, 128, 160}: 018's own
`gap1c_partB.json` has MM_sec **+6.72e-3 / +6.05e-3 / +5.25e-3** and MM_der
**+0.99 / +0.97 / +0.95** there (plain +1.16 / +0.99 / +0.83 — the
`tune_dil` bisection lands on a different dilution branch at n ≥ 96, mm
0.311 vs 0.359 at n = 64, and the whole geometry flips sign). Reproduced
independently: my spec-rebuilt raw ladder at n = 96 gives MM_sec
+6.719528e-3, MM_der +0.9856 (dil 251.49), against −3.400008e-3 at n = 64
(part S3 — both matching 018's checkpoint to all printed digits, so the
error is in the prose, not the data). **Corrected statement:** at λ = 2,
MM_sec and MM_der are negative on the θ\* ladder at every n ∈ {16..160} and
on the raw ladder for n ≤ 64; on the raw ladder's n ≥ 96 dilution branch
both are positive. No downstream damage — the *survival* claims only need
MM_abs > 0 (true on every row), and the *kill* is carried by the witness —
but the sentence as written overstates how universally the signed forms
fail on ladders.

### C2. The golden-ratio parenthetical is false (twice)
Location: 018 Approach, "(z_1(x, x) = x² crosses 1/2 exactly at the
(3−√5)/2 barrier — the sensitivity sign boundary is the record threshold's
golden-ratio structure again)"; repeated in Why-it-failed as "The
golden-ratio coincidence (z crosses 1/2 exactly at marginal (3−√5)/2 for
equal margins) means this sign flip is intrinsic to the record regime".

Arithmetic: at element marginal p = (3−√5)/2 = 0.381966, the zero-margin is
x = 1 − p = (√5−1)/2 = 0.618034 and `x² = 0.381966 = p` — the φ identity
`φ² = 1 − φ` makes x² equal **the marginal itself**, not 1/2. The actual
z = 1/2 crossing at ρ = 1 is at x = 2^{−1/2}, i.e. element marginal
1 − 2^{−1/2} ≈ **0.292893**; at ρ = 2^{3.5} it moves to marginal ≈ 0.389113
(λ-dependent; its proximity to 0.38271 at λ = 3.5 is a coincidence of that
λ and is presumably where the confusion came from). Part S4. What survives:
the load-bearing claim that in-regime histories straddle z = 1/2 is true
(marginals < 0.38271 force x > 0.617, so z̃ ranges over both sides of 1/2),
and the witness's surplus rows do sit at z̃ > 1/2 — the *mechanism* stands,
the golden-ratio numerology attached to it is wrong, exactly the genre of
flourish 004 already killed once in 003 (R2's "p_h ≈ φ").

### C3. Part F shipped with a `PENDING-F` placeholder; only 1 of 4 certification jobs had landed
At review time 018's `--cert` pass had completed only job (1) (the witness
MM cert — the one the Outcome claims); jobs (2)–(4) (MM on the θ\* kill
instances at n = 96/160 at t = 4, window points at t = 51/50) were still
running, `gap1c_partF.json` had one row, `gap1c_cert.log` ends at the
witness line, and the record's §F contains a literal "`- PENDING-F`"
bullet. The Outcome section is careful to claim only the witness
certificate, so nothing certified is overstated — but a record should not
ship an unresolved placeholder, and any future reader should treat the
n = 96/160 MM certificates and window certificates as **absent from 018**
unless a follow-up lands them. (Float values at those points are covered by
K3 below.)

### C4. "positive only for λ ≲ 0.7" — the witness MM_sec crossing is at λ ≈ 0.83
018 part B. My profile: MM_sec = +2.00e-4 at λ = 0.75, +4.95e-5 at 0.8,
−3.21e-4 at 0.9 (part S2), so the sign change is between 0.8 and 0.9, not
at ≈ 0.7. Consistent with 018's own part-R grid (+ at 0.5, − at 1.0), which
is too coarse to support the "0.7". Minor: the claim's role (positivity in
the first-order-dominated region) is unaffected.

### C5 (cosmetic). Runtime drift
Tools bullet says parts A–D "~35 min"; the engine docstring says 20–40 min;
`gap1c_run.log` shows 3496 s ≈ 58 min. Recorded per 013 R1 precedent.

## Claims that survive (and what was done to break them)

### 1. The headline kill (K1) — re-certified from scratch; enclosures agree at width 5.9e-30

- **Hand re-derivations.** (i) *Secant identity:* on the coupling range
  `z ∈ (max(0, x+y−1), min(x, y))`,
  `d/dz ln ρ(z) = 1/z + 1/(1−x−y+z) + 1/(x−z) + 1/(y−z) > 0` (all four
  terms positive), so the margins-plus-odds-ratio parametrization of 2×2
  tables is bijective: the realized table's z̃ IS `z_OR(x, y)`, hence
  `h₂(z̃) = h₂(z_OR(x, y))` and, by the mean value theorem in
  λ' = log₂ ρ, `h₂(z̃) − h₂(z_t) = σ(λ̄)·(log₂OR − λ)` for some λ̄ on the
  secant. 018's "assembly-exact secant form" gloss is correct. (ii)
  *Product zero:* a product potential factorizes the kernel
  `Π_j u_j(a_j)u_j(b_j)2^{λ a_j b_j}`; conditioning on any history leaves
  the response factor `∝ [[u₀u₀, u₀u₁],[u₁u₀, u₁u₁2^λ]]`, so OR = 2^λ at
  every nondegenerate history (and Sinkhorn uniqueness carries this to
  product μ); by (i), z̃ = z_t row by row, so MM_sec ≡ 0 **term by term**,
  not just in aggregate. Engine: exact product anchor at w = 3/7, n = 5,
  t = 4 — plain aggregate exactly 0 (rational identity), MM_sec enclosure
  contains 0 at width 9.3e-30 (part S1). (iii) *Convention one-liner:* a
  margin-degenerate history has x ∈ {0, 1} or y ∈ {0, 1}; then z̃ is pinned
  (x = 0 or y = 0 ⟹ z̃ = 0; x = 1 ⟹ z̃ = y; y = 1 ⟹ z̃ = x) and
  `z_ρ` is pinned to the same boundary value for **every** ρ > 0, so the
  secant contribution is exactly 0 — including degenerate histories changes
  only the denominator. Checked as an exact identity on every degenerate
  positive-mass history of the witness (flag `deg_zero_ok`, part S1).
- **The certificate.** Own kernel (scaled integers), own census, own root
  and log₂ enclosures (part S1):

      witness, t = 181/16:
      MM_sec ∈ [−0.022254658205655785176146, −0.022254658205655785176145]
      plain aggregate ∈ [+1.844669005340829561098944, …945]
      max elementwise marginal = 0.3177701854 (exact < 38271/100000)
      dichotomy exact; degenerate rows contribute 0 exactly; 12 rows

  Strictly inside 018's 21-digit enclosure and matching 013's aggregate
  certificate digit-for-digit — three structurally different kits now agree
  on the same rational-input theorem. The **2-digit tidy witness** is also
  certified negative here (−0.0222736748, max marginal 0.3130, in-regime),
  upgrading 018's float −0.0223 to a certificate.
- **Kit validity attacked first** (part S0): log₂ enclosure vs 67 known /
  random values plus additivity (max width 9.1e-30); Plackett Newton root
  vs the float solver and the closed-form round-trip ρ(z_root)/t = 1 to
  1e-9 on 60 random (x, y, t), bracket endpoint sign identities
  `F(xy) = (1−t)xy(1−x)(1−y)`, `F(min) = min·(1−max)` asserted exactly per
  call; h₂ enclosure spot-checks; ladder-spec equality `MU(7,1) ==`
  witness. All pass.

**CONFIRMED, and strengthened (tidy witness now certified).** The
margin-modulated control in its secant reading is dead at n = 7 in-regime;
by the secant identity + MVT the *derivative-at-target* reading cannot be
salvaged as "the" chain-rule reading either (float MM_der −0.8955
reproduced; not separately certified here, matching 018's own scoping).

### 2. The definitional fork (K2) — fair, and the untested variants change nothing

007 lead 2's text ("weighted by the sensitivity of h(z_ρ(x̃,ỹ)) to ρ …
conjecture E[σ(x̃,ỹ)·(log₂OR − λ)] ≥ 0 … if the weighted mean is ≥ λ")
underdetermines where σ is evaluated; 018's three readings cover the target
and secant evaluations, and its "weighted mean" gloss honestly favors the
unsigned reading. The natural fourth family 018 skipped — σ evaluated at
the **realized** OR — was built and tested here (part S2):

    witness (λ = 3.5):  MM_der_rlz −1.458  MM_abs_rlz +2.133  MM_absrho_rlz +1.808
    θ* ladder n = 96, λ = 2:  MM_der_rlz −4.35e-1  MM_abs_rlz +9.34e-1  MM_absrho_rlz +8.80e-1

The signed realized-OR reading dies on the witness even harder than
MM_der; both unsigned realized variants (including the |dh/dρ| weighting,
which at the realized OR is a genuinely different weight, not a constant
multiple) behave like MM_abs — positive on the witness and on the kill
instance. **No surviving signed variant found; the signed/unsigned
dichotomy 018 draws is robust to the readings it did not test.** (The two
extra unsigned survivors inherit MM_abs's defect: not the derivative of any
gain.)

### 3. The float survival numbers (K3) — reproduced on a spec-rebuilt ladder, different seeds

- Witness: MM_abs **+2.354328**, MM_der **−0.895542**, MM_sec −0.022255 —
  all match (part S3).
- θ\* ladder, n = 96, λ = 2, my own `tune_dil` (dil 31.6228): plain
  **−0.000760** (−7.598652e-4), MM_sec **−3.819411e-4**, MM_abs
  **+0.969045**, mm 0.309 — all match 018 to printed precision, on a
  builder written from 016 §1's prose, not imported.
- Window point n = 96, λ_win = 4.847/93 = 0.052118: A = **+3.747072e-3** ✓.
  The a₁/a₂ 2×2 solves re-derive for **all 10** stored part-D rows
  (rel. 1e-9), and the `a₁·n` growth sequence on θ\* re-computes as
  4.73 → 7.19 → 10.37 → 13.50 → 18.13 ✓ (raw rows re-derive too; note the
  raw fit's a₂ changes sign across n — +0.027 at n = 160 — so 018's "same
  shape" for raw is loose, though its a₁·n ≈ 66 claim is right).
- Perturbation battery with different seeds (7 and 424242, 20 draws each,
  3% multiplicative): **0/40 failures**, MM_sec ∈ [−0.0230, −0.0216] —
  018's seed-2024 result is not a seed artifact.
- The free-climb endpoint (from `gap1c_partC.json` `best_free`):
  re-evaluated with my census — MM_abs **+1.007614e-2**, MM_sec
  −8.327e-5, max marginal **0.153275 < 0.38271** (genuinely in-regime),
  dichotomy ok, H = 1.019 bits (nondegenerate, 14 atoms). The tightest
  MM_abs margin is real (part S4).

CONFIRMED (as EVIDENCE, at 018's own scope).

### 4. Scope honesty (K4) — Outcome matches the logs, with the C1/C3 caveats

"Positive on everything tried" for MM_abs: true of every row in
`gap1c_part[B,C,R].json` (checked). The window-candidate claims match
`gap1c_partD.json` including the λ → 0 boundary convergence of the joint
climb (+1.66e-5 at λ = 0.0008). Part E's dilution-branch jump at
n = 288 (dil 80 → 160) is as recorded. One scope note short of a
correction: of the 8 free-support climbs, only the three witness-seeded
endpoints are in-regime (the five random-seeded endpoints have max
marginal 0.384–0.663), so the "8 climbs" battery carries less in-regime
weight than its count suggests — the floor claim itself is in-regime and
correct. The prior-art entry's leak_terms/gaps do name the new findings
(MM_sec/MM_der/MM_abs, secant form, abs-sensitivity-or-control,
0.0222546582); status REFUTED with the two survivors relegated to `gaps`
is the right bookkeeping.

### 5. Novelty (K5) — checked, low risk

Two web searches ran (2026-08-06): "union-closed sets conjecture Plackett
odds ratio entropy coupling sensitivity weighted" and a narrower
odds-ratio/tilt/secant query. Results surface the known entropy/coupling
literature (Gilmer arXiv:2211.09055 line, Sawin, Chase–Lovett, Liu
arXiv:2306.08824, Yu–Cambie dimension-free) — nothing resembling a
sensitivity-weighted / margin-modulated conditional-OR control or its
refutation. The statement is internal to this lab's 007-lead-2 route, as
expected; the kill is a result about the lab's own conjecture, so novelty
risk is essentially nil.

## Verdict

| # | Claim | Verdict |
|---|-------|---------|
| 1 | Certified witness kill of MM_sec at t = 181/16, in-regime, dichotomy exact | **CONFIRMED** — independent certificate, width 5.9e-30, strictly inside 018's enclosure; tidy witness additionally certified (−0.02227) |
| 1b | Secant identity; product-measure zero; degenerate-contribution-0 | **CONFIRMED** — re-derived by hand; product zero and degeneracy checked as exact identities |
| 2 | Three-way fork; only unsigned readings survive | **CONFIRMED and extended** — realized-OR variants tested: signed dies (−1.458 on witness), unsigned survives; no surviving signed variant exists in the natural family |
| 3a | Witness/ladder float numbers; window value; fit arithmetic; a₁·n growth | **CONFIRMED** (spec-rebuilt ladder, own tuner; all 10 fit rows re-derive) |
| 3b | Ladder sign prose at λ = 2 | **CORRECTED (C1)** — raw ladder at n ≥ 96 has MM_sec, MM_der > 0 per 018's own checkpoint; θ\*-only statement survives |
| 3c | Perturbation robustness of the kill | **CONFIRMED** with different seeds (0/40 failures) |
| 4 | Outcome scope; free-climb endpoint in-regime; prior-art bookkeeping | **CONFIRMED** (with the 3-of-8-in-regime scope note; C3 for the PENDING-F placeholder) |
| 5 | Golden-ratio sign-boundary parenthetical | **REFUTED (C2)** — x² = marginal, not 1/2, at the barrier; true ρ=1 crossing at marginal 1 − 2^{−1/2} ≈ 0.293; straddling claim itself survives |

**Net assessment:** 018's headline — the margin-modulated OR control is
REFUTED in both signed readings at n = 7, certified float-free on 007's own
witness — holds and is now double-certified by structurally disjoint kits;
its two survival claims stand at their stated EVIDENCE scope, and the fork
analysis is complete once the realized-OR variants are added (they change
nothing). The corrections are all prose-level: a checkpoint-contradicting
sentence (C1), a numerological flourish (C2, the same failure genre 004
flagged in 003), an unlanded-certification placeholder (C3), and two
imprecisions (C4, C5). Ledger status VERIFIED_WITH_CORRECTIONS;
`abs-sensitivity-or-control` and `lambda-restricted-or-control` are the
correct surviving gap tags.

## Residual risk

- **Shared-definition risk** (as in 013): every engine including this one
  formalizes the census from 007 §1's prose. Exact arithmetic removes
  numerical error, not a conceptual error in that shared definition; the
  +1.844669005 aggregate anchor ties this review to 013's independently
  hand-derived census, which is the main mitigation.
- **Part-F jobs (2)–(4)** of 018 were still running at review close;
  whatever they produce is unreviewed. The MM_der reading remains float
  EVIDENCE (neither 018 nor this review certified it; the secant
  certificate plus the MVT tie is the argument that it cannot be the
  assembly's reading anyway).
- **Search-completeness claims** (θ-opt floors, free climbs, window climbs)
  were spot-checked for reproducibility, not re-run at scale; they remain
  EVIDENCE with 018's own "one adversary genre" caveat. My variant zoo
  covered the natural sensitivity-evaluation points (target, secant,
  realized; λ- vs ρ-derivative) but not arbitrary nonnegative weight
  functions of the margins — a survivor family there is unconstrained by
  this review (and by design would drift further from chain-rule meaning).
- The two unsigned realized-OR variants found positive here were tested at
  two instances only; they are observations, not survival claims.

## References

- Reviewed record: `problems/union-closed/attempts/018-surviving-controls-vs-ladder.md`;
  its engine `explore/uc_gap1_candidates.py`; checkpoints
  `data/gap1c_part[A-F,R].json`, logs `data/gap1c_{run,ext,robust,cert}.log`.
- Context: `attempts/007-averaged-or-control.md` (§1 census, lead 2, part-T
  witness), `attempts/013-skeptic-review-of-007.md` (exact standard,
  +1.844669005 anchor), `attempts/016-*` / `017-*` (MU(n, r) spec, θ\*,
  window law, interval-squaring idea), `attempts/004-*` (review shape;
  φ-numerology precedent).
- This review's tool/data: `explore/uc_gap1c_skeptic.py`;
  `data/gap1csk_s[0-4].json`; `data/gap1csk_run.log`.
- Novelty check (WebSearch, 2026-08-06): arXiv:2211.09055 (Gilmer),
  arXiv:2306.08824 (Liu), arXiv:2212.12500 (Yu), arXiv:2305.19338,
  arXiv:2306.12351 (survey) — no margin-modulated/sensitivity-weighted OR
  control in the literature surfaced.
