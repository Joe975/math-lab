# 020 — Gap 1 closes: every odds-ratio-control candidate refuted, swarm-guided

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-19
- **Mode:** informed
- **Type:** computational search + certified refutation, with attack
  directions drafted by a dual-family external LLM swarm (protocol
  `docs/SWARM.md`)
- **Tools:** `explore/uc_gap1_deep.py` (parts A–F; deterministic, fixed
  seeds; ~25 min total across runs; checkpoints
  `data/gap1deep_part[A-F]*.json`, logs `data/gap1deep_*run.log`).
  Swarm briefs: `explore/swarm020_attack_brief.md` +
  `swarm020_attack_values.txt` (6 tasks × 2 families: gpt-5.6-luna effort
  medium, gemini-3.7-flash effort medium; outputs in
  `$MATHLAB_OUT/swarm/020-attack-{gpt,gem}/`, ≈ $0.15 of the cycle's
  ≈ $0.21 swarm spend). Reuses 013's exact census
  (`uc_or_avg_skeptic.py`), 018's `census_mm`/enclosure kit, 014/016's
  Sinkhorn engines.
- **Sources:** none external beyond the swarm workers' recall of Plackett
  (1965) and Karlin–Rinott TP2 facts, used only as orientation [T].

## Approach

**Context.** After 016/017 (aggregated control dead at large n) and
018/019 (both signed margin-modulated readings dead at n = 7), exactly two
Gap-1 candidates survived: the unsigned |σ|-weighted control `MM_abs`
(018 lead 1 = queue 1a) and the λ-window-restricted plain aggregate
(017 C2 / 018 lead 4 = queue 1b). The queue's replacement program asked
for probe-first attacks on both.

**Why a swarm-guided attack rather than more blind annealing.** 018's
free-support search at n = 8 had ground MM_abs only to +0.0101 — blind
annealing was hitting a floor. Both external model families were given the
precise functional setup (margins, Plackett root, sensitivity σ_λ, the
known facts F1–F7) and asked independently for kill geometries, floor
proofs, a₁-degenerate window families, the assembly's true requirement,
and a per-history-weighting no-go (six briefs each). The director then
implemented, corrected, and verified everything — worker output was
treated as candidate text only, per `docs/SWARM.md`.

**The two decisive worker convergences** (each produced independently by
both families, which is what earned them implementation):

1. *The first-order coefficient of MM_abs's numerator is a kernel-weighted
   Gram form* Σ_{a,b} ν(a)ν(b)·K(x_a, x_b)·D_a·D_b with K(x, y) =
   |σ_λ(x, y)| — so, unlike Theorem C's perfect square, its positivity is
   equivalent to positive-semidefiniteness of the weight kernel on the
   realized margin set. (One worker supplied a wrong closed form for σ₀ —
   dz/dρ at ρ = 1 has denominator 1, not (1−x)(1−y)+xy; the director's
   first construction failed because of it, and the engine's own
   `dh_dlam` was used from then on.)
2. *Prefix-signed cancellation*: Theorem C's per-i coefficient
   ‖Σ_a ν(a)D_a‖² is the norm of a **signed vector sum**, so a₁ can vanish
   by cancellation across prefixes on measures far from products —
   exactly the a₁-degenerate direction 018 lead 4 asked for, with an
   explicit selector-block recipe.

Direct kernel scans then showed **K is not PSD in-regime at any tested λ**
(2×2 minors already fail, worst near margins (0.2, 0.68); mechanism: the
diagonal weight collapses at the kink z_ρ(x,x) = 1/2 while the cross term
stays large), turning candidate-killing into a *design* problem rather
than a search problem.

## What was done

Statements under test (both from the queue, made precise in 018):

- **(I)** MM_abs(μ, λ) = E[|σ_λ|·(log₂OR − λ)] / E[|σ_λ|] ≥ 0 for every μ
  with element marginals ≤ 0.38271 and every λ > 0.
- **(II)** A(μ, λ) = E[log₂OR − λ] ≥ 0 for every in-regime μ on n
  coordinates and every λ ≤ λ_win(n) = 4.847/(n−3).

### A. Deep free-support anneal on MM_abs (queue 1a as written)

n ∈ {10, 12, 14, 16}, λ ∈ {1, 2, 3.5}, 2500 steps, atom-count mutations,
three seed genres (018's `best_free`, the swarm's frustrated-block shape,
random). Result: the healthy-margin floor did not break by blind search
(frustrated genre floors +0.014…+3.0 in-regime), but the `bestfree` genre
repeatedly escaped into a **near-empty corner** (max marginal ~1e-9,
weights spanning 17 orders of magnitude) with float MM_abs as low as
−24.5. **These are float artifacts**: exact rational recomputation of the
worst endpoint (n = 14, t = 11) flips its sign to +23.39 — the second
confirmed instance of 007's shared-IEEE failure mode
(`gap1deep_partA.json`; the artifact and its exact refutation are kept as
a calibration case). No genuine kill came from part A.

### C/E. The kernel-guided kill of (I) — REFUTED, certified

`model_c1(μ, n, λ)` computes the kernel Gram form directly from μ (prefix
masses, conditional zero-margins, future-bias vectors D_a; weights from
the engine's `dh_dlam`). The model was validated exactly against the
census on a 3-coordinate toy (per-history deviation = λ·D_a·D_b to
0.2% at λ = 0.05, and coordinate-1 deviation = λ‖D_∅‖²).

Annealing μ at n ∈ {4, 5} (≤ 14 atoms, weight floor 1e-6 — deliberately
excluding the part-A corner) to minimize `model_c1` under the marginal
cap, then verifying the best designs against the real Sinkhorn coupling:
**six of ten designs give census MM_abs < 0 in-regime**, e.g. −0.110 at
n = 5, λ = 0.5 (max marginal 0.382669), −0.046 at n = 4, λ = 0.5
(7 atoms, max marginal 0.3826). Convention cross-check: extracting the
potential u from the Sinkhorn scalings and running 018's
`census_mm(u, n, λ)` reproduces −0.1100587 digit-for-digit, dichotomy
exact. 3-significant-digit tidy stays negative (−0.1099).

Exact certification (part E, `exact_mm_abs`): certified interval
enclosure of the MM_abs numerator — z_target by dyadic bisection of the
exact Plackett quadratic, h₂′ via directed log₂ enclosures, dz/dρ by
monotone rational interval arithmetic, the global positive constant
t·ln2 dropped as sign-irrelevant. At rational tilts t = 3/2 and 6/5
(λ ≈ 0.585, 0.263 — both inside the window at these n):

    n = 5 witness: num ∈ [−3.403811e-3, −3.403811e-3]   (t = 3/2)
                   num ∈ [−2.922542e-3, −2.922542e-3]   (t = 6/5)
    n = 4 witness: num ∈ [−1.784190e-3, −1.784190e-3]   (t = 3/2)
                   num ∈ [−2.884705e-3, −2.884705e-3]   (t = 6/5)
    all: max marginal < 0.38271 exactly, dichotomy exact, plain
    aggregate simultaneously certified POSITIVE (+0.033…+0.152).

**(I) is REFUTED at n = 4 with 7 atoms, at healthy margins, inside the
λ-window, with the plain aggregate positive on the same tables** — the
deficit the |σ| weight was designed to suppress was moved to margins
where |σ| is large, via the kernel's non-PSD directions.

### The in-window secant kill (assembly requirement)

Both families' task-4 analyses independently identified the assembly's
true aggregate requirement as the secant condition
E[h₂(z̃) − h₂(z_{2^λ})] ≥ 0 — i.e. MM_sec's numerator, which 018
certified negative only at unrestricted λ on the n = 7 witness. On this
record's witnesses, 018's own `exact_mm` certifies, **inside the
window** (λ_win(5) = 2.42, λ_win(4) = 4.85):

    n = 5: MM_sec ∈ [−2.318587e-3, ·] at t = 3/2;  [−4.024146e-3, ·] at t = 2
    n = 4: MM_sec ∈ [−4.780764e-3, ·] at t = 3/2;  [−8.869717e-3, ·] at t = 2
    (all in-regime, dichotomy exact; enclosure widths < 1e-15)

So the chain-rule assembly cannot demand aggregate gain parity even
λ-window-restricted: whatever bridge survives must be λ-integrated or
non-per-history. The λ-integrated variant ∫₀^λ E[σ_s(log₂OR − s)] ds
was float-tested on the same witnesses and is **also negative**
(−4.9e-3 at λ = 1 on the n = 5 witness; float only, no certificate).

### B. The a₁-cancellation kill of (II) — REFUTED, certified

Implementing the swarm's convergent cancellation mechanism: anneal +
damped Gauss–Newton drives the cancellation residuals Σ_a ν(a)D_a to
machine zero (≤ 1e-33 summed squares) on small free measures, with a
non-product score (max single-prefix ‖D_a‖) bounding the family away from
the trivial product solution. Then descend the λ² coefficient **within
the cancellation manifold** (mutate → re-polish → accept if a₂ drops).

Outcome: both signs of a₂ occur at a₁ = 0. The 3-coordinate selector
genre bottoms out at a₂ ≈ +7e-9 (positive floor); the free n = 4 genre
crosses: a 15-atom μ with marginals ≤ 0.2575, cancellation residuals
≤ 1.7e-16, and A/λ² ≈ −2.0e-5 … −0.9e-5 across λ ∈ [0.005, 0.08]
(consistent quadratic, crossing to positive near λ ≈ 0.13). Exact
certification at three rational tilts:

    t = 21/20 (λ = 0.0704): A ∈ [−5.256144e-8, −5.256144e-8]
    t = 16/15 (λ = 0.0931): A ∈ [−6.242487e-8, −6.242487e-8]
    t = 27/25 (λ = 0.1110): A ∈ [−5.508935e-8, −5.508935e-8]
    all: in-regime exactly (max marginal 0.2575), dichotomy exact.

λ_win(4) = 4.847, so all three certified points are deep inside the
window. **(II) is REFUTED at n = 4.** The kill is not a small-n artifact:
the tilt kernel factorizes over disjoint coordinate blocks, and float
evaluation of the 2- and 3-fold tensor confirms in-window violations at
n = 8 (A = −2.9e-8 at λ = 0.05, window 0.97) and n = 12 (−2.8e-8 at
λ = 0.05, window 0.54); since a₂ adds across blocks while the
cancellation residual stays ~0, a violating λ exists inside λ_win(n)
for every n by tensoring (float + SPECULATION at the "every n" level;
certified at n = 4, float-verified at n = 8, 12).

### F. The per-history-weighting no-go — certified witnesses

One family's task-5 return supplied the right reduction: any weight
w(x, y, λ) ≥ 0 that feeds the assembly must equal σ_λ wherever
σ_λ > 0, i.e. on R₊ = {(x, y) : z_{2^λ}(x, y) < 1/2}; on measures whose
histories all lie in R₊ the choice of w elsewhere is inert. So a
σ-weighted violation with **every history certified inside R₊** rules
out every such w at once. Constraining the part-C designer to R₊
(z_target ≤ 0.485 at the design tilt) and certifying at t = 7/5:

    n = 5: num ∈ [−2.0682e-3, ·], every z_target certified < 1/2
    n = 4: num ∈ [−1.7775e-3, ·], every z_target certified < 1/2
    n = 4: num ∈ [−3.6933e-4, ·], every z_target certified < 1/2
    all in-regime, dichotomy exact.

**No per-history weight function compatible with the chain rule is
nonnegative on all in-regime μ** — 018 lead 2's question is answered
negatively (for the "equals σ on R₊" compatibility notion; a weaker
compatibility notion would need its own statement).

## Outcome

**REFUTED (certified, float-free): both surviving Gap-1 candidates.**
(I) the unsigned |σ|-weighted control fails at n = 4 (7 atoms) and n = 5,
at healthy margins, in-window, at exact rational tilts t ∈ {6/5, 7/5,
3/2}; (II) the λ-window-restricted aggregated control fails at n = 4 at
t ∈ {21/20, 16/15, 27/25}, marginals ≤ 0.2575. Additionally the
assembly's aggregate secant requirement is certified violated inside the
window (n = 4, 5; t ∈ {3/2, 2}), and certified R₊-pure witnesses rule
out every per-history chain-rule-compatible weighting. With 005/006
(pointwise), 007/013 (per-i averaged), 016/017 (aggregated, ∀λ),
018/019 (signed margin-modulated), **every stated form of Gap 1 — the
Plackett odds-ratio control — is now dead.**

Scope and non-claims: nothing here touches the coupling interface itself
(003/004), Gap 2 (mutual-information tax, 009/011 — still EVIDENCE-live),
or gap 3 (recipe totality). The certificates are for the specific
rationalized measures and tilts named in `gap1deep_partE.json` /
`partF.json` / the part-B block above; the tensor extension past n = 12
and the λ-integrated variant's negativity are float-level. The route's
model ceiling 0.4315 is unaffected; what died is the entire family of
per-history odds-ratio bridges between the coupling and the entropy
assembly. **Skeptic pass pending** — certificates reuse 013's log₂
enclosures and this record's own dz/dρ interval step; an independent kit
must re-derive them (queue item).

## Why it failed / what survived

**Why (I) died: the weight kernel is not PSD, and the census can realize
its negative directions.** |σ_λ(x, y)| vanishes on the kink surface
z_{2^λ}(x, y) = 1/2 while staying large at mixed margin pairs straddling
it, so 2×2 minors K(a,b)² > K(a,a)K(b,b) exist entirely in-regime; a
selector coordinate whose classes realize those margins with anti-aligned
future biases makes the weighted first-order term negative while
Theorem C's unweighted square stays nonnegative (it is tiny, not zero,
at the design point — the plain aggregate stays positive on every kill
table). 007 lead 2's intuition ("downweight the light slice") was sound
about the light slice but fatally reweights the healthy-margin bulk:
any weight that varies across margins re-opens the Gram form.

**Why (II) died: products are not the only a₁-degenerate measures.**
014/015's "products are the equality set" (forward inclusion) left open
exactly this: signed cancellation Σ_a ν(a)D_a = 0 at every coordinate
with individual D_a large. On that manifold A = a₂λ² + O(λ³), and a₂ has
no sign protection — the free n = 4 genre finds certified-negative a₂
with margins comfortably in-regime. The window law restricts λ but not
the sign of a₂, so first-order protection (018 part D's growing a₁·n on
the ladder genre) was a property of that genre, not of the statement.

**What survives for the route:** the coupling interface, the
mutual-information-tax line (Gap 2), and one genuinely untested shape —
conditions that are *not* per-history weightings of log₂OR deviations
(e.g. coordinate-level budgets with cross-history structure, or the
λ-integrated forms with weights outside the σ-compatibility class). The
no-go's boundary is precise: it kills w = σ-on-R₊ compatibility, nothing
more.

**Methodological.** (i) The two cross-family convergences were both
implementable and both decisive; the one wrong formula a worker shipped
(σ₀'s denominator) was caught because the director validated the model
against the census before using it — the SWARM.md filter earned its cost.
(ii) The part-A corner artifacts re-confirm 007's lesson with a clean
exact-arithmetic kill; every future MM search needs the cell-floor guard.
(iii) Kernel-PSD analysis of weight functions is a reusable design tool —
tag `sensitivity-kernel-design`.

## Leads generated

1. **Skeptic pass on this record** (mandatory): re-derive `exact_mm_abs`
   (the dz/dρ interval step is new) and the part-B/F certificates with an
   independent enclosure kit (017's interval-squaring path shares no code
   with 013's atanh kit used here); re-run the three part-B tilts; verify
   the R₊ certificates' z < 1/2 legs. Per SWARM.md the skeptic should be
   drafted cross-family from this record's gpt-heavy analysis inputs —
   or by a fresh director session with no swarm at all.
2. **Does the tensor kill survive certification at n = 8?** The m = 2
   tensor (225 atoms) is within reach of the exact census with margin-pair
   caching. Certifying it would upgrade "violating λ exists inside the
   window for every n" from SPECULATION to a per-n certificate schedule.
3. **Restate Gap 1 outside the per-history class, or retire it.** The
   surviving shapes after the no-go: (a) coordinate-level budgets with
   cross-history cancellation (012's corrected B1/B2 machinery is the
   natural host); (b) λ-integrated conditions with weights NOT equal to σ
   on R₊ (the no-go's boundary); (c) abandoning the OR bridge entirely
   and pushing Gap 2's chain-rule-assembly value directly (009/011's
   TAX/ST decomposition survives everything found here — check the two
   new witnesses against it as a first falsifiable step).
4. **Check the 020 witnesses against Gap 2.** Cheap and sharp: compute
   009's CR/TAX/ST on the n = 4/5 kill measures. If Gap 2's candidate
   survives the geometries that killed Gap 1, that is real evidence its
   mechanism is different; if it dies, the route needs restructuring at
   the interface level.
5. **The a₁-cancellation manifold deserves a map**: which (n, genre)
   have negative-a₂ points under the marginal cap, and does the blk3
   genre's positive floor (+7e-9) have a mechanism? A closed-form a₂ on
   the cancellation manifold is plausibly derivable and would decide the
   "every n" tensor claim without computation.

## References

- This repo: 018/019 (candidates + kit), 016/017 (ladder + window),
  014/015 (equality set), 013 (exact standard), 007 (witness, Theorem C),
  009/011 (window law), 005/006, 008/012.
- Swarm: `explore/swarm020_attack_brief.md`, `swarm020_attack_values.txt`;
  models gpt-5.6-luna, gemini-3.7-flash, 12 attack briefs total; outputs
  under `$MATHLAB_OUT/swarm/020-attack-*` (gitignored; prompt hashes in
  the `.meta.json` files).
