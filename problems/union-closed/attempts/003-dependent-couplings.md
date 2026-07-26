# 003 — Dependent couplings (Idea C): formalized, adversarially tested, LIVE

- **Problem:** union-closed sets conjecture (Frankl), `problems/union-closed.md`
- **Date:** 2026-07-26
- **Type:** formalization + adversarial evaluation of Idea C from 001
  (family-adaptive / dependent couplings), executing lead 3 of 002. Outcome:
  **the route survives every adversarial family that killed ideas B and the
  KL ladder** — the first candidate interface in this project not refuted by
  the known counterexample genres. No new bound is claimed.
- **Tools:** `tools/uc_couplings.py` (written and run here; deterministic;
  parts A, B, D, E stdlib-only, parts C, F need numpy; full log and JSON
  checkpoints in `attempts/union-closed/data/003_*`), `tools/uc_weighted_kl.py`
  (002's machinery, re-run as baseline, reproduces exactly). Primary sources:
  Liu arXiv:2306.08824 (fetched and transcribed 2026-07-26; details below
  marked [L] are from that machine transcription — orchestrator should spot-
  check against the PDF before relying on clause-level details), plus the
  001/002 source base (Gilmer, AHS, Chase–Lovett, Sawin, Pebody, Yu, Cambie).

Notation as in 001/002: `h` binary entropy (bits), `ψ = (3−√5)/2 ≈ 0.381966`,
`φ = 1−ψ`, `q = 2p−p²`, record ≈ 0.38271 (Liu). `x, y` denote conditional
ZERO-probabilities of coordinate bits; `U = A∪B`.

## Approach

002 closed the "reweight the functional" direction with a no-go: any
functional of `(law(U_iid), μ)` that is smoothing-insensitive (increase
`O(log 1/δ)` under δ-mass smoothing of μ) has ceiling exactly ψ. Its lead 2/3
pointed here: change the *sampler*, not the functional. Plan: (1) pin down
Liu's conditionally-iid coupling (the proven frontier, 0.38271) and verify —
against the no-go's hypotheses, not by assertion — why it evades the no-go;
(2) formalize the next coupling class, overlap-biased couplings with both
marginals prescribed; derive what replaces the two uses of independence in
the c=0 proof; (3) per 002's protocol, adversarially test the candidate
functional FIRST, on (a) Sawin's geometric mixtures, (b) the Chase–Lovett
smoothed-slice family, plus union-closed controls and a fresh adversary hunt;
only record optimism if the functional separates fakes from genuine families.

## What was done

### 1. The statement under test, and why this interface is legal

For a class `C(μ)` of couplings (joint laws of `(A,B)` with BOTH marginals
equal to μ), define

    Gain_C(μ) = sup_{π ∈ C(μ)} H_π(A∪B) − H(μ)

    (S-coup at p, class C): every μ on 2^[n] with H(μ) > 0 and all
    marginals < p has Gain_C(μ) > 0.

**Licensing lemma** (the family-side interface; one line): if `F` is exactly
union-closed and `μ = Unif(F)`, then every coupling supported on `F × F` —
in particular every `π ∈ C(μ)` — has `A∪B ∈ F` a.s., so
`H_π(U) ≤ log|F| = H(μ)`, i.e. `Gain(Unif F) ≤ 0`. Hence (S-coup at p) ⟹
Frankl at p. Two structural points, both verified computationally below:

- If `C` contains the diagonal coupling (`B = A`), then `Gain_C(μ) ≥ 0` for
  every μ, with equality for union-closed uniforms. So the statement needed
  is *strictness* of a supremum that is always ≥ 0 — qualitatively the right
  shape for an extremal characterization of union-closed families.
- Approximate closure licenses **nothing** here: a coupling may put Θ(1)
  mass on the o(1) fraction of escaping pairs. This is exactly why neither
  the Chase–Lovett ψ-barrier (001 (c)) nor 002's no-go applies: both live at
  the interface `Pr_{iid}[A∪B ∈ F] = 1−o(1)`, and the coupling interface
  consumes the full worst-case closure quantifier. The RHS `log|F| = H(μ)`
  also uses that μ is uniform on its support — the second quantifier every
  (S_c)-type statement relaxed away (002, lead 2).

**The 002 ladder embeds at the bottom of this one** (rigorous, small):
for the mixture family `π_δ = (1−δ)·diag ⊕ δ·π₀`, `law_{π_δ}(U) =
(1−δ)μ + δν` with `ν = law_{π₀}(U)`, and (concavity of entropy plus
`d/dδ H((1−δ)μ+δν)|₀ = CrossEnt(ν,μ) − H(μ)`; verified numerically to 7e-7):

    sup_δ H_{π_δ}(U) > H(μ)  ⟺  H(ν) + D(ν‖μ) − H(μ) > 0.

With `π₀ = iid` the right side is exactly Gilmer's Conjecture-1 functional
`f_μ(1)` from 002. So: **the diag⊕iid coupling rung is equivalent to
Gilmer's refuted Conjecture 1, and (via 002 part C: Sawin gadgets have
`f_μ(1) < 0` at marginals → ψ⁺) its ceiling is exactly ψ.** The KL ladder
was the first-order shadow of the coupling ladder in one particular
direction; its death says that direction fails, not the class.

### 2. Liu's conditionally-iid coupling, and its no-go evasion verified

**The construction** [L, consistent with the abstract fetched verbatim]:
Liu's class `C₃(μ)` consists of couplings built coordinate-by-coordinate:
given histories with conditional one-probabilities `s = P(A_i=1|A_{<i})`,
`t = P(B_i=1|B_{<i})`, the pair `(A_i, B_i)` is drawn, conditionally on a
shared auxiliary `u ~ P_u`, as independent `Bern(q_{u,s}) ⊗ Bern(q_{u,t})`
with `E_u[q_{u,s}] = s` (marginal preservation). Equivalently: symmetric
conditional couplings that are mixtures of product (rank-1) components —
the positive-semidefinite cone inside Sawin's symmetric-coupling class
`C₂`. Chain: iid ⊂ {diag⊕iid} ⊂ C₃ ⊂ C₂ ⊂ all couplings. Proven values:
iid = ψ (Pebody: exactly); Sawin/Yu/Cambie mixture scheme = 0.38234
(optimum of its class); Liu C₃-perturbation ≈ 0.38271 [L: "under
numerically verified hypotheses", 9-dimensional optimization; the paper's
Remark notes a maximal-correlation variant gives no further improvement].
The reason conditioning is the right shape: given `u`, BOTH uses of
independence in Gilmer's argument (below) are restored, at the price of
mixture bookkeeping — no `I(A;B)` tax appears because the correlation is
routed through common randomness rather than direct dependence.

**No-go evasion, verified against the hypotheses** (not asserted). 002's
no-go covers functionals with (i) family-side tautology and (ii) the
smoothing bound "Φ increases ≤ O(log 1/δ) under δ-mass smoothing", in
particular any `Φ ≤ H + O(D)`. The coupling functional `Φ_C(μ) =
sup_{π∈C₃} H_π(U)` violates (ii) maximally, and this is now a *computed
fact* on 002's own certificate instances (part D of `uc_couplings.py`,
exact large-n evaluation; coupling used: share Sawin's latent `k`, then
per-coordinate common-randomness mixtures with union marginal `m_k = 1/2`
for `p_k ≤ 1/2`, diagonal for `p_k > 1/2` — a member of C₃, PSD checked):

    Sawin gadget          max marg   gain_iid      D     gain_coupled
    n=2000  ū=.390 θ=.05  0.402269    −60.63    4.35        +66.88
    n=20000 ū=.386 θ=.02  0.390799   −264.96    5.66       +741.48
    n=60000 ū=.3823θ=.001 0.382536    −51.65    9.97      +2418.56

`Φ_C₃ − (H_iid + cD) = Θ(n)` for every fixed c: the coupling functional is
**smoothing-sensitive** — it responds to a planted δ-mass component at
probability scale Θ(δn) (by re-routing that block of the coupling), not at
log-likelihood scale O(log 1/δ). Consistency check on Liu's theorem itself:
his statement covers all μ with marginals < 0.38271, including Sawin
gadgets with ū ∈ (ψ, 0.38271); if his functional obeyed the no-go's
smoothing bound those gadgets would refute his theorem. The third row above
(marginals 0.382536 < record) exhibits the coupling response explicitly on
exactly the instance that killed every rung of the KL ladder.

### 3. Overlap-biased couplings: formalization and the chain-rule anatomy

**The class.** For `λ ∈ R` (log₂ odds tilt, `ρ = 2^λ`):

    π_λ(A,B) = 2^{ g(A) + g(B) + λ·|A∩B| },

with symmetric Schrödinger/Sinkhorn potentials `g` fitted so both marginals
equal μ (existence/convergence: standard for strictly positive kernels;
observed residuals < 1e-9 everywhere). `λ > 0` up-weights high-overlap
pairs — exactly the pairs whose unions the Chase–Lovett family fails to
keep (001 (c) item 1); `λ = 0` is iid; `λ < 0` anti-overlap. For
exchangeable μ, `g` is a function of `|A|` and everything reduces to weight
profiles (engine: part C; validated against brute force at n=8, against an
independent atom-level engine to 4e-12, and against the closed form below
to 1e-10).

**Where independence is used in the c=0 proof — twice — and what replaces
it.** With histories `a = A_{<i}`, `b = B_{<i}`, the chain rule
`H(U) ≥ Σᵢ E_{(a,b)~π}[h(z_i(a,b))]`, `z_i = Pr[A_i=0, B_i=0 | a,b]`,
holds for EVERY coupling. Independence enters afterwards:

1. *Bit level:* `z = x(a)·y(b)` — conditional independence of the two bits
   given histories, plus locality (`Pr[A_i=0|a,b] = Pr[A_i=0|a]`).
   Replacement: `z = z_ρ(x̃,ỹ)`, the Plackett coupling of the conditional
   margins at conditional odds ratio ρ; for π_λ on product μ the odds ratio
   is exactly `2^λ` at every history and locality survives (verified: engine
   matches per-coordinate closed form to 1e-10). For general μ the
   conditional odds ratios drift with history; the needed lemma — λ ≥ 0
   forces all conditional odds ratios into `[1, 2^λ]` (an FKG/four-function-
   style induction over the future potential) — is **unproven SPECULATION**
   and is the precise technical gap between "functional evaluates well" and
   "theorem".
2. *History level:* `E[x(a)h(y(b))] = E[x]·E[h(y)]` — independence of the
   history pair. Under a genuine coupling this fails AND the bookkeeping
   target degrades: `Σᵢ E_π[h(x̃_i)] = H(A) − Σᵢ I(A_i; B_{<i} | A_{<i})`,
   a mutual-information tax ≈ I(A;B). The diagonal pays full tax (zero
   gain possible) — which is why the class needs interior tilts and why all
   gain curves computed below have interior maxima in λ. Liu's conditioning
   is precisely a device that pays the tax with auxiliary-mixture
   bookkeeping instead.

**The per-coordinate calculus** (part A, rigorous, closed forms verified to
1e-13): a product `Bern(p)^n` under per-coordinate Plackett(ρ) has union
marginal `m(ρ) = 1 − z_ρ(1−p, 1−p)`, and `h(m) > h(p)` iff `m < 1−p`, giving:

    A product Bern(p)^n obstructs the tilt class iff ρ ≤ ρ*(p),
        ρ*(p) = p(3p−1)/(1−2p)²   (ρ*(ψ) = 1: Gilmer/AHS endpoint);
    the FIXED-ρ ceiling over products is
        p*(ρ) = (4ρ−1−√(4ρ+1)) / (2(4ρ−3)),  p*(1) = ψ,  p*(∞) = 1/2.

`ρ*(0.38271) = 1.0302`: a 3% odds-ratio tilt already clears the record's
product obstruction (the tilt analogue of 002's `c_needed = 0.007` — and,
per 002's protocol, treated as a mirage until the real adversaries speak).
With ρ adaptive in p, the product obstruction of the class sits at exactly
1/2.

### 4. ADVERSARIAL TESTS (the decisive section)

All numbers are exact finite-n evaluations (float rounding only); Sinkhorn
marginal residuals < 1e-9, total masses = 1 ± 1e-5. Full tables:
`data/003_partB.json`, `003_partC.json`, `003_partC_ext.json`; log:
`data/003_full_run.log`.

**(a) Sawin geometric-mixture gadgets — the family that killed idea B.**
Sinkhorn tilt sweep, gain = `H_π(U) − H(A)` in bits:

    instance (n, ū, θ)        marg    gain_iid   best tilt gain (λ₂)
    (200, 0.40, 0.05)        0.4124     −9.11      +2.61  (2.5)
    (200, 0.42, 0.08)        0.4404    −17.78      +2.07  (6.0)   [ext]
    (300, 0.3823, 0.02)      0.3871     −2.14     +10.36  (2.5)

Every gadget is **separated** — including marginals 0.387 (record-relevant)
and 0.44. The λ-curves rise smoothly from the negative iid value, cross
zero at moderate tilt, and (0.42 instance) show the interior maximum before
the diagonal limit 0. Mechanism, confirmed against the model overlay: at
fixed λ>0 and large n the Sinkhorn coupling block-assorts the latent k
(supermodular kernel `2^{λ n p_a p_b}`), then applies Plackett(ρ) within
blocks — the tilt *finds the latent mixture structure by itself*, no oracle
access to k needed. The assortative model matches the engine within ~1 bit
at every λ ≥ 0.1 tested.

**(b) Chase–Lovett slice and smoothed slice.** Pure slice `Unif(slice pn)`
(part B; marginals automatic, no Sinkhorn needed — `Σ_B ρ^{|A∩B|}` is
slice-constant): tilt makes `|A∩B|` concentrate at chosen j, and U is then
uniform on the slice `2w₀−j` (exchangeability; brute-force-validated at
n=8). Targeting j so that `|U| ≈ n/2` sweeps the deleted middle band:

    (n, p)          gain_iid     best tilt gain     n(1−h(p))
    (1000, 0.390)     −7.73          +39.09            35.2
    (10000, 0.3823)   +1.25         +407.98           403.5
    (10000, 0.42)   −594.69         +124.34           185.5

Smoothed slice (002 part D's regularized-CL shape, `t = 2^{−8}, 2^{−10}`):
(200, w₀=79): iid −0.68 → tilt **+8.95**; (300, w₀=120, marg 0.4002): iid
−4.62 → tilt **+10.85**. Large-n version under the block coupling (slice
pairs at fixed intersection `s = 2w₀ − n/2`, so U ~ Unif(slice n/2) exactly;
Bernoulli component diagonal): (200000, p₀=.3927): **+6696** bits (002's
best KL value: −2543.6); (2000000, p₀=.3835): **+79047** (002: −3956.2).
The fixed-intersection slice coupling is NOT conditionally-iid (slice
marginals have negative coordinate correlations, not mixtures of products)
— first concrete use of the coupling class beyond C₃.

**(c) Union-closed controls** (part F, atom-level engine): four genuinely
union-closed families (n=1 pair family; power set n=3; up-set∪{∅} n=6;
28-member generator closure n=10): max gain over the λ sweep ≤ −1.6e-2 < 0
in every case (must be ≤ 0 by the licensing lemma; approaches 0 only in the
diagonal limit). Mini CL slice+top family at n=10: gain **+0.101** with 11%
of coupled mass escaping F. The functional separates fake from genuine at
every scale tested.

**(d) Fresh adversary hunt against the one-parameter recipe** (part E) —
looking for the *new* extremal instances rather than only re-testing old
ones. Model of the single-λ sweep on 2-block product mixtures (three
exactly-modeled regimes: iid; fixed λ>0 = block-assortative + Plackett
within, engine-validated to ~1 bit; λ→0⁻ = block-anti-assortative + iid
within, engine-validated sign and magnitude on the bimodal instance:
predicted separation, engine gives +55.1 bits at λ₂ = −0.05, n=300):

    cap on marginals   hardest 2-block mixture                J (best gain)
    0.39               (δ_∅-ish, Bern(.998), P_h=.003)        +0.0094  live
    0.42               (Bern(.066), Bern(.51), P_h=.797)      +0.0081  live
    0.43               (Bern(.019), Bern(.51), P_h=.837)      +0.0026  live
    0.44               (δ_∅-ish,   Bern(.51), P_h=.862)       +2e-7    live
    0.45               (δ_∅-ish,   Bern(.51), P_h=.882)       −1e-7    KILLED
    0.49               (δ_∅-ish,   Bern(.62≈φ), P_h=.79)      −7e-6    KILLED

So the **single-λ recipe has a genuine ceiling ≈ 0.445 against 2-block
mixtures** (a necessary-screen model; unmodeled λ-schedules could only
raise it). The boundary killers are "point mass at ∅ + product just above
1/2" — they defeat one shared λ because gaining on the ∅-block requires
block-anti-assortment while the ≥1/2-block tolerates only the diagonal, and
one parameter cannot do both. A two-scale coupling (anti-assort blocks,
diagonal within the heavy block) separates them at Θ(n) — trivially, but
outside the one-parameter family. Two structural notes: the hardest
high-block parameter at cap→1/2 is `p_h ≈ 0.62 ≈ φ` (the golden-ratio
mechanics resurfacing at the recipe's own extremal frontier), and none of
this affects the record-relevant regime p ≤ 0.44.

**(e) Mini-theorem (rigorous, elementary): block-adaptive couplings
separate every finite product mixture with small mixing entropy.** If
`μ = Σ_k P_k · Bern(p⃗_k)^{⊗}` then the coupling that shares k and applies
per-coordinate common-randomness mixtures achieves

    Gain ≥ Σ_k P_k Σ_i [h(min(2p_{k,i} − p_{k,i}², max(p_{k,i}, 1/2))) − h(p_{k,i})] − H(P),

every bracket being > 0 whenever `p_{k,i} < 1/2` (and = 0 otherwise); if
some marginal of μ is < 1/2 then some component coordinate has
`p_{k,i} < 1/2` with `P_k > 0`, so the gain is positive once n·(rates)
exceeds `H(P)`. Consequence: **no mixture of products with O(1)-many
components (the entire genre of known distributional counterexamples:
Ellis, Sawin, Yu-type two-point constructions) can be a counterexample to
(S-coup) for the conditionally-iid class at any threshold below 1/2.**
The `−H(P)` term is essential (μ = mixture of point masses is arbitrary),
so this is a statement about structured mixtures, not all μ.

## Outcome

**LIVE.** The dependent-couplings interface passes the full 002 protocol:

1. Both mandated adversarial families — Sawin geometric mixtures and CL
   smoothed slices — are separated by the overlap-tilt class at Θ(n) or
   better, at marginals from 0.387 to 0.44, including the exact certificate
   instances of 002 (sign flips of +127 to +83000 bits vs the KL ladder's
   values). Union-closed controls stay ≤ 0 as the licensing lemma demands.
2. The one-parameter tilt recipe (fully family-oblivious once λ is swept)
   has its first genuine obstruction at marginal cap ≈ 0.445 (2-block
   model, engine-validated) — far above the 0.38271 record. The obstruction
   is structural (one λ cannot serve two blocks with opposite needs) and is
   lifted by two-scale/block-adaptive members of the class.
3. No counterexample to (S-coup) itself was found in any tested genre;
   product mixtures are ruled out wholesale by the mini-theorem (e).

**What is NOT claimed:** no theorem, no new constant. (S-coup at p) for any
p > 0.38271 remains unproven; this cycle established that its *statement*
survives everything that killed the previous two interfaces, and located the
exact technical gaps (below). Any future claimed constant from this route
must come with the odds-ratio stability lemma and tax-controlled assembly
proved, and would be UNVERIFIED pending orchestrator adversarial review.

## Why it failed / what survived

Nothing failed at the statement level this cycle — first time in this
project. The honest inventory of what blocks a proof:

- **Gap 1 (bit level):** conditional odds-ratio control for Sinkhorn tilts
  on non-product μ (SPECULATION: λ ≥ 0 keeps all conditional odds ratios in
  `[1, 2^λ]`; FKG-flavored, plausibly provable, unproven).
- **Gap 2 (history level):** the mutual-information tax
  `Σᵢ I(A_i; B_{<i} | A_{<i})` must be beaten by the per-coordinate gains;
  no analogue of the (KEY) inequality with tax is known. Liu's conditioning
  trick avoids the tax only inside C₃; the slice separations above genuinely
  leave C₃, so a proof reaching them needs new machinery.
- **Gap 3 (recipe totality):** the evaluated functional involves sup over λ
  (and, for mixtures/slices, structured couplings). A theorem needs a single
  measurable recipe μ ↦ π(μ) with a provable lower bound for ALL μ, not per-
  genre constructions. The Sinkhorn tilt is the best candidate (it found the
  latent block structure by itself in every test), but its behavior on
  arbitrary non-exchangeable μ is unexplored.

Survived / reusable:

- `tools/uc_couplings.py`: five engines (per-coordinate Plackett calculus;
  slice-tilt closed form; exchangeable Sinkhorn tilt with exact U-profiles;
  block-adaptive large-n mixture evaluator; atom-level Sinkhorn for
  arbitrary small families), all cross-validated against each other, brute
  force, and closed forms. Checkpoints in `data/003_part[A-F]*.json`.
- Closed forms `ρ*(p) = p(3p−1)/(1−2p)²`, `p*(ρ) = (4ρ−1−√(4ρ+1))/(2(4ρ−3))`
  — the tilt-ladder analogue of 001's p(c) table, now with the correct
  interpretation attached (product-only ceiling; adaptive class reaches 1/2
  on products).
- The derivative lemma `sup_δ diag⊕π₀ works ⟺ H(ν)+D(ν‖μ) > H(μ)`:
  identifies 002's entire ladder as the first-order iid-direction shadow of
  the coupling ladder, and gives a one-line necessary condition
  (`CrossEnt(law_π(U), μ) > H(μ) for some π`) to screen any future μ.
- The mini-theorem (e) as a permanent shield: adversaries against couplings
  must be non-product-mixture (slice-like/negatively-correlated) or have
  large mixing entropy.
- The part-E killer genre (`δ_∅ ⊕ Bern(1/2+ε)`) as the new hard-instance
  family for one-parameter recipes.

## Leads generated

1. **(best, concrete) Prove the odds-ratio stability lemma** (Gap 1) for
   the Sinkhorn tilt on general μ, then attempt the tax-controlled assembly
   (Gap 2) for the restricted goal p = 0.383: by part A only ρ ≈ 1.03 of
   uniform conditional tilt is needed, and at ρ = 1+ε the tax is O(ε²)
   per coordinate while the diagonal-calibration slack is O(ε) — a
   perturbative regime where the assembly might close. Everything is
   second-order around the proven c=0 argument. First falsifiable step:
   expand `E[h(z_ρ)] − (1/2φ)(x h(y) + y h(x))` and the tax to O(ε²) at the
   AHS equality point and check the sign of the net coefficient.
2. **Characterize `Gain_all(μ) = 0`.** Necessary condition from the
   derivative lemma: `CrossEnt(law_π(U), μ) ≤ H(μ)` for every coupling π.
   Conjecture (SPECULATION): equality cases are exactly uniform
   distributions on union-closed families plus degenerate limits. Even a
   finite-n proof for exchangeable μ would be a genuine structural theorem;
   the slice analysis (every slice below n/2 has Gain > 0, approaching 0
   like (1/2)log n − O(1) near the half slice) is the boundary case data.
3. **Two-parameter tilts** `2^{λ₁|A∩B| + λ₂ f(|A|,|B|)}` (block-scale +
   coordinate-scale): kills the part-E boundary killers inside a still-
   canonical recipe; find ITS model ceiling (expect > 0.445; the question
   is whether any product-mixture ceiling short of 1/2 survives at all).
4. **Probe C₃ vs slices:** is the pure ψ-slice separated by any
   conditionally-iid coupling? The natural attempts (common planted core)
   collapse to iid-like behavior; if C₃ provably fails on slices, Liu's
   class has a ceiling strictly below the full coupling class and the
   non-PSD fixed-intersection couplings become essential — worth knowing
   before investing in C₃-style assemblies. Small-n optimization over
   mixture components is the falsifiable first step.
5. Bookkeeping for 001/002: 001's Idea C first target ("any coupling class
   whose product obstruction exceeds 0.383") is achieved at the model level
   (ρ ≥ 1.031 tilts; ceiling 0.445 for the one-parameter recipe) — but per
   this cycle, "product obstruction" was the wrong hardness measure; the
   binding obstructions are bimodal mixtures (recipe level) and the three
   proof gaps (theorem level). Future cycles should benchmark against the
   part-E genre, not products.

## References

- Liu, arXiv:2306.08824 — conditionally-iid couplings, ≈ 0.38271 (ISIT
  2024). Abstract fetched verbatim 2026-07-26; construction details [L]
  from ar5iv machine transcription.
- Sawin, arXiv:2211.11504 — Prop. 6 gadgets (re-used here as adversaries);
  the symmetric-coupling class C₂ and the beyond-ψ mixture scheme.
- Yu, arXiv:2212.00658; Cambie, arXiv:2212.12500 — 0.38234 evaluation.
- Chase–Lovett, arXiv:2211.11689 — slice family (Example 1.4), approximate
  barrier.
- Gilmer, arXiv:2211.09055; Alweiss–Huang–Sellke, arXiv:2211.11731;
  Pebody, arXiv:2211.13139; Ellis, arXiv:2211.12401.
- Plackett couplings: standard 2×2 fixed-margin odds-ratio family (any
  reference on copulas; used here purely as a parametrization).
- This repo: `attempts/union-closed/001-*.md`, `002-*.md`;
  `tools/uc_weighted_kl.py`, `tools/uc_couplings.py`,
  `attempts/union-closed/data/003_*`.
