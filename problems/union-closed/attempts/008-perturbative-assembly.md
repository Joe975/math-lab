# 008 — Perturbative assembly around product measures: Prop 6 closed, budgets measured, conditional theorem

- **Problem:** union-closed, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-07-31
- **Mode:** informed
- **Type:** proof completion (005 Prop 6's missing step) + execution of 005's
  leads 3 and 4 + computational budget census + conditional assembly theorem,
  on Gap 3(c) of the LIVE dependent-couplings route (perturbative assembly
  around the proven c = 0 argument).
- **Tools:** `explore/uc_pert.py` (written and run here; standard library
  only; fully deterministic — no RNG anywhere; runtime ~20 s; checkpoints
  `data/pert_P0_validation.json`, `pert_P1_ift.json`, `pert_P2_assembly.json`,
  `pert_P2d_constants.json`, `pert_P3_budgets.json`, `pert_P3b_anatomy.json`,
  `pert_P3c_downward.json`, `pert_P4_twoscale.json`; log
  `data/pert_full_run.log`). Command reproducing every number below:
  `python problems/union-closed/explore/uc_pert.py`
  (individual sections: `--part P0|P1|P2|P2D|P3|P3B|P3C|P4`).
- **Sources:** 003/004/005/006 (this repo; 004's and 006's corrections
  adopted where they differ from 003/005). The implicit-function-theorem
  step is proved from scratch below, so nothing rests on an external IFT
  citation. Gilmer / AHS / Plackett facts are used only through 003's
  closed forms, re-derived in the engine. No external fetches this cycle.

Notation as in 003/005: tilt coupling `π_λ(A,B) ∝ u(A)u(B) 2^{λ|A∩B|}`, both
marginals μ; `ρ = 2^λ`; histories `a = A_{<i}`, `b = B_{<i}`;
`OR_i(a,b)` the conditional odds ratio; `x₀ = 1−p`; `ψ = (3−√5)/2`;
`z_ρ(x,y)` the Plackett both-zero probability;
`ρ*(p) = p(3p−1)/(1−2p)²`; per-coordinate margin
`g(ρ,p) = h(z_ρ(x₀,x₀)) − h(p)`, positive iff `ρ > ρ*(p)`.
The perturbative class throughout is the **capped δ-ball**

    C_n(p, δ) = { μ = (1−δ)·Bern(p)^n + δ·ν :  ν a probability measure on
                  2^[n] with all element marginals ≤ p },

so every μ ∈ C_n(p,δ) has full support, all marginals ≤ p, and is inside the
hypothesis class of (S-coup at p). One calibration note up front: the task
line's "ρ ≈ 1.03" is `ρ*(0.38271) = 1.030222`; at p = 0.383 the threshold is
`ρ*(0.383) = 1.042205`, and `g(1.031, 0.383) = −0.0004 < 0`. The runs below
therefore use ρ = 1.06 and ρ = 1.20 at p = 0.383 (margins +0.000649 and
+0.005252 bits/coordinate; part P2a).

## Approach

Four tasks, in order. (T1) Close the one unproved step of 005's Prop 6 —
smoothness of the Sinkhorn potential in μ — which neither 005 nor 006
re-proved. Rather than cite an IFT, prove it from scratch with explicit
constants, because the *size* of those constants (not just their existence)
turns out to be the whole story of Gap 3(c). (T2) Run 005's lead-3
falsifiable check first. Bookkeeping note: 005's lead 3 (perturbative
assembly, sign of the net second-order coefficient) and lead 4 (the
"one-line kernel change": two-scale tilt vs the crash family) are distinct;
the task prompt conflated them, so both were run — lead 3 as the decisive
check (parts P2/P2d), lead 4 as a cheap adjunct (part P4). Lead 3 is run
*end to end on exact instances* rather than by formal ε-expansion: at n ≤ 9
every term of the assembly (per-coordinate nets, taxes, OR deviations,
margin-variances) is exactly computable from the fitted coupling, which
decides the sign question with no expansion-validity caveats. (T3) The
uniformity question for the δ² budget, which 006 S7 showed is
family-dependent: measure every budget constant across n on adversarial
directions and determine what, if anything, is n-uniform. (T4) Assemble the
conditional theorem, with each budget hypothesis labelled by the gap it
depends on.

The alternative — a purely formal second-order expansion at the AHS
equality point — was rejected because 006 S7 proved the asymptotic-onset
window is family-dependent, so an expansion alone cannot certify any
concrete δ; the exact finite-n computation can, and produces the growth
data the expansion would hide.

## What was done

### 1 (T1). The missing step of Prop 6, identified and proved

**What exactly was missing.** 005 Prop 6 proves `log₂OR_i(a,b) = λ + O(δ²)`
along smooth families `μ_δ` with `μ₀` product, *given* that the Sinkhorn
potential `u(μ)` is differentiable in μ (writing `u_δ = u₀(1 + δv + O(δ²))`).
006 re-derived the first-order cancellation but also took the smoothness as
standard. The missing step is precisely:

> existence, uniqueness, and (twice-)differentiability in μ of the symmetric
> potential u with `u(A) Σ_B 2^{λ|A∩B|} u(B) = μ(A)`, on the open simplex of
> positive μ with fixed finite support, with quantitative bounds on the
> first two derivatives.

**Lemma A (existence, uniqueness).** For every strictly positive probability
μ on a finite support S ⊆ 2^[n] and every λ ≥ 0 there is exactly one
positive u on S with `u(A)(Ku)(A) = μ(A)` for all A, where
`K(A,B) = 2^{λ|A∩B|}`; the resulting π is the unique scaling coupling with
both marginals μ.

*Proof.* Let `Ψ(v) = ½ Σ_{A,B∈S} K(A,B) e^{v_A+v_B} − Σ_A μ_A v_A` on R^S.
Critical points of Ψ are exactly the vectors v with `u = e^v` solving the
marginal equations (the second marginal agrees by symmetry of K; total mass
is automatic since Σμ = 1). Hessian: `HessΨ(v) = D(v) + Q(v)` with
`D(v)_{AA} = e^{v_A}(Ke^v)_A > 0` diagonal and
`Q(v) = D_{e^v} K D_{e^v} ⪰ 0`, because
`K = ⊗_{j≤n} [[1,1],[1,2^λ]]` restricted to S×S is a principal submatrix of
a tensor product of PSD factors (each has det `2^λ−1 ≥ 0`; this is 005/006's
W ⪰ 0 argument applied to the full kernel). So `HessΨ ⪰ D(v) ≻ 0`: Ψ is
strictly convex. Coercivity along every ray `v⁰ + tw`, `w ≠ 0`: if some
`w_A > 0` the term `½K(A,A)e^{2v⁰_A}e^{2tw_A} → ∞` beats the linear part;
otherwise all `w_A ≤ 0` with some `w_{A'} < 0` and
`−Σ μ_A v_A ≥ μ_{A'}|w_{A'}|t − O(1) → ∞` (here μ > 0 is used). A finite
convex function tending to +∞ along every ray is coercive, so a minimizer
exists; strict convexity gives uniqueness. ∎

**Lemma B (smoothness — the step itself).** On
`{μ ∈ R^S : μ > 0, Σμ = 1}` the map `μ ↦ v(μ) = ln u(μ)` is real-analytic,
with differential

    v̇ = (D_μ + Π)^{−1} μ̇,

where Π = (π(A,B))_{A,B} is the coupling matrix. Moreover `D_μ + Π ⪰ D_μ`,
so `‖(D_μ+Π)^{−1}‖₂ ≤ 1/min_A μ(A)`.

*Proof.* `F(v,μ) = ∇_v Ψ` is real-analytic in (v,μ) jointly (finite sums of
exponentials; linear in μ). At the solution,
`D(v)_{AA} = u_A(Ku)_A = μ_A` and `Q(v)_{AB} = π(A,B)`, so
`∂F/∂v = HessΨ = D_μ + Π ≻ 0` is invertible; the analytic implicit function
theorem gives a local analytic branch, and Lemma A's global uniqueness glues
the branches. ∎

**Lemma C (explicit second-order bounds along mixture segments).** Fix a
positive μ₀ on S = 2^[n] with `m₀ = min_A μ₀(A)`, and for a probability ν
let `μ_δ = (1−δ)μ₀ + δν`, `δ ∈ [0, ½]`. Then `min_A μ_δ(A) ≥ m₀/2` and,
along the solution path v(δ):

    ‖v̇(δ)‖∞ ≤ ‖v̇(δ)‖₂ ≤ ‖ν−μ₀‖₂/(m₀/2) ≤ 4/m₀,
    ‖v̈(δ)‖∞ ≤ (2/m₀)·4‖v̇‖∞² ≤ 128/m₀³.

*Proof.* Differentiate `F(v(δ), μ_δ) = 0` twice (`μ̈ = 0`, F linear in μ):
`Hess·v̇ = μ̇` and `Hess·v̈ = −(∂²F/∂v²)[v̇,v̇]`. For any h,

    (∂²F/∂v²)[h,h]_A = e^{v_A}(Ke^v)_A h_A² + 2e^{v_A}h_A(K(he^v))_A
                       + e^{v_A}(K(h²e^v))_A,

and each of the three terms is ≤ `‖h‖∞²·μ_δ(A)` (positivity of K and e^v;
middle term ≤ 2‖h‖∞²μ_δ(A)), so
`‖(∂²F/∂v²)[h,h]‖₂ ≤ 4‖h‖∞²‖μ_δ‖₂ ≤ 4‖h‖∞²`. Apply the Lemma B resolvent
bound with `min μ_δ ≥ m₀/2`. ∎

**Theorem P6′ (005 Prop 6, now fully proved, with explicit constants).**
Fix n, p ∈ (0,1), λ ≥ 0, and let `m₀ = min(p,1−p)^n`. For every probability
ν on 2^[n] and every `δ ∈ [0, ½]`, every nondegenerate history of the tilt
coupling of `μ_δ = (1−δ)Bern(p)^n + δν` satisfies

    |log₂ OR_i(a,b) − λ| ≤ c(n,p)·δ²,   c(n,p) ≤ 924·min(p,1−p)^{−3n},

uniformly over ν, i, (a,b).

*Proof.* By 005's reduction (re-verified by 006 S1),
`log₂OR − λ = log₂ R`, R the alternating cross-ratio of the four
`⟨F_α, W G_β⟩`. Write each term as
`L_{αβ}(δ) = ln⟨F_α(δ), W G_β(δ)⟩` with `F_α(δ) = exp(v(δ))` restricted to
the (a,α)-future-slice. The normalized integrand defines a probability
measure `ω_{αβ,δ}` on future pairs (x,y), and

    L′_{αβ}(δ) = E_ω[ v̇(a,α,x) + v̇(b,β,y) ],
    L″_{αβ}(δ) = E_ω[ v̈(a,α,x) + v̈(b,β,y) ] + Var_ω[ v̇⊕v̇ ],

so `|L″| ≤ 2‖v̈‖∞ + 4‖v̇‖∞²` (osc(v̇⊕v̇) ≤ 4‖v̇‖∞, Var ≤ osc²/4). At δ = 0
the potential is a product (the coupling is the product of per-coordinate
Plackett couplings, 003/004), hence the four slices `F⁰_α = c_α·H_a`,
`G⁰_β = d_β·H_b` are proportional and `ω_{αβ,0}` is the *same* measure for
all four (α,β) (the scalars cancel in the normalization). Therefore
`L′_{αβ}(0) = f(α) + g(β)` exactly, and the alternating sum
(11)+(00)−(10)−(01) annihilates it: `(ln R)′(0) = 0`. Also `R(0) = 1`
(product identity, Prop 1/A of 005). Taylor with integral remainder over
[0,δ] and Lemma C's sup-bounds give
`|ln R(δ)| ≤ (δ²/2)·4·(2‖v̈‖∞+4‖v̇‖∞²) ≤ (δ²/2)(1024/m₀³ + 256/m₀²)`,
and dividing by ln 2 and using m₀ ≤ ½ yields the constant. ∎

Two readings, both load-bearing later: (i) this closes Prop 6 — no step is
cited-but-unverified anymore — and simultaneously proves the **fixed-n
uniform budget**: at fixed n the constant is uniform over the whole
direction simplex, not just per-family (that is stronger than what 005/006
had, where 006 S7 explicitly worried about family-dependence); (ii) the
constant obtained is `exp(Θ(n))`. That is not an artifact of sloppy
bounding: part 3 below gives evidence the true worst-case constant does
grow without bound in n.

**Numerical corroboration (part P1).** At n = 5, p = 0.3, λ = 1, direction
ν = Bern(0.05)^5: the first-order formula `v̇ = (D_μ+Π)^{−1}μ̇` (dense
32×32 Gaussian elimination, solve residual 7.8e-16) predicts the fitted
potential with quadratic remainder: `‖v(δ) − v(0) − δv̇‖∞/δ²` =
6.18, 6.03, 5.77, 5.31 at δ = 0.01, 0.02, 0.04, 0.08 — a clean bounded
ratio, as Lemma C demands. Engine validation (part P0): product identity
`OR ≡ 2^λ` to 2.0e-15 with gain matching the closed form
`n(h(z_ρ(x₀,x₀)) − h(p))` to 4.7e-14 and per-coordinate taxes ≤ 5.6e-16;
005's part-G δ-sweep reproduced from this independent implementation to
3.6e-15 against `data/005_partG.json`; crash identity `log₂OR = λ(3−n)`
exact at n = 5, 8.

### 2 (T2). 005's lead-3 check, run end to end — and lead 4's kernel change

**The assembly bookkeeping** (exact, no expansion; this is what the engine
computes). For any μ and its tilt π, chain rule gives
`Gain ≥ Σ_i N_i` with per-coordinate net

    N_i = E_π[h(z_i)] − H(A_i|A_{<i})
        = E_π[ G(x̃_i, ỹ_i, τ_i) ] − T_i,

    G(x,y,τ) = h(z_{2^{λ+τ}}(x,y)) − ½h(x) − ½h(y),
    τ_i = log₂OR_i − λ,
    T_i  = ½ I(A_i;B_{<i}|A_{<i}) + ½ I(B_i;A_{<i}|B_{<i})   (the tax),

using the identity `H(A_i|A_{<i}) = E_π[h(x̃_i)] + I(A_i;B_{<i}|A_{<i})`
(x̃_i, ỹ_i the conditional zero-margins of the (A_i,B_i) table given both
histories). At the product base point every history has
`(x̃,ỹ,τ) = (x₀,x₀,0)`, `T_i = 0`, and `N_i = g(ρ,p) > 0` for ρ > ρ*(p).
Within the capped δ-ball the multiplicative AHS machinery is *unnecessary*:
the additive margin g does the work, because conditional margins cannot
drift Θ(1) in mass. This is why the perturbative regime is tractable at all.

**The AHS-point coefficient** (003 lead 1's "first falsifiable step",
part P2a): `dg/dρ` at (ρ, p) = (1, ψ) equals `h'(ψ)·φ⁶ = 0.038689` in
closed form (`dz_ρ/dρ|_{ρ=1} = xy(1−x)(1−y)`, re-derived; engine numeric
0.038647 at step 1e-6). Positive — the tilt margin switches on linearly in
(ρ−1) at the AHS equality point, tax and OR-deviations exactly zero at the
product. So the sign question is decided *off* the product, by the δ-terms,
which the end-to-end runs measure.

**End-to-end runs (parts P2b/P2c).** Grid: p = 0.383; ρ = 1.20
(g = 5.25e-3) at n = 5..8 and ρ = 1.06 (g = 6.5e-4) at n = 6; directions
ν ∈ {crash family (marginals ≤ 0.37), δ_∅, Bern(0.05)^n} (P2b) and the
cap-stress set {½Bern(0.083)+½Bern(0.683) ("pm03", marginals exactly p),
Unif(slice round(pn)), zero-prefix embedded crash (below)} (P2c);
δ ∈ {0.02..0.2}. Results, all 78 runs:

- **Every in-class run is positive**: min per-coordinate net ≥ base margin
  in all capped runs (the marginal-slack first-order term
  `G_x·(p − p̄_i) ≥ 0`, `G_x(center) = +0.73`, helps), chain-rule bound
  `Σ N_i > 0`, true gain > gain_lb > 0, no OR below 1 anywhere
  (`mass_or_below_1 = 0` in every run), taxes ≤ 2.5e-4 bits/coordinate even
  at δ = 0.2.
- **The only negative runs are out of class**: Unif(slice 3) at n = 7 has
  marginals 3/7 = 0.4286 > p; there `gain_lb = −0.0066` (δ=0.1) and
  `−0.0436` (δ=0.2), and the loss is quantitatively the first-order term:
  excess marginal 0.0456·δ times `G_x+G_y ≈ 1.47` ≈ the observed net drop.
  The marginal cap is load-bearing, exactly as the c = 0 argument's
  structure predicts — and (S-coup at p)'s hypothesis provides it for free.

**Verdict on lead 3: the check PASSES at every tested (n, δ) in the capped
class.** The assembly does not die at second order; nothing in the measured
range even gets close (worst in-class net = +0.0003, slice n=5 δ=0.2,
still positive). The route restated by 005 survives its assembly test.

**Lead 4's one-line kernel change (part P4).** Two-scale kernel
`2^{λ(|A∩B| − |A||B|/n)}` ("centered overlap", the f(|A|,|B|) = |A||B|/n
member of 003 lead 3's family with λ₂ = −λ₁). The crash history's OR is
kernel-only (006's cancellation argument is kernel-independent), giving the
closed form

    log₂ OR_crash = λ·[(3−n) − (4n−3−n²)/n] = λ(3−n)/n  ∈ [−λ, 0),

engine-exact to 2.8e-16 (n = 5: −0.400; n = 8: −0.625). **The centered
kernel suppresses the crash from `2^{λ(3−n)}` (unbounded) to
`≥ 2^{−λ}` (bounded), by construction, for every n.** The mirror family's
top OR drops from `2^{λ(n−1)}` to `2^{0.8λ}` (n=5). On products the ORs are
no longer exactly `2^λ` but sit in a narrow band (`[0.827, 0.833]` at n=6,
`[0.870, 0.875]` at n=8, in log₂ units — an effective tilt
`λ_eff ≈ λ(1 − p²·(1+o(1)))`) with positive gains (+0.146 at n=8). A
definite outcome for lead 4: the two-parameter family *does* tame the crash
genre — see lead 2 below.

### 3 (T3). Uniformity of the δ² budget across n — the decisive census

Part P3 measures, for directions ν ∈ {crash, zero-prefix embedded crash
(k = 2), δ_∅, Bern(0.7)^n} at (p, λ) = (0.30, 1.0), n = 4..9,
δ ∈ {0.0125, 0.025, 0.05, 0.1}: the pointwise deviation
`max_hist |log₂OR − λ|`, its mass-weighted per-coordinate mean, the
per-coordinate tax, and `E[(x̃−x₀)²]`. The zero-prefix embedded crash
(atoms avoiding coordinates 1..k, crash structure on the tail) keeps all
marginals ≤ 0.37, so it is a *legal* direction for the capped class while
planting structure below a zero-history. Findings:

**(a) The two-sided pointwise budget is NOT n-uniform.** Crash direction,
`max|log₂OR−λ|` at fixed δ = 0.1: 0.019, 0.099, 0.251, 0.525, 0.818, 1.137
at n = 4..9 — increments 0.080, 0.152, 0.274, 0.293, 0.320, i.e. tending to
**≈ +0.3 bits per added coordinate at fixed δ** (and ≈ +0.34 at δ = 0.05
from part P3b). Any budget claim `|log₂OR−λ| ≤ c·δ² for all δ ≤ δ₀, all n`
with (c, δ₀) independent of n is violated once `a·n > c·δ₀²`. EVIDENCE at
n ≤ 9; the continuation of the linear trend beyond n = 9 is SPECULATION,
but it is consistent in both regimes and bounded above by 005 Prop 4's cap
`λ(1+m)`, which itself grows linearly. So Theorem P6′'s exp(Θ(n)) constant
is qualitatively honest: c(n) → ∞, and the fixed-n statement is the most
that is true of the naive budget. Mechanism (SPECULATION, consistent with
part P3b's anatomy): the worst histories are *diagonal* prefixes of the
planted heavy atom (a = b = prefix of S₂ = {3..n}; top history at n = 9 is
`a = b = 000011100`, log₂OR = +1.97, mass 3.5e-3) — Prop 2's
Cauchy–Schwarz excess fed by slice non-proportionality that deepens with
the atom's remaining length; the δ²-coefficient inherits that depth.

**(b) The blow-up is entirely on the harmless side.** Part P3b: in every
capped run, zero mass has OR < 1, and the worst *downward* deviation is
minuscule (worst below-λ history at n = 9, δ = 0.05: log₂OR = +0.985,
i.e. 0.015 below λ, at the i = 2 crash history, mass 0.179). All the
growth in (a) is *upward* (OR ≫ 2^λ) — and upward deviations do not hurt
the assembly until τ exceeds `τ_half ≈ 3.0` (where the coordinate turns
effectively diagonal and the margin, not positivity, is lost); no measured
history comes near that. This is 005 part E's "deviation budget sits above
2^λ" and 006 S8's averaged-control survival, now confirmed *inside the
perturbative ball on adversarial directions*.

**(c) What the assembly actually needs, measured (part P3c).** The
downward one-sided mass-weighted budget
`E_π[(λ − log₂OR_i)_+]` per coordinate, crash direction, δ = 0.05:

    max over i (always i = 2):   /δ² = 0.166, 0.248, 0.331, 0.413, 0.496,
        0.579 at n = 4..9  — ≈ 0.0827·(n−2): LINEAR in n (ρ=1.2 regime);
    average over i:              /δ² = 0.041, 0.050, 0.055, 0.059, 0.063,
        0.065 — increments 0.0083, 0.0056, 0.0041, 0.0032, 0.0026,
        decreasing by ~25% per step: consistent with a bounded limit
        (≈ 0.07–0.08) but not provably convergent from six points.

Same shapes at (p, λ) = (0.30, 1.0) with larger values (max linear
≈ 0.36·(n−2.5)·δ², average flattening toward ≈ 0.3·δ²). So: the
**max-coordinate downward budget is also not n-uniform** (it concentrates
at the crash coordinate i = 2 and grows linearly), but the **per-coordinate
average** — which is what the summed assembly consumes — shows decreasing
increments, consistent with an n-uniform constant. That average bound is
exactly the perturbative form of the restated Gap 1
(`averaged-odds-ratio-control`): EVIDENCE for it at n ≤ 9, SPECULATION
beyond.

**(d) The other budgets.** Taxes per coordinate at δ = 0.1 (crash, λ=1):
3.3e-4 → 1.9e-3 across n = 4..9 — growing (×5.9 over a ×2.25 range of n,
i.e. superlinear in n at these sizes, though decelerating), but three
orders of magnitude below the per-coordinate entropy scale; in the
assembly regime (ρ = 1.2) taxes stay ≤ 2.1e-4 even at δ = 0.2 (part P2c);
tax/δ² ≲ 0.01 throughout. Margin-variance `E[(x̃−x₀)²] ≈ 2δ²`,
n-saturating. Neither is the binding constraint in the measured range; the
tax's growth shape is the part to re-measure at larger n (gap b's
perturbative form).

### 4 (T4). The conditional assembly theorem

**Theorem (conditional, explicit).** Fix p ∈ (ψ, ½), ρ = 2^λ > ρ*(p),
margin g = g(ρ,p) > 0, and s₀ > 0 with [x₀−s₀, x₀+s₀] ⊂ (½,1). Define, for
the class C_n(p,δ) and its tilt couplings, the budget hypotheses (all "per
coordinate, averaged over i = 1..n"):

    (B1) avg_i E_π[(λ − log₂OR_i)_+]      ≤ c₁ δ²    [gap a′: averaged
         odds-ratio control, perturbative downward form]
    (B2) avg_i T_i                        ≤ c₂ δ²    [gap b: tax control,
         perturbative form]
    (B3) avg_i E[(x̃_i−x₀)² + (ỹ_i−x₀)²]  ≤ c₃ δ²    [margin-variance]
    (B4) avg_i P_π(τ_i > τ_half)          ≤ c₄ δ²    [upward tail beyond
         τ_half, the τ where z crosses ½ at the box center]

If (B1)–(B4) hold with constants independent of n, then for every n and
every `δ ≤ δ₀ = sqrt( g / 2(M_τ c₁ + c₂ + (M₂ + ‖G‖_box/s₀²) c₃ +
(g + 2L₂s₀) c₄) )`, every μ ∈ C_n(p,δ) has

    Gain(μ) ≥ Σ_i N_i ≥ n·g/2 > 0,

i.e. (S-coup at p) holds on the δ₀-ball uniformly in n. Constants:
`M_τ = sup|∂G/∂τ|`, `M₂ =` sup negative curvature of (x,y) ↦ G at τ = 0,
`L₂ = sup|∂G/∂x|` (all over the box × [0, τ_half]), `‖G‖_box` the sup of
|G| there; at (p, ρ, s₀) = (0.383, 1.20, 0.10) these evaluate (part P2d,
dense finite differences) to

    g = 0.005252, τ_half = 3.017, M_τ = 0.059, M₂ = 2.252, L₂ = 1.022,
    worst huge-τ dip on the box = 0.0071.

*Proof outline* (each step elementary; the one place needing care is
flagged): split each E_π[G] over the good set
Γ_i = {|x̃−x₀| ≤ s₀, |ỹ−x₀| ≤ s₀, τ_i ≤ τ_half} and its complement. On Γ:
second-order Taylor of G in (x,y) at (x₀,x₀,0) — the first-order term
`(G_x+G_y)·E[x̃−x₀]` is ≥ 0 because `E_π[x̃_i] = 1−p̄_i ≥ x₀` under the
marginal cap and `G_x = G_y > 0` (part P2d) — plus the τ-loss
`M_τ E[(τ)_−]` (G is nondecreasing in τ while z < ½, i.e. up to τ_half);
quadratic remainder ≤ M₂·(B3). Off Γ: for the margin-tail part, Markov
with (B3) at cost `‖G‖_box/s₀²·c₃δ²`; for the τ-tail part, G stays ≥
−(dip) ≥ −2L₂s₀-ish on the box even as τ → ∞ (h is unimodal, so along the
τ-path `h(z) ≥ min(h(z_0), h(min(x,y)))`, and the huge-τ limit is the
diagonal value `h(min(x,y)) − ½h(x) − ½h(y) ≥ −L₂|x−y|`), so the loss per
unit mass is ≤ g + dip, charged to (B4). Taxes subtract directly: (B2).
Summing over i converts every avg-budget to a total, and δ ≤ δ₀ leaves
n·g/2. ∎ [The flagged step: "G nondecreasing in τ up to τ_half" uses
z(x,y,τ) < ½ on the whole box for τ ≤ τ_half, which holds for s₀ = 0.10 at
this (p,ρ) by the monotonicity of z in (x,y,τ) — checked numerically on
the P2d grid; a purist should re-derive the box-uniform τ_half analytically.]

**Instantiation with the measured constants** (c₁, c₂, c₃, c₄) =
(0.07, 0.01, 2.5, 0) from parts P3c/P2c/P3: loss rate 5.64·δ², hence
`δ₀ ≈ 0.022` (part P2d) — dominated by the crude M₂c₃ variance term; the
end-to-end runs (which see the true, partly-positive quadratic term rather
than its absolute-value bound) stay positive to δ = 0.2, an order of
magnitude beyond. So the theorem's constants are pessimistic but the
structure is right.

**Dependency ledger (nothing silently assumed):**
- (B1) with n-uniform c₁ is open: it is the perturbative shadow of gap a′
  (`averaged-odds-ratio-control`). Supported by P3c's flattening average
  (EVIDENCE, n ≤ 9); its max-coordinate version is FALSE with n-uniform
  constant (linear growth, part 3(c)) — so any proof of (B1) must exploit
  the average over coordinates, not per-coordinate uniformity.
- (B2) with n-uniform c₂ is open: perturbative shadow of gap b
  (`mutual-information-tax`). Measured values tiny but growing
  superlinearly in n at small n; shape undetermined.
- (B3), (B4) are new, purely perturbative budgets; both look easiest —
  (B4) is vacuous (c₄ = 0) in the entire measured range — but neither is
  proved. SPECULATION: (B3) should follow from a Sinkhorn TV-stability
  estimate `‖π(μ_δ) − π(μ₀)‖_TV ≤ Cδ` with C n-uniform; no such estimate
  is known for this kernel (Hilbert-metric contraction gives C growing
  with the kernel diameter λn).
- Gap 3 (recipe totality) is untouched: everything here is for the fixed-λ
  Sinkhorn tilt recipe on a δ-ball, which is a *class*, not all μ.

## Outcome

**LIVE, with one sub-result at proof level and one refutation.** Scope:

1. **VERIFIED (proof supplied, pending skeptic pass):** the missing
   smoothness step of 005 Prop 6 — existence, uniqueness, and real-analytic
   μ-dependence of the symmetric tilt potential, with the explicit
   first/second-derivative bounds of Lemmas B/C — and hence Theorem P6′:
   `|log₂OR − λ| ≤ 924·min(p,1−p)^{−3n}·δ²` on the full mixture ball at
   fixed n, uniform over directions. Per repo rules this is a result only
   after an independent skeptic pass; the proof is written to be checked
   line by line, and part P1 corroborates the derivative formula
   numerically.
2. **EVIDENCE (n ≤ 9, exact computation):** 005's lead-3 check passes —
   the perturbative assembly is positive, end to end, on every tested
   in-class instance (78 runs, three ordinary and three adversarial
   directions, δ up to 0.2, ρ ∈ {1.06, 1.2}), with zero mass at OR < 1.
3. **REFUTED (as a route ingredient, by EVIDENCE + the Prop-4 mechanism):**
   the *naive* n-uniform two-sided pointwise budget `|log₂OR−λ| ≤ cδ²`.
   Deviation at fixed δ grows ≈ linearly in n (crash direction; +0.3
   bits/coordinate at δ = 0.1). Its max-coordinate one-sided downward
   version dies too (≈ 0.083(n−2)δ²). What survives observationally is the
   coordinate-averaged downward budget (flattening toward ≈ 0.07δ² in the
   assembly regime).
4. **Conditional theorem** (T4 above): budgets (B1)–(B4) with n-uniform
   constants imply (S-coup at p) on an explicit n-independent δ₀-ball
   around every Bern(p)^n, p up to p*(ρ); with measured constants,
   δ₀ ≈ 0.022 at (p,ρ) = (0.383, 1.20).

**Not claimed:** no unconditional bound and no new constant for Frankl; no
claim that (B1)–(B4) hold for general n (each is labelled, with its gap);
no claim about μ outside the mixture class C_n(p,δ) (in particular nothing
about balls in other metrics, non-product base points, or the recipe's sup
over λ); the linear-growth readings of part 3 are trends at n ≤ 9, not
limits; Theorem P6′'s constant is an upper bound believed wildly loose at
fixed n (measured c ≈ 0.6 vs bound 6.4e10 at n = 5, p = 0.3); the P2d
constants rest on dense finite differences, not interval arithmetic; and
none of this has had its independent skeptic pass yet.

## Why it failed / what survived

Nothing failed at the level this cycle attacked — the assembly survives its
own falsifiable test — but the cycle located exactly where the perturbative
route must next be attacked, and one route ingredient died:

- **The step that breaks (for any pointwise route):** the δ²-coefficient of
  the OR deviation is not uniform in n. The quantity that goes the wrong
  way is the Cauchy–Schwarz excess on diagonal prefixes of a planted heavy
  atom — it grows with the atom's remaining depth (+0.3 bits/coordinate at
  δ = 0.1), so `sup_hist |log₂OR−λ|/δ² → ∞` with n even inside the capped
  ball. Prop 6 cannot be used pointwise-uniformly; every future assembly
  must consume OR control through a mass-weighted, coordinate-averaged,
  **one-sided (downward)** functional. This narrows gap a′ usefully: prove
  `avg_i E[(λ−log₂OR_i)_+] ≤ c δ²` with c independent of n — the upward
  side provably cannot matter below τ_half, and the downward side is where
  all remaining danger lives.
- **Why the fixed-n statement is trivial-adjacent and the n-uniform one is
  the whole game:** Gain is continuous in μ and Gain(Bern(p)^n) = n·g > 0,
  so *some* positive-Gain ball exists at each n by pure compactness. The
  route's value is a ball radius independent of n; that is exactly what
  (B1)–(B4) buy and what nothing weaker buys.

Survived / reusable:

- Lemmas A/B/C + Theorem P6′: permanent structure theory (strict convexity
  of the scaling functional via K ⪰ 0; the resolvent bound
  `‖(D_μ+Π)^{−1}‖ ≤ 1/min μ`; the probability-measure trick that makes
  L′, L″ weighted means/variances of v̇, v̈ — reusable for any kernel of
  the form ⊗ PSD factors, including the centered kernel of part P4).
- The per-coordinate net decomposition `N_i = E[G(x̃,ỹ,τ)] − T_i` with
  G's calculus (τ_half, M_τ, M₂ at any (p,ρ,s₀) via part P2d) — the
  assembly's bookkeeping in a form where every gap plugs in as a budget.
- The zero-prefix embedded crash family: an in-class (marginal-capped)
  adversarial direction genre for any future perturbative claim.
- The centered-overlap kernel `2^{λ(|A∩B| − |A||B|/n)}` with the closed
  form `log₂OR_crash = λ(3−n)/n`: the first recipe member known to bound
  the crash genre's downside uniformly in n.
- `explore/uc_pert.py`: validated against 005's checkpoints, closed forms,
  and the crash identity; computes every assembly term exactly for
  arbitrary μ at n ≤ 9 in seconds.

## Leads generated

1. **(sharpest) Prove budget (B1) in average form.** Target:
   `avg_i E_π[(λ − log₂OR_i)_+] ≤ c(p,λ)·δ²` on C_n(p,δ), c independent of
   n. Falsifiable first step: extend part P3c to n = 10..12 (engine
   handles n = 10 in ~1 min; 11–12 need a profile-space or sparse-support
   variant) and fit the increments — if they stop decreasing, (B1) as
   stated dies and the theorem needs a log n factor (check whether
   δ₀ ~ 1/√log n still beats the compactness-trivial statement; it does).
   Proof direction (SPECULATION): sum 005 Prop 2 over the diagonal (where
   `(λ−log₂OR)_+ = 0` identically) and control the off-diagonal downward
   mass by the crash-history mechanism, which P3c shows is concentrated at
   the single coordinate where the planted structure separates — the
   average then dilutes the one bad coordinate by 1/n, matching the
   observed `max_i ≈ 0.083(n−2)δ²` vs `avg_i → const·δ²` split.
2. **Centered-kernel odds-ratio floor.** Part P4: the kernel
   `2^{λ(|A∩B|−|A||B|/n)}` has crash-OR `2^{λ(3−n)/n} ≥ 2^{−λ}` and
   band-concentrated ORs on products. Conjecture (SPECULATION):
   `OR ≥ 2^{−λ·C}` pointwise for the centered tilt on ALL μ, some absolute
   C. Falsifiable: re-run 005's part-D vertex-extremal analysis for the
   centered kernel — the Prop-4 extremal configurations are point-mass
   slices, so the sharp range is `λ·extrema of [⟨x₁−x₀,y₁−y₀⟩ −
   (|A|-size combination)/n]` — compute whether the size term caps the
   overlap term for 4-atom supports (a finite optimization). If yes, the
   restated Gap 1 should be attacked in the centered class, where a
   pointwise (not just averaged) floor may actually be true.
3. **Sinkhorn TV-stability (feeds B3).** Prove or refute:
   `‖π_λ(μ) − π_λ(μ₀)‖_TV ≤ C(p,λ)·‖μ−μ₀‖_TV` on C_n(p,δ) with C
   independent of n. Falsifiable numerically (compute the TV ratio across
   n for the crash direction; the engine's pi_rows make this a ten-line
   addition). A positive answer gives (B3) and probably (B2); Hilbert-metric
   arguments give only exp(λn) — the question is whether marginal-capped
   mixtures do better.
4. **Tax shape (feeds B2).** Part P3's per-coordinate tax grows
   superlinearly at n ≤ 9 while staying ≲ 0.01δ². Measure `avg_i T_i`
   (not max) at n = 6..12 and fit: if avg tax/δ² flattens like the
   downward OR budget, (B2) is plausible and 009's Gap-2 work applies
   directly here; if it keeps growing, the perturbative assembly needs the
   conditioning trick (Liu's C₃ bookkeeping) even inside the ball —
   a definite fork either way.
5. **Beyond product base points.** The whole analysis perturbs around
   Bern(p)^n. The published extremal families (Sawin mixtures) sit at
   product *mixtures*; redo P2/P3 around a 2-block mixture base (the
   diag⊕block coupling of 003 part D as base coupling) — the first-order
   cancellation of Theorem P6′ fails at non-product bases (slices are no
   longer proportional), so measure how fast the deviation onset moves
   from δ² to δ: that number decides whether "perturbative" can ever mean
   "perturbative around the known extremals" rather than around products.

## References

- This repo: `problems/union-closed/attempts/003-dependent-couplings.md`
  (interface, recipe, Gap 3; leads 1, 3);
  `004-skeptic-review-of-003.md` (corrections used: 0.431496 ceiling,
  restated mini-theorem, half-mixing coupling);
  `005-odds-ratio-control-refuted.md` (Props 1–6, crash family, leads 3–4);
  `006-skeptic-review-of-005.md` (S1 reduction re-check, S7
  family-dependent onset — the constraint this cycle's T3 was designed
  around, S8).
- `explore/uc_pert.py`; checkpoints `data/pert_P0_validation.json`,
  `pert_P1_ift.json`, `pert_P2_assembly.json`, `pert_P2d_constants.json`,
  `pert_P3_budgets.json`, `pert_P3b_anatomy.json`, `pert_P3c_downward.json`,
  `pert_P4_twoscale.json`; log `data/pert_full_run.log`; prior data
  consumed: `data/005_partG.json` (reproduction target), `005_partE.json`
  (context).
- Gilmer arXiv:2211.09055; Alweiss–Huang–Sellke arXiv:2211.11731 (the c = 0
  argument and equality point, used via 003's closed forms); Liu
  arXiv:2306.08824 (record constant, context only). No external fetches
  this cycle; the implicit-function-theorem material is proved inline
  rather than cited.
