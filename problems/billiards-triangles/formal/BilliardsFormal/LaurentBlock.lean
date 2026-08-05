/-
The Laurent-identity block of attempt record 005
(`problems/billiards-triangles/attempts/005-complete-death-law-theorem.md`,
sections (ii)-A Steps 2-4), reviewed and confirmed in 008 (sections 1-2):
the closed-form composition of the half word of W(a,b) and the identities
I1, I2, I3, I4, D1, D2 plus the glide facts, with the family parameters
(a, b) carried by FORMAL variables, so that one kernel-checked proof covers
every integer pair (a, b) at once (005's formal-parameter specialization).

## The encoding (the formal statement IS this model — read before trusting)

005 works in the Laurent ring Q(i)[x^±, y^±, P^±, Q^±] where, on the torus,
x = e^{iα}, y = e^{iβ}, P = e^{iaα}, Q = e^{ibβ}.  Here that ring is
represented by an arbitrary field `K` with distinguished elements
`x y P Q i : K` subject only to `x ≠ 0`, `y ≠ 0`, `P ≠ 0`, `Q ≠ 0`,
`i ≠ 0`, `(2 : K) ≠ 0`.  Note `i` needs no relation `i² = -1`: every
identity in the block turned out to be homogeneous in `i` (each side has
the same total `1/(2i)`-degree), so it holds for ANY nonzero `i` — a
strict generalization of 005's Q(i) setting, discovered because `ring`
closed every goal without ever consuming `i² = -1`.  Every theorem is
universally quantified over `K` and these generators, so instantiating
`K := FractionField` of the Laurent ring recovers the formal identity, and
instantiating `K := ℂ`, `x := exp (I*α)`, `y := exp (I*β)`,
`P := exp (I*a*α)`, `Q := exp (I*b*β)`, `i := Complex.I` recovers the
geometric specialization.  Consequently:

* a monomial `x^ex * y^ey * P^eP * Q^eQ` stands for the unit
  `e^{i(ex·α + ey·β + eP·aα + eQ·bβ)}`; the (a,b)-dependence lives in the
  formal generators `P`, `Q` only (005's exponent-lattice move);
* formal sine and cosine of the angle whose unit is `z` are
  `sinF z i = (z - z⁻¹) / (2*i)` and `cosF z = (z + z⁻¹) / 2`;
* complex conjugation is the substitution
  `(x, y, P, Q, i) ↦ (x⁻¹, y⁻¹, P⁻¹, Q⁻¹, -i)` — on the torus all four
  generators are units of modulus 1, so inversion is conjugation
  (008 §2's star involution).  A geometric point is carried as a PAIR
  `(z, conj z) : K × K`, and each map on pairs has as second component the
  conjugate substitution of its first (an invariant that is syntactically
  visible in the definitions below, and proved for the composite maps);
* `Im ζ = (ζ - conj ζ) / (2*i)`, so the glide-axis projection is
  `proj … v vc = (conj δ * v - δ * vc) / (2*i)` with `vc = conj v`.

The definitions `B, C, A1, mu, delta, w, tau, m, v2` below transcribe
`closed_forms()` of `explore/plaw_general.py` (005's proof object,
independently re-derived in 008) monomial by monomial.

## What is proved (tiers as in the queue-15 task)

* Tier A — `pair_collapse_A`, `pair_collapse_B`, `half_word_closed_form`:
  R₁R₂ = Rot_A(2α), R₀R₂ = Rot_B(−2β) (005 Step 2), and
  R₀ ∘ Rot_A ∘ Rot_B = (z ↦ μ·conj z + w) with μ, w the closed forms
  (005 Step 3).  `iterate_collapse_A/B` and `half_word_letters` extend this
  to the genuine letter word R₀ (R₁R₂)^a (R₀R₂)^b for every a b : ℕ, with
  P = x^a, Q = y^b appearing as actual powers.
* Tier B — `mu_eq_delta_sq`, `tau_parallel_axis`, `mu_mul_conj_mu`,
  `half_word_sq_is_translation`, `glide_action`, and `identity_I1 … I4`,
  `identity_D1`, `identity_D2`: 005's Step 4 identity list, verbatim.
* Tier C — the `*_spec` corollaries: specialization to arbitrary integer
  exponents `P := x^a`, `Q := y^b` (`a b : ℤ`).  In this encoding 005's
  ring-homomorphism specialization is instantiation of universally
  quantified variables, so each corollary is a one-line application.

## What is NOT proved (the bridge)

The bridge (geometry = ring model) is explicitly OUT of scope — 005/008
carry it.  Concretely, the following remain informal, cross-checked in
005 (23-member exact specialization, 400-point float check) and 008
(hand re-derivation + exact off-torus recomposition, 40 members):
the circumdiameter-1 normalization A = 0, B = sin(α+β), C = sin β·e^{iα};
that `reflAB`, `reflCA`, `reflBC` are the reflections in the sides AB, CA,
BC of that triangle; and that the unfolding corridor quantities of the
harness are the `proj`/`m` data below.  A compiling proof here says nothing
about billiards until that mapping is reviewed against 005 §Step 1-3.

Zero `sorry`; no axioms beyond mathlib (audit in `lake-build.log`).
Toolchain pinned in `lean-toolchain` (leanprover/lean4:v4.32.2), mathlib
pinned in `lakefile.toml` / `lake-manifest.json` (tag v4.32.2).
-/
import Mathlib

namespace BilliardsFormal
namespace Laurent

variable {K : Type*} [Field K]

/-! ### The formal trigonometric dictionary -/

/-- Formal sine of the angle whose unit is `z`: on the torus
`z = e^{iθ}` this is `(e^{iθ} - e^{-iθ})/(2i) = sin θ`. -/
def sinF (z i : K) : K := (z - z⁻¹) / (2 * i)

/-- Formal cosine of the angle whose unit is `z`. -/
def cosF (z : K) : K := (z + z⁻¹) / 2

/-! ### Vertices and glide data (005 Step 3 / `plaw_general.closed_forms`) -/

/-- Vertex `B = sin(α+β)` (circumdiameter-1 normalization; real). -/
def B (x y i : K) : K := sinF (x * y) i

/-- Vertex `C = sin β · e^{iα}`. -/
def C (x y i : K) : K := sinF y i * x

/-- `A₁ = R₀(A) = B · (1 - y⁻²)`, the A-fan pivot. -/
def A1 (x y i : K) : K := B x y i * (1 - y⁻¹ ^ 2)

/-- Antilinear multiplier of the half word: `μ = e^{-2iaα + 2i(b-1)β}
= y⁻² P⁻² Q²`. -/
def mu (y P Q : K) : K := y⁻¹ ^ 2 * P⁻¹ ^ 2 * Q ^ 2

/-- Glide-axis direction `δ = e^{i(-aα + (b-1)β)} = y⁻¹ P⁻¹ Q`;
`δ² = μ`. -/
def delta (y P Q : K) : K := y⁻¹ * P⁻¹ * Q

/-- Translation part of the half word:
`w = B · (1 - y⁻² + y⁻²P⁻² - μ)`. -/
def w (x y P Q i : K) : K :=
  B x y i * (1 - y⁻¹ ^ 2 + y⁻¹ ^ 2 * P⁻¹ ^ 2 - mu y P Q)

/-- Glide translation `τ = w + μ · conj w` (conjugation = substitution). -/
def tau (x y P Q i : K) : K :=
  w x y P Q i + mu y P Q * w x⁻¹ y⁻¹ P⁻¹ Q⁻¹ (-i)

/-- Projection onto the normal of the glide axis:
`p(v) = Im(conj δ · v) = (conj δ · v - δ · conj v)/(2i)`, with the
conjugate `vc = conj v` supplied explicitly (pair convention). -/
def proj (y P Q i v vc : K) : K :=
  (delta y⁻¹ P⁻¹ Q⁻¹ * v - delta y P Q * vc) / (2 * i)

/-- Axis offset `m = ½ · Im(conj δ · (w - τ/2))`. -/
def m (x y P Q i : K) : K :=
  proj y P Q i (w x y P Q i - tau x y P Q i / 2)
    (w x⁻¹ y⁻¹ P⁻¹ Q⁻¹ (-i) - tau x⁻¹ y⁻¹ P⁻¹ Q⁻¹ (-i) / 2) / 2

/-- The C-image endpoint of gate `2a+2`:
`v₂ = R₀(Rot_A(2aα) C) = B + y⁻²(P⁻² · conj C - B)` (005's I4 anchor). -/
def v2 (x y P i : K) : K :=
  B x y i + y⁻¹ ^ 2 * (P⁻¹ ^ 2 * C x⁻¹ y⁻¹ (-i) - B x y i)

/-! ### The maps, on pairs `(z, conj z)`

Each second component is the conjugate substitution of the first —
generators inverted, `i` negated, components swapped (for the
orientation-reversing reflections) or kept (for the rotations). -/

/-- `R₂`: reflection in side AB (the real axis): `z ↦ conj z`. -/
def reflAB (p : K × K) : K × K := (p.2, p.1)

/-- `R₁`: reflection in side CA: `z ↦ x² · conj z` (`x² = e^{2iα}`). -/
def reflCA (x : K) (p : K × K) : K × K := (x ^ 2 * p.2, x⁻¹ ^ 2 * p.1)

/-- `R₀`: reflection in side BC: `z ↦ B + y⁻² (conj z - B)`
(`y⁻² = e^{-2iβ}`, and `conj B = B`). -/
def reflBC (x y i : K) (p : K × K) : K × K :=
  (B x y i + y⁻¹ ^ 2 * (p.2 - B x y i),
   B x⁻¹ y⁻¹ (-i) + (y⁻¹)⁻¹ ^ 2 * (p.1 - B x⁻¹ y⁻¹ (-i)))

/-- `Rot_A(2aα)`: rotation about `A = 0` by `2aα`: `z ↦ P² z`. -/
def rotA (P : K) (p : K × K) : K × K := (P ^ 2 * p.1, P⁻¹ ^ 2 * p.2)

/-- `Rot_B(-2bβ)`: rotation about `B` by `-2bβ`: `z ↦ B + Q⁻²(z - B)`. -/
def rotB (x y Q i : K) (p : K × K) : K × K :=
  (B x y i + Q⁻¹ ^ 2 * (p.1 - B x y i),
   B x⁻¹ y⁻¹ (-i) + (Q⁻¹)⁻¹ ^ 2 * (p.2 - B x⁻¹ y⁻¹ (-i)))

/-! ### Tier A: pair collapse and the closed-form composition (005 Steps 2-3) -/

/-- 005 Step 2, first half: `R₁ ∘ R₂ = Rot_A(2α)` (definitional). -/
theorem pair_collapse_A (x : K) (p : K × K) :
    reflCA x (reflAB p) = rotA x p := rfl

/-- 005 Step 2, second half: `R₀ ∘ R₂ = Rot_B(-2β)` (definitional:
`rotB` at `Q := y`). -/
theorem pair_collapse_B (x y i : K) (p : K × K) :
    reflBC x y i (reflAB p) = rotB x y y i p := rfl

/-- 005 Step 3, the closed form of the half word: as maps on pairs,
`R₀ ∘ Rot_A(2aα) ∘ Rot_B(-2bβ) = (z ↦ μ · conj z + w)`, with `μ`, `w`
the closed forms and the parameters entering only through `P`, `Q`. -/
theorem half_word_closed_form (x y P Q i : K)
    (hx : x ≠ 0) (hy : y ≠ 0) (hP : P ≠ 0) (hQ : Q ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) (p : K × K) :
    reflBC x y i (rotA P (rotB x y Q i p)) =
      (mu y P Q * p.2 + w x y P Q i,
       mu y⁻¹ P⁻¹ Q⁻¹ * p.1 + w x⁻¹ y⁻¹ P⁻¹ Q⁻¹ (-i)) := by
  rw [Prod.ext_iff]
  constructor <;>
  · simp only [reflBC, rotA, rotB, mu, w, B, sinF, inv_inv]
    field_simp
    ring_nf

/-- The letter word `(R₁R₂)^a` collapses to `Rot_A(2aα)` with `P = x^a`
as an actual power (the power step of 005 Step 1/2, formal). -/
theorem iterate_collapse_A (x : K) (a : ℕ) (p : K × K) :
    (fun q => reflCA x (reflAB q))^[a] p = rotA (x ^ a) p := by
  induction a generalizing p with
  | zero => simp [rotA]
  | succ n ih =>
    rw [Function.iterate_succ_apply', ih, pair_collapse_A]
    rw [Prod.ext_iff]
    constructor <;>
    · simp only [rotA, mul_inv, pow_succ]
      ring

/-- The letter word `(R₀R₂)^b` collapses to `Rot_B(-2bβ)` with `Q = y^b`
as an actual power. -/
theorem iterate_collapse_B (x y i : K) (b : ℕ) (p : K × K) :
    (fun q => reflBC x y i (reflAB q))^[b] p = rotB x y (y ^ b) i p := by
  induction b generalizing p with
  | zero => simp [rotB]
  | succ n ih =>
    rw [Function.iterate_succ_apply', ih, pair_collapse_B]
    rw [Prod.ext_iff]
    constructor <;>
    · simp only [rotB, mul_inv, pow_succ, inv_inv]
      ring

/-- The full letter half word `R₀ (R₁R₂)^a (R₀R₂)^b` in closed form, for
every `a b : ℕ`: `u(z) = μ conj z + w` with `P = x^a`, `Q = y^b`. -/
theorem half_word_letters (x y i : K) (a b : ℕ)
    (hx : x ≠ 0) (hy : y ≠ 0) (hI : i ≠ 0) (h2 : (2 : K) ≠ 0)
    (p : K × K) :
    reflBC x y i
        ((fun q => reflCA x (reflAB q))^[a]
          ((fun q => reflBC x y i (reflAB q))^[b] p)) =
      (mu y (x ^ a) (y ^ b) * p.2 + w x y (x ^ a) (y ^ b) i,
       mu y⁻¹ (x ^ a)⁻¹ (y ^ b)⁻¹ * p.1 +
         w x⁻¹ y⁻¹ (x ^ a)⁻¹ (y ^ b)⁻¹ (-i)) := by
  rw [iterate_collapse_B, iterate_collapse_A,
    half_word_closed_form x y (x ^ a) (y ^ b) i hx hy
      (pow_ne_zero a hx) (pow_ne_zero b hy) hI h2 p]

/-! ### Tier B: the glide facts and I1-I4, D1, D2 (005 Step 4) -/

/-- Glide fact: `μ = δ²`. -/
theorem mu_eq_delta_sq (y P Q : K) : mu y P Q = delta y P Q ^ 2 := by
  simp only [mu, delta]; ring

/-- Glide fact: `τ` is parallel to the glide axis, `Im(conj δ · τ) = 0`. -/
theorem tau_parallel_axis (x y P Q i : K)
    (hx : x ≠ 0) (hy : y ≠ 0) (hP : P ≠ 0) (hQ : Q ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) :
    proj y P Q i (tau x y P Q i) (tau x⁻¹ y⁻¹ P⁻¹ Q⁻¹ (-i)) = 0 := by
  simp only [proj, tau, mu, w, delta, B, sinF, inv_inv, neg_neg]
  field_simp
  ring_nf

/-- `μ` is a unit with `μ · conj μ = 1` (the half word is an isometry). -/
theorem mu_mul_conj_mu (y P Q : K)
    (hy : y ≠ 0) (hP : P ≠ 0) (hQ : Q ≠ 0) :
    mu y P Q * mu y⁻¹ P⁻¹ Q⁻¹ = 1 := by
  simp only [mu, inv_inv]
  field_simp

/-- The square of the half word is the translation by `τ`:
`u(u(z)) = z + τ` (005 Step 3's "u² is automatically a translation"). -/
theorem half_word_sq_is_translation (x y P Q i : K)
    (hx : x ≠ 0) (hy : y ≠ 0) (hP : P ≠ 0) (hQ : Q ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) (z : K) :
    mu y P Q * (mu y⁻¹ P⁻¹ Q⁻¹ * z + w x⁻¹ y⁻¹ P⁻¹ Q⁻¹ (-i)) +
        w x y P Q i = z + tau x y P Q i := by
  simp only [tau, mu, w, B, sinF, inv_inv]
  field_simp
  ring

/-- The glide action on the axis-normal projection: `p ∘ u = 2m - p`
(the identity behind 005's corridor reduction, 008 §1's `p(u(z)) + p(z)
- 2m = 0`). -/
theorem glide_action (x y P Q i : K)
    (hx : x ≠ 0) (hy : y ≠ 0) (hP : P ≠ 0) (hQ : Q ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) (z zc : K) :
    proj y P Q i (mu y P Q * zc + w x y P Q i)
        (mu y⁻¹ P⁻¹ Q⁻¹ * z + w x⁻¹ y⁻¹ P⁻¹ Q⁻¹ (-i)) =
      2 * m x y P Q i - proj y P Q i z zc := by
  simp only [proj, m, tau, mu, w, delta, B, sinF, inv_inv, neg_neg]
  field_simp
  ring_nf

/-- **I1** (005 (ii)-A Step 4): `m - p(C) = -[cos(aα) sin α sin(bβ)
+ sin(aα) sin β cos(α + (b+1)β)]`. -/
theorem identity_I1 (x y P Q i : K)
    (hx : x ≠ 0) (hy : y ≠ 0) (hP : P ≠ 0) (hQ : Q ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) :
    m x y P Q i - proj y P Q i (C x y i) (C x⁻¹ y⁻¹ (-i)) =
      -(cosF P * sinF x i * sinF Q i +
        sinF P i * sinF y i * cosF (x * y * Q)) := by
  simp only [proj, m, tau, mu, w, delta, B, C, sinF, cosF, inv_inv,
    neg_neg, mul_inv]
  field_simp
  ring_nf

/-- **I2**: `p(A₁) - m = cos(aα) sin((b+1)β) sin(α+β)`. -/
theorem identity_I2 (x y P Q i : K)
    (hx : x ≠ 0) (hy : y ≠ 0) (hP : P ≠ 0) (hQ : Q ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) :
    proj y P Q i (A1 x y i) (A1 x⁻¹ y⁻¹ (-i)) - m x y P Q i =
      cosF P * sinF (y * Q) i * sinF (x * y) i := by
  simp only [proj, m, tau, mu, w, delta, A1, B, sinF, cosF, inv_inv,
    neg_neg, mul_inv]
  field_simp
  ring_nf

/-- **I3**: `p(B) - m = cos((b+1)β) sin(aα) sin(α+β)`. -/
theorem identity_I3 (x y P Q i : K)
    (hx : x ≠ 0) (hy : y ≠ 0) (hP : P ≠ 0) (hQ : Q ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) :
    proj y P Q i (B x y i) (B x⁻¹ y⁻¹ (-i)) - m x y P Q i =
      cosF (y * Q) * sinF P i * sinF (x * y) i := by
  simp only [proj, m, tau, mu, w, delta, B, sinF, cosF, inv_inv,
    neg_neg, mul_inv]
  field_simp
  ring_nf

/-- **I4** (new in 005; the gate-(2a+2) C-image endpoint):
`p(v₂) - m = cos(aα) sin α sin(bβ) - sin(aα) sin β cos(α + (b+1)β)`.
Note I1 `= -(S+T)` and I4 `= S-T` for the same two products. -/
theorem identity_I4 (x y P Q i : K)
    (hx : x ≠ 0) (hy : y ≠ 0) (hP : P ≠ 0) (hQ : Q ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) :
    proj y P Q i (v2 x y P i) (v2 x⁻¹ y⁻¹ P⁻¹ (-i)) - m x y P Q i =
      cosF P * sinF x i * sinF Q i -
        sinF P i * sinF y i * cosF (x * y * Q) := by
  simp only [proj, m, tau, mu, w, delta, v2, B, C, sinF, cosF, inv_inv,
    neg_neg, mul_inv]
  field_simp
  ring_nf

/-- **D1** (collapsing-gate locus): `p(B) - p(C) = sin α sin(aα - bβ)`. -/
theorem identity_D1 (x y P Q i : K)
    (hx : x ≠ 0) (hy : y ≠ 0) (hP : P ≠ 0) (hQ : Q ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) :
    proj y P Q i (B x y i) (B x⁻¹ y⁻¹ (-i)) -
        proj y P Q i (C x y i) (C x⁻¹ y⁻¹ (-i)) =
      sinF x i * sinF (P * Q⁻¹) i := by
  simp only [proj, delta, B, C, sinF, inv_inv, mul_inv]
  field_simp
  ring_nf

/-- **D2** (collapsing-gate locus):
`p(A₁) - p(C) = -sin β sin((a-1)α - (b+1)β)`. -/
theorem identity_D2 (x y P Q i : K)
    (hx : x ≠ 0) (hy : y ≠ 0) (hP : P ≠ 0) (hQ : Q ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) :
    proj y P Q i (A1 x y i) (A1 x⁻¹ y⁻¹ (-i)) -
        proj y P Q i (C x y i) (C x⁻¹ y⁻¹ (-i)) =
      -(sinF y i * sinF (x⁻¹ * y⁻¹ * P * Q⁻¹) i) := by
  simp only [proj, delta, A1, B, C, sinF, inv_inv, mul_inv]
  field_simp
  ring_nf

/-! ### Tier C: specialization to integer exponents (005's ring
homomorphism `P → x^a`, `Q → y^b`)

In this encoding the specialization is instantiation of the universally
quantified generators, so each corollary is a one-line application; the
`a b : ℤ` genericity is inherited, never re-proved.  (`zpow_ne_zero`:
in a field `x ≠ 0 → x^(a:ℤ) ≠ 0`.) -/

section Specialize

variable (x y i : K) (a b : ℤ)

/-- The closed-form composition at `P = x^a`, `Q = y^b`, `a b : ℤ`. -/
theorem half_word_closed_form_spec (hx : x ≠ 0) (hy : y ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) (p : K × K) :
    reflBC x y i (rotA (x ^ a) (rotB x y (y ^ b) i p)) =
      (mu y (x ^ a) (y ^ b) * p.2 + w x y (x ^ a) (y ^ b) i,
       mu y⁻¹ (x ^ a)⁻¹ (y ^ b)⁻¹ * p.1 +
         w x⁻¹ y⁻¹ (x ^ a)⁻¹ (y ^ b)⁻¹ (-i)) :=
  half_word_closed_form x y (x ^ a) (y ^ b) i hx hy
    (zpow_ne_zero a hx) (zpow_ne_zero b hy) hI h2 p

/-- I1 for every pair of integer exponents. -/
theorem identity_I1_spec (hx : x ≠ 0) (hy : y ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) :
    m x y (x ^ a) (y ^ b) i -
        proj y (x ^ a) (y ^ b) i (C x y i) (C x⁻¹ y⁻¹ (-i)) =
      -(cosF (x ^ a) * sinF x i * sinF (y ^ b) i +
        sinF (x ^ a) i * sinF y i * cosF (x * y * y ^ b)) :=
  identity_I1 x y (x ^ a) (y ^ b) i hx hy
    (zpow_ne_zero a hx) (zpow_ne_zero b hy) hI h2

/-- I2 for every pair of integer exponents. -/
theorem identity_I2_spec (hx : x ≠ 0) (hy : y ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) :
    proj y (x ^ a) (y ^ b) i (A1 x y i) (A1 x⁻¹ y⁻¹ (-i)) -
        m x y (x ^ a) (y ^ b) i =
      cosF (x ^ a) * sinF (y * y ^ b) i * sinF (x * y) i :=
  identity_I2 x y (x ^ a) (y ^ b) i hx hy
    (zpow_ne_zero a hx) (zpow_ne_zero b hy) hI h2

/-- I3 for every pair of integer exponents. -/
theorem identity_I3_spec (hx : x ≠ 0) (hy : y ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) :
    proj y (x ^ a) (y ^ b) i (B x y i) (B x⁻¹ y⁻¹ (-i)) -
        m x y (x ^ a) (y ^ b) i =
      cosF (y * y ^ b) * sinF (x ^ a) i * sinF (x * y) i :=
  identity_I3 x y (x ^ a) (y ^ b) i hx hy
    (zpow_ne_zero a hx) (zpow_ne_zero b hy) hI h2

/-- I4 for every pair of integer exponents. -/
theorem identity_I4_spec (hx : x ≠ 0) (hy : y ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) :
    proj y (x ^ a) (y ^ b) i (v2 x y (x ^ a) i)
        (v2 x⁻¹ y⁻¹ (x ^ a)⁻¹ (-i)) -
        m x y (x ^ a) (y ^ b) i =
      cosF (x ^ a) * sinF x i * sinF (y ^ b) i -
        sinF (x ^ a) i * sinF y i * cosF (x * y * y ^ b) :=
  identity_I4 x y (x ^ a) (y ^ b) i hx hy
    (zpow_ne_zero a hx) (zpow_ne_zero b hy) hI h2

/-- The glide facts for every pair of integer exponents. -/
theorem glide_facts_spec (hx : x ≠ 0) (hy : y ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) :
    mu y (x ^ a) (y ^ b) = delta y (x ^ a) (y ^ b) ^ 2 ∧
    proj y (x ^ a) (y ^ b) i (tau x y (x ^ a) (y ^ b) i)
      (tau x⁻¹ y⁻¹ (x ^ a)⁻¹ (y ^ b)⁻¹ (-i)) = 0 :=
  ⟨mu_eq_delta_sq y (x ^ a) (y ^ b),
   tau_parallel_axis x y (x ^ a) (y ^ b) i hx hy
     (zpow_ne_zero a hx) (zpow_ne_zero b hy) hI h2⟩

end Specialize

end Laurent
end BilliardsFormal
