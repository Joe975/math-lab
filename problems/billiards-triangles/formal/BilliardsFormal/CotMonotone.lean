/-
Lemma L1 of attempt record 005 (`problems/billiards-triangles/attempts/
005-complete-death-law-theorem.md`, section (ii)-B), reviewed and confirmed in
008 (section 3):

  **Lemma L1.** For fixed s ∈ (0, π/2], the map c ↦ c·cot(c·s) is strictly
  decreasing on (0, 1].

This is the elementary monotonicity lemma that retires Lemmas C and D of the
W(a,b) death-law theorem (005 uses it three times: the two brackets of
Lemma C's log-derivative, and Lemma D).

Formal statement proved here: `cot_lemma_L1` below —
  `StrictAntiOn (fun c : ℝ => c * Real.cot (c * s)) (Set.Ioc 0 1)`
for every real `s` with `0 < s` and `s ≤ π / 2`.

The proof follows 005's own two-line argument: the derivative of
`u ↦ u · cot u` is `(sin u · cos u − u) / sin² u < 0` (since `sin x < x` for
`x > 0`), giving strict decrease of `u · cot u` on `(0, π/2]`
(`strictAntiOn_id_mul_cot`); Lemma L1 follows by the substitution `u = c·s`,
which rescales the value by the positive constant `1/s`.

Zero `sorry`; no axioms beyond mathlib. Toolchain pinned in `lean-toolchain`
(leanprover/lean4:v4.32.2), mathlib pinned in `lakefile.toml` /
`lake-manifest.json` (tag v4.32.2).
-/
import Mathlib

open Real Set

namespace BilliardsFormal

/-- The derivative of `u ↦ u * cos u / sin u` at a point where `sin x ≠ 0`. -/
lemma hasDerivAt_id_mul_cot {x : ℝ} (hsin : Real.sin x ≠ 0) :
    HasDerivAt (fun u : ℝ => u * Real.cos u / Real.sin u)
      (((1 * Real.cos x + x * -Real.sin x) * Real.sin x
          - x * Real.cos x * Real.cos x) / Real.sin x ^ 2) x := by
  exact ((hasDerivAt_id' (x := x)).mul (Real.hasDerivAt_cos x)).div
    (Real.hasDerivAt_sin x) hsin

/-- Core monotonicity: `u ↦ u * cot u` is strictly decreasing on `(0, π/2]`.
    (005's Lemma L1 after the substitution `u = c·s`.) -/
theorem strictAntiOn_id_mul_cot :
    StrictAntiOn (fun u : ℝ => u * Real.cot u) (Set.Ioc 0 (π / 2)) := by
  have hfun : (fun u : ℝ => u * Real.cot u)
      = fun u : ℝ => u * Real.cos u / Real.sin u := by
    funext u
    rw [Real.cot_eq_cos_div_sin, mul_div_assoc]
  rw [hfun]
  have hpi2 : (π / 2 : ℝ) < π := half_lt_self Real.pi_pos
  have hsinpos : ∀ x ∈ Set.Ioc (0 : ℝ) (π / 2), 0 < Real.sin x := fun x hx =>
    Real.sin_pos_of_pos_of_lt_pi hx.1 (lt_of_le_of_lt hx.2 hpi2)
  apply strictAntiOn_of_deriv_neg (convex_Ioc 0 (π / 2))
  · -- continuity on (0, π/2]
    exact continuousOn_id.mul (Real.continuous_cos.continuousOn) |>.div
      Real.continuous_sin.continuousOn (fun x hx => (hsinpos x hx).ne')
  · -- negative derivative on the interior (0, π/2)
    intro x hx
    rw [interior_Ioc] at hx
    have hx0 : 0 < x := hx.1
    have hsin : 0 < Real.sin x := hsinpos x ⟨hx.1, hx.2.le⟩
    rw [(hasDerivAt_id_mul_cot hsin.ne').deriv]
    have heq : ((1 * Real.cos x + x * -Real.sin x) * Real.sin x
          - x * Real.cos x * Real.cos x) / Real.sin x ^ 2
        = (Real.sin x * Real.cos x - x) / Real.sin x ^ 2 := by
      linear_combination (-x / Real.sin x ^ 2) * Real.sin_sq_add_cos_sq x
    rw [heq]
    apply div_neg_of_neg_of_pos
    · -- sin x · cos x < x, since sin x · cos x ≤ sin x < x
      have h1 : Real.sin x * Real.cos x ≤ Real.sin x :=
        mul_le_of_le_one_right hsin.le (Real.cos_le_one x)
      have h2 : Real.sin x < x := Real.sin_lt hx0
      linarith
    · exact pow_pos hsin 2

/-- **Lemma L1** (record 005, section (ii)-B), stated exactly as written there:
    for fixed `s ∈ (0, π/2]`, the map `c ↦ c·cot(c·s)` is strictly decreasing
    on `(0, 1]`. `StrictAntiOn f I` is mathlib's "f is strictly decreasing on
    I": for all `c₁, c₂ ∈ I` with `c₁ < c₂`, `f c₂ < f c₁`. -/
theorem cot_lemma_L1 {s : ℝ} (hs0 : 0 < s) (hs : s ≤ π / 2) :
    StrictAntiOn (fun c : ℝ => c * Real.cot (c * s)) (Set.Ioc 0 1) := by
  intro c₁ hc₁ c₂ hc₂ hlt
  have mem : ∀ c : ℝ, c ∈ Set.Ioc (0 : ℝ) 1 → c * s ∈ Set.Ioc (0 : ℝ) (π / 2) :=
    fun c hc => ⟨mul_pos hc.1 hs0,
      le_trans (mul_le_of_le_one_left hs0.le hc.2) hs⟩
  have h := strictAntiOn_id_mul_cot (mem c₁ hc₁) (mem c₂ hc₂)
    (mul_lt_mul_of_pos_right hlt hs0)
  -- h : (c₂·s)·cot(c₂·s) < (c₁·s)·cot(c₁·s); divide by s > 0.
  have h' : s * (c₂ * Real.cot (c₂ * s)) < s * (c₁ * Real.cot (c₁ * s)) := by
    calc s * (c₂ * Real.cot (c₂ * s))
        = (c₂ * s) * Real.cot (c₂ * s) := by ring
      _ < (c₁ * s) * Real.cot (c₁ * s) := h
      _ = s * (c₁ * Real.cot (c₁ * s)) := by ring
  exact lt_of_mul_lt_mul_left h' hs0.le

end BilliardsFormal
