# 011 — Lean 4 certificate for 005's Laurent-identity block (I1–I4, glide facts, closed-form composition)

- **Problem:** billiards-triangles, `problems/billiards-triangles/PROBLEM.md`
- **Date:** 2026-08-04
- **Mode:** informed (read `prior-art.json`, 005 in full, 008 §§1–2,
  `docs/FORMALIZE.md`, `STATUS.md` queue item 15 and the 2026-07-31
  formal-parameter-specialization insight, 009 in full, the 009/010 index
  entries, and `explore/plaw_general.py`'s ring layer — the proof object
  being formalized)
- **Type:** formalization — review-shape attempt against
  `005-complete-death-law-theorem.md` (verifies 005 §(ii)-A Steps 2–4: the
  closed-form composition of the half word and the identity block I1–I4,
  D1, D2 + glide facts, with (a,b) as formal parameters, now
  kernel-checked). Second run of the `docs/FORMALIZE.md` lane, following
  the completed L1 pilot (009/010); queue item 15.
- **Tools:** the existing Lean project at
  `problems/billiards-triangles/formal/` — toolchain UNCHANGED from
  009/010: `leanprover/lean4:v4.32.2`, mathlib tag `v4.32.2`
  (rev `905b95818eb32af7874a58b427f50c1711a5e96c`), prebuilt cache reused,
  nothing re-downloaded. New source:
  `formal/BilliardsFormal/LaurentBlock.lean` (~420 lines, 23 theorems),
  imported from `formal/BilliardsFormal.lean`; `formal/AxiomCheck.lean`
  (the audit driver, now on disk — 010's nit). `CotMonotone.lean`
  untouched. Build is deterministic; the new module elaborates in ~24 s;
  build + axiom audit appended to `formal/lake-build.log`. A throwaway
  numeric transcription check (random nonzero complex points, session
  scratchpad, not kept) was used to pin the encoding before proving; the
  Lean proofs supersede it.
- **Sources:** none external beyond mathlib. Repo: 005 (the record
  verified), 008 (skeptic confirmation of the same block, incl. the pair
  "off-torus" algebra this encoding follows), 009/010 (lane precedent),
  `docs/FORMALIZE.md`.

**Forbidden-tactic note:** no `decide`, no `native_decide` anywhere in the
file; the proofs are `field_simp` + `ring`/`ring_nf` (plus two `rfl`s and
two inductions). The axiom audit below confirms nothing beyond the
standard three.

## Claims attacked

1. **005 §(ii)-A Steps 2–4 as ring-model statements**: the pair collapse
   (R₁R₂ = Rot_A(2α), R₀R₂ = Rot_B(−2β)), the closed form
   u = R₀ ∘ Rot_A(2aα) ∘ Rot_B(−2bβ) = (z ↦ μ·conj z + w), and the eight
   formal identities (glide facts μ = δ², Im(conj δ·τ) = 0; I1, I2, I3,
   I4; D1, D2), with the (a,b)-dependence carried by formal variables P, Q
   in the exponent role — one check for all integers (a,b). Attack =
   make Lean's kernel accept them. A wrong sign or a mis-transcribed
   closed form would have refused to close under `ring`.
2. **The informal→formal statement mapping** for that block (the
   FORMALIZE.md bridge caveat): that the Lean statements are 005's
   statements and not a narrowing. Attacked by the encoding table below —
   which, per the lane, must be re-checked by a non-formalizer, so it is
   flagged **PENDING INDEPENDENT STATEMENT REVIEW**.

Not attacked, deliberately: **the BRIDGE (geometry = ring model) is
explicitly OUT of scope — 005/008 carry it.** Concretely: that the
circumdiameter-1 normalization and the three base reflections are what
005 Steps 1–2 say they are geometrically, and that `proj`/`m` are the
harness's corridor quantities, remains informal — cross-checked in 005
(23-member exact specialization + 400-point float check) and 008 (hand
re-derivation + exact off-torus recomposition on 40 members), and stated
verbatim in the module docstring. Also not attacked: Lemmas C/D (L1 was
009), the mod-360 case tree, the necessity theorem, the sufficiency
segments — the analytic layer stays at `VERIFIED`.

## Refutations found

None. Every identity of 005's block closed in the kernel with no
hypothesis strengthened and no case dropped. One **generalization**
(a delta in the safe direction, not a narrowing) was found by the
formalization itself:

- **The block is homogeneous in i.** Every identity holds over an
  arbitrary field for ANY nonzero `i`, with no relation `i² = −1`:
  both sides of each identity carry the same total 1/(2i)-degree, so `i`
  scales out. Discovered because `ring_nf` closed every goal without
  ever consuming `i² = −1`; the hypothesis was then removed from the
  statements. 005's Q(i)-coefficient ring is the instance
  `i := Complex.I` — every 005 use case is covered, and the certificate
  is strictly stronger. (Also confirmed numerically at random points
  with i = 0.7 − 1.3j before the hypothesis was dropped.)

A tooling non-finding worth the ledger (queue item 15 flagged this as the
risk): **no Laurent-polynomial or `AddMonoidAlgebra` infrastructure was
needed.** Universally quantifying over a field K and nonzero generators
x, y, P, Q subsumes the Laurent ring (instantiate K at the fraction field
of ℚ[x,y,P,Q] to recover the formal identity; at ℂ on the torus to
recover the geometry), and `field_simp` + `ring_nf` close every identity
in seconds. The anticipated mathlib API fight never happened.

## Claims that survive

### Tier reached: A + B + C, all green

- **Tier A** (closed-form composition): `pair_collapse_A`,
  `pair_collapse_B`, `half_word_closed_form` — plus, beyond the asked
  scope, `iterate_collapse_A`, `iterate_collapse_B`, `half_word_letters`:
  the genuine letter word R₀ (R₁R₂)^a (R₀R₂)^b, as an a- and b-fold
  function iterate, equals the closed form with P = x^a, Q = y^b as
  actual powers, for every a b : ℕ — this formalizes the power-collapse
  step of 005 Steps 1–2 *inside the ring model*, shrinking the informal
  bridge to the base-geometry mapping only.
- **Tier B**: `mu_eq_delta_sq`, `tau_parallel_axis`, `mu_mul_conj_mu`,
  `half_word_sq_is_translation` (u² = translation by τ), `glide_action`
  (p ∘ u = 2m − p, the corridor-reduction identity 008 §1 verified), and
  `identity_I1`, `identity_I2`, `identity_I3`, `identity_I4`,
  `identity_D1`, `identity_D2`.
- **Tier C** (specialization): `half_word_closed_form_spec`,
  `identity_I1_spec` … `identity_I4_spec`, `glide_facts_spec` — P := x^a,
  Q := y^b for arbitrary a b : ℤ (zpow). In this encoding 005's
  ring-homomorphism specialization is *instantiation of universally
  quantified variables*, so each corollary is a one-line application and
  the "formal zero ⟹ geometric zero" direction is inherited from the
  ∀-statement rather than re-argued.

### The encoding (formal statement ↔ 005's objects) — **PENDING INDEPENDENT STATEMENT REVIEW**

All in `formal/BilliardsFormal/LaurentBlock.lean`, namespace
`BilliardsFormal.Laurent`, over `{K : Type*} [Field K]` with hypotheses
`x ≠ 0`, `y ≠ 0`, `P ≠ 0`, `Q ≠ 0`, `i ≠ 0`, `(2 : K) ≠ 0`:

| 005 / `plaw_general.py` object | formal counterpart | delta |
|---|---|---|
| ring Q(i)[x^±, y^±, P^±, Q^±], x = e^{iα}, y = e^{iβ}, P = e^{iaα}, Q = e^{ibβ} | arbitrary field K, generators x y P Q i ≠ 0, statements ∀-quantified | generalization: any K, any nonzero i (no i² = −1; see above); Laurent ring recovered by instantiation |
| monomial (ex, ey, eP, eQ) | `x^ex * y^ey * P^eP * Q^eQ` | none; exponent lattice = actual exponents |
| sin/cos of ex·α + ey·β + eP·aα + eQ·bβ | `sinF z i = (z - z⁻¹)/(2*i)`, `cosF z = (z + z⁻¹)/2` at the corresponding monomial z | none (Euler form; `sinF` means sine only at i = imaginary unit — bridge note in docstring) |
| conj (exponent negation + coefficient conjugation) | substitution `(x,y,P,Q,i) ↦ (x⁻¹,y⁻¹,P⁻¹,Q⁻¹,−i)` | none on the torus (008 §2's star involution); a point is carried as a pair `(z, conj z)`, conj supplied as the second component / the substituted expression |
| B = sin(α+β), C = sin β·e^{iα}, A₁ = B(1−y⁻²) | `B x y i`, `C x y i`, `A1 x y i` | none (defs transcribe `closed_forms()` monomial by monomial) |
| μ = y⁻²P⁻²Q², δ = y⁻¹P⁻¹Q, w = B(1−y⁻²+y⁻²P⁻²−μ), τ = w + μ·conj w | `mu y P Q`, `delta y P Q`, `w x y P Q i`, `tau x y P Q i` | none |
| p(v) = Im(conj δ·v); m = ½·Im(conj δ·(w−τ/2)) | `proj y P Q i v vc` with vc = conj v; `m x y P Q i` | none (Im ζ = (ζ − conj ζ)/(2i), pair convention) |
| v₂ = B + y⁻²(P⁻²·conj C − B) | `v2 x y P i` | Q dropped from the signature (v₂ does not depend on Q) |
| R₂ = conj; R₁ = e^{2iα}·conj; R₀ = B + e^{−2iβ}·conj(·−B); Rot_A(2aα); Rot_B(−2bβ) | `reflAB`, `reflCA`, `reflBC`, `rotA`, `rotB` on pairs K × K | each map's second component is the conjugate substitution of its first, syntactically |
| specialization P → x^a, Q → y^b by ring homomorphism | instantiation of the ∀-quantified P, Q at `x^(a:ℤ)`, `y^(b:ℤ)` | mechanism differs (instantiation vs ring hom); conclusion identical, and the letters theorems additionally prove the a b : ℕ power-collapse |

### The formal statements (verbatim, main theorems)

```lean
theorem half_word_closed_form (x y P Q i : K)
    (hx : x ≠ 0) (hy : y ≠ 0) (hP : P ≠ 0) (hQ : Q ≠ 0)
    (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) (p : K × K) :
    reflBC x y i (rotA P (rotB x y Q i p)) =
      (mu y P Q * p.2 + w x y P Q i,
       mu y⁻¹ P⁻¹ Q⁻¹ * p.1 + w x⁻¹ y⁻¹ P⁻¹ Q⁻¹ (-i))

theorem half_word_letters (x y i : K) (a b : ℕ)
    (hx : x ≠ 0) (hy : y ≠ 0) (hI : i ≠ 0) (h2 : (2 : K) ≠ 0)
    (p : K × K) :
    reflBC x y i
        ((fun q => reflCA x (reflAB q))^[a]
          ((fun q => reflBC x y i (reflAB q))^[b] p)) =
      (mu y (x ^ a) (y ^ b) * p.2 + w x y (x ^ a) (y ^ b) i,
       mu y⁻¹ (x ^ a)⁻¹ (y ^ b)⁻¹ * p.1 +
         w x⁻¹ y⁻¹ (x ^ a)⁻¹ (y ^ b)⁻¹ (-i))

theorem identity_I1 (x y P Q i : K) (hx : x ≠ 0) (hy : y ≠ 0)
    (hP : P ≠ 0) (hQ : Q ≠ 0) (hI : i ≠ 0) (h2 : (2 : K) ≠ 0) :
    m x y P Q i - proj y P Q i (C x y i) (C x⁻¹ y⁻¹ (-i)) =
      -(cosF P * sinF x i * sinF Q i +
        sinF P i * sinF y i * cosF (x * y * Q))

theorem identity_I2 ⋯ :
    proj y P Q i (A1 x y i) (A1 x⁻¹ y⁻¹ (-i)) - m x y P Q i =
      cosF P * sinF (y * Q) i * sinF (x * y) i

theorem identity_I3 ⋯ :
    proj y P Q i (B x y i) (B x⁻¹ y⁻¹ (-i)) - m x y P Q i =
      cosF (y * Q) * sinF P i * sinF (x * y) i

theorem identity_I4 ⋯ :
    proj y P Q i (v2 x y P i) (v2 x⁻¹ y⁻¹ P⁻¹ (-i)) - m x y P Q i =
      cosF P * sinF x i * sinF Q i -
        sinF P i * sinF y i * cosF (x * y * Q)

theorem identity_D1 ⋯ :
    proj y P Q i (B x y i) (B x⁻¹ y⁻¹ (-i)) -
        proj y P Q i (C x y i) (C x⁻¹ y⁻¹ (-i)) =
      sinF x i * sinF (P * Q⁻¹) i

theorem identity_D2 ⋯ :
    proj y P Q i (A1 x y i) (A1 x⁻¹ y⁻¹ (-i)) -
        proj y P Q i (C x y i) (C x⁻¹ y⁻¹ (-i)) =
      -(sinF y i * sinF (x⁻¹ * y⁻¹ * P * Q⁻¹) i)

theorem tau_parallel_axis ⋯ :
    proj y P Q i (tau x y P Q i) (tau x⁻¹ y⁻¹ P⁻¹ Q⁻¹ (-i)) = 0

theorem mu_eq_delta_sq (y P Q : K) : mu y P Q = delta y P Q ^ 2

theorem glide_action ⋯ (z zc : K) :
    proj y P Q i (mu y P Q * zc + w x y P Q i)
        (mu y⁻¹ P⁻¹ Q⁻¹ * z + w x⁻¹ y⁻¹ P⁻¹ Q⁻¹ (-i)) =
      2 * m x y P Q i - proj y P Q i z zc
```

(elided `⋯` = the same six hypotheses as `identity_I1`; full text of all
23 theorems, including the Tier C `*_spec` corollaries at a b : ℤ, in the
source file.)

### Build and axiom audit

- `lake build`: `Build completed successfully (8658 jobs)`;
  `BilliardsFormal.LaurentBlock` elaborates in ~24 s; zero `sorry`, zero
  warnings; the root module imports both `CotMonotone` and `LaurentBlock`.
- `lake env lean AxiomCheck.lean`: all 23 `Laurent` theorems (and 009's
  two) depend on exactly `[propext, Classical.choice, Quot.sound]` — no
  `sorryAx`, no `Lean.ofReduceBool` (i.e. no `native_decide`).
- Both outputs appended to `formal/lake-build.log` (timestamped block).
- Toolchain pins unchanged: `lean-toolchain` = `leanprover/lean4:v4.32.2`,
  mathlib tag `v4.32.2` = rev `905b95818eb32af7874a58b427f50c1711a5e96c`
  in `lakefile.toml` / `lake-manifest.json`.

### What the certificate does and does not cover

Covered: exactly the formal statements above — 005's identity block and
closed-form composition in the ring model, for all formal P, Q (hence all
integer (a,b) by Tier C), over any field with the six nonzeroness
hypotheses. Not covered (unchanged from 005/008's own scoping, and per
the lane): **the bridge (geometry = ring model) is explicitly OUT of
scope — 005/008 carry it**; the necessity case tree, Lemmas C/D, the
sufficiency segments, and every numerical claim of 005 remain `VERIFIED`,
not `FORMALIZED`. Per `docs/FORMALIZE.md`, this pass is complete only
after (a) an independent statement review of the encoding table above and
(b) an independent rebuild — neither is done by this record.

## Leads generated

1. **Complete the FORMALIZED pass** (the 010 pattern): an independent
   agent (not this one) (a) re-runs `lake build` + `AxiomCheck.lean` from
   the committed source in a fresh container (definite outcome: green or
   not), and (b) reviews the encoding table hypothesis-by-hypothesis
   against 005 §(ii)-A and `plaw_general.closed_forms()` — the
   high-value checks are the conj-as-substitution convention, the pair
   second components, and the sinF/cosF monomial dictionary (definite
   outcome: sign-off or a filed delta).
2. **The i-homogeneity is a checkable statement about 005's own ring:**
   every coefficient in `plaw_general.py`'s eight zero-polynomial checks
   should live in ℚ·i^k with k determined by the term's sin-depth, i.e.
   the Q(i) coefficients are never genuinely mixed. One afternoon Python
   check of the coefficient dictionaries; if true, 005's ring can be
   restated over ℚ with i a formal unit — worth knowing before any
   future port (definite outcome: true/false).
3. **Formalize Lemma D and the L1 → C/D chain** (009 lead 2, still
   open): with L1 and now the block formal, Lemmas C/D are the remaining
   informal layer between the certificates and 005's necessity theorem.
4. **Shrink the bridge with a torus-evaluation lemma:** over ℂ,
   `sinF (Complex.exp (θ * I)) I = Complex.sin θ` (Euler) would connect
   `sinF`/`cosF` to actual trigonometry, moving the "sinF means sine"
   line of the mapping table from prose into the kernel. Small, definite.
5. **Port the skeleton to the design-track families** (005 lead 5): the
   file contains nothing billiards-specific except the three base
   reflections; any family (w₀B₁^aB₂^b)² with vertex-pivoting blocks
   gets the same `half_word_letters`-style theorem by editing the three
   `refl*` definitions (definite outcome per family: builds or surfaces
   a structural mismatch).

## References

- `problems/billiards-triangles/attempts/005-complete-death-law-theorem.md`
  §(ii)-A (the block formalized), §Approach (the formal-ring method).
- `problems/billiards-triangles/attempts/008-skeptic-review-of-005.md`
  §§1–2 (skeptic confirmation; the off-torus pair algebra and the
  specialization-soundness argument this encoding leans on).
- `problems/billiards-triangles/attempts/009-lean-pilot-cot-lemma.md` and
  `010-statement-review-of-009.md` (lane precedent; toolchain; the
  AxiomCheck nit fixed here).
- `problems/billiards-triangles/explore/plaw_general.py` (`closed_forms`,
  `cmd_formal` — the proof object transcribed).
- `docs/FORMALIZE.md` (the lane: qualification, bridge caveat, three-part
  pass).
- mathlib4, tag v4.32.2, github.com/leanprover-community/mathlib4
  (rev `905b95818eb32af7874a58b427f50c1711a5e96c`). No papers consulted;
  no machine transcriptions.
