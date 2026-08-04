# 012 — Statement review + independent rebuild of 011 (Lean Laurent block)

- **Problem:** billiards-triangles, `problems/billiards-triangles/PROBLEM.md`
- **Date:** 2026-08-04
- **Mode:** informed (read `prior-art.json`, 011 in full, 005 in full —
  the source of truth for the statements — 008 §§1–2 (the second witness:
  star involution, off-torus pair algebra, v₂'s geometric identification),
  `docs/FORMALIZE.md`, 010 for the pass shape,
  `explore/plaw_general.py` line by line (`closed_forms`, `cmd_formal` —
  the proof object 011 transcribed), and every line of the five files
  under `formal/` outside `.lake/`)
- **Type:** skeptic review of `011-lean-laurent-block.md` (default stance:
  REFUTE). Independent completion of the three-part `FORMALIZED` pass of
  `docs/FORMALIZE.md` for 005's Laurent-identity block: (1) statement
  review by a non-formalizer — the main event here, since the encoding is
  nontrivial and a wrong mapping is invisible to the kernel — (2) clean
  build, (3) reproduction via my own `lake build`. I wrote none of the
  Lean code under review.
- **Outcome in one line:** 011 survives on every attacked front. The 23
  Lean statements are 005's statements under the encoding table, checked
  object by object and sign by sign; every theorem re-proven from scratch
  in my own 7-variable symbolic ring (i formal — independently confirming
  the i-homogeneity claim); the statements at the torus instance match
  005's prose quantities to ~4e-15 at 200 random (a, b, α, β); my own
  fresh `lake build` is green and my own axiom audit of all 23 theorems
  shows exactly `[propext, Classical.choice, Quot.sound]`. Zero
  corrections, three cosmetic nits. The FORMALIZED pass on the Laurent
  block is now **complete** (bridge permanently out of scope, as designed).
- **Tools:** new `explore/lbsk_review.py` (stdlib-only, deterministic,
  seed 20260804, ~5 s; output `data/lbsk_review.json`): Layer 1 is a
  from-scratch 7-variable Laurent-polynomial ring over ℚ — variables
  x, y, P, Q, **i**, z1, z2, so `i` is a formal unit exactly as in the
  Lean statements — re-proving every identity by exact polynomial
  subtraction, with the definitions transcribed from `LaurentBlock.lean`
  (NOT from `plaw_general.py`, which it never imports); Layer 2
  instantiates the Lean expressions at x = e^{iα}, y = e^{iβ},
  P = e^{iaα}, Q = e^{ibβ}, i = 1j and compares against 005's own prose
  formulas. Rebuild: elan / Lean `leanprover/lean4:v4.32.2` via
  `$HOME/.elan/bin`, `rm -rf .lake/build && lake build` (packages cache
  kept, per the lane), then my own scratch `SkepticAxiomCheck011.lean`
  (deleted after the run; contents = `#print axioms` on all 23 theorems,
  enumerated by me from the source file, not copied from
  `AxiomCheck.lean`).
- **Sources:** repo + the pinned mathlib checkout only; no external
  papers.

No file of 011 (or any prior attempt) was modified.

Reproduce (repo root):

```
python3 problems/billiards-triangles/explore/lbsk_review.py --out problems/billiards-triangles/data/lbsk_review.json
cd problems/billiards-triangles/formal && rm -rf .lake/build && lake build   # PATH needs $HOME/.elan/bin
# axiom audit: scratch file importing BilliardsFormal with #print axioms on all 23 Laurent theorems
```

## Claims attacked

The assigned list (a)–(h) plus two of my own (i)–(j):

1. **(a) Statement mapping, theorem by theorem** — that each of the 23
   Lean statements says what 005 §(ii)-A says: variable correspondence,
   the sinF/cosF dictionary at the right angles (including the
   factor-of-2 conventions in Rot_A(2aα), Rot_B(−2bβ)), sign
   conventions, and that `identity_I1..I4`, `D1`, `D2` are 005's I1–I4,
   D1, D2 and not renamed/permuted variants.
2. **(b) The conjugation encoding** — conj as the substitution
   (x,y,P,Q,i) ↦ (x⁻¹,y⁻¹,P⁻¹,Q⁻¹,−i), points as pairs with
   syntactically-substituted second components: does anything in the
   block need conj as a ring involution on arbitrary elements; does the
   pair encoding smuggle in "second component = star of the first"; is
   this 008's star involution?
3. **(c) The dropped i² = −1** — is the generalization direction safe,
   i.e. do the statements at the intended instance i = Complex.I say
   005's statements, with sinF/cosF reducing to actual sin/cos there?
4. **(d) Genericity** — (a,b) enter only as exponents of ∀-quantified
   generators; the ℕ-iterate + ℤ-specialization combination covers 005's
   integer parameters a ≥ b ≥ 1 with no silent narrowing; v₂'s dropped Q
   argument is genuinely absent in 005's v₂ too.
5. **(e) Cheat scan** — sorry/admit/native_decide/axiom/unsafe/partial/
   implemented_by/notation shadowing/vacuous hypotheses over all project
   .lean files; that the root module imports `LaurentBlock` so the build
   elaborates it.
6. **(f) Independent rebuild + axiom audit** — my own `lake build` after
   deleting the project build output; my own `#print axioms` on all 23
   theorems from a scratch file I wrote, not the committed
   `AxiomCheck.lean`.
7. **(g) Record audit** — does 011 or its index entry overclaim; is the
   bridge genuinely excluded in prose; is the numeric i-check labelled a
   sanity check; are range/gaps/one_line honest.
8. **(h) Spot re-derivation** — re-derive identities independently of
   both 005's prose and the Lean proofs. Done for ALL of them
   (Layer 1's from-scratch exact symbolic computation is exactly the
   "Python Fractions over a rational function field at symbolic
   exponents" route), plus a full by-hand product-to-sum derivation of
   D1 as a second, computer-free witness.
9. **(i, own) Coefficient-rationality of the conj convention** — the
   substitution is complex conjugation only if every constant in the
   definitions is rational (conjugating c(i) ∈ ℚ(i) is i ↦ −i); a single
   non-rational literal would silently break the encoding.
10. **(j, own) Toolchain-pin consistency** — `lean-toolchain`,
   `lakefile.toml`, `lake-manifest.json`, and the log agree on
   Lean 4.32.2 / mathlib rev `905b958…`.

## Refutations found

**None.** No statement delta, no sign error, no narrowing, no cheat, no
axiom beyond the standard three, no overclaim. Three cosmetic nits (none
is a correction; nothing downstream changes):

1. 011's Tools line says `LaurentBlock.lean` is "~420 lines"; the file
   is 423 lines. Rounding, noted only for completeness.
2. 011 reports the module elaborating in ~24 s; my rebuild took 41 s for
   the same module (33 s for CotMonotone vs 010's 27 s). Load-dependent,
   not a reproducibility problem — both builds are green with identical
   job counts (8658).
3. The conj-as-substitution convention is sound *because* every literal
   in the definitions is rational (attack (i) below). This invariant is
   stated nowhere in the .lean file itself (the record's mapping table
   implies it). A future edit inserting a genuinely complex coefficient
   would silently decouple the substitution from conjugation — worth a
   one-line comment if the file is ever extended. Maintenance caveat,
   not an error.

## Claims that survive

### (a) Statement mapping — CONFIRMED, object by object, no deltas

Every definition in `LaurentBlock.lean` lined up against 005 §(ii)-A
Steps 2–4 and `plaw_general.closed_forms()` (read line by line):

| Lean | 005 / plaw_general | verdict |
|---|---|---|
| `sinF z i = (z − z⁻¹)/(2i)`, `cosF z = (z + z⁻¹)/2` | `fsin`/`fcos` (Euler form, coefficient ±1/(2i), 1/2) | match |
| `B = sinF (x*y) i` | B = sin(α+β) = fsin(1,1,0,0) | match |
| `C = sinF y i * x` | C = sin β · e^{iα} | match |
| `A1 = B(1 − y⁻²)` | A₁ = R₀(A) = B(1−y⁻²) | match |
| `mu = y⁻²P⁻²Q²` | μ = e^{−2iaα+2i(b−1)β} = mono(0,−2,−2,2) | match |
| `delta = y⁻¹P⁻¹Q` | δ = e^{i(−aα+(b−1)β)}, δ² = μ | match |
| `w = B(1 − y⁻² + y⁻²P⁻² − μ)` | w = B(1 + P⁻²y⁻² − y⁻² − μ) | match (commuted) |
| `tau = w + mu·w[subst]` | τ = w + μ conj w | match |
| `proj y P Q i v vc = (delta y⁻¹P⁻¹Q⁻¹·v − delta·vc)/(2i)` | p(v) = Im(conj δ·v); note delta y⁻¹ P⁻¹ Q⁻¹ = yPQ⁻¹ = conj δ | match |
| `m = proj (w − τ/2) (…subst)/2` | m = Im(conj δ(w − τ/2))/2 | match |
| `v2 x y P i = B + y⁻²(P⁻²·C[subst] − B)` | v₂ = R₀(Rot_A(2aα)C) = B + y⁻²(P⁻² conj C − B); no Q anywhere (see (d)) | match |
| `reflAB p = (p.2, p.1)` | R₂ = conj | match |
| `reflCA x p = (x²p.2, x⁻²p.1)` | R₁ = e^{2iα} conj | match |
| `reflBC x y i p` | R₀ = B + e^{−2iβ} conj(·−B) | match |
| `rotA P p = (P²p.1, P⁻²p.2)` | Rot_A(2aα) about A = 0 | match (angle 2aα via P², factor-of-2 correct) |
| `rotB x y Q i p` | Rot_B(−2bβ) about B via Q⁻² | match (sign −2bβ correct) |

Composition order checked against 005 Step 1's word-ordered product:
u = R₀ ∘ (R₁R₂)^a ∘ (R₀R₂)^b applies (R₀R₂)^b first — the Lean nesting
`reflBC (rotA P (rotB … p))` and the iterate nesting in
`half_word_letters` both match. `pair_collapse_A/B` are 005 Step 2
exactly (R₁R₂ = Rot_A(2α) = `rotA` at P := x; R₀R₂ = Rot_B(−2β) = `rotB`
at Q := y). The RHS monomials of the identities decode to exactly 005's
angles: `cosF P` = cos(aα), `sinF Q i` = sin(bβ), `cosF (x*y*Q)` =
cos(α+(b+1)β), `sinF (y*Q) i` = sin((b+1)β), `sinF (P*Q⁻¹) i` =
sin(aα−bβ), `sinF (x⁻¹y⁻¹PQ⁻¹) i` = sin((a−1)α−(b+1)β) — each row also
float-verified at 200 random points (Layer 2, worst error 1.5e-15).
I1–I4/D1/D2 are 005's identities verbatim, same LHS orientation
(m − p(C) for I1; p(·) − m for I2–I4), same signs, no permutation. The
I1 = −(S+T), I4 = S−T structure is preserved. `glide_action` is 008 §1's
p∘u = 2m − p; `half_word_sq_is_translation` is Step 3's u² = z + τ with
conj(u(z)) expanded per the pair convention — checked to be the correct
expansion. The record's verbatim theorem block matches the file
character-for-character where quoted.

**Attacks that found nothing:** hunted for a factor-of-2 slip in the
angle conventions (P² vs P in `rotA` — P² is right: Rot_A(2aα) with
P = e^{iaα}); for a sign flip in `rotB` (Q⁻² encodes −2bβ — right); for
an I2/I3 swap (cos(aα)·sin((b+1)β) vs cos((b+1)β)·sin(aα) — not
swapped); for an LHS orientation flip in I1 (it is m − p(C), the others
p(·) − m — matches 005); for a δ vs conj δ swap in `proj` (the first
factor is `delta y⁻¹ P⁻¹ Q⁻¹` = conj δ, correctly multiplying v, not
vc). Each of these would have made Layer 2's instance comparison against
005's prose formulas fail at 1e-15 precision; none did.

### (b) The conjugation encoding — CONFIRMED; nothing smuggled

- The substitution (x,y,P,Q,i) ↦ (x⁻¹,y⁻¹,P⁻¹,Q⁻¹,−i) equals complex
  conjugation on the torus instance because (1) each generator has unit
  modulus there, so inversion = conjugation, and (2) every literal
  coefficient in the file is rational (attack (i): read all defs — the
  constants are 1 and 2 only), and conjugating c(i) ∈ ℚ(i) is i ↦ −i.
  This is exactly 008 §1's star involution ("x→1/x, y→1/y, i→−i is a
  ring automorphism equal to complex conjugation on the torus"),
  extended by P, Q inversion — which 008 had implicitly, since it
  specialized P → x₀^a before starring. Same operation.
- Nothing in the block needs conj as an involution on *arbitrary*
  elements: every conjugate the theorems consume is either the
  substituted form of an explicit generator expression (I1–I4, D1, D2,
  τ, m) or a second pair component that is *universally quantified*
  (`glide_action` takes arbitrary z, zc; `half_word_closed_form` takes
  an arbitrary pair p). That is the safe direction: the Lean theorems
  are componentwise affine identities in the pair, so they are
  *stronger* than the conjugate-pair reading — the intended semantics is
  recovered by instantiating zc := conj z, with no hypothesis anywhere
  demanding p.2 = conj p.1. Nothing is smuggled: the certificate never
  assumes the second component is genuinely a conjugate; it proves
  identities that hold for any second component, of which the conjugate
  is one instance.
- The pair invariant itself (second components of the composite maps =
  substituted first components) is proved by the theorems' RHS shapes,
  and I verified it independently: symbolically in Layer 1 (the second
  components of my recomposed maps match the substituted closed forms
  exactly) and numerically in Layer 2 (`pair C`, `pair w`, `pair v2`
  rows: substituted expression = complex conjugate at the instance,
  ~2e-15).

### (c) The dropped i² = −1 — SAFE; the homogeneity claim independently proven

Direction check: the theorems are ∀-quantified over i ≠ 0, so the
intended instance i := Complex.I is an *instantiation* — the certificate
at the instance can only say more, never less, than 005's statement.
What remained to check is that at i = I the statements ARE 005's
statements: Layer 2 confirms sinF/cosF at the torus instance reduce to
the actual sines and cosines of exactly 005's angles (11 dictionary
rows, worst 1.5e-15), so nothing weaker is being said at the instance —
no identity's two sides differ by an i²-collapsible factor. Stronger:
Layer 1 re-proves every identity with i as a *formal Laurent variable*
in my own ring — an independent proof of 011's homogeneity finding
(each identity is an exact zero polynomial in ℚ[x^±,y^±,P^±,Q^±,i^±];
had any side needed i² = −1, my formal-i check would have failed, since
my ring imposes no relation on i). 011's "discovered because ring_nf
never consumed i² = −1" is thus confirmed by a second, non-Lean method.

### (d) Genericity — CONFIRMED; no silent narrowing

- (a,b) enter the Tier A/B theorems only through the ∀-quantified
  generators P, Q (grep: no natural-number parameter appears outside
  `iterate_collapse_*`, `half_word_letters`, and the `Specialize`
  section). One identity check covers all (a,b) — 005's move, faithfully.
- ℕ vs ℤ: 005 needs integer parameters a ≥ b ≥ 1. `half_word_letters`
  (a b : ℕ) covers the genuine letter words for all a, b ≥ 1 (and adds
  the degenerate a = 0 / b = 0 cases 005 never claims — extra scope in
  the safe direction). The Tier C `*_spec` corollaries at a b : ℤ (zpow)
  cover every integer pair, a superset of 005's range. The specialization
  mechanism differs from 005's (instantiation of ∀-variables vs ring
  homomorphism), as 011's table says; the conclusion is identical, and
  the dangerous direction (formal zero ⟹ geometric zero) is inherited
  from the ∀-statement — 008 §2's soundness argument is not even needed
  in this encoding.
- v₂'s dropped Q: confirmed unused in 005 too. `plaw_general`'s v2 is
  built from monomials with eQ = 0 throughout, and geometrically
  v₂ = R₀(Rot_A(2aα)C) involves no b. I4's RHS still carries Q (through
  p and m) — consistent in both.

### (e) Cheat scan — CLEAN

`grep -inE "sorry|admit|native_decide|decide|axiom|unsafe|partial|implemented_by|notation|macro|elab|opaque|extern"`
over all four project .lean files: only hits are the word "sorry" inside
two doc comments. No local notation, no instance declarations, no
shadowing (all defs live in the fresh `BilliardsFormal.Laurent`
namespace). Hypotheses x, y, P, Q, i ≠ 0 and (2 : K) ≠ 0 are jointly
satisfiable (ℚ, ℂ), satisfied at *every* torus instance (units are
nonzero; char 0), so nothing is vacuous and the char-2 exclusion costs
no intended instance. `BilliardsFormal.lean` imports both `CotMonotone`
and `LaurentBlock`, and my own build log shows all three modules
elaborated — the build cannot pass without elaborating the certificate.
`AxiomCheck.lean` is on disk (010's nit fixed) and enumerates all 23 +
009's 2 theorems; I did not rely on it (see (f)).

### (f) Independent rebuild + axiom audit — GREEN

From `problems/billiards-triangles/formal/`, PATH prefixed with
`$HOME/.elan/bin`:

```
$ rm -rf .lake/build          # project output only; .lake/packages kept
$ lake build
✔ [8655/8658] Built BilliardsFormal.CotMonotone (33s)
✔ [8656/8658] Built BilliardsFormal.LaurentBlock (41s)
✔ [8657/8658] Built BilliardsFormal (7.1s)
Build completed successfully (8658 jobs).
```

Then my own scratch `SkepticAxiomCheck011.lean` (`import BilliardsFormal`
+ `#print axioms` on each of the 23 `Laurent` theorems, the list
enumerated by me from the source): all 23 report exactly
`[propext, Classical.choice, Quot.sound]` — no `sorryAx`, no
`Lean.ofReduceBool` (so no `native_decide` anywhere in the dependency
cone), no project-local axiom. Scratch file deleted after the run.
Independence scope as in 010: fresh session, different agent, project
build products deleted and re-elaborated from committed source; mathlib
oleans from the pinned on-disk cache (the trust base every mathlib user
accepts). Toolchain pins consistent across `lean-toolchain`,
`lakefile.toml`, `lake-manifest.json` (mathlib rev
`905b95818eb32af7874a58b427f50c1711a5e96c` = tag v4.32.2), and both
build logs (attack (j)).

### (g) Record audit — HONEST

- The bridge exclusion is genuine and prominent: 011's "Not attacked,
  deliberately" paragraph, the module docstring's "What is NOT proved
  (the bridge)" section, and the index entry's
  `geometry-bridge-not-formalized` gap all scope it out identically, and
  nothing in the record's Outcome-equivalent sections claims geometric
  content for the certificate.
- The numeric i = 0.7 − 1.3j check is labelled a pre-proof sanity check
  ("the Lean proofs supersede it") — correctly not load-bearing.
- prior-art.json 011 entry: `range` names exactly the formal statements
  and the encoding; `one_line` matches the record (23 theorems, tiers
  A+B+C, i-homogeneity, bridge out of scope, review pending); `gaps`
  listed the two pending pass parts plus the permanent bridge gap —
  the correct under-claiming direction. The two pending gaps are closed
  by this record (index updated; 011's .md untouched); the bridge gap
  stays, as permanent scope.
- Numbers spot-checked: 23 theorems (counted in source: 6 + 11 + 6),
  8658 jobs (my rebuild agrees), zero sorries (kernel + grep agree),
  "~420 lines" vs 423 (nit 1).

### (h) Spot re-derivation — ALL identities re-proven from scratch; D1 also by hand

- **Layer 1** (`lbsk_review.py`, `data/lbsk_review.json`): my own
  7-variable Laurent ring over ℚ (exponent-tuple dict, Fraction
  coefficients, division only by unit monomials — written for this
  review, no code or design shared with `plaw_general.py`, which
  represents coefficients in ℚ(i) with i NON-formal, or with 008's
  `psk_review.py`). The Lean definitions were transcribed from the .lean
  file and every theorem re-checked by exact polynomial subtraction:
  **31/31 exact zero polynomials** — I1–I4, D1, D2, μ = δ², τ ∥ axis,
  μ·conj μ = 1, u² = z + τ (z formal), glide action (z, zc formal),
  both components of the closed form (pair formal), both pair collapses,
  and `half_word_letters` re-composed by genuine symbolic iteration for
  all 0 ≤ a, b ≤ 3 (16 cases, both components). This includes at least
  one of I1–I4 as required — in fact all four.
- **By hand** (computer-free second witness), D1: p(B) − p(C) =
  Im(conj δ·B) − Im(conj δ·C) with conj δ = e^{i(aα−(b−1)β)} gives
  sin(α+β) sin(aα−(b−1)β) − sin β sin((a+1)α−(b−1)β); product-to-sum on
  both, the (b−2)-terms cancel, leaving
  ½[cos((a−1)α−bβ) − cos((a+1)α−bβ)] = sin α · sin(aα−bβ). This is
  005's D1 and Lean's `sinF x i * sinF (P*Q⁻¹) i` exactly.
- **Layer 2** closes the mapping triangle: Lean statement (proven true
  by Layer 1 and the kernel) ↔ 005's prose formulas (I1–I4/D1/D2 trig
  RHS, Step-2 reflections, Step-3 closed forms μ, δ, w, τ, m, the
  p(v) = Im(conj δ v) projection, v₂ = R₀(Rot_A C), and u composed
  directly from 005's geometric maps) — 30 mapped quantities, 200 random
  (a, b, α, β) with a ≤ 12, worst error 4.4e-15.

### Verdict on the three-part FORMALIZED pass (docs/FORMALIZE.md)

1. **Statement review** (the point): DONE by this record, by a
   non-formalizer, definition by definition and theorem by theorem — no
   narrowing, no deltas; the encoding table of 011 is signed off
   unchanged.
2. **Clean build**: CONFIRMED — zero sorries, axioms exactly the
   standard three for all 23 theorems, full source in the repo, verified
   by my own audit from my own theorem list.
3. **Recorded toolchain**: CONFIRMED — pins consistent, and my own
   independent `lake build` (fresh session, build output deleted) is
   green.

**The FORMALIZED status of 011 (= the machine-checked certificate for
005's Laurent-identity block) is complete.** Scope unchanged from 011's
own statement: the ring-model block only. Not claimed by this review:
the geometry bridge (permanently out of the certificate's scope — 005's
Steps 1–3 prose + 008's off-torus recomposition carry it, and Layer 2's
float agreement with 005's geometric maps is a mapping check at 200
points, not a proof); Lemmas C/D, the case tree, and the sufficiency
segments (VERIFIED via 008, not formal); mathlib internals; the Lean
kernel.

## Leads generated

Nothing new beyond 011's own leads, which this review leaves as the
queue state (Lemma D / C-chain formalization; the torus-evaluation
lemma `sinF (exp (θI)) I = sin θ`, which would move Layer 2's float
mapping rows into the kernel; the i-homogeneity Python check of 005's
coefficient dictionaries — note Layer 1 here already proves the
homogeneity *of the identity block*; 011's lead 2 asks the finer
per-coefficient question about `plaw_general`'s ring). One addition:

1. **Comment the rationality invariant** (nit 3): a one-line note in
   `LaurentBlock.lean` that the conj-as-substitution convention requires
   all literal coefficients rational. Definite outcome: the invariant is
   stated where a future editor will see it, or deliberately left to the
   record.

## References

- `problems/billiards-triangles/attempts/011-lean-laurent-block.md` (the
  record under review).
- `problems/billiards-triangles/attempts/005-complete-death-law-theorem.md`
  §(ii)-A Steps 1–4 (the statements of truth), §(i) Reduction (glide
  action use site), Conventions block.
- `problems/billiards-triangles/attempts/008-skeptic-review-of-005.md`
  §§1–2 (star involution, off-torus pair algebra, v₂ identification —
  the second witness).
- `problems/billiards-triangles/attempts/009-…`/`010-…` (lane precedent;
  pass shape).
- `docs/FORMALIZE.md` (the three-part pass this record completes).
- `problems/billiards-triangles/explore/plaw_general.py`
  (`closed_forms`, `cmd_formal` — read line by line; never imported by
  my tool).
- `problems/billiards-triangles/formal/` — `BilliardsFormal/LaurentBlock.lean`
  (every line), `BilliardsFormal.lean`, `AxiomCheck.lean`,
  `lakefile.toml`, `lean-toolchain`, `lake-manifest.json`,
  `lake-build.log`.
- New code and data: `explore/lbsk_review.py`, `data/lbsk_review.json`.
- mathlib4 tag v4.32.2 (rev `905b95818eb32af7874a58b427f50c1711a5e96c`),
  pinned checkout. No papers; no machine transcriptions.
