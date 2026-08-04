# 009 — Lean 4 formalization pilot: Lemma L1 of 005 machine-checked

- **Problem:** billiards-triangles, `problems/billiards-triangles/PROBLEM.md`
- **Date:** 2026-08-04
- **Mode:** informed (read `prior-art.json`, 005 §(ii)-B and front matter,
  008 §3, `docs/FORMALIZE.md`, `STATUS.md` queue item 16, `AGENTS.md`,
  `CONTRIBUTING.md`)
- **Type:** formalization — review-shape attempt against
  `005-complete-death-law-theorem.md` (verifies 005's Lemma L1 with a
  kernel-checked Lean 4 certificate). This is the lab's first run of the
  `docs/FORMALIZE.md` lane (queue item 16, first half; the Laurent-identity
  block I1–I4 is NOT attempted here).
- **Tools:** elan 4.2.3; Lean 4 toolchain `leanprover/lean4:v4.32.2`;
  mathlib pinned at tag `v4.32.2`
  (rev `905b95818eb32af7874a58b427f50c1711a5e96c`); project at
  `problems/billiards-triangles/formal/` (`lakefile.toml`, `lean-toolchain`,
  `lake-manifest.json`, source `BilliardsFormal/CotMonotone.lean`).
  Prebuilt mathlib via `lake exe cache get` (nothing of mathlib rebuilt from
  source). Build is deterministic; the pilot's own module compiles in ~8 s;
  build + axiom-audit output captured in `formal/lake-build.log`.
- **Sources:** none external beyond mathlib itself. Repo: 005 (the record
  whose lemma is formalized), 008 (its skeptic review, §3 CONFIRMED),
  `docs/FORMALIZE.md` (the lane definition, including the bridge caveat this
  record's mapping section exists for).

## Claims attacked

1. **005 Lemma L1** (§(ii)-B, radians): "For fixed s ∈ (0, π/2], the map
   c ↦ c·cot(cs) is strictly decreasing on (0, 1]." Already
   skeptic-confirmed by hand in 008 §3; the attack here is the stronger
   one the FORMALIZE lane defines: state it formally and make Lean's kernel
   accept a proof. A false or mis-stated lemma would have refused to close.
2. The **informal→formal bridge** for that statement: that the formal
   theorem proved is the record's sentence and not a narrowing (the failure
   surface `docs/FORMALIZE.md` warns about). Attacked by a
   hypothesis-by-hypothesis mapping, written out below — but per the lane,
   the mapping must be re-checked by an independent reviewer, so it is
   flagged **PENDING INDEPENDENT STATEMENT REVIEW**.

Not attacked: Lemmas C and D themselves (L1's consumers), the I1–I4
Laurent block, and everything else in 005. No claim of 005 is re-derived
informally here; 008 already did that.

## Refutations found

None. The lemma formalized exactly as stated in the record: no hypothesis
had to be strengthened, no case dropped, no strictness weakened. (One
non-finding worth recording: mathlib at v4.32.2 has `Real.cot` and
`Real.cot_eq_cos_div_sin` but no derivative API for `Real.cot`, so the
quotient-rule step is done inline — see "what was needed" below. That is a
tooling observation, not a gap in 005.)

## Claims that survive

### The formal statement proved (verbatim from `formal/BilliardsFormal/CotMonotone.lean`)

```lean
theorem cot_lemma_L1 {s : ℝ} (hs0 : 0 < s) (hs : s ≤ π / 2) :
    StrictAntiOn (fun c : ℝ => c * Real.cot (c * s)) (Set.Ioc 0 1)
```

proved with zero `sorry`, via the intermediate theorem (same file)

```lean
theorem strictAntiOn_id_mul_cot :
    StrictAntiOn (fun u : ℝ => u * Real.cot u) (Set.Ioc 0 (π / 2))
```

Axiom audit (in `formal/lake-build.log`): both theorems depend on exactly
`[propext, Classical.choice, Quot.sound]` — the three standard Lean/mathlib
axioms, nothing else, no `sorryAx`.

### Mapping to 005's informal statement — **PENDING INDEPENDENT STATEMENT REVIEW**

005 §(ii)-B: "**Lemma L1.** For fixed s ∈ (0, π/2], the map c ↦ c·cot(cs)
is strictly decreasing on (0, 1]." Hypothesis by hypothesis:

| 005's informal statement | formal counterpart | delta |
|---|---|---|
| "fixed s ∈ (0, π/2]" | `{s : ℝ} (hs0 : 0 < s) (hs : s ≤ π / 2)` | none: (0, π/2] written as the two inequalities |
| "the map c ↦ c·cot(cs)" | `fun c : ℝ => c * Real.cot (c * s)` | none; `Real.cot` is radian, matching 005's "radians throughout this subsection" |
| "strictly decreasing on (0, 1]" | `StrictAntiOn … (Set.Ioc 0 1)` | none: `StrictAntiOn f I` is mathlib's "∀ c₁ c₂ ∈ I, c₁ < c₂ → f c₂ < f c₁", and `Set.Ioc 0 1` is (0, 1] |
| cot well-defined on the range | no side condition needed | for c ∈ (0,1], cs ∈ (0, π/2] where sin > 0; `Real.cot` is total (junk value where sin = 0) but the proof never evaluates it there |

Deliberate proof-shape delta (statement unchanged): the Lean proof first
proves the substituted form `strictAntiOn_id_mul_cot` — u·cot u strictly
decreasing on (0, π/2] — and derives L1 by u = c·s (rescaling by the
positive constant s). That is 005's own proof reparametrized: 005
differentiates in c and gets d/dc[c·cot(cs)] = [sin(2cs) − 2cs]/(2sin²(cs));
the Lean proof differentiates in u = cs and gets
(sin u·cos u − u)/sin²u = [sin(2u) − 2u]/(2 sin²u), the identical
expression. The inequality driving it is the same `sin x < x for x > 0`
(mathlib `Real.sin_lt`). The theorem exported as `cot_lemma_L1` is the
record's statement, not the substituted one.

What L1's certificate does **not** cover (unchanged from 005/008's own
scoping): Lemmas C and D are informal consumers of L1 — their derivations
(log-derivative split, endpoint identity, factorization; Lemma D's
b·cot(bx) − cot x < 0) remain at `VERIFIED` (008 §3 hand re-derivation),
not `FORMALIZED`. Formalizing them is the natural next increment; see
leads.

### Toolchain and reproduction

- `formal/lean-toolchain`: `leanprover/lean4:v4.32.2`
- mathlib: tag `v4.32.2` = rev `905b95818eb32af7874a58b427f50c1711a5e96c`
  (pinned in `formal/lakefile.toml`, locked in `formal/lake-manifest.json`)
- Commands, from `problems/billiards-triangles/formal/`:
  `lake exe cache get && lake build` — result:
  `Build completed successfully (8657 jobs)`, log in `formal/lake-build.log`
  together with the `#print axioms` audit of both theorems.
- Per `docs/FORMALIZE.md`, CI does not build Lean; the reproducibility
  contract is a skeptic's own independent `lake build` on another machine.
  That, plus the statement review above, completes the FORMALIZED pass —
  neither is done by this record.

### Pilot verdict on the lane itself (queue item 16 asked for this)

The obstruction-vs-tooling question has an answer: **no obstruction, either
mathematical or tooling, for a lemma of this class.** Toolchain install to
green build was well under the timebox: elan install ~1 min, mathlib clone +
`lake exe cache get` ~25 min (≈7 GB in `.lake/`, git-ignored), the proof
itself compiled on the first `lake build`. The one friction point: mathlib
v4.32.2 has no `HasDerivAt` lemma for `Real.cot`, so the derivative is
computed inline by the quotient rule on cos/sin
(`hasDerivAt_id_mul_cot` in the file) — three lines, not a blocker. The
useful mathlib surface: `strictAntiOn_of_deriv_neg`,
`Real.hasDerivAt_sin/cos`, `HasDerivAt.mul/div`, `Real.sin_lt`,
`Real.cot_eq_cos_div_sin`, `linear_combination` for the Pythagorean
rewrite.

## Leads generated

1. **Complete the FORMALIZED pass on this record.** An independent agent
   (not this one): (a) re-runs `lake build` in a fresh container from the
   committed source (definite outcome: builds green or does not), and
   (b) performs the statement review of the mapping table above against
   005 §(ii)-B (definite outcome: signs off or files a delta). Until both,
   the certificate is a formalizer-claimed one.
2. **Formalize Lemma D** (sin(bx)/sin x strictly decreasing on (0, π/(2b)],
   b ≥ 2): it is a two-line consumer of L1 (`b·cot(bx) − cot x < 0` via L1
   at s = bx, c = 1/b vs c = 1) plus a log-derivative wrapper. Definite
   outcome: compiles against `cot_lemma_L1` or exposes a domain subtlety in
   the b·x substitution. Cheap; would make the first *chain* L1 → D formal.
3. **Queue item 16, second half — the Laurent-identity block I1–I4.**
   005 proves these by exact polynomial subtraction in
   ℤ[x, x⁻¹, y, y⁻¹, S, T]; the Lean analogue is an identity of Laurent
   polynomials, decidable by `decide`/`ring` over ℤ. The geometry bridge
   (u = R₀·Rot_A(2aα)·Rot_B(−2bβ)) is exactly what the statement-review
   step must scrutinize there — the formal check cannot see a wrong bridge
   (STATUS insight 2026-07-31). Definite outcome per identity: `ring`
   closes it or the formalization surfaces a mismatch with
   `plaw_general.py`'s ring.
4. **Lemma C** is the hard one of the pair (endpoint identity
   π/2 = (b+1)v₀, factorization with sign analysis of each factor); attempt
   only after 2 establishes the substitution pattern. Definite outcome:
   `H2 > 0` on (0, v₀) formal, or a pinned list of which of 005's five
   inequality steps resists mathlib's API.

## References

- `problems/billiards-triangles/attempts/005-complete-death-law-theorem.md`
  §(ii)-B (Lemma L1, statement and informal proof; the record verified).
- `problems/billiards-triangles/attempts/008-skeptic-review-of-005.md` §3
  (hand re-derivation of L1/C/D, 14,183-point adversarial numerics — the
  skeptic-confirmation that `docs/FORMALIZE.md` requires before
  formalizing).
- `docs/FORMALIZE.md` (the lane: qualification, bridge caveat, three-part
  pass).
- mathlib4, tag v4.32.2, github.com/leanprover-community/mathlib4
  (rev `905b95818eb32af7874a58b427f50c1711a5e96c`). No papers consulted;
  no machine transcriptions.
