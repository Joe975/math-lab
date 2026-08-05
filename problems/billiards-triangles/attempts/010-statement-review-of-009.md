# 010 — Statement review + independent rebuild of 009 (Lean pilot, Lemma L1)

- **Problem:** billiards-triangles, `problems/billiards-triangles/PROBLEM.md`
- **Date:** 2026-08-04
- **Mode:** informed (read `prior-art.json`, 009 in full, 005 §(ii)-B and its
  L1 use sites, 008 §front matter for review conventions, `docs/FORMALIZE.md`,
  the four `formal/` project files, `formal/lake-build.log`, and — for the
  definition audit — the vendored mathlib source under
  `formal/.lake/packages/mathlib`)
- **Type:** skeptic review of `009-lean-pilot-cot-lemma.md` (default stance:
  REFUTE). This record is the independent completion of the three-part
  `FORMALIZED` pass that `docs/FORMALIZE.md` defines: (1) statement review by
  a non-formalizer, (2) clean-build confirmation, (3) toolchain/reproduction
  check via my own `lake build`. I wrote none of the Lean code under review.
- **Outcome in one line:** 009 survives on every attacked front. The formal
  statement is 005's Lemma L1 hypothesis-for-hypothesis with no narrowing;
  the source is cheat-free and non-vacuous; my own fresh `lake build`
  (project build output deleted first) is green and my own axiom audit shows
  exactly `[propext, Classical.choice, Quot.sound]`; all three of 005's uses
  of L1 instantiate inside the formalized hypotheses. Two cosmetic nits, no
  corrections. The FORMALIZED pass on 005's Lemma L1 is now **complete**.
- **Tools:** grep over `formal/` and the vendored mathlib (definition audit,
  cheat scan); elan 4.2.3 / Lean `leanprover/lean4:v4.32.2` via
  `$HOME/.elan/bin` for the rebuild (`rm -rf .lake/build && lake build`,
  deterministic, CotMonotone rebuilt from source in 27 s) and the axiom audit
  (a temporary `SkepticAxiomCheck.lean` in `formal/`, deleted after the run;
  exact commands and output below). No new committed code.
- **Sources:** repo + the vendored mathlib source only; no external papers.

No file of 009 (or any prior attempt) was modified.

## Claims attacked

The attack list, (a)–(f) as assigned plus (g)–(h) of my own:

1. **(a) Statement fidelity.** That
   `theorem cot_lemma_L1 {s : ℝ} (hs0 : 0 < s) (hs : s ≤ π / 2) :
   StrictAntiOn (fun c : ℝ => c * Real.cot (c * s)) (Set.Ioc 0 1)`
   says exactly what 005 §(ii)-B's L1 says ("For fixed s ∈ (0, π/2], the map
   c ↦ c·cot(cs) is strictly decreasing on (0, 1]") — variable roles,
   interval endpoints and open/closed-ness, strictness direction, s-range.
2. **(b) Definition audit.** That mathlib's `Real.cot`, `Set.Ioc`,
   `StrictAntiOn` mean what 009's mapping table asserts, checked against the
   vendored mathlib *source* (not docs), and that no junk value of `Real.cot`
   at a pole of sin is relied on inside the claimed domain.
3. **(c) Vacuity / cheat scan.** Every line of `CotMonotone.lean` and
   `BilliardsFormal.lean` read for `sorry`/`admit`/`native_decide`/
   `axiom`/`unsafe`/`partial`/`@[implemented_by]`/custom `notation`/`macro`/
   shadowed `cot`/`π`; that the hypotheses are satisfiable and the domain
   nonempty; that the root module actually imports the theorem file (a
   `lake build` that never elaborates it would pass trivially).
4. **(d) Independent rebuild + axiom audit** (FORMALIZE.md parts 2–3): my own
   `lake build` after deleting the project's own build output, and my own
   `#print axioms` run, not trusting `lake-build.log`.
5. **(e) Record audit.** Whether 009's prose or its prior-art.json entry
   claims more than the certificate shows (in particular anything about
   "retiring Lemmas C and D", which are informal consumers 009 did not
   formalize), and whether `range` is scoped honestly.
6. **(f) Usage check.** That 005 actually leans on L1 *as formalized* — each
   use site instantiates (s, c-pair) inside `0 < s ≤ π/2`, `c ∈ (0, 1]` —
   so the certificate covers a statement 005 really consumes, not a variant.
7. **(g, own) Proof-mathematics hand-check.** The quotient-rule expression,
   the Pythagorean simplification, and the final inequality re-derived by
   hand (the kernel guarantees the formal steps; this checks that what the
   kernel accepted is the argument 005 and 009 describe, i.e. no smuggled
   alternative proof with different domain needs).
8. **(h, own) Toolchain-pin consistency + 009's tooling claims.** That
   `lean-toolchain`, `lakefile.toml`, `lake-manifest.json`, and the log's
   `lean --version` line all pin the same versions 009 reports, and that
   009's "mathlib v4.32.2 has no derivative API for `Real.cot`" is true.

## Refutations found

**None.** No hypothesis mismatch, no narrowing, no cheat, no axiom beyond
the standard three, no vacuity, no overclaim in the record. Two cosmetic
nits, recorded here because a skeptic pass should say what it noticed
(neither is a correction; nothing downstream changes):

1. `formal/lake-build.log` line 8 shows the formalizer's axiom audit run as
   `lake env lean AxiomCheck.lean`, but no `AxiomCheck.lean` is committed —
   it was evidently a scratch file. Reproduction therefore requires
   recreating it (two lines; contents inferable from the log). I re-ran the
   audit myself from a fresh scratch file, so the certificate does not rest
   on that log line.
2. The module docstring of `CotMonotone.lean` (lines 9–11) says the lemma
   "retires Lemmas C and D". As a description of 005's usage this is
   accurate, but read in isolation from the .lean file it could be taken as
   a claim about what is *formalized*. The record 009 itself scopes this
   correctly and explicitly ("Lemmas C and D … remain at `VERIFIED` …, not
   `FORMALIZED`"), so this is a docstring-tone nit only.

## Claims that survive

All of them. Item by item, with what I did to try to break each:

### (a) Statement fidelity — CONFIRMED, no deltas

Checked against 005 §(ii)-B lines 157–158 verbatim:

| 005's L1 | formal | verdict |
|---|---|---|
| "fixed s ∈ (0, π/2]" | `{s : ℝ} (hs0 : 0 < s) (hs : s ≤ π / 2)` — universally quantified over exactly (0, π/2] | match ("fixed s" = ∀ s; two inequalities = the half-open interval, closed at π/2 as 005 needs for Lemma D's endpoint s = bx = π/2) |
| variable roles: c is the moving variable, s the parameter | the lambda binds `c`; `s` is fixed outside `StrictAntiOn` | match — no c/s swap |
| "the map c ↦ c·cot(cs)" | `fun c : ℝ => c * Real.cot (c * s)` | match (radians: `Real.cot`; 005 says "radians throughout this subsection") |
| "strictly decreasing on (0, 1]" | `StrictAntiOn … (Set.Ioc 0 1)` | match — see (b) for both definitions read from source |

No hypothesis is stronger than 005's, no case dropped, no strictness
weakened; the interval is (0, 1] closed at 1 exactly as 005 needs (c = 1 is
one side of all three of 005's comparisons). 009's mapping table is signed
off unchanged.

### (b) Definition audit — CONFIRMED from vendored mathlib source

Read at `formal/.lake/packages/mathlib` (rev `905b958…`, tag v4.32.2):

- `Real.cot` (`Mathlib/Analysis/Complex/Trigonometric.lean:88`) is
  `(Complex.cot x).re` with `Complex.cot z = cos z / sin z` (line 44), and
  `Real.cot_eq_cos_div_sin : cot x = cos x / sin x` (line 641). Standard
  convention, radian.
- Poles: Lean's `x / 0 = 0`, so `Real.cot` is total with junk value 0 where
  sin = 0. On the claimed domain nothing touches a pole: c ∈ (0,1],
  s ∈ (0,π/2] gives cs ∈ (0, π/2] ⊂ (0, π), where the proof itself
  establishes sin > 0 (`Real.sin_pos_of_pos_of_lt_pi`, line 51–52 of
  `CotMonotone.lean`). The one global step, the `funext` rewrite
  `cot u = cos u / sin u`, is a definitional identity valid at every real
  including poles (both sides are the same junk there), and monotonicity is
  only ever asserted on the pole-free `Ioc`. No junk value is load-bearing.
- `StrictAntiOn f s` (`Mathlib/Order/Monotone/Defs.lean:106`) is
  `∀ ⦃a⦄, a ∈ s → ∀ ⦃b⦄, b ∈ s → a < b → f b < f a` — i.e. strictly
  decreasing on s, the direction 005 claims. (Attack considered: a
  StrictAnti/StrictAntiOn or Anti/AntiOn confusion, or ≤ in the conclusion.
  Neither is present.)
- `Set.Ioc a b` (`Mathlib/Order/Interval/Set/Defs.lean:63`) is
  `{x | a < x ∧ x ≤ b}` = (a, b]. So `Ioc 0 1` = (0, 1] and
  `Ioc 0 (π/2)` = (0, π/2], as 009 states.

### (c) Vacuity / cheat scan — CLEAN

- `grep -inE "sorry|admit|native_decide|axiom|unsafe|partial|implemented_by|notation|macro|elab"`
  over both project .lean files: the only hit is the word "sorry" inside
  `CotMonotone.lean`'s doc comment. No `local notation`, no custom `cot` or
  `π`: the file's `open Real Set` resolves `π` to `Real.pi` and every `cot`
  is written qualified as `Real.cot`.
- Vacuity: the hypotheses are jointly satisfiable (e.g. s = 1) and
  `Set.Ioc 0 1` is nonempty, so the theorem has content; `hs : s ≤ π / 2`
  cannot empty anything.
- Elaboration reachability: `BilliardsFormal.lean` is exactly
  `import BilliardsFormal.CotMonotone`; `lakefile.toml`'s default target is
  the `BilliardsFormal` lib rooted there; my own build log (below) shows
  both `BilliardsFormal.CotMonotone` and `BilliardsFormal` compiled. These
  two files are the only .lean sources outside `.lake/` — nothing else to
  audit.

### (d) Independent rebuild + axiom audit — GREEN (FORMALIZE.md parts 2–3)

Exact commands, from `problems/billiards-triangles/formal/`, PATH prefixed
with `$HOME/.elan/bin` (elan 4.2.3):

```
$ rm -rf .lake/build          # project's own output only; .lake/packages kept
$ lake build
✔ [8655/8657] Built BilliardsFormal.CotMonotone (27s)
✔ [8656/8657] Built BilliardsFormal (7.1s)
Build completed successfully (8657 jobs).

$ printf 'import BilliardsFormal\n#print axioms BilliardsFormal.cot_lemma_L1\n#print axioms BilliardsFormal.strictAntiOn_id_mul_cot\n' > SkepticAxiomCheck.lean
$ lake env lean SkepticAxiomCheck.lean
'BilliardsFormal.cot_lemma_L1' depends on axioms: [propext, Classical.choice, Quot.sound]
'BilliardsFormal.strictAntiOn_id_mul_cot' depends on axioms: [propext, Classical.choice, Quot.sound]
$ rm SkepticAxiomCheck.lean   # scratch file, not committed
```

Exactly the three standard axioms, no `sorryAx`, no project-local axiom.
Scope of independence, stated honestly: this is a different agent and a
fresh container session from 009's, the project's own build products were
deleted and re-elaborated from the committed source by the pinned toolchain,
but the mathlib *dependency* oleans came from the same on-disk
`.lake/packages` cache (per the review brief: mathlib was not re-downloaded
or rebuilt). The trust base is therefore the pinned mathlib cache + Lean
kernel, which is the trust base every mathlib user accepts.

### (e) Record audit — HONEST

- 009 flags its own mapping table `PENDING INDEPENDENT STATEMENT REVIEW`,
  states outright that Lemmas C and D "remain at `VERIFIED` …, not
  `FORMALIZED`", and its final section says the pass is incomplete until an
  independent rebuild + statement review — i.e. it under-claims, which is
  the correct direction. (This record is that completion.)
- prior-art.json 009 entry: `range` names exactly the formal statement and
  says C/D and I1–I4 are NOT formalized; `gaps` lists
  independent-statement-review-pending and independent-rebuild-pending.
  Honest; both gaps are closed by this record and removed from 009's index
  entry (index entries are living state; 009's .md is untouched).
- Numbers spot-checked against the log: "compiles in ~8 s" vs log's 7.5 s,
  8657 jobs, timestamps — consistent.
- The one place a hasty reader could over-read is the .lean docstring's
  "retires Lemmas C and D" (nit 2 above); the record itself never makes
  that mistake.

### (f) Usage check — 005 consumes L1 exactly as formalized

All three use sites re-read in 005 §(ii)-B (lines 174–177, 187):

1. Lemma C, first bracket: s = v ∈ (0, v₀], v₀ = π/(2(b+1)) ≤ π/4 < π/2;
   compares c = b/a vs c = 1, with 1 ≤ b ≤ a so b/a ∈ (0, 1]. In-domain.
   (005 writes the bracket "≤ 0, equality iff a = b": when b/a = 1 the two
   c's coincide and the bracket is identically 0 — the strict formal
   statement is exactly what the b/a < 1 case needs, nothing weaker
   suffices, nothing stronger is used.)
2. Lemma C, second bracket: s = π/2 − v ∈ (0, π/2); compares c = 1/a ≤ 1/2
   vs c = 1 (a ≥ 2). In-domain, strict as 005 claims.
3. Lemma D: s = bx with x ∈ (0, π/(2b)], so s ∈ (0, π/2] — this is the use
   that needs the closed endpoint s = π/2; compares c = 1/b vs c = 1
   (b ≥ 2). In-domain.

No other form of L1 (different constants, different interval) appears
anywhere in 005; lines 379 and 447–448 refer back to the same statement. The
formalized statement is the used statement; the certificate is not vacuous.

### (g) Proof-mathematics hand-check — the accepted proof is 005's argument

- `hasDerivAt_id_mul_cot`: quotient rule on (u·cos u)/sin u with
  f′ = cos x − x·sin x, g′ = cos x — the stated derivative expression is the
  textbook (f′g − fg′)/g².
- The `linear_combination`-justified simplification:
  (cos x − x sin x)·sin x − x·cos²x = sin x cos x − x·(sin²x + cos²x)
  = sin x cos x − x. Hand-checked; matches 005's
  [sin(2cs) − 2cs]/(2 sin²(cs)) after the double-angle rewrite, as 009's
  "proof-shape delta" paragraph claims.
- Negativity: sin x·cos x ≤ sin x (sin ≥ 0, cos ≤ 1) and
  `Real.sin_lt : 0 < x → sin x < x` (read from source,
  `Mathlib/Analysis/SpecialFunctions/Trigonometric/Bounds.lean:42`) — the
  same "sin x < x" driver 005 cites.
- The u = c·s transfer in `cot_lemma_L1`: c₁ < c₂ gives c₁s < c₂s
  (`mul_lt_mul_of_pos_right`, s > 0), the core theorem gives
  (c₂s)cot(c₂s) < (c₁s)cot(c₁s), and dividing by s > 0
  (`lt_of_mul_lt_mul_left`) gives the goal. Direction-checked by hand: this
  really is f(c₂) < f(c₁), the StrictAntiOn conclusion, not its reverse.

### (h) Toolchain pins — CONSISTENT; tooling claims TRUE

`lean-toolchain` = `leanprover/lean4:v4.32.2`; `lakefile.toml` requires
mathlib at `rev = "v4.32.2"`; `lake-manifest.json` locks mathlib at
`905b95818eb32af7874a58b427f50c1711a5e96c` with `inputRev v4.32.2`; the
log's `lean --version` says 4.32.2. All four agree with 009's front matter.
009's claim that this mathlib has no derivative API for `Real.cot`:
confirmed by grep over the vendored `Mathlib/Analysis/` — no
`hasDerivAt_cot`/`deriv_cot` for `Real` (only `arccot`/`coth` hits), so the
inline quotient-rule lemma was genuinely necessary.

### Verdict on the three-part FORMALIZED pass (docs/FORMALIZE.md)

1. **Statement review** (the point): DONE by this record, by a
   non-formalizer, hypothesis by hypothesis — no narrowing, no deltas.
2. **Clean build**: CONFIRMED — zero `sorry`, axioms exactly
   `[propext, Classical.choice, Quot.sound]`, full source in the repo,
   verified by my own audit, not the formalizer's log.
3. **Recorded toolchain**: CONFIRMED — pins consistent across all four
   files, and my own independent `lake build` (fresh session, project build
   output deleted) is green.

**The FORMALIZED status of 009 (= the machine-checked certificate for 005's
Lemma L1) is complete.** Scope unchanged from 009's own statement: exactly
`cot_lemma_L1` and its core `strictAntiOn_id_mul_cot`; Lemmas C, D and the
I1–I4 Laurent block remain informal (`VERIFIED` via 008).

What is *not* claimed by this review: no re-verification of 005's Lemmas
C/D or anything else in 005 (008 did that); no audit of mathlib itself
beyond reading the definitions used; no rebuild of mathlib from source; no
check of Lean-kernel soundness.

## Leads generated

Unchanged from 009's leads 2–4 (Lemma D next, then I1–I4, then Lemma C) —
this review adds no new ones and found nothing that spawns any. One
process lead:

1. **Commit the axiom-audit scratch file pattern.** Add the two-line
   `AxiomCheck.lean` to `formal/` (it is part of the reproduction recipe
   the log references but the repo does not contain). Definite outcome:
   `lake env lean AxiomCheck.lean` works from a clean checkout with no
   file-recreation step, or the file is deliberately left out and the log
   line should carry the file's contents instead.

## References

- `problems/billiards-triangles/attempts/009-lean-pilot-cot-lemma.md` (the
  record under review).
- `problems/billiards-triangles/attempts/005-complete-death-law-theorem.md`
  §(ii)-B (Lemma L1 statement, lines 157–161; use sites lines 174–177, 187).
- `problems/billiards-triangles/attempts/008-skeptic-review-of-005.md`
  (prior hand verification of L1; review-shape conventions).
- `docs/FORMALIZE.md` (the three-part pass this record completes).
- `problems/billiards-triangles/formal/` — `BilliardsFormal/CotMonotone.lean`,
  `BilliardsFormal.lean`, `lakefile.toml`, `lean-toolchain`,
  `lake-manifest.json`, `lake-build.log`.
- mathlib4 tag v4.32.2 (rev `905b95818eb32af7874a58b427f50c1711a5e96c`),
  read from the vendored copy under `formal/.lake/packages/mathlib`:
  `Analysis/Complex/Trigonometric.lean` (cot), `Order/Monotone/Defs.lean`
  (StrictAntiOn), `Order/Interval/Set/Defs.lean` (Ioc),
  `Analysis/SpecialFunctions/Trigonometric/Bounds.lean` (sin_lt). No
  papers; no machine transcriptions.
