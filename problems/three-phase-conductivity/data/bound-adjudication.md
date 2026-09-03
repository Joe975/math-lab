# Adjudication: our laminate vs the transcribed Cherkaev 2009 B2 bound

**STATUS: DONE (best effort). Verdict: (d) — the discrepancy survives every
check attempted here and looks genuine, but is NOT proven, and one specific
step in Cherkaev's derivation (the "localized polyconvexity" structural-
variation argument, Remark 4.1) could not be independently re-derived within
this pass. See "What remains unresolved" for exactly what would close this.**

## The question

At sigma = (1,2,5), m2 = 1/4 (so sqrt(m2) = 1/2 is rational and everything is
exact in Q), laminates built in this repo sit strictly BELOW the transcribed
Cherkaev 2009 Theorem 7.1 "B2" lower bound at m1 just under m11 = 1/8:

| m1 | our structure | HS_lo | transcribed B2 | ours - B2 |
|---|---|---|---|---|
| 31/250 | 3.0053470065 | 3.0053404539 | 3.0053523724 | -5.4e-06 |
| 3/25 | 3.0270098722 | 3.0268456376 | 3.0271504085 | -1.4e-04 |
| 11/100 | 3.0831692377 | 3.0816326531 | 3.0845383760 | -1.4e-03 |

A structure below a valid lower bound is impossible, so something has to be
wrong somewhere. Default stance going in: assume our side is wrong.

## Sources fetched and read directly (not from memory, not trusting the two
prior transcriptions blindly)

- A. Cherkaev, "Bounds for effective properties of multimaterial
  two-dimensional conducting composites, and fields in optimal composites"
  (Mech. Mater. 41 (2009) 411-433), preprint PDF fetched fresh from
  `http://www.math.utah.edu/~cherk/publ/newbounds7.pdf`, extracted page-by-page
  with `pypdf` (52 pages, full text). This is the SAME preprint the two
  earlier workers used, but re-fetched and re-extracted independently rather
  than reusing their quoted excerpts, specifically to catch OCR artifacts they
  might have silently corrected without flagging (see finding 0 below — this
  paid off).
- WebSearch pass for an erratum/correction on this bound: **none found.** No
  secondary source (including Cherkaev-Zhang 2011, arXiv:1009.3060, and the
  2025 Fel duality paper already on file) flags a problem with Theorem 7.1.

## Task 1 — independent re-transcription of Theorem 7.1 and eqs 7.5-7.24

**Result: character-for-character match to both prior transcriptions**, with
one raw-OCR artifact worth recording so nobody re-discovers it and panics.

The extracted text of eq. (7.15) literally reads:

```
B(m1, m2) = { B1 if m11 <= m1 <= 1 ; B2 if m11 <= m1 <= m12 ; B3 if 0 <= m1 <= m12 }
```

Read literally, the B2 case is `m11 <= m1 <= m12`, which is backwards from
both prior transcriptions (`m12 <= m1 <= m11`) and would be an empty interval,
since the paper itself proves `m12/m11 <= 1` two paragraphs earlier and the t0
piecewise formula (7.4), extracted from the *same* page without ambiguity,
unambiguously reads `Z1 if m12 <= m1 <= m11`. This is a pdftotext linearization
artifact from the stacked LaTeX `\begin{cases}` block, not a real difference —
both prior transcriptions had it right. Flagging only so a future re-transcription
attempt doesn't get spooked by the raw extraction and doesn't need to redo this
check.

Formulas B1, B2, B3, Z5, Z6, Z7, m11, m12, Z0-Z4 all match the repo's
`data/literature-check.md` transcription exactly, term for term. **No
transcription error in the closed-form B2 formula itself.**

### A second, more consequential OCR artifact — resolved, not a bug

Eq. (7.2) extracts as:

```
where
1
H1(t) = m1/(2k1) + m2/(k2+t) + (m1(k1-t)+2k1m3)^2 / [2k1(2k1m3(k3+t)+m1(k1^2-t^2))]
```

Note the stray `1` on its own line, immediately before `H1(t) = ...`. Cross-checking
against Section 5's general N-material derivation (eq. 5.8: `H1 = 1 / (M^T(R1+Y1 P^T)^-1 M)`,
eq. 5.10: an explicit formula for **`1/H1`**, not `H1`) shows this stray `1` is
the numerator of a stacked fraction `1/H1(t) = ...` that pdftotext split across
lines. So eq. (7.2) actually defines `1/H1(t)`, and eq. (7.1)'s `kL = max_t(-t + H1(t))`
uses the true reciprocal `H1(t) = 1 / [that sum]`, exactly matching the pattern of the
classical two-phase HS formula `HS_lo = -k1 + 1/(m1/(2k1)+m2/(k1+k2))` (harmonic-mean-of-resistances
form). **Neither prior transcription in this repo made this mistake** — both
correctly used `B2 = k2 + (1-sqrt(m2))^2 * Z5/Z6` (Z5/Z6, not the raw sum), so this
was never live in the repo's actual comparison; it is flagged here only because
it is exactly the kind of silent artifact that would corrupt a from-scratch
re-implementation of eq. (7.1)-(7.2) done by literally trusting the extracted
label. Verified numerically (Task 2 below) that the correct reciprocal
interpretation is what reproduces Cherkaev's own closed-form B1/B2/B3 exactly;
the literal (non-inverted) reading does not.

## Task 2 — independent maximization of the translation family, bypassing the
closed-form split entirely

Implemented eq. (7.1)-(7.2) fresh in `sympy`/`Fraction`, with `sqrt(m2)` exact
since m2 = 1/4, and maximized `-t + 1/H1(t)` over `t` by solving `d/dt = 0`
symbolically and filtering real roots in `[k1, k2]` (plus checking both
endpoints), independent of the m11/m12 region-split machinery entirely:

| m1 | max_t(-t+1/H1(t)) | transcribed closed-form B2 |
|---|---|---|
| 31/250 | 3.0053523724414615 | 3.0053523724 |
| 3/25   | 3.0271504084574725 | 3.0271504085 |
| 11/100 | 3.084538375973304  | 3.0845383760 |

Exact agreement (to the 10 digits reported in the original table) at all three
points. **This rules out a Maple-simplification error in going from (7.1)-(7.2)
to the closed form B2** — the closed form is not just correctly transcribed,
it is correctly *derived* from the stated translation family, independent of
any algebra Cherkaev's Maple did.

Sanity check on the optimization method itself: tried a naive full-grid
maximization (not filtering by the sign of the denominator
`2k1m3(k3+t)+m1(k1^2-t^2)`) on a near-two-phase test case (m3 -> 0), and it
finds a *spurious* larger "maximum" near t=k2 where that denominator goes
negative — an artifact of a domain violation, not a real critical point. The
critical-point-solve method used for the three real points above lands well
inside the region where the denominator is positive (checked explicitly, e.g.
at m1=0.12, t*=1.1388: denominator ≈ 15.4 > 0), so this artifact does not
affect the three reported values, but it's a live risk for anyone
re-implementing this by grid search rather than solving the first-order
condition.

## Task 3 — from-scratch audit of our own structure, sharing no code with
`harness/three-phase-conductivity/laminate.py` or `verify_laminate.py`

Re-derived the 2D laminate homogenization formula from physical first
principles (not read from either harness file): in the frame aligned with
layering normal `n` and tangent `t`, continuity of tangential field and
normal current across the interface gives, for two sub-tensors
`K_i = [[a_i,b_i],[b_i,c_i]]` in the `(n,t)` basis with fractions `f_i`,

```
J_n(E_n,E_t) solves  f1(J_n - b1 Et)/a1 + f2(J_n - b2 Et)/a2 = En
E_{n,i} = (J_n - b_i Et)/a_i
J_t = f1(b1 E_{n,1} + c1 Et) + f2(b2 E_{n,2} + c2 Et)
```

applied recursively bottom-up through the tree (all normals in the saved
record are axis-aligned, `(1,0)` or `(0,1)`, so no rotation/sqrt is needed —
everything stays in exact `Fraction` arithmetic throughout). Sanity-checked
the formula against the textbook 2-phase axis-aligned laminate identities
(series/harmonic mean normal to layers, parallel/arithmetic mean along
layers) before trusting it on the real tree.

Ran this against `data/attained/below-m11-f1-0.12.json` (the only saved,
`verify_laminate.py`-passing structure among the three table rows; the other
two rows in the adjudication table have no saved record to audit — see "Not
checked" below):

- **Volume fractions**: recomputed `{p1: 3/25, p2: 1/4, p3: 63/100}` — bit-for-bit
  identical to the claimed `fractions` field.
- **Tensor**: recomputed `(a, b, c)` — `a` and `c` are bit-for-bit identical to
  the claimed `tensor` field (both are the same huge exact rationals down to
  the last digit), `b = 0` exactly.
- **Isotropy**: `(a-c)^2 + b^2 ≈ 1.57e-133` — even smaller than the record's
  claimed "< 1e-61", i.e. `|a-c|` is on the order of `1e-66`. This is not
  exact isotropy (`a != c` exactly, as fractions), but the deviation is so far
  below any of the reported discrepancies (`5.4e-6` to `1.4e-3`) that it
  cannot be the explanation (see Task 4).

**Verdict on Task 3: the structure's tensor and fractions check out under a
genuinely independent implementation.** No computational bug found.

**Not checked**: the m1 = 31/250 and m1 = 11/100 rows in the table have no
saved JSON record in `data/attained/`, and `explore/tp_attain.py` (the only
generator script found) builds structures at exactly `f1 = m11`, not below it
— so whatever produced those two rows is not currently reproducible from the
repo as it stands. Flag this as a genuine gap: **only the middle row (m1 =
3/25 = 0.12) has been independently audited end to end.** The other two rows
should be treated as unverified pending either locating/regenerating their
trees or re-deriving them via the same family.

## Task 4 — scope check: does Theorem 7.1 require exact isotropy or another
hypothesis our instance might violate?

Read Section 2.2 (definition of K*) and Section 5 (the general N-material
bound Theorem 5.1) closely:

- **K* is defined via `J(e0,χ) = ½ Tr[K*(χ) e0 e0^T]`** (eq. 2.16) using two
  *orthogonal* applied fields simultaneously (eq. 2.13-2.15). Theorem 5.1's
  `kL` is proved as a bound on this same energy functional for **"any
  N-material composite that satisfies (5.18)"** — no isotropy hypothesis
  appears in Theorem 5.1 itself. Taking `e0 = I` (equal orthogonal loadings),
  `J(I,χ) = ½ Tr[K*] = (a+c)/2`, so the underlying general bound applies to
  `(a+c)/2` for *any* admissible composite, isotropic or not. Theorem 7.1
  ("isotropic composite... k* ≥ B(m1,m2)") is the isotropic corollary, where
  `(a+c)/2 = a = c = k*` collapses to a single number — but the general
  (non-isotropic) version, if valid, already forces `(a+c)/2 ≥ kL` regardless.
  Since our structure's `a` and `c` agree to ~66 digits, `(a+c)/2` and the
  reported single value are identical to far more precision than the
  discrepancy — **the "not exactly isotropic" reading does not rescue the
  bound**: even the maximally charitable trace-based quantity sits below B2 by
  the same ~1.4e-4 at m1 = 0.12.
- **Condition (5.18)** (positive-definiteness of `R_r + Y_r P^T`, required for
  the minimization step to be valid) — per **Remark 5.2**, "for three-material
  mixtures, either (5.18) is satisfied, or Hashin-Shtrikman bound holds."
  Checked the companion sufficient condition (5.5), `m3 >= m1 * 2t(t-k1) /
  [(k1+t)(k3+t)]` at `t=k2` (its stated tightest point), at all three m1
  values: RHS is `~0.021-0.024`, `m3` is `~0.63-0.64` — satisfied with huge
  margin, nowhere near the boundary where (5.18) could fail.
- **k3 = ∞ special case**: not applicable — our k3 = 5 is finite, so Theorem
  7.1 (not the k3=∞ Theorem 7.2) is the right theorem, matching what both
  prior transcriptions used.
- **Phase ordering k1 < k2 < k3**: satisfied (1 < 2 < 5).

**No scope/hypothesis mismatch found in the parts of the derivation that were
checked.** The one thing that could not be ruled out is buried one level
deeper — see below.

## What remains unresolved — the one place a genuine error could still hide

Theorem 5.1's general bound is built in two stages: a *basic* field constraint
`S^2+V^2 >= D^2` (eq. 4.1) that Cherkaev states explicitly holds "for all
structures, whether they are optimal or not" — this is the Nesi-level
constraint and is not in question. B2, however, needs the *tighter* bound from
"localized polyconvexity" (Section 4.1-4.2), which plugs in a sharper
pointwise field constraint `D^2 <= Theta_i(S)` derived in **Remark 4.1** via a
**structural-variation argument**: swap two infinitesimal inclusions of
materials `i` and `j`, argue the energy change must be non-negative "if the
tested configuration is optimal," and extract an inequality on the fields from
that non-negativity.

This is a standard technique in this literature (first-order necessary
conditions of a minimizing sequence, used to derive a lower bound valid for
*all* admissible structures via "the global minimizer is in particular
locally optimal, so it satisfies this condition; therefore its energy — which
lower-bounds every other structure's energy — satisfies kL"). The logic, if
executed correctly, does produce a bound valid for every microstructure, not
just optimal ones. But it is a real argument with real hypotheses (about
which infinitesimal swaps are admissible, connectedness/accessibility of each
phase at every point, etc.), and **it was not independently re-derived here**
— doing so would mean redoing Cherkaev's Section 4 calculation of `Theta_i(S)`
from the structural-variation argument itself, which is a multi-page
derivation, not something a numeric check can validate.

This is the one candidate explanation not yet excluded: that the
structural-variation necessary condition, as derived, implicitly assumes
something about microstructure topology (e.g. simple/coated-ellipse-type
local geometry) that a deep, six-level hierarchical laminate — one that was
specifically engineered to carry an extreme, non-generic field split in phase
1 (attempt 004 reports mirror-image fields `(10/3,2/3)` and `(2/3,10/3)` in
phase 1, saturating the field-variance floor from attempt 002/003's own
diagnostics) — does not actually satisfy, even though it is a perfectly valid
admissible composite by the basic physics (harmonic potentials + interface
jump conditions, which is all Task 3 checked).

## Verdict

**(d), tentatively, with real uncertainty remaining** — not (a): the B2
formula is transcribed correctly and independently reproduced by direct
maximization of the stated translation family. Not (b): the structure's
tensor and volume fractions are independently confirmed exact at the one
point (m1 = 3/25) that has a saved, auditable record. Not (c) on every
hypothesis checked: no isotropy requirement in the underlying general bound
that our near-isotropic (residual ~1e-133) structure could plausibly violate
at the observed scale, (5.18)/(5.5) hold with large margin, phase ordering
and finite-k3 are both fine.

What is genuinely NOT settled: whether Cherkaev's "localized polyconvexity"
field constraint (Remark 4.1, the structural-variation step that produces the
tighter-than-Nesi bound) extends validly to arbitrary hierarchical
microstructures, or implicitly assumes something our construction's extreme,
engineered field distribution violates. This is the load-bearing step behind
B2 and B3 (not behind B1 = plain HS, which is unaffected) that this pass could
not re-derive from first principles in the time available.

**To be sure, in either direction, the following would be needed:**

1. **Independently re-derive Theta_i(S) from Remark 4.1's structural-variation
   argument** (Section 4.4, "Extremal constraints," not yet read closely here)
   and check whether its derivation implicitly restricts to a class of
   microstructures (e.g. requires each phase to be locally "swappable" at
   every point of every other phase's domain, or bounds the *local* field
   variation within Ω1 in a way a rank-6 laminate's phase-1 sublaminate could
   violate) that this repo's hierarchical laminates fall outside of.
2. **Audit the remaining two table rows** (m1 = 31/250 and m1 = 11/100) with
   the same from-scratch method as Task 3 — right now only the middle point is
   independently confirmed; if either of the other two turns out to be a
   genuine computational error, that would change the picture (though it
   would not explain the middle point, which is solid).
3. **A second, structurally different valid lower bound** at the same three
   points — e.g. a correctly-transcribed Nesi (1995) three-phase formula for
   general (non-symmetric) fractions, which neither this pass nor the two
   priors could locate in accessible form — to see whether it is also
   violated (would strengthen "genuine discrepancy") or sits above our
   structure's value but below B2 (would sharpen exactly where the error, if
   any, is localized).
4. Ideally, contact with someone who has independently implemented Cherkaev
   2009's translation bound (e.g. via the wheel-assemblage or Cherkaev-Zhang
   2011 follow-up work) to compare numerically at these exact rational
   points — the fastest real external check, not attempted here.

**Until (1)-(3) are done, per PROBLEM.md's verification contract nothing
should be recorded as killed against B2 in either direction**, and any new
attempt record building on the below-m11 structures should cite this file and
explicitly flag the comparison against B2 as unresolved rather than as either
"beats the published bound" or "invalid."
