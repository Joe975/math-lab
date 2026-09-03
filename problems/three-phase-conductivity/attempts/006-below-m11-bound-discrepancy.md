# 006 — Laminates below Cherkaev's B2 bound: a systematic discrepancy, REQUIRING ESCALATION

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (follows 004/005)
- **Type:** computational finding + primary-source check; **unresolved**
- **Tools:** `explore/tp_cherkaev_bound.py` (the transcribed bound and its
  consistency tests, `--selftest`), `explore/tp_fields.py`, `explore/tp_attain.py`,
  `explore/tp_below_m11.py` (rebuilds every tabulated row, `--selftest`,
  `--write`), the harness. Records: `data/attained/below-m11-m1-*.json`, all
  three passing `verify_laminate.py`; `data/cherkaev-primary-excerpts.md`;
  `data/bound-adjudication.md`.
- **Sources:** A. Cherkaev, "Bounds for effective properties of multimaterial
  two-dimensional conducting composites, and fields in optimal composites",
  preprint dated 2008-05-06 of Mech. Mater. 41 (2009) 411–433, obtained from
  `math.utah.edu/~cherk/publ/newbounds7.pdf` and **text-extracted and read
  here** — Theorem 7.1, eqs (7.1)–(7.18), Section 4.1. Marked [T] because the
  extraction is mechanical, but this is the primary document, not a secondary
  quotation.

## Approach

004 built a laminate attaining HS_lo at f₁ = m₁₁. The obvious next question was
whether anything attains *below* m₁₁. Cherkaev's Theorem 7.1 answers it: below
m₁₁ his bound B2 is strictly greater than HS, so HS is not attainable there.
Checking our own structures against B2 was meant to be a routine consistency
test of the whole picture. It failed, and the failure did not go away under
every check applied to it.

## What was done

### 1. The discrepancy

Structures from 004's family, with f₁ pushed just below m₁₁, sit **strictly
below B2**. At σ = (1,2,5), m₂ = 1/4 (so √m₂ = 1/2 and everything is exact in ℚ):

| m₁ | our structure | HS_lo | B2 | ours − B2 |
|---|---|---|---|---|
| 31/250 | 3.0053470065 | 3.0053404539 | 3.0053523724 | −5.4e−06 |
| 3/25 | 3.0270098722 | 3.0268456376 | 3.0271504085 | −1.4e−04 |
| 11/100 | 3.0831692377 | 3.0816326531 | 3.0845383760 | −1.4e−03 |

**It is systematic, not a one-off.** Across four conductivity triples the same
family lands consistently between HS and B2, capturing roughly half the claimed
improvement, with (B2−HS)/(ours−HS) settling near 1.89:

| σ | m₁₁ | m₁ | ours − HS | B2 − HS | ratio |
|---|---|---|---|---|---|
| (1,2,5) | 0.125 | 0.100 | 4.524e−03 | 8.566e−03 | 1.893 |
| (1,3,7) | 0.0833 | 0.0667 | 6.664e−03 | 1.168e−02 | 1.753 |
| (1,4,9) | 0.0625 | 0.0500 | 8.719e−03 | 1.476e−02 | 1.693 |
| (2,5,11) | 0.0952 | 0.0762 | 8.713e−03 | 1.646e−02 | 1.889 |

A structure below a valid lower bound is impossible, so something is wrong.

### 2. Everything checked on our side

- Exact rational volume fractions hitting the target exactly.
- Off-diagonal **exactly 0**; diagonal difference below 1e−61, so both
  eigenvalues coincide. Genuinely isotropic, not "isotropic in the trace".
  (This was my leading suspicion — that I was comparing the mean of an
  anisotropic tensor's eigenvalues against an isotropic bound. It is not that.)
- Both harness algebra routes agree (projection form vs rotated-frame interface
  averaging), and `verify_laminate.py` passes from the command line on the
  saved record.
- Keller–Dykhne duality holds.
- An **explicit admissible-field certificate**: `tp_fields` produces the
  piecewise-constant field, and its average equals E₀, its flux average equals
  σ\*E₀, and the energy identity holds, all exactly. That is the definition of
  the effective tensor, established by a route independent of the lamination
  algebra.

### 3. Everything checked on the bound's side

- **Provenance corrected.** An earlier worker report in this directory
  attributes its transcription to "Cherkaev2009-preprint arXiv:1009.3060". That
  arXiv id is a *different paper* — Cherkaev & Zhang, "Optimal anisotropic
  three-phase conducting composites: Plane problem" — whose abstract assumes
  "the conductivity of one of the materials is infinite". Our instances violate
  that, so this looked like the resolution. It is not: the correct paper was
  located, fetched and text-extracted here.
- **The transcription is correct.** Against the primary text, B1 (7.16),
  B2 (7.17), B3 (7.18), Z5, Z6, Z7 and m₁₁ (7.5) all match what we used.
- **The hypotheses are satisfied.** Theorem 7.1 [T]: "The effective conductivity
  k\* of a two-dimensional isotropic composite of three isotropic materials with
  conductivities k1<k2<k3 taken in the fractions m1, m2 and m3 … is bounded from
  below by the bound kL = B(m1,m2)". No further restriction. The k₃ = ∞ case is
  his separate Theorem 7.2.
- **The closed form is faithful to his own variational statement.** Implementing
  (7.1)–(7.2) directly, kL = max over t ∈ [k1,k2] of (−t + H₁(t)), and
  maximising over 400 001 values of t, reproduces the closed-form B2 at all
  three points, and at t = k1 reduces to HS_lo exactly.
- **Three internal consistency tests pass**: B2(m₁₁) = HS_lo(m₁₁) exactly in ℚ;
  B2(m₁₂) − B3(m₁₂) ≈ 1e−34; B2 ≥ HS_lo throughout its own region.
- **The field-ordering condition (4.2) does not exclude us.** The paper requires
  of optimal structures that field norms be ordered by material. Our below-m₁₁
  structure satisfies it cleanly: p1 ∈ [2.438891, 2.439933], p2 ∈ [1.334483,
  1.334483], p3 ∈ [0.669500, 0.673043], ordered and non-overlapping.

### 4. The adjudication pass

A dedicated adjudicator, briefed to assume our side was wrong, re-did the work
independently (`data/bound-adjudication.md`). Its verdict: **tentative genuine
discrepancy, explicitly not proven.** It ruled out, with its own tooling:

- *Transcription error.* Re-fetched and re-extracted the preprint with `pypdf`
  and re-transcribed Theorem 7.1 and eqs 7.5–7.24 from scratch; matches both
  earlier transcriptions. It flags two OCR hazards for future re-implementers
  (a case-ordering artefact in eq 7.15, and a stacked `1/H1(t)=` that
  `pdftotext` splits into a stray `1`), noting neither was actually made here.
- *Maple-simplification error.* Independently maximised the translation family
  by solving the first-order condition symbolically rather than trusting the
  closed form; matches B2 to 10 digits at all three points. It also caught a
  methodology trap worth recording: **naive grid maximisation over t finds a
  spurious larger maximum near t = k2** in a near-two-phase limit, from a sign
  flip in a denominator outside the formula's valid domain. My own grid scan
  (§3) was of that naive kind, so its agreement was partly luck; the
  critical-point solve is the trustworthy route.
- *A bug in our structure.* Re-derived 2D axis-aligned laminate homogenisation
  from the interface continuity conditions, sharing no code with either harness
  route, and reproduced our tensor and volume fractions **bit for bit**.
  Isotropy residual ~1e-66, far too small to explain a 1.4e−4 gap.
- *Scope.* Read §2.2 and §5: the underlying Theorem 5.1 bounds Tr[K\*]/2 for any
  composite satisfying (5.18), not only isotropic ones, so even the most
  charitable trace reading does not rescue B2. Condition (5.5) is satisfied at
  all three points with large margin. Phase ordering and finite k₃ both fine.

**What it could not settle,** and what is now the most likely home for an error
in either direction: B2's improvement over plain translation rests on a
localized-polyconvexity field constraint derived in Remark 4.1 by a
structural-variation argument (infinitesimal inclusion swaps, energy change
argued non-negative "if the tested configuration is optimal"). That is a
multi-page derivation it did not reproduce. A hidden restriction on
microstructure topology could live there, and our engineered hierarchical
laminate — whose phase-1 field is a deliberate mirror-image extreme pair — is
exactly the sort of object that might fall outside it invisibly.

It also flagged a real reproducibility gap: only the middle row had a saved
record. **Fixed here** — `explore/tp_below_m11.py` now rebuilds all three rows
and writes `data/attained/below-m11-m1-*.json`, each passing
`verify_laminate.py`.

### 5. The remaining lead

The bound improves on the plain translation bound by adding constraints the
paper repeatedly qualifies as holding **in optimal structures**. Section 4.1 is
titled "Boundedness of the Fields in Optimal Structures" and states [T] "The
fields in optimal microstructures satisfy certain additional local optimality
conditions that pointwise restrain the ranges ω_i of the fields in optimal
composites. These conditions, implemented into the polyconvex envelope
procedure, result in better bounds." If those conditions are assumed of the
minimiser rather than established without loss of generality, the resulting
function need not lower-bound every structure. Condition (4.2) is one such
condition and our structure satisfies it; the constraints in Sections 4.2, 4.4
and 5 have **not** been enumerated or tested, and that is where the
adjudication is now working.

## Outcome

`EVIDENCE` about our structures, and an **unresolved conflict with a published
theorem, reported as REQUIRING ESCALATION** rather than as a refutation. The
adjudication pass leans the same way — "tentative genuine discrepancy, not
proven" — and its recommendation, which the ledger follows, is that this be
recorded as *unresolved and under active adjudication* rather than as either a
win or a kill.

**Explicitly not claimed.** We do **not** claim Cherkaev's Theorem 7.1 is
wrong. The live hypotheses, in the order I would bet on them, are: (i) a
constraint in the derivation, not yet enumerated, that excludes structures like
ours — in which case the theorem is fine and its scope is narrower than its
statement reads; (ii) an error in the published bound; (iii) an error on our
side that five independent checks have not caught. Nothing about attainability
below m₁₁ is claimed in either direction, and 004's verified attainment *at*
m₁₁ is untouched by any of this.

**Scope of the computation:** the family of 004, at the (σ, m₁, m₂) listed, in
exact rational arithmetic. It is a search within one laminate family, so it
gives upper estimates of what laminates can do — which is the direction that
matters here, since it is a structure beating a bound.

## Why it failed / what survived

Nothing of ours failed; the conflict is the finding. What makes it hard to
dismiss is the accumulation: the discrepancy is exact, systematic across four
conductivity triples, robust in shape (ratio near 1.89), reproduced from the
primary source rather than a secondary quotation, and survives every hypothesis
tested so far — an anisotropy artefact, a transcription slip, a wrong paper, a
mis-simplified closed form, and the one field constraint we could identify.

Method note worth keeping: the provenance error was real and nearly resolved
the whole thing in the wrong direction. Two separate workers cited "Cherkaev
2009" while one of them was quoting a paper with an incompatible hypothesis.
Fetching the actual PDF and extracting its text cost little and changed the
conclusion. **Prefer the primary document over any number of agreeing secondary
quotations.**

Reusable: `tp_cherkaev_bound.py` (bound, region split, and the four consistency
tests, with the discrepancy pinned as an assertion so it cannot silently
vanish); the primary-source excerpts; the admissible-field certificate as a way
to verify an effective tensor independently of the lamination algebra.

## Leads generated

1. **Enumerate every constraint in Sections 4.2, 4.4 and 5** and test each
   against our structure's exact fields. This is the highest-value step and is
   what the adjudication was asked to do. A constraint we violate resolves
   everything and narrows the theorem's scope; if we satisfy all of them, the
   conflict sharpens considerably.
2. **A fourth, fully independent implementation of the effective tensor** —
   ideally not laminate algebra at all, e.g. a direct numerical solve on an
   explicit finite-scale realisation. Our discrepancy is ~5e−5 relative, so this
   needs care to be conclusive, but even three digits would be informative.
3. **Search for errata or later corrections** to the 2009 paper, and check
   whether any citing work reproduces the B2 branch independently.
4. **Find the true minimum below m₁₁** over a wider laminate class than 004's
   family. If a much lower structure exists, the conflict is larger and easier
   to adjudicate; if 004's family is near-optimal, the ratio near 1.89 may
   itself be a clue to the derivation's slack.
5. **Contact-the-author path**, recorded as an option not taken: this is the
   kind of discrepancy that a short note to the author would settle faster than
   any amount of re-derivation. Flagging for the human, not for an agent.

## References

- A. Cherkaev, preprint `math.utah.edu/~cherk/publ/newbounds7.pdf` (2008-05-06),
  of Mech. Mater. 41 (2009) 411–433. Theorem 7.1, eqs (7.1)–(7.18), §4.1. [T]
- A. Cherkaev and Y. Zhang, arXiv:1009.3060 — the *anisotropic, k₃ = ∞* paper
  that an earlier report conflated with the above. [T]
- This repo: `attempts/004`, `attempts/005`, `data/cherkaev-primary-excerpts.md`,
  `data/bound-adjudication.md`.
