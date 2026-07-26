# Union-Closed — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`. Record the mode
> in your attempt front-matter either way.

Machine-readable index: `prior-art.json`. Full records: `attempts/`.
Route-specific tooling: `explore/`.

## Editorial view of the attack surface

This is *this lab's* framing, not published consensus:

- Understand precisely why (3−√5)/2 is the barrier for the entropy method and
  what structural property of approximate-union-closed families breaks it;
  look for an added constraint (exact union-closure) the argument discards.
- Exhaustive/randomized search over small ground sets for families minimizing
  max element frequency; characterize extremal families.
- Lattice reformulation (union-closed families ↔ join-semilattices); test
  whether known equivalent forms (graph version, lattice version) admit
  sharper small-case analysis.

Assessed as the best "serious problem with genuine traction" candidate in the
problem set — active frontier with live technique.

## Attempts

### 001 — Entropy barrier map · `MAP`

Mapped the entropy-method frontier: the barrier statement made precise (any
*average-case* use of closure stops at ψ = (3−√5)/2; the Chase–Lovett
slice+top family is the universal blocker), a re-derivable writeup of the
sharp argument, and four candidate ideas (A–D) for injecting *exact*
closure. Established the n ≤ 4 exhaustive baseline (minimum max frequency is
exactly 1/2 over all union-closed families on ground sets of size ≤ 4).

**Key structural finding.** The entropy method's only use of closure is
H(A∪B) ≤ log|F| for iid uniform A, B — an average-case fact, tight at ψ.
Any advance must use worst-case closure of atypical/overlapping pairs,
dependent couplings, or counting structure.

**Dead end recorded:** k-wise unions strictly worsen the constant
(0.382 → 0.318 → 0.276 for k = 2, 3, 4).

### 002 — Weighted-KL ladder (idea B) · `REFUTED`

Formalized and closed idea B. The family-level version is vacuous for c ≤ 1
and false for c > 1; the distributional version is killed for **every** c ≥ 0
by Sawin's geometric-mixture family, with an exact finite-n certificate below
the 0.38271 record.

**Root obstruction.** KL charges escaping union-mass by log-likelihood
(log(1/δ) for planted mass δ) against a Θ(n) entropy drop. This generalizes
to a **no-go covering all smoothing-insensitive functionals** Φ(law(U), μ) —
read 002 before attempting *any* entropy-side strengthening.

### 003 — Dependent couplings (idea C) · `LIVE` (superseded by 004's corrections)

Formalized overlap-biased (Sinkhorn/Plackett-tilt) couplings and adversarially
tested them. The route separates every adversarial family that killed idea B
and the KL ladder — the first interface in this project not refuted by the
known counterexample genres. Verified that the functional genuinely evades
002's no-go (it is smoothing-*sensitive* by computation on 002's own
certificate instances).

**No bound is claimed.** Three labeled proof gaps remain: Plackett
odds-ratio control on non-product μ, the mutual-information tax at history
level, and recipe totality.

### 004 — Skeptic review of 003 · `VERIFIED with corrections`

Independent re-implementation (no shared code with 003) as an adversarial
verification pass, default stance refute. The interface and the headline
separations are **real**. Corrections applied:

- Mini-theorem (e) bracket positivity in 003 is **false as stated**; restated,
  and the `{0} ∪ [½,1)` genre it fails on is covered by a verified
  half-mixing conditionally-iid coupling.
- The single-parameter recipe ceiling is **0.431496** in closed form, not the
  ≈0.445 claimed in 003 (grid-floor artifact — the part-E hunt's grid excluded
  the extremal adversary). The corrected value sits *below* 0.44, contradicting
  003's claim that the result "does not affect the record-relevant regime".
- One non-reproducible data file, several overstated sentences.

**Use 004's statements, not 003's**, wherever they differ.

## Current standing

Route status `LIVE` — ceiling 0.4315 against the current published record
0.38271, with three labeled proof gaps as the open work. No new constant is
claimed by this lab.
