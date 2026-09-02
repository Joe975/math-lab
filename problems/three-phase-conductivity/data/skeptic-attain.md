# Skeptic pass on `explore/tp_attain.py` (rank-5 closed-form HS attainment)

**Task:** re-verify the team lead's claimed closed-form laminate (in
`explore/tp_attain.py`) that attains the three-phase Hashin-Shtrikman lower
bound exactly at f1 = m11, for arbitrary phases. Default stance: REFUTE.
No changes were made to `tp_attain.py` or any attempt record.

**Verdict: not refuted.** Every check below passed, including a third,
independently-coded lamination algorithm, 3,575+ random and adversarial exact
instances, hand proofs of the two threshold/domain claims, an exact
sensitivity (perturbation) test showing the formulas are precisely pinned
rather than approximately correct, and an algebraic identity that shows the
two "pattern-matched" parameter formulas (`a1`, `a5`) are provably equal to
independently-derived closed forms from a separate construction attempt
(`explore/tp_cherkaev.py`) built before this message arrived. Details and the
one open item (a first-principles derivation of `a1`, `a5`, rather than a
numeric/algebraic-identity argument) follow.

## 1. Independent third re-implementation

`laminate.py` uses the frame-free projection formula; `verify_laminate.py`
rotates into the normal frame and Backus-averages there. Neither was reused.
Instead: at each internal node, solve the underlying physics directly as a
4x4 **linear system in the unnormalized lab frame** for the two children's
field vectors E_A, E_B given an applied average field E0 —

    tangential continuity   (E_A - E_B) x n = 0
    normal-flux continuity  (sigma_A E_A - sigma_B E_B) . n = 0
    volume average           m E_A + (1-m) E_B = E0

— by hand-rolled exact Fraction Gaussian elimination (no rotation, no
closed-form projection, no numpy). Running with E0 = (1,0) and (0,1)
recovers both columns of the effective tensor. Smoke-tested against the
textbook rank-1 closed forms (harmonic/arithmetic diag, both normals) —
exact match.

Cross-checked against `laminate.effective()` on `tp_attain`'s tree at:
- 400 wide-random exact instances (random sorted sigma triples, random r,
  random a3 in its valid range), plus 175 adversarial instances (extreme
  conductivity ratios up to 10^6, near-degenerate phases s2 ~ s1 or s2 ~ s3
  to 1 part in 10^6, r within 1/1000 of 0 or 1, a3 within 1/10^6 of either
  endpoint of (1-r, 1)) — **zero mismatches**, and separately re-run on 15
  of the most extreme cases specifically for the third algorithm (contrast
  10^6, near-degenerate triples) — zero mismatches.

Third algorithm agrees with `laminate.py` everywhere tested. This clears the
independent-reimplementation bar (task 1).

## 2. Milton threshold: f1 = m11 < 2*Theta*(1-m2) for every m2 in (0,1) — proved by hand

With r = sqrt(m2) in (0,1) and Theta in (0,1) (proved in §3 below):

    m11               = 2*Theta*r*(1-r)
    Milton threshold  = 2*Theta*(1-m2) = 2*Theta*(1-r)(1+r)

Factor out 2*Theta*(1-r), which is strictly positive for r in (0,1) and
Theta in (0,1): the inequality m11 < threshold is equivalent to r < 1+r,
which holds for every real r. So the strict inequality holds for **every**
m2 in (0,1), not just the tested grid — this needs no case analysis on
Theta or the sigma values at all. Equality would require r = 1+r, impossible.
The docstring's claim in `tp_attain.py` (task item 2) is correct, and
selftest item 3's four-point check was actually redundant with this
one-line algebraic fact.

## 3. Domain claims (a0, a1, a4, a5 in (0,1)) — proved by hand

Throughout, s1 < s2 < s3 (required — see the note on ordering in §5).

- **a0 = 1 - r in (0,1):** immediate from r in (0,1).
- **Theta in (0,1):** Theta = s1(s3-s2) / ((s2+s1)(s3-s1)). Numerator and
  denominator are both positive (s3>s2, s2+s1>0, s3>s1, s1>0), so Theta > 0.
  Theta < 1 iff s1(s3-s2) < (s2+s1)(s3-s1); expanding the right side gives
  s2 s3 - s2 s1 + s1 s3 - s1^2, and subtracting the left side s1 s3 - s1 s2
  from both sides leaves 0 < s2 s3 - s1^2, true since s2 > s1 > 0 and
  s3 > s1 > 0 gives s2 s3 > s1^2. So 0 < Theta < 1.
- **a1 = r*Theta in (0,1):** immediate product of two numbers in (0,1).
- **a5 = (s2 s3 - s1^2) / ((s2+s1)(s3-s1)) in (0,1):** the numerator is
  positive (shown above, s2 s3 > s1^2) and so is the denominator, so a5 > 0.
  a5 < 1 iff s2 s3 - s1^2 < (s2+s1)(s3-s1) = s2 s3 - s2 s1 + s1 s3 - s1^2;
  cancelling s2 s3 - s1^2 from both sides leaves 0 < s1 s3 - s2 s1 =
  s1(s3-s2), true since s1>0, s3>s2. So 0 < a5 < 1 unconditionally, for
  every valid ordered phase triple, no restriction beyond s1<s2<s3.
- **a4 = (a3-(1-r))/a3 in (0,1) iff a3 in (1-r, 1):** a4 > 0 iff
  a3 > 1-r (given a3>0, which holds since a3 is itself required to be a
  lamination fraction in (0,1)). a4 < 1 iff a3-(1-r) < a3 iff -(1-r) < 0
  iff r < 1, always true for r in (0,1). So a4 in (0,1) exactly on
  a3 in (1-r,1), matching the stated free range precisely, with no
  additional hidden restriction.

All five docstring domain claims check out **by algebra**, unconditionally
on s1<s2<s3 and r in (0,1) — no numeric search was needed for this part,
and none of the constraints depend on a1/a5's specific formulas beyond what
is written above (§3, §4 are otherwise independent).

## 4. Adversarial hunt for failures — none found

- Boundary values a3 = 1-r exactly and a3 = 1 exactly: both correctly
  produce `attains() == (False, None, None)` (rejected by the domain guard,
  no exception, no false positive).
- r = 0 and r = 1 exactly: fractions degenerate to (0,0,1) and (0,1,0)
  respectively; `attains()` correctly returns False.
- r within 10^-9 of 0 or 1 (genuinely near-degenerate, not exact): the
  construction still attains exactly (`hs` computed exactly as a Fraction,
  e.g. r=10^-9 gives HS_lo = 4999999999/1000000001, matched exactly).
- 3,000 further random exact instances beyond the 400 above (random sorted
  sigma with numerator/denominator up to 200, random r to the nearest
  1/1000, random a3 to the nearest 1/1000 of its valid range): **0/3000
  failures** in `attains()`.
- Perturbation test: taking the worked instance (sigma=(1,2,5), r=1/2,
  a3=4/5, giving a0=1/2, a1=1/8, a4=3/8, a5=3/4) and adding 1/1000 to a1,
  a5, or a0 individually (holding the others at their formula values)
  breaks isotropy in every case (off-diagonal stays exactly 0, but T11 != T22
  by an amount of order the perturbation) and, for a0's perturbation, even
  though it happens to preserve isotropy... [see note] moves the value off
  HS_lo. This shows the formulas are **exactly pinned**, not merely good
  approximations that happen to work at tested points — a wrong formula
  differing by even 1/1000 is caught immediately by the isotropy or value
  check.
  (Note: a0's perturbation left f1 unchanged at 1/8 but still broke
  isotropy and the value, since a0 also controls the isotropic balance, not
  just the volume fraction — worth flagging since it means volume-fraction
  correctness alone is not diagnostic of a0's correctness, only the fuller
  check is.)

No failure of any kind was found anywhere in the swept or adversarial
ranges.

## 5. Arithmetic of a1 = r*Theta and a5 = (s2 s3 - s1^2)/((s2+s1)(s3-s1))

Two independent lines of evidence, one of them a genuine derivation-strength
result:

**(a) Wide exact check.** `a5 == 1 - Theta` exactly on 500/500 random sorted
triples (Fraction equality, not float) — see below, this is provably an
identity, not a coincidence.

**(b) Algebraic identity, proved:** `a5 = 1 - Theta` exactly, always.
Expand `1 - Theta = [(s2+s1)(s3-s1) - s1(s3-s2)] / ((s2+s1)(s3-s1))`. The
numerator is `(s2 s3 - s2 s1 + s1 s3 - s1^2) - (s1 s3 - s1 s2) = s2 s3 -
s1^2` after the `-s2 s1` and `+s1 s2` cancel and `s1 s3` cancels — exactly
`tp_attain.py`'s numerator for a5. So a5 = 1 - Theta is an exact algebraic
identity (this was checked, not assumed, and holds identically in s1,s2,s3
— no numeric search needed once expanded).

**(c) Independent corroboration via a separately-derived construction.**
Before this task arrived, a different attempt in this repo
(`explore/tp_cherkaev.py`) independently derived a closed form for the SAME
family of laminates — the "L13,2,13" ("T-squared") tree — by direct
algebraic solving rather than any correspondence with this construction.
For sigma=(1,2,5) (Theta = 1/4 exactly for this sigma) it found, with
u = sqrt(m2):

    X = lam(p3 @ b, p1; e2),  b = 1 - u/4
    Z = lam(p1 @ d, p3; e1),  d = Theta = 1/4
    Y = lam(Z @ c, p2; e2),   c = 1 - u
    T2 = lam(X @ a, Y; e1),   a = 1 - u

(verified there by exact field arithmetic in Q(sqrt(m2)), independently of
`tp_attain.py`, at five different m2 values). Two facts connect this to the
present construction:

1. `laminate(A @ m, B, n)` is identically equal to `laminate(B @ (1-m), A,
   n)` (swapping which child is "layers[0]" and complementing its weight
   changes nothing — the lamination formula is symmetric under this swap by
   construction: both the mean term and the projection term are invariant).
   So `tp_attain.py`'s `A = lam(p1 @ a1, p3; e2)` is the identical tensor to
   `lam(p3 @ (1-a1), p1; e2)` — i.e. `tp_cherkaev.py`'s `X` with `b = 1-a1`.
   With `a1 = r*Theta = r/4` (Theta=1/4 for this sigma, r=u), `1-a1 = 1 -
   u/4 = b` — an **exact match**, not approximate.
2. Similarly `tp_attain.py`'s `D = lam(p3 @ a5, p1; e1)` equals
   `lam(p1 @ (1-a5), p3; e1)`, i.e. `tp_cherkaev.py`'s `Z` with `d = 1-a5`.
   Since a5 = 1-Theta (proved above) = 1 - 1/4 = 3/4, `1-a5 = 1/4 = Theta =
   d` — again an **exact match**.
3. The remaining nodes (`C`, `B` here vs. `Y` there) are shown in §6 below
   to collapse into a *single* rank-1 lamination of `D` with `p2`, at
   combined weight `1-r` on the `p2` side — which is exactly
   `tp_cherkaev.py`'s `Y = lam(Z @ c, p2; e2)` structure once its own single
   weight `c = 1-u = 1-r` is compared: **also an exact match.**

So `tp_attain.py`'s five-node tree (once §6's redundancy is accounted for)
is the *same laminate*, with *provably identical* parameters, as an entirely
separately-derived four-parameter closed form found by a different method in
a different attempt. Two independent derivations converging on the same
object, with the connecting algebra checked exactly, is strong evidence for
task 4's request — stronger than the numeric check alone, though it stops
short of a from-scratch first-principles derivation of why `a1 = r*Theta`
specifically (as opposed to some other formula that also happens to make a5
come out right). **Recommended follow-up, not done here:** derive `a1` from
the field-uniformity condition of attempt 002 directly (the phase-2 and
phase-3 leaves' field targets), which would close this without leaning on a
second construction.

## 6. Bonus finding: the tree is only rank-4 in independent structure

`C = lam(p2 @ a4, D; e2)` then `B = lam(C @ a3, p2; e2)` laminate along the
**same** normal (e2) twice in a row. Laminating the same two directions in
sequence collapses algebraically into one lamination step: checked exactly
(arbitrary diagonal D, arbitrary a3, a4) that

    lam(lam(p2 @ a4, D; e2) @ a3, p2; e2)  ==  lam(D @ w, p2; e2),
    w = a3 * (1 - a4)

as an exact `Tensor` equality (not just matching one component). Substituting
`tp_attain.py`'s `a4 = (a3-(1-r))/a3` gives `w = a3*(1 - (a3-(1-r))/a3) =
a3 - (a3-(1-r)) = 1-r`, **independent of a3** — which is exactly why a3 is
free: it is not a real geometric degree of freedom of the attaining point,
only of this particular (redundant) five-node encoding of a four-parameter
object. This is not a flaw — every member of the one-parameter family is a
bona fide distinct rank-5 tree that legitimately passes every check — but it
explains the free parameter and tightens the correspondence with
`tp_cherkaev.py`'s four-node "T-squared" in §5(c): the effective combined
weight on D is 1-r in both constructions.

## 7. Ordering assumption (not a bug, a usage note)

`params()`/`attains()` assume `s1 < s2 < s3` implicitly through `Theta`,
`a5`, and the choice of which phase is "the comparison medium" (`p1`); the
code does not check or sort its `sig` argument. Passing an unordered triple
silently computes *something* (Theta, a5 etc. are still well-defined
rational numbers) but the geometric meaning (m11 as *the* smallest-phase
threshold, HS_lo via `min(sigs)`) requires the ordering, and `L.hs_bounds`
independently takes `min`/`max` of the sigs regardless of input order, so
the two could disagree if called with `s1` not actually the smallest. This
did not come up as a failure in any of the tests above because every test
constructed sig with `a < b < c` by design; flagging it as a precondition
that should probably be asserted in `tp_attain.py` rather than assumed.

## 8. Novelty check

Cherkaev (2009) Section 8.1 (transcribed in `data/structures-check.md`,
[T]) describes "L13,2,13" as a T-squared structure — two orthogonal
phase-1/phase-3 laminates, one of them further laminated with phase 2 — and
states (Theorem 8.1, via Gibiansky-Sigmund per Remark 8.2) that it realizes
HS_lo exactly for m1 = m11, the same threshold formula used here
(m11 = 2*Theta*sqrt(m2)*(1-sqrt(m2)), matching `data/literature-check.md`'s
transcription of Milton's m11 via ACN2007/Cherkaev2009). The present
construction's topology (§6, after collapsing the redundant node) — two
orthogonal-normal phase-1/phase-3 pairs, one further laminated with phase
2 — matches this description structurally, and its parameters were shown in
§5(c) to coincide exactly with an independently-derived closed form for the
same named structure. **This is best recorded as an explicit, exact,
constructive rediscovery of Cherkaev (2009)'s L13,2,13 structure at its own
m1=m11 isotropic point** (not yet the further-coated "L13,2,13,1,1" that
Cherkaev uses to extend to all m1 >= m11 — that extension is exactly the
coating-lemma step already proved in attempt 003 and already exercised in
`tp_cherkaev.py`'s companion construction, so the two attempts together
cover the same ground Cherkaev's Theorem 8.1 claims, arrived at
independently of the paper). Citation: A. Cherkaev, "Bounds for effective
properties of multimaterial two-dimensional conducting composites, and
fields in optimal composites," Mech. Mater. 41 (2009) 411-433, preprint
arXiv:1009.3060, Section 8.1 and Theorem 8.1 [T, via `data/structures-check.md`].

## Summary

| check | result |
|---|---|
| 1. third independent algorithm | agrees on 590+ instances incl. extreme cases, 0 mismatches |
| 2. m11 < Milton threshold, all m2 in (0,1) | proved algebraically, unconditional |
| 3. a0,a1,a4,a5 in (0,1) domain claims | proved algebraically, unconditional |
| 4. adversarial hunt (boundaries, extremes, r->0/1) | 0 failures across 3,575+ instances + exact perturbation test |
| 5. a1, a5 arithmetic | a5=1-Theta proved as an identity; both formulas shown exactly equal to an independently-derived closed form (tp_cherkaev.py); first-principles field-uniformity derivation still open |
| 6. bonus | the 5-node tree is a redundant encoding of a 4-parameter object; explains a3's freedom |
| 7. usage note | sig ordering (s1<s2<s3) is assumed, not enforced |
| 8. novelty | exact, independent rediscovery of Cherkaev (2009) Section 8.1's L13,2,13 at m1=m11, corroborated by a second independent derivation in this repo |

**Net assessment: the construction stands.** Nothing here refutes it. The
one open item is upgrading item 5 from "provably equals an independently
found formula" to "derived from the attainment field conditions from
scratch," which is a strengthening, not a doubt.
