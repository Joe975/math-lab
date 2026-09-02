# Structures check: is HS_lo attained on m1 in [m11, 2*Theta*(1-m2))?

Mode: informed side, literature-only task (no attempt records read or written).
Method: downloaded and full-text-extracted (pypdf) the three primary PDFs below;
quoted passages are machine transcriptions of the extracted text, marked [T] with
page/section. One additional very recent paper was located and checked but found
off-target (see "Sources not usable" below).

Sources read in full text:
- N. Albin, A. Cherkaev, V. Nesi, "Multiphase laminates of extremal effective
  conductivity in two dimensions," J. Mech. Phys. Solids 55 (2007) 1513-1553.
  Free PDF: https://www.math.k-state.edu/~albin/pubs/albin_cherkaev_nesi_2007_mle.pdf
  (cited **ACN2007**; page numbers below are the paper's own printed page numbers,
  e.g. "p.19" = the page footer reading "MULTIPHASE LAMINATES... 19")
- A. Cherkaev, "Bounds for effective properties of multimaterial two-dimensional
  conducting composites, and fields in optimal composites," Mechanics of
  Materials 41 (2009) 411-433. Preprint PDF:
  http://www.math.utah.edu/~cherk/publ/newbounds7.pdf (cited **Cherkaev2009**;
  page numbers below are the preprint's own internal page numbers, e.g. "p.42")
- A. Cherkaev and Y. Zhang, "Optimal anisotropic three-phase conducting
  composites: Plane problem," arXiv:1009.3060v3 [math-ph], 22 May 2011,
  published Int. J. Solids Struct. 48(20) (2011) 2800-2813. Fetched directly
  from arXiv: https://arxiv.org/pdf/1009.3060 (cited **CZ2011**)

Sources located but not usable for the numbered questions:
- Milton (1981), Appl. Phys. A 26, 125-130 — paywalled, no free copy found.
  Everything about it below is via ACN2007's and Cherkaev2009's quotations of
  it, as already flagged in the prior `literature-check.md`.
- Nesi (1995), Proc. Roy. Soc. Edinburgh A 125, 1219-1239 — paywalled, no free
  copy found (same as prior report).
- L. G. Fel, "Isotropic conductivity of two-dimensional three- and four-phase
  symmetric composites: duality and universal bounds," arXiv:2512.20401v1
  [cond-mat.dis-nn], 23 Dec 2025. Fetched and read
  (https://arxiv.org/html/2512.20401). Does transcribe an explicit closed-form
  three-phase Nesi bound (its eq. 3.24, a pair of cubics in the elementary
  symmetric functions of σ1,σ2,σ3) — but **only for the cyclically symmetric
  case**, i.e. m1=m2=m3=1/3. [T, Fel2025, Sec. "Nesi's bounds (three-phase
  composite)"]: "new bounds σ_Ne+ and σ_Ne− were derived for σe(σ1,σ2,σ3) of
  the isotropic 2D three-phase composite made of **cyclically symmetric**
  isotropic constituents." Not usable for task 3's requested points
  f=(1/8,1/8,3/4) and f=(1/4,1/4,1/2), which are not symmetric. Not chased
  further.

---

## 1. ACN2007's general laminate family and named structures — VERDICT: found

ACN2007 Section 4 builds up a hierarchy of laminates, all special cases of a
6-parameter "orthogonal laminate of high rank" (their Figure 6, p.15):

- **T-structure**, `L(13,2)`: laminate K1+K3 (normal n1), then laminate that
  with K2 (normal n2, orthogonal). [T, p.18-19, Fig. 8 caption region]
- **Milton-Kohn / coated-T structures**, region `L(KT)`: laminate the
  T-structure with more K1 ("coating" it). [T, p.19] "This construction
  proves the optimality of the bound (12) in a region of anisotropic points
  for a smaller value of m1 than was previously known possible: m1 ≥
  Θ(1 − m2) rather than m1 ≥ 2Θ(1 − m2)." Note explicitly [T, p.19]: "the
  only optimal **isotropic** structures found in L(KT) were already known to
  be optimal" — i.e. the coated-T family alone reaches a new isotropic point
  only at Milton's own threshold m1 = 2Θ(1−m2); it doesn't extend isotropic
  attainability below it.
- **T²-structure**, `L(13,2,13)` composed further: laminate the T-structure
  with a *second*, independently-parametrized K1+K3 laminate (orthogonal
  direction). Defined by eq. (26), p.20; parametrized by ω1, ω2 with
  ω1+ω2 = (1/Θ)(m1+2Θm2), ω1ω2 = m2 [T, eq. (28), p.22].

**Direct answer to "does any of these attain HS_lo on [m11, 2Θ(1−m2))":**
Yes, explicitly claimed, twice, in ACN2007 itself:

> [T, ACN2007 Theorem 2, p.14, restating Gibiansky & Sigmund 2000] "Let
> m=(m1,m2,m3) ... such that 2Θ(√m2 − m2) ≤ m1 ≤ 2Θ(1 − m2) ... Then there
> exists an isotropic structure with the given volume fractions and optimal
> for the bound (12) [the translation/HS lower bound]."

2Θ(√m2−m2) = 2Θ√m2(1−√m2) is *exactly* Cherkaev's m11 (matches the brief's
formula and the previously-reported numeric value m11=0.114277 at m2=1/8).
So Theorem 2, taken at face value, already covers the whole team interval
[m11, 2Θ(1−m2)) and beyond, up to m1=2Θ(1−m2) itself.

ACN2007 immediately complicates this, however — right after stating Theorem 2:

> [T, ACN2007 Remark 5, p.14] "The results of Gibiansky and Sigmund raise an
> interesting question. If the volume fractions satisfy the inequalities
> 2Θ(√m2−m2) ≤ m1 < 2Θ(1−m2), is the isotropic point on the translation bound
> attainable? In particular, is it possible to find an explicit formula for
> the partition suggested by the computer-generated structure...?"

So ACN2007 itself flags Theorem 2 (as stated by Gibiansky-Sigmund, discovered
by numerical topology-optimization, not a closed-form proof) as **not yet
rigorously constructive** on exactly the team's interval, and poses closing
that gap as an open question at the time of writing (2007).

ACN2007's own answer attempt is **Theorem 6** (p.24), which explicitly claims
to close this:

> [T, ACN2007 Theorem 6(ii), p.24] "If 2Θ(√m2−m2) ≤ m1 ≤ Θ(1−m2), then (12)
> is optimal. There exists a set of optimal points on the bound which
> includes the isotropic point and whose most anisotropic member is that
> given by the optimal T²-structure..."
> [T, ACN2007 Theorem 6(i), p.23] "If m1 > Θ(1−m2), then (12) is optimal.
> There exists a set of optimal points on the bound which includes the
> isotropic point and whose most anisotropic member is [the] coated
> T-structure..."
> [T, ACN2007 p.25, remark right after Theorem 6] "the corresponding figure
> for Theorem 2 coincides with the right side of Figure 11. This is because
> Theorem 2 and Theorem 6 are identical for isotropic structures."

**Caution — I could not independently verify Theorem 6(ii) algebraically at
the team's own numeric point.** Plugging m1=m2=1/8, k=(1,2,5), Θ=0.25 into
the T²-structure's own defining equations [T, eq. (28), p.22]
(ω1+ω2 = (1/Θ)(m1+2Θm2) = 0.75, ω1ω2 = m2 = 0.125) gives ω1=0.5, ω2=0.25 —
**not equal**, and ACN2007 states elsewhere [T, p.22] that this T²-structure
is isotropic *only* when ω1=ω2, which by their own algebra only occurs
exactly at m1 = 2Θ(√m2−m2) = m11. So the single T²-structure obtained by
directly substituting the team's (m1,m2) is anisotropic, not isotropic — a
literal reading of eq. (28) alone does not hand you the isotropic point at
m1=1/8 for free. Theorem 6's "set of optimal points ... includes the
isotropic point" evidently refers to some larger family (varying additional
freedom in ν, ν' beyond the "most anisotropic member" values ACN2007 gives
explicitly), which the excerpted text does not spell out in a form I could
directly re-derive. **Flag this loudly**: Theorem 6(ii) is stated in ACN2007
as an established theorem, but I was not able to independently reconstruct
the isotropic member of its structure family from the formulas given, and
recommend the team's own exact-rational harness attempt this construction
directly rather than trust my re-derivation or the bare theorem statement.

## 2. Cherkaev (2009) §7-8 and Cherkaev-Zhang (2011) — VERDICT: found, resolves the above

Cherkaev (2009) is a *later* paper than ACN2007 and gives what reads as the
actual closing of the gap that ACN2007's Remark 5 left open, via a
self-contained isotropic-by-construction structure (not the plain T²-structure):

> [T, Cherkaev2009 p.42, Section 8.1] "Finally, the T²-structures are
> sequentially laminated by the two orthogonal layers of k1, forming the
> structure L13,2,13,1,1... The above-listed conditions for the fields in
> laminates form a system of equations for the unknown volume fractions of
> layers. If the system has a solution, the optimal structure is found. The
> solvability conditions restrict the range of volume fraction m1 as
> m1 ≥ m11, see [4] [=ACN2007]. The described structure realizes
> Hashin-Shtrikman bound because sufficient conditions (3.27), (3.28), and
> (4.26) are satisfied everywhere."
>
> [T, Cherkaev2009 p.42, Section 8.1, describing the B1-structures shown in
> Fig. 7 left] "They consist of inclusions sequentially laminated by two
> orthogonal layers of the amount m1 − m11 of k1. The inclusions are
> T²-structures L13,2,13 in which the amount of the k1-material is equal to
> m11."

This construction is self-consistent with my algebra above: taking the
T²-structure *at exactly* m1=m11 (where, per ACN2007's own condition ω1=ω2,
it genuinely is isotropic), then adding the *extra* k1 material (fraction
m1 − m11) as two further equal-fraction orthogonal laminations, preserves
isotropy (laminating an isotropic tensor with equal orthogonal fractions of
an isotropic material stays isotropic) while raising the effective
conductivity from the m11 value up toward k2. This resolves the tension I
flagged in §1 without needing Theorem 6(ii)'s T²-structure formula directly.

> [T, Cherkaev2009 Theorem 8.1, p.43] "The bound (7.19)-(7.24) is exact in
> each point: There exist laminates of a finite rank that realize the
> bounds."

And independently, citing Gibiansky-Sigmund's *own* proof (not just the
numerical discovery ACN2007 called into question):

> [T, Cherkaev2009 Remark 8.2, p.42] "Gibiansky and Sigmund [15] proved
> optimality of this construction: It realizes Hashin-Shtrikman bound in the
> interval m1 ∈ [m11, 1]."

(Caveat on this last quote: ref [15] in Cherkaev2009's own bibliography is
Gibiansky & Sigmund 2000, JMPS 48(3):461-498, "Multiphase composites with
extremal **bulk modulus**" — an *elasticity* paper, not conductivity. Remark
8.2 is Cherkaev's own claim that the analogous construction, reinterpreted
for conductivity, attains HS_lo on all of [m11,1]; it is not a literal
restatement of a conductivity theorem inside ref [15] itself. ACN2007 itself
makes the same "reinterpreting their results" move for the conductivity
problem — see ACN2007 footnote 3, p.14 [T]: "The paper focused mainly on the
problem of bulk moduli, but the results easily apply to the conductivity
problem as they describe in section 5.3.")

**Region-B1 solvability condition m1 ≥ m11 is exactly Cherkaev's own
attainability threshold for the piece of his bound B(m1,m2) that equals HS**
(the previously-reported region split: B1=HS for m1≥m11). Note this is a
**different bound from "Nesi's bound"** referenced elsewhere in the same
paper — the prior report's quoted Remark 4.4 ("Nesi bound is generally not
achievable by a structure") is about Nesi's *original* 1995 formula, a
distinct (and, per Cherkaev, less refined) bound that Cherkaev's own
Theorem 7.1 improves on. Theorem 8.1's attainability claim is specifically
about Cherkaev's own bound B(m1,m2), not Nesi's — there is no contradiction
between "Nesi bound generally not achievable" and "Cherkaev's B(m1,m2), which
equals HS for m1≥m11, is exact everywhere." This should resolve the ambiguity
the prior worker flagged.

**At the team's numeric point** (k=(1,2,5), m2=1/8, m1=1/8=0.125):
m11=0.114277 < 0.125, so m1 ≥ m11 — squarely inside the region where
Cherkaev's B1-structure (Section 8.1, above) is claimed to realize HS_lo
exactly, per Theorem 8.1.

### Cherkaev-Zhang (2011) — gap region is for k3=∞ only, and doesn't touch the isotropic line

CZ2011 treats a *different, more restrictive* problem: it assumes k3=∞ (a
superconductor phase) and works out the *anisotropic* G-closure boundary in
an (r, m1) parameter plane (r = degree of anisotropy of the imposed field;
r=1 is the isotropic case). It reports genuine, quantified gaps:

> [T, CZ2011 abstract] "The found structures match the bounds in all but one
> region of parameters; we discuss the reason for the gap and numerically
> estimate it."
> [T, CZ2011 §6.4, p.21-22] "The bound (74) in region E is not attainable...
> We have not found optimal structures in the complementary region D2... and
> conjecture that the bound probably is not exact there... The part of region
> D (region D2) close to E is probably not attainable as well."
> [T, CZ2011 §6.4, p.24, numerical result] "the relative differences are
> rather small, [on] the order of 10⁻⁴ and even ... 10⁻⁷" between the guessed
> structure's energy and the unattained bound BE, in region E.

But CZ2011 explicitly says the isotropic line (r=1) is **not** part of this
gap:

> [T, CZ2011 §4, p.11] "The upper interval r = 1 correspond[s] to isotropy of
> K∗, the isotropic bound have been derived in (Nesi 1995), and optimal
> structures have been found in (Cherkaev, 2009). Other bounds and structures
> are new."

So per CZ2011's own framing, the isotropic case is claimed *solved* (bound
and structure meet) via Cherkaev (2009) — consistent with §2's finding
above — and the paper's actual unresolved gap (regions D2, E) is confined to
genuinely anisotropic effective tensors, and only in the k3=∞ special case,
which does not directly apply to the repo's finite-k3=(1,2,5) problem.
**Caveat**: I did not locate a finite-k3 analogue of CZ2011's gap analysis;
if one exists in Cherkaev2009 §6 (anisotropic bounds, general finite k3, only
skimmed here) or elsewhere, it was not found. Flag as **not fully checked**
for finite k3, not "confirmed absent."

## 3. Nesi (1995) explicit three-phase formula — VERDICT: not found (for general fractions)

Confirmed the prior report's "not found" for the primary source. The one
explicit closed-form three-phase Nesi bound I located (Fel 2025, arXiv:2512.20401,
eq. 3.24 — a pair of cubics in I1=Σσi, I2=Σσiσj, I3=σ1σ2σ3) is **restricted to
the cyclically symmetric case m1=m2=m3=1/3** [T, quoted in full under "Sources
not usable" above], so it cannot be evaluated at the team's requested
f=(1/8,1/8,3/4) or f=(1/4,1/4,1/2) — those aren't symmetric fractions. No
numeric comparison against HS_lo or the rank-4 value was possible. This
remains a genuine gap, not a resolved "absent."

## 4. Does the literature make an explicit attainability claim for [m11, 2Θ(1−m2))? — VERDICT: yes, attained (claimed)

Three independent-but-connected published claims, all pointing the same way:

1. ACN2007 Theorem 2 (restating Gibiansky-Sigmund 2000): isotropic HS-optimal
   structure exists on 2Θ(√m2−m2)=m11 ≤ m1 ≤ 2Θ(1−m2) — but ACN2007's own
   Remark 5 flags this as not yet rigorously constructive over that whole
   range at the time of writing (2007).
2. ACN2007 Theorem 6: claims to supply that missing construction (coated-T
   structures for m1≥Θ(1−m2); T²-structures for m11≤m1≤Θ(1−m2)) — but I could
   not independently reconstruct the isotropic member of the T²-structure
   family from the given formulas at the team's numeric point (§1, flagged
   loudly).
3. Cherkaev (2009) Theorem 8.1 + Section 8.1: a *different*, self-consistent,
   isotropic-by-construction laminate (L13,2,13,1,1 — a T²-structure fixed at
   m1=m11, coated with orthogonal layers of the extra k1) that provably stays
   isotropic and realizes HS_lo for all m1 ≥ m11, i.e. all of [m11,1] ⊃
   [m11, 2Θ(1−m2)). This is the cleanest and most self-consistent of the
   three, and it is the one I'd trust first.

**So: the literature's own position is "attained," not "left open" or "known
not attainable."** No source found states or implies non-attainability on
this interval; CZ2011's real, quantified gap is a different problem (k3=∞,
genuinely anisotropic tensors) that does not reach the isotropic line. The
one substantive uncertainty is that I could not personally re-derive
Cherkaev's B1-structure's exact laminate tree (volume fractions, directions)
from the formulas transcribed above closely enough to hand the team a
ready-to-verify exact-rational construction — Section 8.1 and the Appendix
(Section 9, "Calculation of parameters of optimal laminates," Cherkaev2009
p.45-52) contain more detail (explicit field formulas e11 etc., eq. 9.1-9.5)
that I did not fully work through into a lamination-tree the harness could
consume; recommend a follow-up pass specifically on Cherkaev2009's Appendix
if the team wants to build this structure into `laminate.py` and check it
against HS_lo exactly, rather than relying on the rank-4 5.46e-4 approximation.
