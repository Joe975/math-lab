# Closed-form analysis of two screened laminate shapes (rank 3 and rank 4)

sigma = (1, 2, 5), f = (f1, f2, f3) = (1, 1, 6)/8 unless stated otherwise.
All numbers below are produced and checked exactly (Fraction) by
`explore/tp_shapes.py --selftest` / `--report`; nothing here is a float claim.
Mode: informed (per team-lead brief). This note is a working record for the
requester, not a `problems/*/attempts/*` record.

## 0. Correction to the shape as first described

The rank-3 shape as verbally summarized to me used normal **e2** for the
outer lamination and **e1** for both L12 and L13. That shape is real (I
derived and checked it) but its unique isotropic point is `sigma* ~
3.454362`, nowhere near the screen's `3.3848002`.

Re-deriving the actual node order `tp_search.py` uses (`internal_nodes()` is
`[top, left, right]` in DFS-preorder, but the **printed** `normals` tuple is
indexed by that *same* list, not by the order the fractions are consumed)
shows the real winning shape has **outer normal e1**, **L12 normal e1**,
**L13 normal e2** — i.e. the outer split and L12 share a normal, L13 uses the
orthogonal one. With that correction the closed form matches the screen
exactly (`3.384800169...`, see selftest). Treat the shape below as the
corrected/verified one.

## 1. Rank-3 shape: closed form

```
L12 = laminate(p1 (frac m_L), p2; normal e1)
L13 = laminate(p1 (frac m_R), p3; normal e2)
result = laminate(L12 (frac m_top), L13; normal e1)
```

Elementary fact used twice: laminating two axis-aligned diagonal tensors
along normal e_k gives a diagonal result whose e_k-component is the
weighted **harmonic** mean of the two inputs' e_k-components, and whose
transverse component is the weighted **arithmetic** mean of the transverse
components (proved directly from `laminate()`'s lamination formula for
`n = e_k`, `harness/three-phase-conductivity/laminate.py`, and checked
numerically and via `L.laminate` against hand algebra during derivation).

With `h(m;p,q) = 1/(m/p + (1-m)/q)`, `a(m;p,q) = m p + (1-m) q`:

```
h12 = h(m_L; s1, s2)   a12 = a(m_L; s1, s2)     (L12, normal e1: T11=h12, T22=a12)
a13 = a(m_R; s1, s3)   h13 = h(m_R; s1, s3)     (L13, normal e2: T11=a13, T22=h13)

sigma11 = h(m_top; h12, a13)  =  h12*a13 / ((1-m_top) h12 + m_top a13)
sigma22 = a(m_top; a12, h13)  =  m_top a12 + (1-m_top) h13
```

Fraction constraints (2 independent, since f1+f2+f3=1 is automatic):

```
f2 = m_top (1 - m_L)     f3 = (1 - m_top)(1 - m_R)
=>  m_L = 1 - f2/m_top    m_R = 1 - f3/(1-m_top)
```

leaving **m_top as the single free parameter**, on the domain
`(f2, 1 - f3)` (both endpoints are where m_L or m_R hit 0 or 1).

### Isotropy polynomial

Clearing denominators in `sigma11(m_top) = sigma22(m_top)` gives a cubic
whose constant term vanishes — `m_top = 0` is an extraneous root introduced
by clearing `m_top` itself out of a denominator, not a physical solution
(`m_top=0` is outside the domain). Dividing it out leaves, at
`f = (1,1,6)/8`, `sigma = (1,2,5)`:

```
2917 x^2 - 3374 x + 520 = 0
```

**Degree 2**, not higher. Discriminant `= 3374^2 - 4*2917*520 = 5316516 =
2^2 * 3^4 * 61 * 269`; `61*269 = 16409` is squarefree and not a perfect
square, so **the root is genuinely irrational** and this quadratic is its
exact minimal polynomial over Q (irreducible: no rational root by the
rational root test over divisors of 520/2917, and a degree-2 primitive
polynomial with no rational root is automatically irreducible over Q).

```
m_top* = (1687 - 9*sqrt(16409)) / 2917  ~  0.18310642055899...
sigma*  ~  3.3848001699                  (gap to HS_lo = 37/11: 0.02116381)
```

Certified: `rank3_isotropic_point` bisects this exactly in Fraction
arithmetic to width `< 1e-60`; the corresponding tree (built with those
rational m_L, m_R, m_top) reproduces the tensor via `laminate.effective()`
bit-for-bit (`selftest`).

## 2. Why this beats the coated (Milton-type) construction

Propagating the imposed macroscopic field `E0 = (E, E)` (any direction works
by isotropy) down the tree using tangential-field / normal-flux continuity
at each interface gives the exact per-sublayer fields:

```
p1 in L12: (sigma11 E / s1,  E)
p2 in L12: (sigma11 E / s2,  E)
p1 in L13: (sigma11 E / a13, h13 E / s1)
p3 in L13: (sigma11 E / a13, h13 E / s3)
```

**Phase 2 and phase 3 each see a single, spatially uniform field** (they
appear in only one branch). **Phase 1 sees two distinct field values** —
`(sigma11 E/s1, E)` in the L12 branch and `(sigma11 E/a13, h13 E/s1)` in the
L13 branch — because it is split across both. Since `s1` is the smallest
conductivity and `f1` the smallest fraction here, phase 1 is exactly the
phase the classical Milton coated-assemblage construction needs to coat the
*other two* phases with a spatially uniform field in each grain to hit
HS_lo; that construction is feasible only when the coating-fraction
`rho` it requires stays in `[0,1]`, which fails once `HS_lo > sigma2`
(see `explore/tp_coated.py::iso_coat`, and `PROBLEM.md`'s citation of
Milton 1981 / Nesi 1995). The rank-3 shape above sidesteps that requirement
entirely: it never asks phase 1 to see a uniform field, only that the
*aggregate* (branch-averaged) response come out isotropic. That extra
freedom is exactly why it attains a value below the coated construction's
(here: infeasible) target, without needing the failed uniform-field
geometry.

## 3. Leading order as f1 -> 0 (f2:f3 fixed at 1:6)

```
f1        m_top*        sigma*        gap = sigma*-HS_lo   gap/f1
1/8       0.18310642    3.38480017    2.116381e-02          0.169310
1/16      0.15330065    3.80416417    3.820672e-02          0.611308
1/32      0.13908522    4.04688139    5.059886e-02          1.619164
1/64..    --            --            NO ISOTROPIC POINT (shape fails)
```

**This exact shape stops attaining isotropy below a finite threshold in
f1.** Bisecting for where the sign change in `sigma11-sigma22` across the
domain `(f2, 1-f3)` disappears gives

```
f1* in (0.030077570569119, 0.030077570583671)   (f2:f3 fixed at 1:6)
```

Below `f1*` the gap between `sigma11` and `sigma22` has one sign
throughout the whole domain — no member of *this* 3-parameter family is
isotropic at all, so "gap to HS_lo" is not even defined there for this
shape (a different rank-3 topology/normal choice would be needed). So the
"leading order in f1" the task asked for is not a clean power law over the
requested points `f1 = 1/8, 1/32, 1/128, 1/512`: the shape is only valid at
the first two of those (`1/8`, `1/32`; `1/16` also computed above), and
*fails outright* at `1/128` and `1/512`. Over the valid range the gap is
growing, not shrinking (`gap/f1` roughly 0.17 -> 0.61 -> 1.6), i.e. this
particular rank-3 shape is getting *worse* relative to HS_lo as f1 shrinks,
right up until it stops working. This is consistent with (not a
re-derivation of) Nesi's 1995 result that the true sharp bound pulls away
from HS_lo as the smallest phase's fraction shrinks — see part 4.

## 4. Rank-4 shape

```
X = laminate(p3 (frac b), p1; normal e2)
Z = laminate(p1 (frac d), p3; normal e1)
Y = laminate(Z (frac c), p2;  normal e2)
result = laminate(X (frac a), Y; normal e1)
```

In words: X is a plain rank-1 mixture of the best and worst phases (p3
dominant, small p1 fraction) laminated in e2; Z is the complementary rank-1
mixture of the same two phases (p1 dominant) laminated in the orthogonal
direction e1; Y layers that Z-mixture with the middle phase p2 in e2; the
whole thing is finally laminated with X in e1. This is **not** a member of
the "coat p3 by p1 in e2, then laminate with a p2-coated (p1|p3)" family in
the Milton sense — a true coating needs the *same* core laminated in **two**
orthogonal normals with the SAME core tensor (`tp_coated.py::iso_coat`); here
X and Z are two *different* rank-1 mixtures of p1/p3 (different fraction, b
vs d) built once each, not one core coated twice. It's better described as a
generic second-order (rank-2-of-rank-2) nested laminate that happens to put
p1 in both of its two top-level branches (via X directly, and via Z inside
Y), the rank-4 analogue of the rank-3 shape's "split the worst phase across
both branches" trick.

Closed form (same harmonic/arithmetic-of-diagonal-tensor rule as Section 1,
applied twice more):

```
T11_X = a(b; s3,s1)     T22_X = h(b; s3,s1)
T11_Y = c h(d; s1,s3) + (1-c) s2         T22_Y = 1/(c/a(d;s1,s3) + (1-c)/s2)
sigma11 = 1/(a/T11_X + (1-a)/T11_Y)      sigma22 = a T22_X + (1-a) T22_Y
```

Fraction constraints: `f2 = (1-a)(1-c)` fixes `c(a) = 1 - f2/(1-a)`; and
`f1 = a(1-b) + (1-a) c d` fixes `b(a,d) = 1 - (f1-(1-a) c d)/a`, leaving
**two** free parameters `(a, d)`; imposing isotropy removes one more,
leaving a genuine 1-parameter family to optimize over — exactly what
`rank4_optimize` does (golden-section over `a`, with `d` pinned at each `a`
by exact bisection to the isotropic point, all Fraction-valued).

Optimizing the exact closed form (not the float screen) gives

```
a* ~ 0.63864142   d* ~ 0.26597457   (b*~0.90271, c*~0.65408)
sigma* ~ 3.3641819563
gap to HS_lo (37/11 = 3.3636363636...):  5.455927e-04
```

versus the screen's rounded point (`a,b,c,d ~ 0.631443, 0.902048, 0.660840,
0.259276`), which the harness reproduces exactly as `sigma* ~
3.36422216 / 3.36422442` (tiny residual anisotropy from the rounded
decimals) — gap `~5.8e-4`. **The gap is real, not an optimization
artifact**: exact optimization over the shape's one remaining degree of
freedom only shaves the gap from `5.8e-4` to `5.46e-4`, it does not close
it. This is a small, bounded-rank census result (`EVIDENCE` about this rank
and topology only), and it is exactly the qualitative picture Nesi (1995)
predicts: `f1 = 1/8` is small enough that the sharp bound should sit
strictly above HS_lo.

## 5. SPECULATION: does rank -> infinity approach HS_lo, or a strictly larger limit?

The rank-3 -> rank-4 progression (`gap` 0.0212 -> 0.00055, roughly a
40x drop for one more rank at fixed f1=1/8) is consistent with either

- (a) laminates of this general "split the worst phase across branches"
  family converging to HS_lo as rank -> infinity (rate not yet estimated —
  two points don't establish a rate), with the residual gap at any finite
  rank being purely a finite-rank effect and not the Nesi obstruction, or
- (b) laminates converging to some rank-independent value strictly above
  HS_lo, with the Nesi bound (or an even tighter bound at this small f1) as
  the true finite floor, and the rank-3 -> rank-4 improvement just closing
  the gap between "unstructured guess" and that floor.

**SPECULATION, labelled:** given Section 3's finding that the *specific*
rank-3 topology used here has a hard f1 threshold below which it isn't even
isotropic, and Nesi's theorem is specifically about the smallest-phase-small
regime this problem lives in, (b) seems more likely than (a) at f1=1/8 --
but this is not evidenced by anything computed here beyond two data points.

**Falsifiable next test:** run the same "worst-phase-split-across-branches"
family at rank 5 and 6 (screen, then exact-optimize the winning topology the
way Section 4 does) at the same f=(1,1,6)/8. If the gap keeps shrinking by a
roughly constant *ratio* per added rank (geometric decay, e.g. ~40x per
rank continuing), that is evidence for (a) / convergence to HS_lo. If the
gap's decay slows and appears to be leveling off toward a nonzero floor (or
if, as in Section 3, the family becomes infeasible again at some rank before
converging), that is evidence for (b). A second, independent falsifiable
check: evaluate Nesi's or Cherkaev's published tighter bound (properly
transcribed per CONTRIBUTING rule 6, not yet done in this repo) at
f=(1,1,6)/8 and see whether `3.3641819563` (this rank-4 exact optimum) sits
above it, at it, or (impossible if the bound is correctly transcribed)
below it.

## Reproduce

```
python problems/three-phase-conductivity/explore/tp_shapes.py --selftest
python problems/three-phase-conductivity/explore/tp_shapes.py --report
```
