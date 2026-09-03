# 017 — The below-m₁₁ curve is a ratio of curvatures at m₁₁, in closed form

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (executes 015 lead 1)
- **Type:** exact computation + closed-form derivation
- **Tools:** `explore/tp_curvature.py` (new; jet arithmetic over ℚ, the closed
  forms, `--selftest --table --m2-scan --verify --closed-form`),
  `explore/tp_ratrecon.py` (new; how the closed forms were reconstructed,
  `--derive`). Tests: `tests/test_three_phase_curvature.py`. Deterministic,
  standard library only; the whole module runs in seconds.
- **Sources:** this repo's 004, 005, 006, 009, 013, 014, 015;
  `explore/tp_cherkaev_bound.py` for the transcribed B2 [T].

## Approach

015 measured the position of the 006 family's minimum between the two bounds,

    ratio(m₁) = (family min − HS_lo)/(B2 − HS_lo)  →  c   as m₁ → m₁₁⁻,

got c ≈ 0.53–0.69 on five conductivity triples, and asked for a closed form. It
attacked that by refining a numerical minimisation at finitely many m₁ — the
obvious route, and the wrong one, because it treats a limit as something to be
approached rather than something to be evaluated.

The observation that changes the cost. **Both differences vanish at m₁₁**: 004's
structure attains HS_lo exactly there, and the transcribed B2 meets HS_lo there
(`tp_cherkaev_bound --selftest`, test A). And both vanish to *second* order —
the first derivatives agree too, checked exactly below. So c is not a limit to
be approached numerically at all; it is a ratio of second derivatives at one
point,

    c = γ/β,   β = ½(B2 − HS_lo)″(m₁₁),   γ = ½(family min − HS_lo)″(m₁₁),

and every ingredient is an exactly computable rational number.

γ is a **Schur complement**. Write g(m₁, a₀) = value − HS_lo for the family
(a₅ fixed by isotropy, a₁ by the phase-1 volume constraint, a₃ redundant by
005). Then g ≥ 0 everywhere — HS is a valid bound — and g = 0 at the interior
base point (m₁₁, 1−r), so the gradient vanishes there automatically and

    g = A₂₀u² + A₁₁uv + A₀₂v² + O(3),  u = m₁ − m₁₁, v = a₀ − (1−r),
    γ = min over v = A₂₀ − A₁₁²/(4A₀₂).

## What was done

### 1. Machinery

`tp_curvature.py` implements truncated Taylor ("jet") arithmetic over
`Fraction`, re-implements the lamination formula over that ring (the harness's
`laminate()` coerces with `Fraction()` and asserts `0 < m < 1`, so it cannot
take jets), and solves the isotropy condition for a₅ as a series by
chord-Newton. The constant term is cross-checked against `harness/laminate.py`
on every run, so the module doubles as a **second independent implementation**
of the family's effective tensor.

```
python tp_curvature.py --selftest      # 5 triples, all assertions exact in Q
```

The selftest pins, exactly in ℚ at each triple: the base point attains HS_lo;
the gap and **both** its partial derivatives vanish there; B2 − HS_lo and its
first derivative vanish there; the Hessian is PSD (as it must be); a₃ is
redundant (two different a₃ give identical jets); and 0 < c < 1.

### 2. c, exactly

| σ | m₁₁ | β | γ | **c exact** | 015's measured c |
|---|---|---|---|---|---|
| (1,2,5) | 1/8 | 320/27 | 448/71 | **189/355** = 0.53239437 | 0.532369 |
| (1,3,7) | 1/12 | 256/7 | 3520/167 | **385/668** = 0.57634731 | 0.576304 |
| (1,4,9) | 1/16 | 14080/171 | 44800/909 | **665/1111** = 0.59855986 | 0.598509 |
| (1,2,9) | 7/48 | 200448/9295 | 24088320/1618513 | **191675/277733** = 0.69014125 | 0.690125 |
| (2,5,11) | 2/21 | 101871/2560 | 349982451/16509200 | **1209312/2270015** = 0.53273304 | 0.532696 |

Reproduce with `python tp_curvature.py --table`.

**Correction to 015, recorded here rather than by editing it.** All five of
015's values are low, by 2.5e−5 to 4.6e−5. Its table of ratios *increases*
toward m₁₁ (0.532192 at h = 1e−3, 0.532374 at h = 1e−4) yet its quoted "limit"
is 0.532369 — below both. The values were read off at a finite m₁ rather than
extrapolated. The exact limit is 189/355 = 0.5323944, consistent with 015's own
sequence.

### 3. Closed forms

β, γ and c are rational functions of (x, y, r) = (σ₂/σ₁, σ₃/σ₁, √m₂) — the
problem is homogeneous in σ, so only the ratios enter. They were reconstructed
by exact rational-function interpolation in r at fixed σ, factoring the
denominators by exact deflation, then fitting the shape numbers against x and y
(`python tp_ratrecon.py --derive` re-runs the whole reconstruction). With

    p  = (1+x)/(y−x),          q  = (y−1)/(x−1),
    A  = p(y−1),               n₁ = (1+x)(y−1)/((x−1)(y+1)),
    K  = σ₁·p³(y−1)²(y+1)/4,
    C₂ = (1+x)³y(y−1) / (2(y−x)²(xy−1)),
    B₂ = (1+x)[2x(y³−1) − (x²+2x+3)y² + (3x²+2x+1)y] / (2(y−x)²(xy−1)),
    q₂(r) = C₂ + B₂r − r²,

the three quantities are

    β  =  K (r + n₁) / [(1−r)(r+p)² (r+pq)]
    γ  =  K r (A − r) / [(1−r)(r+p)² q₂(r)]
    c  =  r (A − r)(r + pq) / [(r + n₁) q₂(r)]

with the K and the (1−r)(r+p)² cancelling between the two curvatures. At
σ = (1,2,5), m₂ = 1/4: p = 1, q = 4, A = 4, n₁ = 2, q₂ = 10/3 + 17r/3 − r², and
c = 3r(16−r²)/[(2+r)(10+17r−3r²)] = 189/355 at r = 1/2.

**Verified as identities, not fitted:** `python tp_curvature.py --closed-form`
checks the closed forms against the jet expansion at **2610** (σ, r) instances —
199 triples including non-integer σ and σ₁ ≠ 1, ten values of r — with zero
mismatches.

### 4. Two facts the closed form exposes

**(a) c depends on m₂, not on σ alone.** 015 wrote c(σ) after sampling only
m₂ = 1/4. At σ = (1,2,5) the exact values are

| r = √m₂ | 1/4 | 1/3 | 1/2 | 2/3 | 3/4 |
|---|---|---|---|---|---|
| c | 17/45 = 0.3778 | 143/322 = 0.4441 | 189/355 = 0.5324 | 7/12 = 0.5833 | 2223/3707 = 0.5997 |

a swing of 0.22 — larger than the whole σ-spread 015 reported. So the constant
is c(σ, m₂), and 015's "σ-dependent limit" understates what it varies with.
(`python tp_curvature.py --m2-scan`.)

**(b) 0 < c < 1 at every admissible instance tested, and c → 1 only at infinite
contrast.** A closed-form sweep over 76 038 admissible (σ, m₂) points finds no
exception, with c ranging over 0.0377 to 0.9538. Pushing the parameters:
c → 1 as y = σ₃/σ₁ → ∞ (0.9791 at y = 100, 0.99979 at y = 10⁴, 0.9999979 at
y = 10⁶, always strictly below 1), and c → 0 as r → 0.

That limit is the striking part. σ₃/σ₁ → ∞ is exactly the case Cherkaev treats
*separately* — Theorem 7.2 is the k₃ = ∞ statement, and the Cherkaev–Zhang paper
(arXiv:1009.3060) assumes one phase has infinite conductivity. So B2's curvature
at m₁₁ is asymptotically correct in the infinite-contrast regime his optimal
structures were built for, and increasingly over-optimistic as contrast falls.

### 5. Verification

- **Independent implementation.** The jet ring's algebra is a rewrite; its
  constant term is compared to `harness/laminate.py` (both `effective` and the
  Keller–Dykhne identity) at every base point in the selftest.
- **Independent numeric route.** `python tp_curvature.py --verify` ignores the
  expansion entirely: it minimises the family over a₀ directly at
  h = 10⁻², 10⁻³, 10⁻⁴, 10⁻⁵ with a₅ from bisection, and reports gap/h². At
  σ = (1,2,5): 6.6501, 6.3429, 6.3132, 6.3102 against γ = 6.3098591549, with
  the measured ratio 0.53239214 against c = 0.53239437 — O(h) convergence to
  the predicted values on all five triples. The optimal a₀ offset also tracks
  the predicted slope −A₁₁/(2A₀₂) to three digits.
- **One fully exact point.** At h = 10⁻³ on the predicted optimal path, the
  structure is rebuilt in exact ℚ, handed to `harness/laminate.py`, and its
  volume fractions and anisotropy (< 10⁻⁵⁵) checked there before the gap is
  taken.
- **Out-of-sample.** The closed forms were reconstructed from six triples at
  m₂ = 1/4-ish samples and then tested on 2610 instances they were not fitted
  to, including σ₁ ≠ 1 and non-integer σ.

## Outcome

`VERIFIED` for the exact values, scoped to: the 006 family (the rank-4
axis-normal topology of 004/005), the second-order coefficient at m₁ = m₁₁, and
the transcribed B2 [T].

- c = γ/β is an exactly rational number at every rational (σ, m₂), given in
  closed form above and verified as an identity at 2610 instances.
- The five values 015 reported are each 2.5e−5 to 4.6e−5 too low; the exact
  ones are tabulated in §2.
- c varies with m₂ as much as with σ.

`EVIDENCE`, scoped to the 76 038-point sweep: 0 < c < 1 at every admissible
instance, with the supremum approached only as σ₃/σ₁ → ∞.

**Not claimed.**

- That γ is the true sharp curvature over *all* microstructures. It is the
  family's. 018 tests one large enrichment; higher rank and other topologies
  remain open.
- That c < 1 is *proved*. It is a polynomial inequality in (x, y, r) that the
  closed form now makes attackable, but 76 038 exact points are `EVIDENCE`, not
  a proof.
- Nothing about m₁ far below m₁₁. The whole record is a statement about the
  second-order behaviour at m₁₁; 015's observation that the ratio stays near c
  down to m₁₁ − 0.027 is a separate, still-empirical fact.

## Why it failed / what survived

Nothing failed. What this replaces is a *method*: three records (013, 014, 015)
spent more than 23 000 structure evaluations establishing that the family sits
below B2 and measuring by how much. The same content, sharper and exact, is one
Hessian at one point — and it generalises to every σ and m₂ at once instead of
to the handful a search can afford.

The method lesson is specific enough to reuse. **When a bound is met exactly at
a boundary point of the region where it is in question, the disagreement just
inside that region is a derivative comparison at the point, not a search over
the region.** Here the entire below-m₁₁ escalation, in the limit, collapses to
one inequality between two rational numbers: γ < β. That is why the exact
values matter — 015's numbers were close enough to be right about the sign and
still wrong in the fourth decimal, which is exactly the resolution at which a
curvature comparison lives.

The second lesson is about tangency order. It was not obvious in advance that
B2 meets HS_lo to *first* order at m₁₁ and not merely continuously. If the
meeting had been transversal, the ratio would tend to 0 and 015's whole curve
would have been an artefact of finite h. The order had to be checked before any
of this was worth doing, and checking it took one exact computation.

Reusable: `tp_curvature.py` (jet arithmetic over ℚ, and a second independent
implementation of the family's tensor); the closed forms; `tp_ratrecon.py` (the
reconstruction pipeline, applicable to any exactly-sampled rational function in
this problem).

## Leads generated

1. **Prove 0 < c < 1.** With the closed form it is a polynomial inequality in
   (x, y, r) over the admissible region — the numerator/denominator difference
   is an explicit polynomial whose sign is the whole claim. A `FORMALIZED`
   candidate if it holds, and a very sharp refutation if it does not.
2. **Ask what c → 1 as σ₃/σ₁ → ∞ means for the escalation.** If B2's curvature
   is exactly right at infinite contrast, the natural conjecture is that B2 is
   the m₁₁-curvature of the *true* bound only in that limit, and that the
   correct finite-contrast bound has curvature γ. Testing that needs a second
   family that also attains at m₁₁ — see 018.
3. **Do the same expansion at m₁₂**, the other breakpoint of the transcribed
   piecewise bound, where B2 meets B3. Same machinery, different point;
   m₁₂ is irrational so the arithmetic moves to ℚ(√Z₂).
4. **Expand to third order.** The jet module already carries degree 3; the u³
   coefficient would say whether the family's curve is convex or concave
   relative to B2 just below m₁₁ and how far the quadratic model reaches
   (015's data holds to about m₁₁ − 0.027).
5. Carried, unchanged from 015: re-derive Remark 4.6, and the note to the
   author. Now with a much sharper question to ask — is the curvature at m₁₁
   supposed to be β?

## References

- This repo: `attempts/004`, `attempts/005`, `attempts/006`, `attempts/009`,
  `attempts/013`, `attempts/014`, `attempts/015`; `explore/tp_curvature.py`,
  `explore/tp_ratrecon.py`, `explore/tp_cherkaev_bound.py`.
- A. Cherkaev, *Bounds for effective properties of multimaterial
  two-dimensional conducting composites, and fields in optimal composites*,
  Mech. Mater. 41 (2009) 411–433, Thm 7.1 and 7.2 [T] — see
  `data/cherkaev-primary-excerpts.md`.
