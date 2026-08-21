# 048 — Branch L1 in scale-free form: the infimum is exactly 1, approached only logarithmically, and c\* is essential

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** proof probe (047 lead 1).
- **Tools:** `explore/uc_hu_L1.py` (new; the ratio formulation, the
  boundary identity, the binding arm of t, the c\*-free test, and the
  fixed-q minimiser map; deterministic, seeds 4901–4904; checkpoint
  `data/hu_L1.json`); `explore/uc_hu_L1_skeptic.py` (new; terms
  rebuilt in nats from this record's prose, the ratio floor attacked
  by a descent that actively minimises it, the boundary identity and
  the two-regime classification re-derived independently; exit 0).
  Reproduce: run the two in that order.
- **Sources:** none.

## Approach

047 isolated branch **L1**: where the first coordinate carries the
maximum marginal (q = f₀ = 1−x, so the slack A vanishes identically),
the pair interaction alone must cover the conditional deficit,

    C := t·Δ ≥ −B,   t = min(x−1/2, 1−x),
    Δ = 2ψ(p₀+p₁) − ψ(2p₀) − ψ(2p₁),   ψ(s) = h(min(1/2,s)),
    −B = x·h(p₀)[c*(q) − c*(p₀)] + (1−x)·h(p₁)[c*(q) − c*(p₁)],

subject to f₁ = x p₀ + (1−x) p₁ ≤ q.

047 advised expanding at x → 1, on the grounds that the margin C + B
gets small there. That advice was about the *absolute* margin, which
vanishes at x → 1 only because both terms do — the same absolute-vs-
scale-free trap 046 §G hit at the q → 1/2 boundary. This record works
in the ratio C/(−B) throughout.

## What was done

**P1. The ratio has no tight interior point.** Over 600,000 branch
samples (145,103 with a deficit) the minimum ratio is **1.1492**, and
the skeptic's descent — which actively minimises the ratio rather than
sampling it — bottoms out at **1.0881**, at q = 1e-7, p₀ = 0,
p₁ = 1/2. Nothing in the interior approaches 1.

**P2. The boundary identity.** At the configuration (p₀,p₁) = (0,1/2)
the ratio is **exactly 1/c\*(q)** — verified to 2.2e-16 independently:

    q       ratio      1/c*(q)
    1e-1    1.85423    1.85423
    1e-2    1.33218    1.33218
    1e-3    1.21278    1.21278
    1e-4    1.15712    1.15712
    1e-5    1.12459    1.12459
    1e-6    1.10323    1.10323

  Since c\*(q) < 1 for every q ∈ (0,1/2), this configuration always
  satisfies L1 strictly; and since c\*(q) → 1 as q → 0, the branch
  infimum is **exactly 1** — approached only as q → 0, and only
  **logarithmically** (c\*(q) = h(2q)/h(q) − 1 with
  h(2q)/h(q) → 2 like 1 + 1/log₂(1/q)). So L1 is asymptotically tight
  with **no expansion to take** at the tight corner: there is no
  leading order to match, only a log.

**P3. Which arm of t binds.** The branch samples split evenly between
t = 1−x and t = x−1/2, but of the configurations with ratio < 1.5,
**every single one has t = 1−x = q** (375 versus 0). So at every tight
configuration x ≥ 3/4 and L1 reads simply

    q·Δ ≥ −B.

**P4. The c\*-free simplification is useless.** The crude bound
−B ≤ x h(p₀) + (1−x) h(p₁) is true on all 600,000 branch samples (it
just drops the −c\*(p) terms), but t·Δ dominates that bound on
**0 of 600,000** — it fails on 100% of the branch. Dropping c\*
discards the entire margin: **c\* is essential to L1**, and any proof
must carry it.

**P5. The minimiser has two regimes.** At fixed q, the minimising
(p₀,p₁) is the corner (0,1/2) only for small q; above q ≈ 0.03 it
moves into the interior and beats 1/c\*(q):

    q       min ratio   at (p₀,p₁)          1/c*(q)    regime
    0.3     2.35102     (0.42857, 0.00000)  7.42395    interior
    0.2     1.21754     (0.18661, 0.25357)  2.89905    interior
    0.1     1.39645     (0.07937, 0.28571)  1.85423    interior
    0.05    1.47058     (0.03515, 0.33214)  1.56845    interior
    0.03    1.46101     (0.00000, 0.50000)  1.46101    corner
    0.01    1.33218     (0.00000, 0.50000)  1.33218    corner
    0.001   1.21278     (0.00000, 0.50000)  1.21278    corner

  (Classification independently reproduced by the skeptic at every q.)

## Outcome

- **NOT PROVED: L1.** What this record delivers is the sharp shape of
  the statement, not the statement.
- **VERIFIED (identity): ratio(q, 0, 1/2) = 1/c\*(q) exactly**, so the
  branch infimum is exactly 1 and is attained nowhere — approached
  only as q → 0.
- **EVIDENCE (strong, descent-attacked): C/(−B) ≥ 1 on the whole
  branch**, floor 1.0881 under a minimising descent.
- **PROVED (negative): the c\*-free simplification cannot work** —
  t·Δ dominates the c\*-free bound on 0 of 600,000 branch samples.
- **CORRECTION to 047's lead 1:** its advice to "expand in (1−x) at
  x → 1" is misdirected. In the scale-free ratio there is no
  polynomial corner to expand: the approach to the infimum is
  logarithmic in q, and separately, at every tight configuration
  t = 1−x = q, which is the simplification actually worth having.
  (047's underlying observation — that the *absolute* margin vanishes
  at x → 1 — is correct but not useful, the same trap 046 §G recorded
  at the other boundary.)
- **Not claimed:** any proof; anything outside branch L1; that the
  interior minimiser of P5 is the true minimum (it is a grid
  minimum on a 140² mesh per q).

## Why it failed / what survived

The proof did not close, and the reason is now specific rather than
vague: L1's infimum is exactly 1, so **no argument with any slack to
spare can work** — every bound in a proof must be asymptotically sharp
as q → 0. That rules out the whole family of crude-bound arguments
(P4 kills the most natural one outright) and explains why the
concavity route of 046 was never going to be enough: concavity gives a
bound tight only at p₀ = p₁, whereas L1's extremal is
(p₀,p₁) = (0,1/2), maximally *un*equal.

What survives, and is worth carrying: the tight configuration is
completely explicit — one conditional at 0, the other at 1/2, with
t = q — and along it the inequality is the single scalar statement
c\*(q) ≤ 1, i.e. h(2q) ≤ 2h(q). That is elementary and true. The open
part is everything off that configuration, where the ratio is larger
but the current tools do not see why.

## Leads generated

1. **Prove L1 as a perturbation of its extremal.** The extremal
   configuration and its value are now exact: at (0,1/2), ratio =
   1/c\*(q) ⟺ h(2q) ≤ 2h(q). Show the ratio is minimised there for
   q ≤ 0.03 (P5's corner regime) by a monotonicity argument in p₁,
   and handle q > 0.03 separately, where P5 says the minimum is
   comfortably above 1 (≥ 1.21 on the mesh).
2. **h(2q) ≤ 2h(q) is the shadow of the whole conjecture.** The
   extremal statement of L1 is exactly subadditivity of h at the
   doubling point. Worth checking whether the n ≥ 3 analogue of L1
   reduces to h(kq) ≤ k·h(q) — if so, the (HU-TAX) family has an
   elementary skeleton and the difficulty is entirely in the
   off-extremal region.
3. **Certify L1 on a q-grid** with the exact-rational kit, using P3's
   simplification (t = q at every tight configuration) to drop a
   variable. Cheaper than the 046 box, and it would upgrade P1 from
   sampled to certified.

## References

- This repo: 047 (branch L1 and the dichotomy), 046 (the
  decomposition, N2-CONC, the corrected c\*, and §G's absolute-vs-
  scale-free lesson), 045/044 (the promoted conjecture and the
  own-constant standard), 031 (the averaging obstruction).
  `data/hu_L1.json`.
- No external sources.
