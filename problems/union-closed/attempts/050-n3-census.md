# 050 — n = 3, the size the line skipped: best-order survives, and the order quantifier is vacuous only at n = 2

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** census + adversarial descent (045 lead 1 at an untested
  size; 047 lead 3).
- **Tools:** `explore/uc_hu_n3_census.py` (new; own-constant census at
  n = 3 with full 6-order enumeration, n = 4 contrast, the
  order-quantifier measurement, best-order descents; deterministic,
  seeds 5000–5300; checkpoint `data/hu_n3_census.json`);
  `explore/uc_hu_n3_certify.py` (new; the order-failure witness
  certified under each kit alone; checkpoint
  `data/hu_n3_certify.json`); `explore/uc_hu_n3_census_skeptic.py` (new; HU/CR/c* rebuilt in nats
  by an explicit history recursion sharing no code with the 037/046
  stack, every census block regenerated and recomputed, plus its own
  60-start descent attacking n = 3; exit 0). Reproduce: run the two in
  that order.
- **Sources:** none.

## Approach

Every measurement on this route has been at n = 2 (046–048) or n ≥ 4
(033/035/037/038/040/044/045). **n = 3 was never examined at all** —
and it is the first size with a second interaction layer, so it is
exactly where the n = 2 structure either generalises or breaks. Two
questions, both falsifiable:

1. does best-order positivity — 045's promoted conjecture — survive at
   n = 3, judged against c\*(own max marginal) per 044's standard?
2. is the order quantifier vacuous at n = 3? 046 proved it vacuous at
   n = 2 (identity always suffices); 033/035/044 proved it essential
   at n ≥ 4. n = 3 decides whether n = 2 was special or whether the
   vacuity persists.

## What was done

**A/B. Own-constant census** (random in-regime measures, full order
enumeration, margins against c\*(own fmax)):

    block          instances   min margin: best / rollout / identity   violations (b/r/i)
    n=3 cap 0.45      3000      +0.01732 / +0.01732 / +0.00437              0 / 0 / 0
    n=3 cap 0.49      3000      +0.00124 / +0.00124 / −0.00696              0 / 0 / 2
    n=4 cap 0.45      1500      +0.03447 / +0.03447 / +0.01544              0 / 0 / 0
    n=4 cap 0.49      1500      +0.01403 / +0.01403 / −0.00515              0 / 0 / 1

  At cap 0.49, 21 of the 3000 n = 3 instances have *some* order with
  CR < 0 (13 of 1500 at n = 4) — the order-dependence 033 found at
  n = 4 is already present at n = 3.

**C. The order quantifier is NOT vacuous at n = 3.** The identity
order fails the own-constant bound on 1 of 2000 fresh n = 3 instances
at cap 0.49 (and 2 of 3000 in the census block above), while rollout
meets it on 2000/2000. So **046's vacuity result is special to
n = 2**: at n = 3 the coupling already needs a genuine order choice,
and n = 3 is not a boundary case — it behaves like n ≥ 4. The finding
has content: at n = 3 the identity order coincides with the best order
on only 54 of 300 instances.

**D. Best-order descents at n = 3** (8 starts per cap, own-constant
flagged): zero violations at caps 0.45 and 0.49, with the minimum
own-constant margin descending to **+0.000000** — i.e. the adversary
saturates the equality bound rather than crossing it, which is the
same behaviour 045 measured at n ≤ 5 and consistent with 042's
block-tensor equality family. The skeptic's independent 60-start
descent reaches −9.6e-13, i.e. zero to float precision, with no
crossing.

**E. The failure, certified.** The sharpest identity-order failure is
stronger than "below c\*" — its CR is **negative**. The witness is a
4-atom n = 3 measure with max marginal exactly 0.485:

    mu = {∅: 0.515, {0,1}: 0.046779089, {0,2}: 0.317092622,
          {0,1,2}: 0.121128289}

    identity order (0,1,2):  CR_HU certified in [·, −1.005300143e-2]
    best order (1,2,0):      CR_HU certified in [+2.775518457e-1, ·]

  both under **each kit alone** (029 standard), with the exact
  marginal < 1/2 and H > 1/2 checked in rationals. So the statement
  "the order quantifier is essential at n = 3" is certified, not
  merely sampled.

## Outcome

- **EVIDENCE: best-order positivity survives at n = 3** — 0 violations
  across 4,500 census instances and 16 descents at caps 0.45/0.49,
  with descents saturating equality rather than crossing it. 045's
  promoted conjecture now has support at every size the line has
  examined (n = 2 through 7).
- **EVIDENCE: rollout also survives at n = 3** (0 violations,
  2000/2000 in part C), consistent with 044's finding that rollout is
  safe below cap ≈ 0.495.
- **REFUTES (certified) the natural extension of 046's vacuity claim:
  the order quantifier is essential from n = 3 onward.** The identity
  order fails the own-constant bound at n = 3 (2 of 3000 at cap 0.49,
  min margin −0.00696), and at the sharpest witness its CR is
  certified **negative** (−1.005300143e-2, each kit alone) while the
  best order is certified positive (+2.775518457e-1). n = 2 is the
  only size where the quantifier is empty.
- **Not claimed:** certificates for the positive results (the census
  and descents are float-level; only the n = 3 order-failure is
  certified); anything about caps above 0.49 at n = 3; that 16
  descents exhaust the n = 3 adversary space.

## Why it failed / what survived

The interesting negative is the vacuity: 046 established that at n = 2
the deficit cells contribute exactly 0 and one order always suffices,
and it was tempting to read that as the beginning of a pattern. It is
not — one extra coordinate is enough to restore the order-dependence,
which fits 046's own mechanism: the n = 2 proof of N2-ONE-BAD used
that *both* conditionals being out of regime forces f₁ > 1/2, and that
argument has no analogue once there are three coordinates and a
history of length two. So the n = 2 structure is a genuine accident of
the smallest case, and a proof strategy built on it (046/047/048's
lemmas) will need new ingredients at n = 3, not just more algebra.

What survives and grows stronger: best-order positivity, now measured
at every size from 2 to 7, always with the equality family binding and
never crossed.

## Leads generated

1. ~~Certify the n = 3 identity-order failures~~ **DONE in part E**
   (each kit alone; CR certified negative, best-order control
   certified positive). Remaining: certify one of the positive
   best-order floors at n = 3 too, so both directions are exact.
2. **What breaks N2-ONE-BAD at n = 3?** The lemma's proof is a
   one-liner about two conditionals; the n = 3 analogue would be about
   the conditional structure after a length-2 history. Write down what
   the correct statement is (or find the counterexample instance in
   this record's checkpoint) before attempting any n = 3 proof.
3. **Descend harder at n = 3**: 16 descents is thin next to the 17-start
   campaigns of 033/035. Port the 040 anneal to n = 3 best-order and
   check whether the equality saturation is a stall or the true floor.

## References

- This repo: 046 (the n = 2 vacuity this bounds, and N2-ONE-BAD),
  047/048 (the n = 2 proof structure whose generality this tests),
  045 (the promoted conjecture), 044 (the own-constant standard),
  042 (the equality family the descents saturate), 038 (the census
  generator). `data/hu_n3_census.json`.
- No external sources.
