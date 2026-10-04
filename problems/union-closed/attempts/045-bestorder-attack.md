# 045 — Best-order HU survives its direct attack at caps 0.495–0.499; the oracle's floors saturate the equality family

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-20
- **Mode:** informed
- **Type:** adversarial attack (044 lead 1 — the line's last unrefuted
  positivity statement, attacked head-on where both fixed rules died).
- **Tools:** `explore/uc_hu_bestorder.py` (new; weight-descents on the
  max-over-all-orders CR/H at caps 0.495/0.497/0.499, n ≤ 5 with full
  enumeration inside the objective; 044's audit lesson baked in — every
  endpoint flagged both as KILL (ratio < 0) and SHARP (ratio <
  c*(own fmax)); deterministic; checkpoint `data/hu_bestorder.json`);
  `explore/uc_hu_bestorder_skeptic.py` (new; every endpoint's
  best-order ratio, own fmax and own margin re-derived by full
  enumeration through the independent 037-skeptic stack; exit 0).
  Reproduce: run the two in that order.
- **Sources:** none.

## Approach

After 044 the ladder reads: canonical order dies at cap 0.49 (035),
rollout dies at 0.495 (044), and best-order positivity — "for every
in-regime μ SOME revelation order has CR_HU ≥ 0" — is the only
for-all-μ statement left standing (038's sandwich top). This attempt
attacks it directly in the caps where the fixed rules died, with
hostile seeding: the 044 kill witness itself, the sharpest
D_bestorder endpoints from 037, equality-family members perturbed
off-family (043's finding that the adversary stalls near the family —
so start it there, off balance), and the standing starts.

## What was done

Descents at three caps (21/23/23 starts; every candidate scored by
full 24/120-order enumeration; endpoints flagged against their OWN
constant per 044):

    cap 0.495: global floor +2.190730e-4,  0 kills, 0 sharp violations   [c*(cap) = 7.21e-5]
    cap 0.497: global floor +2.190730e-4,  0 kills, 0 sharp violations   [c*(cap) = 2.60e-5]
    cap 0.499: global floor +1.295407e-4,  0 kills, 0 sharp violations   [c*(cap) = 2.89e-6]

Three structural observations:

1. **The binding endpoints saturate their own constants**: several
   descents (notably those seeded at 037's D endpoints) stall at
   own-margin +2.7e-9 — the attacked best-order floor equals
   c*(own fmax) to float noise. Under the oracle, what binds is the
   042 equality family: exactly the (HU-TAX, best-order) shape with
   the family as equality set.
2. **Off-family perturbations descend BACK onto the family**: the
   (d,4) family point perturbed by ±35% multiplicative noise returns
   to own-margin exactly 0; the (d,2)⊗(d,2) perturbation stalls just
   above (+1.3e-4) — the family is an attractor for the oracle
   adversary, not just a stall shelf (contrast 043's correction about
   the ROLLOUT polish, which stalls near but off the family).
3. **The 044 kill witness is comfortably positive under best-order**
   (+3.07e-3 after further descent from it) — the rollout kill does
   not seed an oracle kill.

Skeptic: every endpoint's ratio, own fmax and own margin re-derived by
full enumeration through the independent evaluator stack; global
floors, kill lists and sharp-violation lists all consistent; exit 0.

## Outcome

- **EVIDENCE: best-order HU positivity survives its direct attack**
  at caps 0.495/0.497/0.499 (n ≤ 5, 67 descents, hostile-seeded),
  and so does the SHARP form — zero endpoints below c*(own fmax).
  The attacked floors saturate the equality family's constant from
  above, which is what a true (HU-TAX, best-order) conjecture would
  make them do.
- **038 lead 2 is hereby promoted to the line's main conjecture,
  SPECULATION as always:** for every μ with all marginals ≤ p̄ < 1/2,
  max over orders of CR_HU ≥ c*(p̄)·H(μ), equality exactly on the 042
  block-tensor family (which needs p̄ ≥ 1/4 per 034 for the constant
  as written). If proved it gives the recipe an existence license
  with rollout as the constructive finder below 0.49.
- **Not claimed:** anything at n ≥ 6 (the objective's enumeration cost
  confines the descents to n ≤ 5); that the descent adversary is
  strong at these caps (it is 035's move set; the anneal adversary of
  040 was not ported to the oracle objective — queued); no
  certificates this round (the floors are stall points except the
  family-saturated ones, whose family values 042 already certified).

## Why it failed / what survived

Nothing failed — and after a window in which two fixed order rules
were refuted exactly where their predecessors' campaigns had looked
safe, a survival claim earns its skepticism: the honest statement is
that the ORACLE form resisted the same families and seeds that killed
canonical and rollout, including transfer seeding, with the equality
family visibly binding. The structural gap to a proof is unchanged
(the per-cell surplus/deficit averaging of 031, now with an order
quantifier); what this record adds is that the target's measured
shape — floors pinned to c*(own fmax) on the family — is exactly the
conjectured one at the caps where everything else died.

## Leads generated

1. **Port the 040 anneal to the oracle objective** (n = 4 only, 24
   orders per candidate — affordable): the descent stalls are the weak
   point of this record's negative result.
2. **n = 6 oracle attack via sampled orders** (branch-and-bound or
   best-of-720-sampled): the n ≤ 5 confinement is the other weak
   point.
3. **Proof probe for the 2-block case first**: is
   max-order CR_HU ≥ c*(p̄)H on n = 2 blocks provable by the 041
   count-pair DP structure? The DP collapse makes n = 2 an
   honest-sized hand calculation and would be the first order-
   quantified positive theorem on the line.

## References

- This repo: 044 (the lead and the kill witness), 042/041 (the
  equality family and DP), 040 (DIAG ceiling; anneal to port), 038
  (the sandwich), 037 (D endpoints), 036/043 (audit standards),
  035/034/031/030 (the line). `data/hu_bestorder.json`.
- No external sources.
