# 044 — The own-constant audit finds the kill the weaker flag hid: rollout-HU is REFUTED at cap 0.497, certified

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-20
- **Mode:** informed
- **Type:** audit sweep + refutation + exact certification (043's
  standing audit caveat, executed window-wide).
- **Tools:** `explore/uc_hu_ownconst.py` (new; every stored endpoint of
  the window's checkpoints plus the regenerated 1200-instance 038
  census, re-scored under rollout against c*(own max marginal);
  deterministic; checkpoint `data/hu_ownconst.json`);
  `explore/uc_hu_ownconst_certify.py` (new; the kill certified under
  each kit alone with the 60-digit rollout fixed point and a certified
  positive best-order control; checkpoint
  `data/hu_ownconst_certify.json`). The kill was confirmed by BOTH
  float stacks (engine and the independent 037-skeptic stack, agreeing
  to 1e-12) before certification. Reproduce: run the two in that
  order.
- **Sources:** none.

## Approach

043 flagged that every engine of this window tested violations as
CR/H < 0, weaker than the conjectured CR/H ≥ c*(fmax) with fmax the
instance's OWN maximum marginal (034's lesson). This attempt closes
that audit gap: re-score all 212 stored endpoint measures across
`hu_order2` / `hu_rollcensus` / `hu_roll_anneal` / `hu_blocks` /
`hu_canon` (the last under rollout — 035's canonical endpoints are the
hardest instances of record), plus the regenerated census, against
their own constants under the rollout order.

## What was done

**A. 212 stored endpoints.** 211 have own-margin ≥ 0 (the
family-saturated anneal endpoints sit at exactly 0, as they must —
they ARE the equality family). **One violation:** the
`D_bestorder_cap_0.497` endpoint bred from `floor:windowkill` — a
measure optimized to minimize the BEST-ORDER ratio, stored in 037's
checkpoint, and never before scored under rollout.

**B. The kill, dissected and certified.** The witness is a 9-atom
n = 4 measure, max marginal exactly 0.49500 < 0.497, H = 3.137:

    CR under the rollout order (1,0,2,3):  −9.596936657e-4   CERTIFIED < 0 (each kit alone)
    18 of 24 orders negative; rollout ranks 13/24
    best order (0,3,2,1):                  +1.610189393e-2   CERTIFIED > 0 (each kit alone)

  The rollout order of the rationalized measure is confirmed as a
  60-digit fixed point (037/043 pattern), so the certificate is about
  rollout-HU, no order-provenance caveat. Both float stacks agreed to
  1e-12 before certification.

**C. The census, own-constant.** All 1200 regenerated 038 census
instances: zero violations, minimum own-margin +4.49e-2 — random
in-regime measures sit far above their own constants; the kill genre
is, as always on this line, a descent product.

## Outcome

- **REFUTED (certified): rollout-order HU positivity above cap
  ≈ 0.495.** CR_roll < 0 in-regime at max marginal 0.495. Since
  0.495 < 1/2, **HU with the rollout order cannot by itself deliver
  Frankl** — the same verdict 035 delivered for the canonical order at
  0.49, one rung higher. The progression is now: canonical rescues
  0.45 and dies at 0.49 (035); rollout rescues 0.49 and dies at 0.495
  (here). 037's survival story at 0.497/0.499 was an artifact of
  attacking rollout directly — the killing measure was bred against
  the best-order oracle and transferred.
- **EVIDENCE: the best-order form survives this witness** (certified
  positive control, ratio +0.00513 ≥ c*(0.495) = 7.2e-5) and every
  other endpoint of record. The sandwich's top — best-order
  positivity, 038 lead 2 — is now the ONLY unrefuted for-all-μ
  positivity statement on the HU line, exactly as 038's framing
  anticipated.
- **VERIFIED (audit): the window's other 211 stored endpoints and all
  1200 census instances clear their own constants** under rollout.
- **Method lesson, adopted:** violation flags in future engines must
  test CR/H < c*(own fmax), not CR/H < 0 — the weaker flag sat on
  this kill for four records; and endpoints bred against one
  objective must be re-scored under every standing rule before a
  survival claim is recorded (transfer attacks are free).
- **Not claimed:** that 0.495 is rollout's exact threshold (unrefuted
  at 0.49; the crossing is bracketed in (0.49, 0.495] only by these
  families); anything against best-order HU.

## Why it failed / what survived

"It" here is rollout-HU's candidacy for the whole (0, 1/2) range, and
it failed the way every fixed order rule so far has failed: a
descent-bred measure near the cap where the rule's one-step lookahead
picks a poisoned prefix (rollout ranks 13/24 on its own kill — worse
than the middle). What survives is exactly the structure 038 and 040
predicted: the DIAG ceiling says no rule can guarantee more than
c*(p̄), and this record says fixed rules do not even keep positivity
to 1/2 — so the proof-shaped object is the EXISTENCE statement
(best-order positivity, oracle-measured floors still clean), with any
fixed rule certified per-instance where it happens to work. Recipe v3
should carry rollout-HU to a defended cap of 0.49 on current
evidence, and the interval [0.49, 1/2) currently belongs to no rule.

## Leads generated

1. **Attack best-order HU at caps 0.495–0.499 directly** (n ≤ 5 full
   enumeration inside the descent, seeded at this witness and the D
   endpoints): the one remaining positivity statement, and 037 part
   D's floors were never own-constant audited at these caps beyond
   this sweep. A best-order kill below 1/2 would end the
   order-quantifier program entirely; survival sharpens 038 lead 2
   into the line's main conjecture.
2. **Transfer-attack protocol:** every stored endpoint of every future
   campaign gets scored under all standing rules (the sweep here is
   the template and is cheap); adopt as a standing skeptic step.
3. **Where exactly does rollout cross?** Bracket in (0.49, 0.495] by
   descents against rollout seeded at THIS witness's family at
   intermediate caps — cheap, and calibrates how much cap each rung of
   lookahead buys (canonical 0.44→0.49, rollout 0.49→~0.495 …
   diminishing).

## References

- This repo: 043 (the audit caveat this executes), 042/041/040 (the
  equality-family context), 038 (the sandwich framing this confirms),
  037 (the checkpoint the kill was hiding in), 035 (the canonical
  precedent), 034 (own-marginal lesson), 036 (fixed-point pattern via
  its committed module), 031/030 (HU). `data/hu_ownconst.json`,
  `data/hu_ownconst_certify.json`.
- No external sources.
