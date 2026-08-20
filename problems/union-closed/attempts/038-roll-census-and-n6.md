# 038 — Roll census + n=6/7: rollout is no best-order proxy, but it provably dominates canonical and never goes negative

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-20
- **Mode:** informed
- **Type:** census + adversarial attack + one transcribed mini-theorem
  (037 leads 3 and 1).
- **Tools:** `explore/uc_hu_rollcensus.py` (new; the census with full
  order enumeration, n = 6 and n = 7 rollout descents; deterministic,
  fixed seeds 938000+/941000+; checkpoint `data/hu_rollcensus.json`);
  `explore/uc_hu_rollcensus_skeptic.py` (new; summaries recomputed from
  rows, census spot re-evaluations and every descent endpoint through
  the 037-skeptic independent stack; exit 0). Reproduce: run the two in
  that order.
- **Sources:** rollout improvement property: D. Bertsekas, rollout
  algorithms for combinatorial optimization (classical; [T] — argument
  transcribed and re-derived below, not looked up verbatim).

## Approach

037 left two edges open: is rollout actually a best-order proxy (it
found the exact best order on all three record witnesses), and does it
survive above n = 5, where all of 037's descents lived? Both are
directly falsifiable: a census family where roll ranks low is the
counterexample seed for rule design, and an n = 6 kill would bound
rollout's reach the way 035 bounded canonical's.

## What was done

**P1. Census, 300 random in-regime instances per block** (n ∈ {4, 5},
caps ∈ {0.45, 0.49}, full order enumeration alongside, canonical
ranked for contrast):

    block          roll=best   worst rank   worst (best−roll)/H   roll<0   canon>roll
    n=4 cap 0.45     60.3%        7/24           0.101              0         0
    n=4 cap 0.49     57.0%        9/24           0.207              0         0
    n=5 cap 0.45     39.0%       61/120          0.153              0         0
    n=5 cap 0.49     52.7%       41/120          0.091              0         0

  Three readings. (i) **Rollout is not a best-order proxy**: it is
  exactly best only 39–60% of the time and can leave 0.2·H on the
  table — the three 037 witnesses (all exactly best) were descent
  products, not typical measures. (ii) **It never went negative**:
  0 of 1200, including the 5 instances where some order IS negative.
  (iii) **Canonical is strictly better than roll on 0 of 1200**
  (equal on 52/52/24/28, strictly worse on the rest) — which turned
  out to be a theorem, below.

**P2. n = 6 rollout descents** — the three starts 037 skipped plus
four fresh hostile seeds, caps 0.49 and 0.497: **zero violations**,
global floors +0.030430 and +0.053740. (The floors sit far above the
n ≤ 5 floors — at n = 6 this descent stalls early; the honest reading
is "no kill found", not "the floor is sharp".)

**P3. n = 7 spot descents** (3 starts, reduced rounds, cap 0.49):
zero violations, floor +0.099692. Same caveat.

**P4. Mini-theorem (ROLL-DOM), a transcription:** for every μ and
every n, CR_roll(μ) ≥ CR_canon(μ). Proof: the canonical rule is a
deterministic function f of (μ, revealed set), so its completion
satisfies canon-completion(P) = [f(P)] + canon-completion(P + f(P)).
Let V_k be the rollout's step-k score max_i CR(P_k + i +
canon-completion). The sequence realizing V_k is available verbatim as
the candidate i = f(P_{k+1}) at step k+1, so V_{k+1} ≥ V_k; and V_0 ≥
CR(canon) because i = f(∅) is a step-0 candidate. Since V_{n−1} =
CR(roll), done. This is the classical rollout improvement property
specialized to HU orders; the only fact it needs is that the
completion policy is a deterministic function of the revealed SET,
which the canonical rule is by construction (035). **Consequences:**
rollout inherits every canonical positive of record wholesale — the
030 HU-mix sub-genre theorem's floor, the 0.38271 and 0.45 descents of
035 — and the census's canon>roll = 0 column is explained, not
coincidence. Numerically: canon_rank ≤ roll_rank on all 1200 rows.

## Outcome

- **EVIDENCE: rollout-HU survives at n = 6 (caps 0.49, 0.497) and
  n = 7 (cap 0.49)** — zero violations; floors are stall points, not
  sharp.
- **EVIDENCE (negative, informative): rollout is not a best-order
  proxy** — worst census rank 7/24, worst gap 0.207·H. 037 lead 3's
  optimistic reading is refuted; the gap between roll and the oracle
  is real on typical measures, just never (in 1200 instances)
  sign-changing.
- **PROVED (transcribed, same-session — reviewer re-derivation queued
  with the 036/037 batch): (ROLL-DOM)** CR_roll ≥ CR_canon for every
  μ. Rollout can only improve on canonical, so every canonical
  positive result transfers.
- **Not claimed:** sharpness of any n ≥ 6 floor (the descent visibly
  stalls); anything about caps above 0.497 at n ≥ 6; that the census
  generator's instance class is adversarial (it is the 035-style
  random-support class; descent products are known to be harder).

## Why it failed / what survived

The census kills the tempting conflation: rollout's value is NOT that
it approximates the best order (it often doesn't) — it is that its
one-step-lookahead-to-the-end never seems to land below zero, and
provably never lands below canonical. That reframes 037 lead 2: the
proof target worth having is not "roll ≈ best" but the sandwich
CR_best ≥ CR_roll ≥ CR_canon, where the left member is what part D
measures and the right member is now proved. A proof of best-order
positivity plus an effective bound on the roll-vs-best gap would make
roll a certified implementation of an existence statement.

## Leads generated

1. **(ROLL-DOM) generalizes:** the same argument gives CR_rollout ≥
   CR_π for ANY deterministic set-function completion policy π — e.g.
   rollout-on-rollout, or rollout completing with the identity order.
   Is there a completion policy whose rollout provably exceeds 0 on a
   genre where canonical does not? The [1/4,1/2]-component genre is
   free (HU-mix is order-free); the first non-trivial target is the
   {0}∪[1/2,1) genre of 004/025.
2. **The worst-gap census instances** (0.1–0.2·H below best, stored in
   the checkpoint rows by index) are the seeds for any future rule
   improvement — a rule should be tested on exactly those before any
   new descent campaign.
3. **n = 6+ descents need a better adversary** before the caps story
   is believed above n = 5: the current move set stalls 40× above the
   n ≤ 5 floors. Port 033's atom add/drop + order-pool moves to the
   rollout objective, or anneal.
4. Unchanged: 037 leads 1 (n = 6–7 at larger supports feeds into 3
   above), 2 (best-order positivity as proof target — now with the
   sandwich framing), 4 (near-product attractor).

## References

- This repo: 037 (the rules and floors under test), 035 (kill
  witnesses, descent move set), 036 (independent stack used by the
  skeptic), 031/030 (HU, HU-mix). `data/hu_rollcensus.json`.
- External: Bertsekas' rollout improvement property [T], transcribed
  in P4; no other sources.
