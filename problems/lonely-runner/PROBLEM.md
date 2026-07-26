# Lonely Runner Conjecture

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** k runners on a unit circular track start together and run at
pairwise distinct constant speeds. The conjecture: each runner is at some time
"lonely" — at distance ≥ 1/k from every other runner.

## Published status

Proved up to k = 7 (k = 7 by Barajas–Serra). A September 2025 preprint by
Rosenfeld (arXiv:2509.14111) claims a proof of k = 8; treat k = 8 as likely
settled pending peer review, and k ≥ 9 as the open frontier.

Standard reduction: speeds may be assumed integers with one runner stationary,
and the problem is scale-invariant, so primitive (gcd = 1) tuples suffice. Tao
proved it suffices to check speeds up to roughly k^{O(k²)}, making each k a
(huge) finite problem.

The conjectured rigidity "tight ⇔ speeds {1, …, k−1}" is **false as stated**:
Goddyn–Wong constructed additional tight instances (accelerations and
sporadics). Fan–Sun have an amended spectrum conjecture.

## Verification contract

- ML(v) = sup_t min_i ‖v_i·t‖ must be computed in **exact rational
  arithmetic**. Floating-point ML values are not acceptable for a tightness
  claim, since the entire question is about exact comparison to 1/k.
- A **scan claim** must state the speed bound V and the primitivity
  convention, and must show that scaled copies of a primitive tuple were
  recognized as such rather than counted separately.
- A scan is `EVIDENCE` bounded by its V. Say the V.

## Harness (tier 0)

- `harness/lonely-runner/lonely_runner.py` — exact ML computation and
  tuple scanning; the reference implementation.
- `harness/lonely-runner/lonely_runner_analyze.py` — independent
  re-verification of scan output. Run it on any scan you intend to claim.
