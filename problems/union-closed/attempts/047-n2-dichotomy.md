# 047 — The n = 2 target splits into a dichotomy: every deficit is covered by ONE term, and the q = f₀ branch is always the interaction

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** structural analysis + one proved lemma (046 lead 1).
- **Tools:** `explore/uc_hu_n2_dichotomy.py` (new; the A/B/C
  decomposition at scale, the dichotomy census, the split by dominant
  marginal, tightness of each branch; deterministic, seeds 4701–4704;
  checkpoint `data/hu_n2_dichotomy.json`);
  `explore/uc_hu_n2_dichotomy_skeptic.py` (new; terms rebuilt from
  this record's prose in nats, and — the check that matters — the
  decomposition re-verified against a from-scratch history-recursion
  HU evaluator sharing no code with the 037/046 stack, plus an
  adversarial descent that maximises the shortfall; exit 0).
  Reproduce: run the two in that order.
- **Sources:** none.

## Approach

046 reduced the n = 2 case of the promoted best-order conjecture to an
explicit scalar inequality and showed the pair-interaction term is
load-bearing. Write that inequality as **A + B + C ≥ 0**:

    A = (c*(f₀) − c*(q))·h(f₀)                    first-coordinate slack
    B = x(σ(p₀) − c*(q)h(p₀)) + (1−x)(σ(p₁) − c*(q)h(p₁))   diagonal ledger
    C = t·Δ,  t = min(x−1/2, 1−x),
        Δ = 2ψ(p₀+p₁) − ψ(2p₀) − ψ(2p₁)           pair interaction

with ψ(s) = h(min(1/2,s)), σ(p) = c*(p)h(p) = h(z*(p)) − h(p). A ≥ 0
because c* is decreasing and q ≥ f₀; C ≥ 0 is 046's Lemma N2-CONC
(ψ concave). **Only B can be negative.** The question this record asks
is the obvious next one: when B < 0, what pays for it?

## What was done

**A. Lemma N2-ONE-ABOVE (proved).** B < 0 requires a conditional
marginal above q, and **at most one can be**: the x-weighted mean of
the two conditionals is exactly f₁ ≤ q, so they cannot both exceed q.
(0 counterexamples in 200,000 in-regime samples, and 0 under 400,000
draws targeted at the region where both are large.) Together with
046's N2-ONE-BAD this pins the deficit structure completely: one
conditional above q, contributing a deficit x·h(p)·[c*(q) − c*(p)],
and everything else nonnegative.

**B. The dichotomy.** Over 200,000 in-regime samples, 42,575 have
B < 0. Of those:

    covered by the slack A alone (A + B ≥ 0):        31,982
    covered by the interaction C alone (B + C ≥ 0):  41,530
    needing BOTH terms:                                   0
    failing outright:                                     0

**So the n = 2 Case-A theorem is equivalent to: B < 0 ⟹
max(A, C) ≥ −B.** The two mechanisms are not additive partners; each
deficit is paid by one of them alone.

**C. The split is by which coordinate carries the max marginal.**

    q = f₀ (first coordinate dominant, so A = 0):
        7,435 deficit cases; the interaction alone fails on 0
    q = f₁ (second coordinate dominant):
        34,650 deficit cases; interaction alone fails on 1,017,
        slack alone fails on 2,899 — and never both

  The q = f₀ branch is the clean one: there A vanishes identically, so
  the statement is **L1: if q = f₀ then C ≥ −B** — the pair
  interaction alone covers the conditional deficit whenever the
  *first* coordinate carries the maximum marginal. That is a
  two-variable inequality (x and the pair p₀,p₁) with no competing
  terms, and it is the smallest unproved statement on the whole HU
  line.

**D. Tightness.** Minimum of B + C on the L1 branch is **+1.17e-04**
(at x = 0.9992, p₀ = 0.0002, p₁ = 0.0712) and the minimum of
max(A,C) + B over all deficit cases is **+1.66e-05**; the skeptic's
descent, which actively maximises the shortfall, gets no lower than
**+4.3e-07** (at x = 0.98447, p₀ = 0.01555, p₁ = 0.01466). Both
branches approach equality only as x → 1 — the n = 1 degeneration,
where t → 0, C → 0 and B → 0 together. Away from that corner the
margins are not small, so **a proof of L1 can afford crude bounds
except near x = 1**, where the correct move is to expand in
(1−x) rather than to bound.

## Outcome

- **PROVED (hand, skeptic-attacked): Lemma N2-ONE-ABOVE** — at most
  one conditional marginal exceeds q, because their x-weighted mean is
  f₁ ≤ q. With 046's N2-ONE-BAD, the deficit structure at n = 2 is now
  completely determined.
- **EVIDENCE (strong, skeptic-attacked): the dichotomy** — in 42,575
  deficit cases, every one is covered by A alone or by C alone, none
  by their sum, none failing; and an adversarial descent on the
  shortfall bottoms out at +4.3e-07.
- **EVIDENCE: branch L1** — when the first coordinate carries the max
  marginal, the interaction alone always covers (0 failures in 7,435
  cases). This is the record's proposed next proof target and it is
  now stated in two variables.
- **VERIFIED: the decomposition** — A + B + C equals the margin to
  **2.9e-15** against an independent history-recursion evaluator, so
  the dichotomy is a statement about the right quantity.
- **Not claimed:** a proof of L1 or L2; anything at n ≥ 3; that the
  dichotomy extends beyond n = 2 (it is stated for the n = 2
  decomposition only, where "the" interaction term is a single
  scalar).

## Why it failed / what survived

Nothing failed; the line advanced from "the interaction is
load-bearing" (046) to "the interaction is load-bearing *exactly when
the first coordinate dominates*, and then it suffices by itself". The
mechanism reading: when q = f₀ the revealed coordinate is the one
carrying the extremal marginal, so the ledger's first-coordinate slack
A is exactly zero by construction — there is nothing to donate — and
the surplus must come from the coupling's own pair structure. When
q = f₁ the first coordinate is *sub*-extremal, and the difference
c*(f₀) − c*(q) is a genuine donation that usually suffices.

That is a satisfying shape, and it also explains 031's averaging
problem from the other side: the schemes that failed there were
averaging schemes, and averaging destroys precisely the term that the
q = f₀ branch depends on entirely.

The honest limit: this is a census plus two proved sign lemmas, not a
theorem. Both branches remain unproved inequalities, and the
degenerate corner x → 1 is where both approach equality.

## Leads generated

1. **Prove L1: q = f₀ ⟹ t·Δ ≥ −B.** Two variables after eliminating
   x via q = 1−x. Δ has the closed form
   2ψ(p₀+p₁) − ψ(2p₀) − ψ(2p₁); −B is
   c*(q)(x h(p₀)+(1−x)h(p₁)) − xσ(p₀) − (1−x)σ(p₁). Expand at x → 1
   (the only tight corner) and bound crudely elsewhere.
2. **Prove L2 by the same split**: when q = f₁ and A < −B, show
   C ≥ −B. The census says this sub-case has 1,017 instances in
   200,000 draws — a thin region worth characterising before
   attacking.
3. **Does the dichotomy survive at n = 3?** There the interaction is
   no longer one scalar; the falsifiable version is "each deficit is
   covered by a single term of the n = 3 decomposition". If it fails
   at n = 3, the n = 2 proof will not generalise and that is worth
   knowing before investing in it.

## References

- This repo: 046 (the decomposition, N2-ONE-BAD, N2-CONC, the
  corrected c*), 045 (the promoted conjecture), 044 (own-constant
  margins), 042/041 (the equality families), 037 (the cell ledger),
  031 (the averaging obstruction this explains from the other side).
  `data/hu_n2_dichotomy.json`.
- No external sources.
