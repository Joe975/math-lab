# 051 — Correction to 050 part D: the n = 3 descents saturated by degenerating, and the honest n = 3 floor is certified

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** self-correction + adversarial descent + certification
  (050's own part D, reviewed same day; the repo forbids editing 050,
  so the correction lives here).
- **Tools:** `explore/uc_hu_n3_essential.py` (new; the degeneracy
  audit, essential-support descents, and the certificate;
  deterministic, seeds 6100/6149; checkpoint
  `data/hu_n3_essential.json`);
  `explore/uc_hu_n3_essential_skeptic.py` (new; every value recomputed
  with the 050 skeptic's independent nats history-recursion evaluator,
  the degeneracy re-derived from 050's committed checkpoint, and the
  constraint itself tested for vacuity; exit 0). Reproduce: run the
  two in that order.
- **Sources:** none.

## Approach

050 part D reported that best-order descents at n = 3 reach a minimum
own-constant margin of **+0.000000** with zero violations, and read
that saturation as the equality family binding — "consistent with
042's block-tensor equality family". Inspecting the endpoints
immediately after committing showed that reading is wrong, and the
correction is worth more than the original claim.

## What was done

**A. The degeneracy.** 050's two sharpest part-D endpoints are:

    cap 0.45 (n3seed1): marginals [0.4462, 0.0, 0.0], live atoms 2
    cap 0.49 (n3seed6): marginals [0.4895, 0.4895, 0.0], live atoms 2

  Both have an atom driven to weight 0 and at most two coordinates
  carrying mass: they are **n = 1 instances embedded in n = 3**, where
  the bound holds with equality automatically (a single coordinate at
  its own extremal — one of 046's own equality families). The descent
  escaped to a degenerate face of the simplex rather than finding an
  n = 3 near-extremal. Re-derived independently by the skeptic from
  050's committed checkpoint.

**B. The honest n = 3 floor.** Re-running the descent under an
**essentiality constraint** — every coordinate marginal ≥ 0.05, at
least 4 atoms with weight ≥ 1e-3 — confines the adversary to genuinely
three-dimensional instances (the constraint is not vacuous: it rejects
659 of 3000 random in-regime draws):

    cap 0.45: 10 descents, min own-constant margin  +0.000018, 0 violations
    cap 0.49: 10 descents, min own-constant margin  +0.000249, 0 violations

  Three orders of magnitude above the degenerate 2.6e-09, and still
  strictly positive.

**C. Certified.** The sharpest essential endpoint (max marginal
0.49952 — the descent is regime-bounded, not cap-bounded, and pushed
the marginal above its starting cap) is certified positive under
**each kit alone**: CR_HU ∈ [+3.691112270e-05, ·] (kit A) and
[+3.691112263e-05, +3.691112270e-05] (kit B), with the exact marginal
< 1/2 and H > 1/2 checked in rationals.

## Outcome

- **CORRECTION to 050 part D:** its "+0.000000, saturating the
  equality bound" is a **degeneration artifact** — the descent left
  the three-dimensional problem. 050's part-D sentence "consistent
  with 042's block-tensor equality family" should read "the descent
  escapes to the n = 1 face, where equality is automatic". Per the
  repo rule 050 is left as written.
- **EVIDENCE (certified at the sharpest point): the n = 3 best-order
  floor over genuinely 3-dimensional instances is positive** —
  +0.000018 at cap 0.45 and +0.000249 at cap 0.49 over 20 descents,
  with the sharpest certified > 0 under each kit alone.
- **050's headline results are untouched:** the census (parts A–C),
  the order-quantifier finding and its certified negative witness do
  not involve the part-D descents at all.
- **Method lesson (third in this window, and the same shape as the
  other two):** an adversary minimising a margin will find the
  *cheapest* way to make it vanish, and degenerate faces are usually
  cheapest. 046 §G caught the q → 1/2 face, 046 §G′ caught depth-first
  non-coverage, and this catches the low-dimensional face. Descent
  campaigns on this route should carry an explicit non-degeneracy
  constraint, and report what fraction of draws it rejects.
- **Not claimed:** that 20 essential descents exhaust the n = 3
  adversary space; anything about caps above 0.49; that 0.05 / 4 atoms
  is the canonical essentiality threshold (it is a choice, stated so
  the numbers are reproducible).

## Why it failed / what survived

What failed is a reading, not a computation: every number 050 reported
in part D is correct, and the skeptic reproduced them — they simply
describe a degenerate family. The lesson generalises past this route:
three times in one window, an adversary reported a vanishing margin
that came from a degeneracy (a vanishing constant, a vanishing
coverage, a vanishing dimension) rather than from the inequality being
nearly violated. The fix each time was to make the objective
scale-free or the region essential.

What survives is stronger than before: with degeneracies excluded, the
n = 3 best-order floor is positive by a real margin and certified.

## Leads generated

1. **Re-audit the earlier descent campaigns for degeneracy.** 037/038/
   040/045 descents were never checked for support collapse; the
   essentiality diagnostic here is three lines and would say whether
   any recorded "floor" is a degenerate face. Cheap, and it would
   either confirm those floors or correct several at once.
2. **Adopt essentiality as a standing descent constraint** for this
   route, alongside 044's own-constant flag and 046's scale-free
   margin.

## References

- This repo: 050 (the record corrected here), 046 (§G and §G′, the two
  earlier degeneracy lessons, and the n = 1 equality family), 044 (the
  own-constant standard), 042 (the equality family 050 invoked).
  `data/hu_n3_essential.json`.
- No external sources.
