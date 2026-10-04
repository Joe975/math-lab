# 002 — Independent audit of the eight-body certificates

- **Problem:** central-configurations-8
- **Date:** 2026-09-06
- **Mode:** informed
- **Type:** skeptic review
- **Reviewer:** independent_checker agent; record compiled from its written audit.
- **Reviews:** 001-three-route-certified-pilot
- **Tools:** independently written lifted-polynomial Fraction checker, exact distance-order bounds, adversarial tests.
- **Sources:** 001, data/skeptic-notes.md, and data/verification-summary.json.

## Claims attacked

1. Direct-force certificate existence and uniqueness, including the removed gauge equation.
2. Physical meaning of roots of a positive inverse-distance polynomial lifting.
3. All 39 saved polynomial certificates (19 discovery, 20 catalogue) represent n=8 unit-mass roots.
4. The 20 catalogue roots are pairwise inequivalent.
5. Complete uniqueness within the aligned, distinct concentric-square family.
6. Whether numerical searches and sampled continuation are being promoted to a proof of completeness.

## Refutations found

No mathematical refutation of the scoped claims was found. One API boundary
needs care: distinct_by_distances does not authenticate its inputs; callers
must verify certificates first. The current callers do.

The catalogue matching is numerical, not a proved isometry between exact roots.
The exact conclusion is existence of 20 inequivalent local roots. Neither
sampled trajectories nor the rectangular partition proves a complete catalogue.
The review does not independently certify continuous branch tracking, novelty,
or the coverage pilot's untried exclusion strength.

## Claims that survive

**VERIFIED, local scope:** the independent 43-variable polynomial implementation
uses no direct-kernel or interval code from central.py. Every s_ij interval is
positive, so s_ij^2*r_ij^2=1 selects physical inverse distance. The independent
Jacobian and rational contraction/self-map checks accept all 39 certificates.
Adversarial tests reject negative inverse distances, zero anchors, zero inverses,
wrong masses/dimensions, collisions, and shifted centers. Exact secant tests
exercise Jacobian enclosures. An eight-gon verifies with NumPy/SciPy imports
prohibited in the verifier.

The gauge proof survives: pairwise torque cancellation and x0!=0 recover the
omitted force equation, and summation recovers center of mass zero. The
contraction condition makes the proposed preconditioner nonsingular, so a fixed
point is a true zero. This does not prove global gauge uniqueness.

All 190 catalogue pairs have strictly separated squared-distance order
statistics. These were independently enclosed from polynomial coordinate boxes,
without importing central.py; exact separation margins are saved in
verification-summary.json. A separated multiset proves inequivalence; a matching
multiset would not prove equivalence.

The two-square radial equation and the strict decrease of H(t) were independently
re-derived, including cancellation of perpendicular cross terms, positivity of
the common lambda, endpoint limits, and the exclusion of t=1 collisions. The
argument proves precisely the aligned concentric-square family. The literature
check in 001 identifies it as a rediscovery.

The continuation code explicitly labels its curves numerical and sets complete
false. Branch jumps or missed branches affect search coverage, not validity of
independently certified endpoints. No completeness claim survives because none
was established.

Reproduce the saved-certificate audit:

```powershell
& C:/Repos/math-lab/.venv/Scripts/python.exe C:/Repos/math-lab/problems/central-configurations-8/explore/audit_certificates.py
& C:/Repos/math-lab/.venv/Scripts/python.exe -m pytest C:/Repos/math-lab/tests/test_cc8_polynomial.py -q -p no:cacheprovider
```

## References

- 001-three-route-certified-pilot.md.
- Full independent derivations and review: data/skeptic-notes.md.
- Machine-readable certificate and pairwise-separation audit: data/verification-summary.json.
- Polynomial-checking background uses the equations stated in PROBLEM.md; no new external theorem is claimed by this review.
