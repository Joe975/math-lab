# Eight equal masses: planar central configurations

Classify the collision-free configurations of eight equal positive point masses
in the Newtonian planar n-body problem, up to similarity and permutation.
Collinear configurations are included in the planar problem. A central
configuration supports homothetic motion and, in the plane, rigid rotation.

With masses and gravitational constant set to one, choose the scale lambda=1:

    q_i = sum_{j != i} (q_i-q_j)/|q_i-q_j|^3,    q_i in R^2.

The equations imply center of mass zero. In moment-of-inertia normalization,
divide positions by sqrt(sum_i |q_i|^2); lambda then changes, so do not impose
both scale normalizations simultaneously.

## Published frontier

Moczurad and Zgliczynski (2019) give complete computer-assisted classifications
for equal planar masses through n=7 and certify asymmetric examples at n=8,9,10.
They explain that global exclusion, rather than local root certification, is the
obstacle to extending their complete enumeration. Their Theorems 6 and 11 give
explicit separation and radius bounds. Doicu, Zhao and Doicu (2022) provide a
larger numerical catalogue with local checks, explicitly without a completeness
guarantee. A complete equal-mass planar n=8 classification was not found in the
literature check dated 2026-09-06.

## What would settle a claim

This is a classification problem, not a conjectured numeric bound. A claimed
complete list is refuted by one certified configuration proved inequivalent to
every member, or by a certified positive-dimensional family outside that list.
A finite list of verified roots alone cannot establish completeness.

An individual certificate must specify exact rational isolating boxes, the
normalization and gauge, collision exclusion, a rigorous existence/uniqueness
test, and why the reduced equations imply the full force equations. A small
floating-point residual is insufficient. Independent checking must use a
separately derived implementation. Certification concerns the nearby exact
root, not the rounded center of its box.

Counts modulo symmetry require a sound equivalence argument. Unequal distance
multisets prove inequivalence; matching sorted distances alone does not prove
equivalence. Full coverage requires a justified compact domain and a verifiable
partition into exclusions and isolating regions, including all gauge charts
and collision boundary cases. Report every unclassified region explicitly.

## References

- Moczurad and Zgliczynski, *Central configurations in planar n-body problem
  with equal masses for n=5,6,7* (2019),
  https://doi.org/10.1007/s10569-019-9920-6 .
- Doicu, Zhao and Doicu, *A stochastic optimization algorithm for analyzing
  planar central and balanced configurations in the n-body problem* (2022),
  https://doi.org/10.1007/s10569-022-10075-7 .
- Published numerical coordinates:
  https://github.com/AlexandruDoicu/Balanced-and-Central-Configurations .
