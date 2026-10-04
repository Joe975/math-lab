# Independent skeptical verification

Mode: informed. Scope: local certificates and the aligned two-square reduction;
no global enumeration or novelty claim.

## Direct checker review

Independently differentiated each pair contribution. For d=q_i-q_j,
D(d/|d|^3)=I/|d|^3-3dd^T/|d|^5. The signs in the two diagonal
and two off-diagonal blocks of `central.py` agree with f=q-force.
The rational square-root bracket is outward: floor(sqrt(floor(u S^2)))/S
is a lower bound, and adding 1/S gives a strict upper bound. Interval
division correctly reverses reciprocal endpoints and rejects zero crossings.
Collision rejection is conservative and safe. I found no soundness defect.

For T(x)=x-Cf(x), the row bounds enclose DT throughout the convex box.
Strict row-sum contraction and strict image inclusion imply a unique fixed
point by Banach. The contraction bound also makes C J(x) invertible at every
box point (Neumann series), so square C is invertible. Thus a fixed point
really is a zero of f; no independent determinant check is needed.

The omitted y equation is recovered from exact pairwise torque cancellation:
sum_i (x_i f_iy-y_i f_ix)=0. At a reduced zero with y0=0 this is x0*f0y=0.
The anchor interval excludes zero, so the omitted equation vanishes. Summing
all equations then gives center of mass zero. This removes rotational freedom
only locally: the certificate does not prove global uniqueness, orientation
uniqueness, or uniqueness up to relabeling. Physical central configurations
have positive lambda, so scaling to lambda=1 is legitimate.

Distance order-statistic separation is sound even when individual intervals
overlap: sort lower and upper endpoint lists independently. It proves only
inequivalence, never equivalence. Callers of `distinct_by_distances` must first
verify the input certificates; that helper itself does not authenticate them.

## Independent implementation and attacks

`verify_polynomial.py` uses different equations, polynomial lifting, exact
Fraction arithmetic, and an independently written sparse interval Jacobian.
Positive variables s_ij obey s_ij^2 |q_i-q_j|^2=1 and the force uses s_ij^3.
This establishes positive physical inverse distances and excludes collisions.
Thirteen deterministic tests pass, including a lambda-normalized eight-gon,
negative inverse distance, zero anchor, perturbed center, zero preconditioner,
collision, wrong mass, malformed dimension, and zero radius. Exact rational
secants test the Jacobian enclosure. Verification of the eight-gon succeeds
while imports of numpy/scipy/mpmath are explicitly prohibited.

Reproduce:
`C:/Repos/math-lab/.venv/Scripts/python.exe -m pytest C:/Repos/math-lab/tests/test_cc8_polynomial.py -q -p no:cacheprovider`

## Aligned distinct concentric squares

Independent derivation: put the inner square at (±1,0),(0,±1), outer at
(±t,0),(0,±t), t>1. A single radius-r square has radial coefficient
c4/r^3, where c4=1/4+1/sqrt(2)=(1+2sqrt(2))/4. Cross-square perpendicular
vertices contribute the same 2/(1+t^2)^(3/2) coefficient at both radii,
so they cancel on subtraction. The remaining coefficient difference is

    D(t)=c4(1-t^-3)-(1+1/t)/(t-1)^2+(1-1/t)/(t+1)^2.

Thus D(t)=0 is equivalent, with all divisions positive, to

    c4=H(t)=2t^2(3t^2+1)/[(t^3-1)(t^2-1)^2].

For numerator N and denominator P,

    N'/N=2/t+6t/(3t^2+1)<4/t,
    P'/P=3t^2/(t^3-1)+4t/(t^2-1)>7/t.

Consequently H'/H<-3/t<0. H tends to infinity as t decreases to 1 and
to zero as t tends to infinity. Exactly one ratio solves the equation.
All outer-square force terms have positive radial projection, so the common
lambda is positive, allowing a unique scaling to lambda=1. Tangential
components cancel by reflection. The case t=1 is a collision and excluded;
t<1 is exchanged radii. This proves existence and uniqueness in precisely
the aligned, distinct-radius, concentric two-square family, modulo scale,
orthogonal transformations, and relabeling. It does not classify general D4
configurations, twisted squares, or all eight-body configurations. No gap
was found in this restricted-family proof. Novelty requires a separate
literature check; this note makes no novelty assertion.

## Full saved-certificate audit

Independently reran the polynomial verifier on all 39 saved polynomial
certificates: 19 discovery certificates and 20 published-catalogue certificates.
Every certificate has n=8 and passed. Separately implemented exact squared
distance order-statistic enclosures from the polynomial coordinate boxes,
without importing `central.py`. All 190 unordered pairs of the 20 catalogue
boxes have at least one strictly separated order statistic. Exact positive
separation margins and certificate paths are in `verification-summary.json`.
These establish 20 inequivalent local configurations, not completeness.

Reproduce the independent audit:
`C:/Repos/math-lab/.venv/Scripts/python.exe C:/Repos/math-lab/problems/central-configurations-8/explore/audit_certificates.py`

The independent checker/audit test suite now has 14 passing tests, including
an exact extremal-box test of the independently implemented distance bounds.

Reviewed `search.py` and `continuation.py`. Both explicitly describe numerical
discovery; continuation reports `complete: false`, and equal-mass endpoints
must pass both certificate builders before being saved as new objects.
Sampled singular values, SVD tangents, and predictor-corrector steps are not
treated as certified continuation. Numerical graph isomorphism is explicitly
labeled discovery deduplication, not exact equivalence. No false proof claim
was found. Search can miss branches, jump branches, or merge close numerical
candidates; these affect coverage only and do not invalidate independently
certified endpoints. The claimed matching of 19 discoveries to a catalogue
remains numerical evidence, not an exact equivalence theorem.
