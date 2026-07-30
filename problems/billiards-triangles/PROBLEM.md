# Periodic Orbits in Triangular Billiards

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** Does every triangle admit a periodic billiard orbit? A billiard
orbit is the path of a point mass moving in straight lines inside the triangle
and reflecting off each side with angle of incidence equal to angle of
reflection. It is periodic if it returns to its initial position *and*
direction after finitely many reflections. Trajectories that hit a vertex are
undefined and do not count.

Only the angles matter: the question is about a point (α, β) in the parameter
triangle {α, β > 0, α + β < π}, since similarity does not change the dynamics.

## Published status

Open in general for irrational-angled obtuse triangles; settled in every case
below, and claimed in full by a 2026 preprint (see the last paragraph of this
section).

- **Acute triangles.** Fagnano (1775): the orthic triangle — the triangle
  joining the three altitude feet — is a periodic orbit.
- **Right triangles.** Holt (1993); the same result appeared independently in
  Gal'perin–Stepin–Vorobets (1991).
- **Isoceles triangles.** Known; surveyed in Tokarsky–Garber–Marinov–Moore
  (below).
- **Rational-angled triangles** (all angles rational multiples of π). Masur
  (1986): infinitely many periodic orbits of distinct periods.
  Boshernitzan–Galperin–Krüger–Troubetzkoy, *Periodic billiard orbits are dense
  in rational polygons*, Trans. AMS **350** (1998) 3523–3535.
- **Obtuse angle ≤ 100°.** R. E. Schwartz, *Obtuse Triangular Billiards II: One
  Hundred Degrees Worth of Periodic Trajectories*, Experimental Mathematics
  **18** (2009) 137–171. Computer-assisted, using the McBilliards program of
  Schwartz and Hooper, and carrying a stability conclusion.
- **Obtuse angle ≤ 112.3°.** J. Garber, B. Marinov, K. Moore, G. Tokarsky,
  *One Hundred and Twelve Point Three Degree Theorem*, arXiv:1808.06667.
  Computer-assisted, without the stability conclusion of Schwartz's result.
  The authors state that the method does not reach past roughly 112.5° and
  that a new idea is needed beyond that.

A. Katok listed the polygonal periodic-orbit problem among his "Five Most
Resistant Problems in Dynamics".

**A 2026 preprint claims the general case.** G. Forni, *Existence of a Periodic
Orbit for Billiards in Polygons*, arXiv:2606.10102 (submitted 8 June 2026),
claims that the billiard flow in any finite polygon has at least one periodic
orbit, which subsumes the triangle question. The argument is by contradiction,
built on the billiard-flow dynamics of Galperin–Krüger–Troubetzkoy, a
one-parameter scaling of the natural Riemannian metric on the unit tangent
bundle, and the topology of the cut locus of those metrics. Treat existence as
likely settled pending peer review. Note what the argument does not do: being
non-constructive, it exhibits no orbit, bounce word or period for any specific
triangle, so the **constructive** frontier is untouched and remains the
112.3° theorem — for an irrational-angled obtuse triangle with largest angle
above 112.3°, no explicit periodic orbit is known.

### The unfolding framework (published)

Rather than reflecting the trajectory in a side, reflect the triangle: a
candidate orbit becomes a straight segment crossing a chain of reflected
copies. A **bounce word** — the sequence of sides struck — therefore determines
a linear-algebraic condition on (α, β), and the condition holds on an open set:
a word that yields a periodic orbit at one triangle yields one throughout a
neighbourhood. Each word thus certifies an open **orbit tile** of the parameter
triangle, and the conjecture becomes a covering problem: does some family of
words tile the whole obtuse region? This is the framework McBilliards is built
on, and the mechanism behind both the 100° and 112.3° theorems. Schwartz
reports that deep searches reveal infinite patterns of tiles that stop short of
covering a neighbourhood of certain parameter points, the (π/6, π/3) triangle
among them.

## Verification contract

Any claim recorded against this problem must meet the bar in `CONTRIBUTING.md`.
Specifically here:

- A claim that a word **W** certifies a region **R** must state R exactly —
  rational vertices, or explicit defining inequalities — and must be
  established in exact arithmetic. A floating-point evaluation of an unfolding
  proposes a candidate; it never certifies one.
- A certificate must rule out degeneracy explicitly: an orbit through a vertex
  is not a periodic orbit, and neither is a segment that leaves the unfolded
  corridor.
- A **coverage** claim must state the parameter region covered, the full
  certificate list, and the arithmetic used. Partial coverage is `EVIDENCE`
  about the covered region only — never about the conjecture.
- A claimed orbit for a specific triangle must be re-derived by an independent
  method (direct simulation under certified enclosures, not a second run of the
  unfolding code) before it is recorded.
- Any statement about where a word search stalls is `EVIDENCE` bounded by the
  word length and the search design, both of which must be stated.

## Harness (tier 0)

- `harness/billiards-triangles/unfold.py` — the reference implementation.
  Decides whether a bounce word certifies an open region of triangles, working
  in apex coordinates (A = (0,0), B = (1,0), C in a rational box) so that no
  trigonometry enters the certificate path, and in first-order affine forms
  over ℚ so that a verdict covers every triangle in the box at once. Verdicts
  are TRUE, FALSE or UNKNOWN; UNKNOWN means subdivide, never "probably". The
  translation condition is decided in integer arithmetic on the word alone.
  Note the scope: an open-region certificate can only ever find orbits that
  survive perturbation, so an unstable orbit is invisible to it by
  construction.
- `harness/billiards-triangles/verify_cover.py` — independent re-verification
  by direct billiard simulation. At a rational apex a trajectory is exactly
  rational, so it re-derives any accepted certificate with no enclosures at
  all, checking the bounce sequence, the closure, and that no bounce lands on
  a vertex. It does not unfold, which is the point: the corridor criterion is
  about the set of gates rather than their order, and only simulation can
  catch that. Run it on any certificate you intend to claim.
