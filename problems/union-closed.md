# Union-Closed Sets Conjecture (Frankl)

**Statement.** If F is a finite union-closed family of sets with F ≠ {∅},
then some element belongs to at least half the sets in F.

**Status.** Open at 1/2. Major recent progress: Gilmer (2022) proved a
constant fraction (~0.01) via an entropy/information-theoretic argument;
follow-ups (Alweiss–Huang–Sellke, Chase–Lovett, Sawin, and others) pushed the
constant to ≈ 0.38 (the (3−√5)/2 barrier). It is known the pure Gilmer-style
argument cannot pass (3−√5)/2 without new ideas (Chase–Lovett constructed
approximate counterexamples to the strengthened entropy statement).

**Attack surface.**
- Understand precisely why (3−√5)/2 is the barrier for the entropy method and
  what structural property of approximate-union-closed families breaks it;
  look for an added constraint (exact union-closure) the argument discards.
- Computational: exhaustive/randomized search over small ground sets for
  families minimizing max element frequency; characterize extremal families.
- Lattice reformulation (union-closed families ↔ join-semilattices); test
  whether known equivalent forms (graph version, lattice version) admit
  sharper small-case analysis.

**Realism.** Active frontier with live technique — best "serious problem with
genuine traction" candidate in the set.
