# Hadwiger–Nelson Problem

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** How many colours are needed to colour every point of the plane
so that no two points at distance exactly 1 share a colour? The value χ(R²) is
unknown.

By De Bruijn–Erdős compactness (which uses the axiom of choice) χ(R²) equals
the largest chromatic number of a *finite* unit-distance graph — a finite set
of points in the plane, with an edge between two points exactly when they are
at distance 1. So the question is finitary: exhibit a finite unit-distance
graph that needs k colours, and χ(R²) ≥ k follows.

## Published status

- 4 ≤ χ(R²) ≤ 7 was the classical range: the Moser spindle forces 4, and
  Isbell's 7-colouring of a hexagonal tiling of diameter slightly under 1
  gives the upper bound.
- de Grey (2018) exhibited a 1581-vertex non-4-colourable unit-distance graph,
  raising the lower bound to 5. So χ(R²) ∈ {5, 6, 7}, and all three remain
  possible.
- Polymath16 (2018–2021) minimised the smallest *known* 5-chromatic
  unit-distance graph. The standing record is 509 vertices and 2442 edges, due
  to Jaan Parts (arXiv:2010.12665). The project's final thread describes the
  effort as concluded, not the problem.
- Whether a 6-chromatic unit-distance graph exists is open. No approach is
  known that would push the lower bound to 6, and no approach is known that
  would bring the upper bound below 7.
- Variants behave differently and are not substitutes: colourings with
  measurable colour classes are known to need at least 5, colourings whose
  classes are bounded by Jordan curves at least 6, and the answer for some
  formulations depends on the set-theoretic axioms assumed (Shelah–Soifer).

**Sources.** The figures above were assembled from encyclopaedia summaries and
paper abstracts rather than the primary papers, so treat them as machine
transcribed `[T]` and check any one you intend to lean on. The vertex counts
(1581, 509) and the bound χ ∈ {5,6,7} are the load-bearing ones.

## Verification contract

The objects here are geometric, and floating point silently destroys the only
property that matters. A claim about this problem must therefore do the
following to be checkable.

- **Coordinates are exact.** Vertices must lie in an explicitly named real
  number field (Q(√3, √11) and its relatives are the usual homes) with exact
  arithmetic in that field. A distance is 1 or it is not; a distance that is
  1.0000000001 is not an edge, and a search that treats it as one has proved
  nothing.
- **Unit-distance claims certify both directions.** Every declared edge has
  squared distance exactly 1, and every declared non-edge does not. Vertices
  claimed distinct must be certified distinct in the field, not by tolerance.
- **Non-k-colourability is a refutation and needs an independent certificate.**
  Ship a DRAT/LRAT proof checked by a checker that did not produce it, or an
  exhaustive search reproduced by a second, independently written
  implementation. A solver reporting UNSAT is not by itself a result.
- **k-colourability ships the colouring.** Verification is linear and must be
  run by a checker that did not construct it.
- **A vertex count is EVIDENCE about a graph, never about the plane.** "5
  colours are forced by this 509-vertex graph" bounds χ(R²) below by 5 and says
  nothing else. State the graph.

## Harness (tier 0)

None yet. A contributor adding one should put the exact-field arithmetic, the
unit-distance certifier and the colouring checker here, since all three verify
the objects themselves rather than any particular route.
