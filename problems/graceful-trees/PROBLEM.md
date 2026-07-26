# Graceful Tree Conjecture (Ringel–Kotzig)

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** Every tree on n vertices admits a graceful labeling: an
injection V → {0, …, n−1} such that the edge labels |f(u) − f(v)| are exactly
{1, …, n−1}.

## Published status

Open in general. Verified computationally for all trees up to at least 35
vertices. Known for many classes (caterpillars, paths, stars, trees of
diameter ≤ 5). Ringel's conjecture itself — which gracefulness would imply —
was proved asymptotically by Montgomery–Pokrovskii–Sudakov (2020), but
gracefulness proper remains open. Lobsters (trees within distance 2 of a path)
are a well-known open subclass.

Anick's database of graceful labeling counts covers trees up to 16 edges.

## Verification contract

- A claimed **labeling** is checked by exact enumeration of its edge-label
  multiset against {1, …, n−1}.
- A claimed **count** must state its convention explicitly: raw labelings, or
  essential counts normalized by |Aut|. The two differ and mixing them silently
  invalidates comparisons.
- Counts overlapping Anick's ≤16-edge database must be reconciled against it.

## Harness (tier 0)

- `harness/graceful-trees/graceful.py` — graceful labeling enumeration and
  counting; the reference implementation.
- `harness/graceful-trees/graceful_core.c` — bitmask backtracking kernel.
  Build: `cc -O2 -o graceful_core graceful_core.c`.
