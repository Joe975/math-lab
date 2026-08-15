You are writing an independent verification program. Work from the definitions
below only. Do not look for or assume any published result about this object.

## Definitions

Let T be the map on positive integers with T(n) = n/2 when n is even and
T(n) = 3n + 1 when n is odd.

Fix an integer k >= 2 and set W = 2^(k+1).

Define a finite directed graph H_k:

- Vertices: the integers 1, 2, ..., W, plus one extra vertex called OUT.
- Edges: every vertex v in {1, ..., W} has exactly one outgoing edge.
  It goes to T(v) if T(v) <= W, and to OUT otherwise.
  OUT has no outgoing edge.

Let A = {1, 2, ..., 2^k} (a subset of the vertices of H_k).

Define S_k as the MINIMUM NUMBER OF EDGES that must be deleted from H_k so
that afterwards there is no directed path from any vertex of A to OUT.

## Your task

Write a single self-contained Python 3 program, standard library only, that
computes |S_k| for every k from 2 to 14 inclusive and prints one line per k
in exactly this format, nothing else on the line:

k=<k> S=<value>

Approach you must use: {{value}}

## Requirements

- The program must be correct rather than fast; k = 14 is small.
- Do not hard-code any table of answers. The program must actually compute.
- Print nothing but the 13 result lines.
- Put the complete program in a single ```python fenced block. Before the
  block, state in at most 5 sentences why your method computes a genuine
  minimum, not merely some cut.
