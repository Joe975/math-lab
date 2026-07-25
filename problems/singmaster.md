# Singmaster's Conjecture

**Statement.** There is a finite bound N such that no entry > 1 appears more
than N times in Pascal's triangle. (Empirically N = 8 may suffice; 3003
appears 8 times and no entry is known to appear more.)

**Status.** Open. Best general bound is O((log n · log log log n)/(log log n)^3)
appearances (Kane). Matomäki–Radziwiłł–Shao–Tao–Teräväinen (2022) proved the
conjecture holds in the "interior" of Pascal's triangle. Infinitely many
entries appear ≥ 6 times (Fibonacci-parameterized family from the Pell-like
equation for C(n, k) = C(n−1, k+1)).

**Attack surface.**
- Computational: extend the search for entries with multiplicity ≥ 8 (needs
  smart algorithms — values are astronomically large; search in (n, k) space
  via collision-finding among binomials with small k).
- Study the Diophantine equations C(n, j) = C(m, k) for small fixed j < k;
  each pair is a curve — which are resolved, which are open? Build the table.
- Boundary vs interior: the MRST-T result leaves near-edge cases; map exactly
  what parameter region remains open.

**Realism.** The multiplicity-8 collision search is a well-defined
computational project with a real (if unlikely) discovery payoff, and the
equation table is solid library content.
