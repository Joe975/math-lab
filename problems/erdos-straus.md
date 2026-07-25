# Erdős–Straus Conjecture

**Statement.** For every integer n ≥ 2 there exist positive integers x, y, z
with 4/n = 1/x + 1/y + 1/z.

**Status.** Open. Verified computationally to very large bounds (well beyond
10^17). It suffices to prove it for primes. Known modular identities settle
all n outside certain residue classes; the hard cases are primes in classes
like n ≡ 1 (mod 24) — specifically squares mod 840 obstruct the classical
identity constructions (Mordell's analysis).

**Attack surface.**
- Extend/complete the table of polynomial identities by residue class; map
  exactly which classes remain uncovered and why (obstruction structure).
- Study the number of representations f(n) statistically — Elsholtz–Tao give
  average bounds; look for structure in primes with few representations.
- Computational: find primes with minimal representation counts and mine
  them for common structure.
- Dead-end watch: pure identity approaches provably cannot cover all n with
  finitely many polynomial identities (quadratic-residue obstruction) — any
  new identity scheme must explain how it evades this.
