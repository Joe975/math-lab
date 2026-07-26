# Erdős–Straus Conjecture

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** For every integer n ≥ 2 there exist positive integers x, y, z
with 4/n = 1/x + 1/y + 1/z.

## Published status

Open. Verified computationally to very large bounds (well beyond 10^17). It
suffices to prove it for primes. Known modular identities settle all n outside
certain residue classes; the hard cases are primes in classes like
n ≡ 1 (mod 24) — specifically, squares mod 840 obstruct the classical identity
constructions (Mordell's analysis). Elsholtz–Tao give average bounds on the
number of representations f(n).

Known limitation of the classical route: finitely many polynomial identities
provably cannot cover all n, because of the quadratic-residue obstruction. Any
new identity scheme has to explain how it evades this.

## Verification contract

- A claimed **solution** (x, y, z) is checked in exact integer arithmetic:
  `4·x·y·z = n·(x·y + y·z + z·x)`. No floating point.
- A claimed **identity family** must be machine-verified over its whole
  qualifying range, not spot-checked. State the range.
- A claimed **representation count** f(p) must be reproducible from the
  harness with a recorded command, and cross-checked against an independent
  implementation on a subrange.

## Harness (tier 0)

- `harness/erdos-straus/es_fcount.c` — parallel C kernel computing exact
  f(p) = #{x ≤ y ≤ z : 4/p = 1/x + 1/y + 1/z}.
  Build: `cc -O2 -fopenmp -o es_fcount es_fcount.c`.
- `harness/erdos-straus/es_fcount_run.sh` — chunked, checkpointed driver.
