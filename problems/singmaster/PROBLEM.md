# Singmaster's Conjecture

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** There is a finite bound N such that no entry > 1 appears more
than N times in Pascal's triangle. (Empirically N = 8 may suffice: 3003
appears 8 times and no entry is known to appear more often.)

## Published status

Open. The best general bound is O((log n · log log log n)/(log log n)³)
appearances (Kane). Matomäki–Radziwiłł–Shao–Tao–Teräväinen (2022) proved the
conjecture holds in the "interior" of Pascal's triangle, leaving near-edge
cases. Infinitely many entries appear ≥ 6 times, via a Fibonacci-parameterized
family from the Pell-like equation for C(n, k) = C(n−1, k+1).

de Weger published a list of binomial coincidences that any search should
reproduce.

## Verification contract

- **State the multiplicity convention explicitly.** Whether symmetric cells
  C(n,k) and C(n,n−k) count separately, and whether the trivial cells C(V,1),
  C(V,V−1) are included, changes every reported number. Claims using different
  conventions are not comparable.
- A **census claim** ("complete to bound B") requires a stated completeness
  argument for why no qualifying value ≤ B was missed, not just a search log.
- Hits must be re-verified independently, and the census must agree with
  de Weger's coincidence list on the overlapping range.

## Harness (tier 0)

- `harness/singmaster/singmaster_search.py` — canonical-cell enumeration and
  collision search, checkpointed.
