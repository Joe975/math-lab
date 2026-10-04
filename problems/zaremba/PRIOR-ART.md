# Zaremba — prior art

> **Tier 1.** What this lab has tried. Index in `prior-art.json`; full records
> in `attempts/`.

## 001 — baseline z(n) census (2026-07-27, informed, VERIFIED)

The problem was added and immediately given its baseline. Two halves, per the
PROBLEM.md contract:

- **Witness scan.** z(n) ≤ 5 for every n ≤ 10^7 with a stored certificate m
  per n, every certificate independently re-verified by the stdlib
  implementation; extended to n ≤ 10^8 as a single-implementation check with
  certificates discarded. Zero failures anywhere.
- **Exact census.** z(n) computed exactly for n ≤ 10^5 by exhaustive scan,
  and re-computed by a structurally different algorithm (forward generation
  over continuants, no Euclid, no gcd) — the two 99,999-line tables are
  identical.

What it found: the exceptional set {n : z(n) ≥ 4} has exactly 38 members
below 10^5 — three with z = 5 (6, 54, 150) and 35 with z = 4 — and its
largest member is 6234, leaving the remaining ~94% of the range exceptional-
free. 37 of the 38 are even (the odd one: 1155 = 3·5·7·11). First witnesses
concentrate tightly at the left edge of the K ≤ 5 Cantor set, (3√5−5)/10;
exploiting that edge is what makes the scan fast (and missing it is a ~60×
slowdown, recorded in the attempt as a methodological note).

Standing speculation (labelled in the record): {z ≥ 4} is finite. Leads: the
exact census to 10^6, a bitmap mark-all scan to 10^9, odd-n restriction,
z = 2 residue structure vs the corrected local–global conjecture, and the
divisibility anatomy of scan-hard n.
