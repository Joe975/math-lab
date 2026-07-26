# Union-Closed Sets Conjecture (Frankl)

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** If F is a finite union-closed family of sets with F ≠ {∅},
then some element belongs to at least half the sets in F.

## Published status

Open at 1/2. Gilmer (2022) proved a constant fraction (~0.01) via an
entropy/information-theoretic argument; follow-ups (Alweiss–Huang–Sellke,
Chase–Lovett, Sawin, Pebody, Yu, Cambie) pushed the constant to
(3−√5)/2 ≈ 0.381966, and Liu (arXiv:2306.08824) to ≈ 0.38271, the current
published record. Chase–Lovett constructed approximate counterexamples to the
strengthened entropy statement, so it is known that the pure Gilmer-style
argument cannot pass (3−√5)/2 without new ideas.

Primary literature: Gilmer arXiv:2211.09055; Alweiss–Huang–Sellke
arXiv:2211.11731; Chase–Lovett arXiv:2211.11689; Sawin arXiv:2211.11504;
Pebody arXiv:2211.13139; Ellis arXiv:2211.12401; Yu arXiv:2212.00658;
Cambie arXiv:2212.12500; Liu arXiv:2306.08824.

## Verification contract

Any claim recorded against this problem must meet the bar in `CONTRIBUTING.md`.
Specifically here:

- A claimed **bound** must state the threshold p and the class of families or
  distributions it covers, and must be demonstrated against the adversarial
  families in the published literature — at minimum Sawin's Proposition 6
  geometric mixtures and the Chase–Lovett slice family — not merely asserted
  to survive them.
- A claimed **counterexample family** must be verified union-closed by exact
  enumeration, not by construction argument alone.
- Computational evidence is labelled `EVIDENCE`, never `VERIFIED`.

## Harness (tier 0)

- `harness/union-closed/uc_search.py` — exhaustive/randomized search over
  small ground sets for union-closed families minimizing max element
  frequency. Reproduce the small-case baseline with this before claiming
  anything about extremal families.
