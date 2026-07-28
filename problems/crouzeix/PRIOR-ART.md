# Crouzeix — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`.

## Attempts

**None.** This problem was onboarded and has not been worked. There is no prior
art to be informed by, so `blind` and `informed` mode are currently equivalent
here — which makes the first attempts worth running blind, since blind costs
nothing while the record is empty.

## Editorial view of the attack surface

Operator theory is a new attack surface for the lab, and this problem is the
riskiest of the three onboarded with it. The risk is not mathematical, it is
methodological: the quantities in the conjecture are eigenvalue and
singular-value problems, so nothing here is exactly rational, and every claim
has to be carried by certified enclosures written from scratch in the standard
library. Budget accordingly — medium, not high — and scope the certification
machinery before designing any census.

What makes it worth having anyway is that a single number decides everything.
The Crouzeix ratio exceeding 2 anywhere refutes the conjecture, so a search
over matrices and polynomials is self-policing in the same way the volume
product is, and the published numerical work gives known answers to validate
tooling against.

Concrete lines, if you want them:

- Self-test: reproduce the 2×2 theorem numerically and the ratios of the known
  near-extremal families. Its purpose is to validate the certification
  machinery against answers that are already known — treat a mismatch as a bug
  in our enclosures, not as a finding.
- Local-maxima census in dimension 3, polynomial degree ≤ 3: do all basins
  terminate at known extremal structure? `EVIDENCE`, scoped by dimension,
  degree and search design, all three of which must appear in the record.

Kill condition, inline: if certified enclosures cost more than roughly 10× the
floating-point computation, the wide census is dead on arrival. Fall back to
structured families where the norm has a closed form and say so in the record —
the fallback is a legitimate outcome, and quietly shrinking the scope while
keeping the original framing is not.

One more caution specific to this problem. The published numerical work already
reports that optimization flows to the known near-extremizers. Reproducing that
is not a result; the value of a census here is in the search design being
independent of theirs, so state precisely how ours differs before running it,
not afterwards.
