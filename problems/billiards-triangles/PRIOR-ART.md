# Triangular billiards — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`.

## Attempts

**None.** This problem was onboarded and has not been worked. There is no prior
art to be informed by, so `blind` and `informed` mode are currently equivalent
here — which makes the first attempts worth running blind, since blind costs
nothing while the record is empty.

## Editorial view of the attack surface

New attack surface for the lab: classical dynamics. The structural fit is the
unfolding framework described in `PROBLEM.md` — a bounce word certifies an open
region of the (α, β) parameter triangle, so the question is a covering problem
over a two-parameter family. That is the same shape as a tuple census: word ↔
tuple, region certificate ↔ an exact quantity computed over ℚ, with the whole
claim standing or falling on exact arithmetic.

Note before spending budget here. Forni's June 2026 preprint claims the general
existence statement, so the *existence* question may be settled. What is not
settled, and what this lab can actually move, is the **constructive** side:
explicit certified orbits past the 112.3° frontier, and a map of where a finite
word census stalls. Frame any attempt against that, not against existence.

Concrete lines, if you want them:

- Self-test: re-derive coverage of the acute, right and isoceles cases with our
  own implementation. Expected `VERIFIED` with scope = the re-derived classes;
  its real purpose is validating the certificate machinery against known
  answers.
- Word census by length L: which short words certify which obtuse regions, and
  where coverage stalls as the largest angle approaches and passes 112.3°. The
  geometry of the uncovered set is a `MAP`-grade deliverable even if no new
  region is covered.

Kill condition to carry into the queue: if the certificate words needed near
the frontier grow super-exponentially in length as the angle increases, the
finite-census route is dead. Record that rather than pushing the search deeper —
a measured growth rate that kills the route is a better outcome than another
week of enumeration.
