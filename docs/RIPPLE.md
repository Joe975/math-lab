# Ripple scan: one result, many problems

The second cross-pollination flow: take **one new result** and scan the other
problems for where it bites. The inverse of `docs/IDEATE.md`.

A "result" here is either:

- **external** — a paper or preprint touching one of our problems or one of
  our mechanism families (the way a preprint can settle a frontier the queue
  still lists as open); or
- **internal** — a route in this lab reaching `LIVE`, `VERIFIED`, or a
  generalized no-go. A strong *negative* result ripples too: a no-go that
  covers a family of functionals may foreclose queued routes on other
  problems, which is worth knowing before a cycle is spent on them.

Run a scan when such a result appears, and as a standing check during the
orient step of a cycle (`docs/CYCLE.md`) when something new has landed since
the last one. Informed-side only, like everything in this layer.

## The procedure

### 1. Characterize the result

Write three lines before scanning:

- **What it establishes**, with its actual scope (range, family, hypothesis).
  For external results: verify the source says what it seems to say —
  transcription errors at this step poison every downstream conclusion, and
  the scan multiplies the blast radius.
- **The mechanism**, as tags — existing ones from `mechanisms.json` where
  they fit, and what field lens it lives in.
- **What it needs to apply**: the structural features of a problem that the
  argument actually uses (a monotone structure, a residue-class split, an
  expansion property...). This list is what the scan matches against.

### 2. Scan the portfolio

For each other problem — all of them; the scan is cheap and the point is
reaching problems nobody was looking at:

```bash
python scripts/mechanisms.py where <tag>      # has this family touched it?
```

then its `prior-art.json` (`route_summary`, `gaps`, dead-end statuses). Ask
three questions, in order of value:

1. **Does it unblock a recorded gap?** Labeled proof gaps and `LIVE`-route
   obstructions are the highest-value targets — they are precise statements
   of what is missing, so matching against them is almost mechanical.
2. **Does it foreclose a queued route?** If the result is a no-go whose scope
   covers a queued lead, the queue item should die now, not after a cycle.
3. **Does it open an untried lens?** If the mechanism's field shows `·` for
   that problem in `mechanisms.py matrix`, the result may be the seed for an
   ideation sweep there.

A miss is recorded in one line; a hit gets a paragraph: which gap or queue
item, why the result's hypotheses plausibly hold there, and the first
concrete step. `SPECULATION` labels apply to transfer hypotheses like
anything else.

### 3. Record and route the hits

- **Small scan, no hits:** a dated note in the `STATUS.md` insights section
  ("scanned <result> across portfolio; no purchase, because <reason>") so the
  next cycle does not redo it.
- **Hits:** queue items in `STATUS.md` in the standard lead format, citing
  the result. A hit on a problem's recorded gap goes at the *top* of that
  problem's queue — unblocking a labeled gap outranks new exploration.
- **A scan that changed the picture** (killed a queue item, or produced a
  serious transfer route): record it as a `MAP` attempt on the problem most
  affected, so the reasoning survives with the usual permanence.

External results enter tier-1 prose only (`PRIOR-ART.md`, attempt records,
`STATUS.md`) unless they are established literature suitable for a
`PROBLEM.md` — and editing tier 0 follows the usual rule that nothing of ours
goes there.

## Calibration

Most scans find nothing, and a recorded nothing is the deliverable — it is
what makes the next scan cheap. The failure mode is the opposite one:
analogy inflation, where a result "morally" applies everywhere and generates
five vague queue items. The test is the hypothesis list from step 1: a hit
must name which recorded gap it meets and why the hypotheses hold, or it is
not a hit.
