# Lonely Runner — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`. Full records: `attempts/`.

## Editorial view of the attack surface

- Verify k = 8 for structured/small speed sets; look for near-violating speed
  tuples and study their arithmetic structure (they tend to be
  arithmetic-progression-like).
- Test tightness rigidity computationally for k = 8, 9.
- View-obstruction / zonotope reformulation: covering radius of certain
  polytopes; try LP/SDP relaxations for lower bounds on the gap.

Assessed as: full k = 8 is likely out of computational reach, but structure
mining of near-tight examples is tractable and could support a partial result.

## Attempts

### 001 — k = 8 near-tight census to V = 72 · `VERIFIED` (range only)

Among all 1,473,109,704 speed 7-tuples with max speed ≤ 72, exactly 3
primitives have ML < 13/100, all with ML = 1/8 exactly:

| primitive tuple | structure |
|---|---|
| (1,2,3,4,5,6,7) | canonical tight instance |
| (1,2,3,4,5,7,12) | Goddyn–Wong acceleration |
| (1,4,5,6,7,11,13) | Goddyn–Wong sporadic |

Both Goddyn–Wong instances were **recovered from scratch**, which is the
run's main correctness check. The 21 raw hits are precisely the scaled copies,
each independently computed with identical ML.

- **No tuple with ML < 1/8** among all 1.47·10⁹ — consistent with the
  conjecture and with Rosenfeld's k = 8 preprint.
- **No tuple with 1/8 < ML < 13/100**: the spectrum has an empirical gap
  above 1/8 for all speeds ≤ 72.
- The ML spectrum below 1/7 (V ≤ 40) is exactly {s/(7s+k), k ∈ {1,2}} —
  consistent with Fan–Sun's amended spectrum conjecture; a new data point
  at n = 7.

Exact rational arithmetic throughout, independently re-verified.

## Open lines

- k = 9 near-tight scan (threshold near 1/9). Do the feasibility analysis
  first — 8-tuples grow fast. Consider restricting to accelerations/near-APs
  of known structures plus a bounded full scan.
- k = 8 is likely settled by Rosenfeld, so further k = 8 search is low value.
