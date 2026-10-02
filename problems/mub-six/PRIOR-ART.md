# MUBs in dimension six — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`. Full records: `attempts/`.

## Editorial view of the attack surface

- Numerical searches for four MU bases are well trodden in the literature;
  a new negative search adds little unless it is calibrated or certified.
- The productive directions are certification (prove the numerical floor is a
  true minimum, ideally a global one within a family) and structured searches
  seeded from the known families of 6×6 complex Hadamard matrices.

## Attempts

### 001 — Multi-start census of the four-basis defect landscape · `EVIDENCE`

22,000 random starts of projected-gradient descent on four bases of C^6: no
start reaches L = 0. 70.7% of the 20,000-start run converge to
L = 0.0512492189962838 (ASD 0.998291692700), which reproduces the published
Raynal–Lü–Englert optimum (rediscovery, not new). At the floor all three
non-standard bases are Hadamard (unbiased to the standard basis). The
residual defect is shared equally by the three Hadamard–Hadamard pairs, and
each row of their overlap matrices is a permutation of {0.124398, 0.151402,
0.18105 ×4}.

Controls under identical settings: d = 2 sits at the proved floor 2/3; known
MU sets are found in d = 3, 4, 5, 7 and for three bases in d = 6. The d = 7
hit rate is only 0.7%, which is the main caveat on reading the d = 6 null.
Leads: certify the floor; closed forms for the three overlap values;
Hadamard-family seeds; equal-sample controls.
