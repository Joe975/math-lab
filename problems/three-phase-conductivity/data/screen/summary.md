# Bounded-rank laminate screen: summary

Generated from every `data/screen/*.log` by the aggregation step recorded in
attempt 001. Full machine-readable table: `summary.csv`.

All values are FLOAT SCREEN results (seed 1, 2 restarts, Nelder-Mead, kept only
when anisotropy and fraction error are both < 1e-6). They are candidates; the
exact values are in `data/certified/`. `milton` is whether the classical
coated assemblage attains the bound at that point (lower side: HS_lo <= sigma2).

## sigma = (1,2,5), lower side

| f (x/8) | milton | HS_lo | rank 3 | gap | rank 4 (axes) | gap |
|---|---|---|---|---|---|---|
| 1,1,6 | no | 3.3636360 | 3.3848002 | 2.12e-02 | 3.3642207 | 5.85e-04 |
| 1,2,5 | no | 3.0000000 | 3.0252805 | 2.53e-02 | 3.0004589 | 4.59e-04 |
| 1,3,4 | no | 2.6923080 | 2.7130975 | 2.08e-02 | 2.6926380 | 3.30e-04 |
| 1,4,3 | no | 2.4285710 | 2.4413911 | 1.28e-02 | - | - |
| 1,5,2 | no | 2.2000000 | 2.2046898 | 4.69e-03 | - | - |
| 2,1,5 | no | 2.6923080 | 2.6978163 | 5.51e-03 | 2.6978163 | 5.51e-03 |
| 2,2,4 | no | 2.4285710 | 2.4333680 | 4.80e-03 | 2.4333680 | 4.80e-03 |
| 2,3,3 | no | 2.2000000 | 2.2018769 | 1.88e-03 | - | - |
| 3,1,4 | no | 2.2000000 | 2.2004689 | 4.69e-04 | - | - |
| 3,2,3 | yes | 2.0000000 | 2.0000000 | 0 (exact) | - | - |
| 6,1,1 | yes | 1.2857140 | 1.2938023 | 8.09e-03 | - | - |

## sigma = (1,2,5), upper side (gap = HS_hi - best)

| f (x/8) | milton | HS_hi | rank 3 | gap |
|---|---|---|---|---|
| 1,1,6 | yes | 3.7958120 | 3.7758347 | 2.00e-02 |
| 2,5,1 | no | 1.9709540 | 1.9708071 | 1.47e-04 |
| 3,4,1 | no | 1.8292680 | 1.8263697 | 2.90e-03 |
| 4,2,2 | no | 1.9421490 | 1.9420303 | 1.19e-04 |
| 4,3,1 | no | 1.6932270 | 1.6876217 | 5.61e-03 |
| 5,1,2 | no | 1.8016190 | 1.8010245 | 5.95e-04 |
| 5,2,1 | no | 1.5625000 | 1.5560484 | 6.45e-03 |
| 6,1,1 | no | 1.4367820 | 1.4320731 | 4.71e-03 |

## sigma = (1,3,20), lower side, rank 3

| f (x/8) | HS_lo | rank 3 | gap |
|---|---|---|---|
| 1,1,6 | 6.7241380 | 6.8862415 | 1.62e-01 |
| 1,2,5 | 5.4615380 | 5.6206731 | 1.59e-01 |
| 1,3,4 | 4.5537190 | 4.6630588 | 1.09e-01 |
| 1,4,3 | 3.8695650 | 3.9249190 | 5.54e-02 |
| 1,5,2 | 3.3354840 | 3.3494575 | 1.40e-02 |
| 2,1,5 | 4.3760000 | 4.3980377 | 2.20e-02 |
| 2,2,4 | 3.7323940 | 3.7467204 | 1.43e-02 |
| 2,3,3 | 3.2264150 | 3.2289042 | 2.49e-03 |
| 3,1,4 | 3.1226990 | 3.1228779 | 1.79e-04 |

## Observations (facts, then labelled interpretation)

1. **Controls behave.** At f = (3,2,3)/8, where HS_lo = sigma2 = 2 exactly
   (Milton's boundary), the rank-3 screen returns exactly 2 and the exact
   certification confirms an exact hit. So the pipeline can reach HS when a
   low-rank structure exists.
2. **Milton-attainable does not mean rank-3-attainable.** At f = (6,1,1)/8
   lower the assemblage attains HS_lo but rank 3 stalls 8.1e-3 above it, and
   on the upper side at (1,1,6)/8 rank 3 stalls 2.0e-2 below HS_hi. The
   assemblage is a higher-rank object.
3. **The gap shrinks monotonically as f1 grows toward Milton's threshold**
   along every ray tried, on both sides and both conductivity triples.
4. **Larger conductivity contrast, larger gap.** At sigma = (1,3,20) the
   rank-3 gaps are roughly an order of magnitude bigger than at (1,2,5) for
   the same fractions.
5. **Rank 4 helps at f1 = 1/8 but returned no improvement at f1 = 2/8.**
   At (1,1,6), (1,2,5), (1,3,4) rank 4 cuts the gap by 40-60x; at (2,1,5) and
   (2,2,4) it reproduces the rank-3 value to every digit printed. Since the
   rank-3 winner's normals are axis normals, the rank-4 axis search space
   contains it, so this is a failure to improve, not a contradiction.
   SPECULATION: this is more likely two-restart optimiser non-convergence on a
   harder landscape than a real feature of the family; the falsifiable test is
   a rank-4 rerun at (2,2,4)/8 with more restarts and a different seed.
