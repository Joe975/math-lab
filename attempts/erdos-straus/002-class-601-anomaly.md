# 002 — Is the class-601 anomaly real? f(p) to 10^6 for p ≡ 1 mod 24

**Problem:** Erdős–Straus conjecture; follow-up to attempt 001 (queue item 2).
**Date:** 2026-07-26.
**Tools:** `tools/es_fcount.c` (new parallel C kernel), `tools/es_fcount_run.sh`
(chunked/checkpointed driver), `tools/es_class_stats.py`,
`tools/es_structure_hunt.py`. Data: `attempts/erdos-straus/data/f1mod24_*.txt`.

## Approach

Attempt 001 found (p < 10^5, bottom 50 by raw f among p > 1000) that low-f
primes concentrate in QR-related classes {1,49,73,97} mod 120, but that the
non-QR class **601 mod 840** held 6 of the bottom 50 — unexplained. Plan:
recompute exact f(p) = #{x ≤ y ≤ z : 4/p = 1/x+1/y+1/z} for **all** primes
p ≡ 1 mod 24 up to 10^6 (after spot-checking that other classes mod 24 can be
excluded), run class-enrichment statistics with binomial p-values mod
840/168/120, and structure-hunt whatever classes are genuinely enriched.

A methodological fix over 001: f(p) grows with p, so a "bottom 0.5% by raw f"
set is dominated by the smallest primes. Class comparisons stay fair under raw
selection (all classes are equidistributed in p), but the raw bottom set wastes
most of its mass re-measuring small primes. We therefore also use a
**size-normalized selection**: rank f within dyadic bands of p
(band = bit-length of p), take the bottom of the percentile ranks. Both
selections are reported.

## What was done

1. **New kernel `es_fcount.c`** — same mathematics as 001's counter
   (x ∈ (p/4, 3p/4], A = 4x−p, B = px, count divisors d | B², dmin ≤ d ≤ B,
   A | B+d, with d = d0 or d0·p, d0 | x²), plus OpenMP over primes,
   prime-range chunking with a residue filter, and a prime-list mode.
   **Validation:** output identical to 001's kernel on all 9592 primes < 10^5,
   identical to the independent pure-Python counter for p < 3000, and
   list-mode vs range-mode agree on all shared primes (130/130).
2. **Exclusion check for classes ≢ 1 mod 24** (130 sampled primes per class
   mod 24 in [10^5, 10^6)): class 1 is systematically lowest — its 5th
   percentile of f (86) sits below every other class's sampled *minimum*
   (112–243); medians 170 (r=1) vs 220–537 (others), matching 001's mod-24
   stratification at 10× the scale. (Caveat: full-population minima of other
   classes lie below their n=130 sample minima, but at 10^5 the full-data
   gap was already 23 vs 119; exclusion is safe for bottom-tail statistics.)
3. **Full run:** f(p) for all **9732** primes p ≡ 1 mod 24, p < 10^6, in 10
   checkpointed chunks (`data/f1mod24_<lo>_<hi>.txt`, atomic writes, restart
   skips finished chunks). ~4 min wall on 4 cores.
4. **Statistics** (`es_class_stats.py`): per-class counts in the bottom 0.5%
   (49 primes) and bottom 2% (195 primes), raw and size-normalized, for the
   24 classes mod 840 that are ≡ 1 mod 24, and projections mod 168 and 120,
   each with exact counts, enrichment ratios, and one-sided binomial
   p-values.
5. **Structure hunt** (`es_structure_hunt.py` + within-class follow-up):
   bottom-100 primes vs size-matched random control from the same
   population — factorizations of p±1, QR memberships (mod 5, 7, 840),
   and N_typeI(p; B) = #{(a,b,c) : a ≤ b ≤ B, c | a+b, 4ab | cp+a+b}, the
   per-prime count of applicable small Elsholtz–Tao Type I identities
   (001, lead 1). Class-level: covering-family count (001's ab | 210 sweep)
   vs class mean f.

## Outcome

### 1. The class-601 anomaly is **noise** — it does not persist at 10^6

- **Size-normalized selection:** class 601 holds **0 of the bottom 49**
  (expected 2.07) and **1 of the bottom 195** (expected 8.26; enrichment
  0.12) — if anything it is *under*-represented.
- **Raw selection** (001's style): 601 holds 4/49 (enrichment 1.93,
  binomial p = 0.15) and 12/195 (enrichment 1.45, p = 0.13) — the direction
  of 001's observation replicates but is far from significance, and with 24
  classes tested a p of 0.13 is unremarkable.
- What 001 saw was the tail of a real but modest **class-level** effect:
  601 has the lowest mean f (149.1 on [10^5,10^6)) of the 18 non-QR
  classes, which lets its small primes leak into a raw bottom set. That
  mean is fully in line with its covering count (see below), not anomalous.

### 2. The genuine signal: the bottom tail lives exactly on Mordell's six QR classes

Size-normalized bottom 2% (195 primes): **191/195 (98%)** lie in the six
obstructed classes {1,121,169,289,361,529} mod 840, whose population share is
2370/9732 = 24.4%. Aggregate binomial p ≈ 10^−110 (bottom 0.5%: 48/49,
p ≈ 10^−28). Rate form: P(bottom 2% | QR class) = 8.1% vs
P(bottom 2% | non-QR class) = 0.054% — a **149× rate ratio**. Per-class
(bottom 2%, normalized): 1 (obs 54, exp 7.9, enr 6.8, p = 1.4e−29),
121 (35, 4.3×, 4.6e−13), 169 (33, 4.2×, 2.7e−12), 361 (31, 3.9×, 9.3e−11),
529 (22, 2.7×, 2.1e−5), 289 (16, 2.1×, 4.9e−3). **No non-QR class is
significantly enriched** (best: 241 with 2/195). Projections localize the
signal the same way at every modulus: mod 120 only {1, 49} (QR) are enriched
(124 and 71 of 195; classes 73, 97 hold 0); mod 168 only {1, 25, 121} (QR)
are enriched (87/53/51; classes 73, 97, 145 hold 2/1/1). The signal is
QR structure, period; it is strongest for the full mod-840 refinement, and
class 1 (p ≡ 1 mod both 5 and 7) is the deepest.

### 3. Characterization found (rate form, falsifiable)

**Low-f primes are the primes to which few small Type I identities apply.**
Concretely, at p < 10^6, size-normalized selection:

- Class level: mean f on [10^5, 10^6) is monotone in the number of
  full-class Type I covering families (001 sweep): 0 families (six QR
  classes) → mean f 118–136; 1 family → 149–177; 2 → 167–205; 3 → 203.
  Pearson corr(class mean f, covering count) = **0.92** across the 24
  classes. The six 0-family classes contain 98% of the bottom 2%
  (191/195, p ≈ 10^−110); predicted false-negative rate ~2%.
- Prime level: N_typeI(p; 30) has median 4 on the bottom-100 vs 11 on
  size-matched controls; threshold N ≤ 4 captures 54/100 of the bottom at a
  5/100 false-positive rate. Spearman corr(f, N_typeI) = 0.57 on 800 random
  primes in [5·10^5, 10^6), and stays +0.26…+0.40 *within* fixed classes
  mod 840 — so identity-applicability predicts f beyond class membership.
- The four non-QR interlopers in the bottom 195 fit the same mechanism:
  411841, 652081 (cls 241), 210601 (cls 601), 140761 (cls 481) all have
  N_typeI ≤ 6 (control median 11) and very smooth p−1
  (e.g. 411840 = 2^6·3^2·5·11·13).
- Secondary covariate: within class 1 mod 840, the number of divisors ≤ 100
  of p−1 correlates *negatively* with f (Spearman −0.27) and negatively
  with N_typeI (−0.24); bottom-100 primes have strikingly smooth p−1
  (median largest prime factor 48 vs 395 for controls). Mechanistic
  reading: a Type I triple applies iff p ≡ −(a+b)/c mod 4ab, and
  p ≡ 1 mod 4ab (forced when p−1 is divisible by many small numbers) never
  satisfies this for small (a,b,c) — being ≡ 1 mod everything small is the
  extreme point of the QR obstruction. This is why class 1 mod 840 is the
  deepest class and why smooth p−1 deepens it further. The correlation
  fades in shallower classes (−0.08 in 529, −0.04 in 601), consistent with
  the covering family already dominating f there.

What this is *not*: a complete characterization. Spearman 0.57 leaves a
large residual; N_typeI ≤ 4 misses half the bottom set. The sharp exact
statement "bottom 2% ⊆ six QR classes" fails at rate 4/195 ≈ 2%.

### 4. Side result: growth of min f (001, lead 4)

Min f per dyadic band, p ≡ 1 mod 24: 34 (2^15–2^16), 37, 46, 52, **67**
(2^19–2^20, at p = 589681 ≡ 1 mod 840). Min/median ratio is stable (~0.35),
and the growth 34 → 67 over log p 11.1 → 13.5 matches (log p)³ scaling
(ratio predicted 1.85, observed 1.97). No sign of a slow-growing hard
subsequence; the minimum is always attained in a QR class (classes of the
band minima: 121, 169, 1, 1, 1). f(p) ≥ 34 for all p ≡ 1 mod 24 in
[2^15, 10^6); overall min for p > 1000 remains f(2521) = 9.

## Why it failed / what survived

Nothing failed procedurally; the anomaly under test was refuted. The 001
observation "601 holds 6 of the bottom 50" was a small-sample fluctuation
sitting on top of a real but unremarkable class-level gradient (601 = lowest
mean f among covering-count-1 classes). The mod-7 "secondary effect"
speculation in 001 (601 ≡ −1 mod 7) is dead: classes 649, 97 mod 840 share
601's non-QR mod-7 structure and show nothing, and 601 itself vanishes from
the bottom tail under size normalization.

Survived / strengthened:
- The QR-obstruction stratification is far sharper at 10^6 than 001 could
  see: 98% of the size-normalized bottom tail is the six QR classes
  (aggregate p ≈ 10^−110), and *no* non-QR class is enriched.
- 001's lead 1 (f tracks identity applicability) is confirmed
  quantitatively at class level (r = 0.92) and per-prime level
  (Spearman 0.57), with a mechanism (p ≡ 1 mod many small moduli blocks
  every small Type I congruence) that also explains the within-class
  p−1-smoothness effect.
- Validated parallel kernel `es_fcount.c` (~4 min for this run; list mode
  enables targeted sampling at larger p) and the checkpointed dataset of
  f(p) for all 9732 primes ≡ 1 mod 24 below 10^6.

## Leads generated

1. **Predict f, not just rank it.** N_typeI with a=b≤30 explains part of f;
   a regression of f against the full applicable-triple count with proper
   weighting (each triple contributes solutions ∝ divisor structure of d)
   might explain most of the variance. If f ≈ g(applicable identities) with
   small residual, the "hard primes" are exactly the identity-poor ones and
   the conjecture's computational frontier can be screened by N_typeI alone
   (cheap: O(B²) congruences per p, no factoring, no divisor sums) —
   relevant for scanning p ≫ 10^7 where exact f is unaffordable.
2. **The ≡ 1 mod M direction.** The deepest primes are ≡ 1 mod many small
   moduli (class 1 mod 840, smooth p−1). Conjecture-shaped question: is
   min f over p ≡ 1 mod lcm(1..k) attained arbitrarily deep, i.e. does
   "p ≡ 1 mod everything small" capture the true extremal family? Test by
   computing f on primes ≡ 1 mod 2520·k for k = 1, 2, 3 up to 10^7 with the
   list-mode kernel.
3. **Where do the residual representations come from?** For the six QR
   classes no small Type I triple applies, yet f ≥ 34 at 10^6-scale. Dump
   full solution lists for the band minima (589681, 471241, 132721, …) and
   classify solutions by Elsholtz–Tao type and by which modulus finally
   admits them; the shape of the *guaranteed* solutions for QR-class primes
   is exactly what any proof must produce.
4. **Min-f growth tracking (log p)³** deserves a wider baseline: with the
   sampling kernel, estimate band-minima up to 10^8 on QR-class primes only
   (they own every band minimum) to test whether min f / (log p)³ converges.

## Reproduction

```
cc -O2 -fopenmp -o tools/es_fcount tools/es_fcount.c
tools/es_fcount_run.sh                         # ~4 min on 4 cores; resumable
python3 tools/es_class_stats.py                # enrichment tables + p-values
python3 tools/es_structure_hunt.py             # bottom-100 vs control, classes
```
