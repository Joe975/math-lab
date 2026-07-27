# 001 — Baseline z(n) census and witness scan

- **Problem:** zaremba, `problems/zaremba/PROBLEM.md`
- **Date:** 2026-07-27
- **Mode:** informed
  (the problem was added this same cycle, so the prior-art index was empty;
  "informed" because the author also wrote PROBLEM.md and read the whole
  repo, so nothing here is independent of the lab's framing)
- **Type:** computational search
- **Tools:** `harness/zaremba/zscan.c` (C kernel, unsigned 64-bit integer
  arithmetic throughout, deterministic; `cc -O2 -o zscan zscan.c`),
  `harness/zaremba/zaremba.py` (stdlib reference implementation),
  `problems/zaremba/explore/crosscheck_tree.py` (independent generation-based
  implementation), `problems/zaremba/explore/verify_witnesses.py` (witness
  file integrity + fixed-seed sampling). Runtimes, single core: exact census
  to 10^5 ~40 s; witness scan to 10^7 ~25 s; check-only scan to 10^8 ~8 min;
  tree cross-check to 10^5 ~2 min.
- **Sources:** Bourgain–Kontorovich arXiv:1107.3776 (statement, history,
  Hensley counterexample); Kontorovich arXiv:1208.5460 [T] and Huang
  arXiv:1310.3772 [T] — both consulted via search-result summaries rather
  than the PDFs themselves.

## Approach

Baseline census, in two halves matching the two claim types the PROBLEM.md
contract defines:

1. **Witness scan**: for every n up to a large bound, find the smallest
   coprime m whose canonical partial quotients are all ≤ 5. Any n with no
   witness would refute Zaremba's conjecture outright; every n found gets a
   checkable certificate.
2. **Exact census**: z(n) exactly — minimum over *all* coprime m — on a
   smaller range, which is what every structural question (how big is
   {z ≥ 4}? where does it stop?) actually needs.

Why this rather than the obvious alternative: the problem had no prior art
in this lab, and censuses are the substrate every later claim will lean on.
Reproducing the density-one machinery (Bourgain–Kontorovich/Huang) would
cost a theory build-out and produce no per-n certificates; a witness-only
scan without the exact half would have missed the exceptional-set structure,
which turned out to be the interesting output.

Two kernel design points, both load-bearing for the runtimes quoted:

- The canonical expansion (last quotient ≥ 2) is exactly what plain Euclid
  emits for reduced m/n with 0 < m < n, so no normalization step exists to
  get wrong.
- The witness scan starts just below m = n·(3√5−5)/10 rather than at m = 1
  or m = n/6. The infimum of the K ≤ 5 fraction set is
  [0; 5,1,5,1,…] = (3√5−5)/10 = 0.1708203932…, so every m below that line is
  wasted work. This was found the slow way: starting at n/6 makes mean tries
  grow *linearly* (≈ 0.00415·n, the measure of the dead zone between 1/6 and
  the infimum — observed 2092 at slice mean n = 5×10^5, predicted 2077),
  which is a ~60× slowdown at 10^7.

## What was done

Statements under test, precisely:

- **S1.** For every n, 1 ≤ n ≤ 10^7, there exists m with gcd(m, n) = 1 whose
  canonical partial quotients are all ≤ 5 (i.e. z(n) ≤ 5), with the witness
  recorded.
- **S1′.** Same for 1 ≤ n ≤ 10^8, witnesses checked in-kernel but not
  retained.
- **S2.** Exact z(n) for 2 ≤ n ≤ 10^5.

Commands (repo root; `out/` is git-ignored scratch — everything quoted below
regenerates from these):

```bash
cc -O2 -o zscan harness/zaremba/zscan.c
./zscan check 5 1 10000000 out/zaremba-witness-1e7.txt   # S1
./zscan check 5 1 100000000                              # S1'
./zscan exact 2 100000 out/zaremba-ztable-1e5.txt        # S2
python problems/zaremba/explore/crosscheck_tree.py 100000 out/zaremba-tree-1e5.txt
diff out/zaremba-ztable-1e5.txt out/zaremba-tree-1e5.txt  # must be empty
python problems/zaremba/explore/verify_witnesses.py out/zaremba-witness-1e7.txt \
    2000 858840 1200870 2648448 3744300 4044090 5454930 6865824 7602132 \
    8604750 9495570   # file structure + per-slice worst cases from the S1 log
python harness/zaremba/zaremba.py verify-file out/zaremba-witness-1e7.txt --max 5
```

Cross-checks performed, against what independent thing:

- **S2, full range:** `crosscheck_tree.py` recomputes the entire z table for
  2 ≤ n ≤ 10^5 by the opposite algorithm — forward generation of every
  canonical word with digits ≤ 5 via the continuant recurrence
  q_i = a_i·q_{i−1} + q_{i−2}, marking each reached denominator with the
  smallest max-digit over all words reaching it. No Euclidean division and no
  gcd occurs anywhere in it; continuants are coprime by construction. Python
  vs C, generation vs per-pair expansion, and the two 99,999-line tables are
  **identical** (`diff` empty), so agreement is evidence about correctness,
  not determinism.
- **S1, complete:** `verify_witnesses.py` confirmed the witness file has
  exactly one line per n = 1..10^7 in order (plus a fixed-seed sample and
  the ten per-slice worst cases); then `zaremba.py verify-file` re-verified
  **all 9,999,999 certificates** with the independent stdlib implementation:
  0 bad, 18 s.
- **Kernel vs reference, overlapping range:** kernel and `zaremba.py` agree
  on z(n) for n ≤ 200 and on first witnesses for n ≤ 300 (identical files).
- z(6) = 5 reproduces the known sharp case by hand: 5/6 = [0; 1, 5].
- S1′ has **no independent pass** — same kernel, larger range, certificates
  discarded. It is listed as weaker evidence deliberately.

No proof steps are claimed anywhere in this record; the only unproven
assertions are in "Why it failed / what survived" and "Leads generated",
labelled there.

## Outcome

`VERIFIED` for the computations, in the ranges stated below; `EVIDENCE` —
never more — for the conjecture itself:

- **S1: z(n) ≤ 5 for all 1 ≤ n ≤ 10^7**, zero failures; every one of the
  9,999,999 witness certificates independently re-verified (file is 152 MB,
  regenerable by the command above; not committed).
- **S1′: z(n) ≤ 5 for all 1 ≤ n ≤ 10^8**, zero failures, single
  implementation, certificates not retained.
- **S2, for 2 ≤ n ≤ 10^5:** z(n) = 5 for exactly three n — **6, 54, 150** —
  and z(n) = 4 for exactly **35** more, the largest being **6234**. Every n
  with 6235 ≤ n ≤ 10^5 has z(n) ≤ 3. Histogram: z=2: 47,049; z=3: 52,912;
  z=4: 35; z=5: 3. Full exceptional set with witnesses and expansions:
  `data/z-exceptional-1e5.txt`.
- **Exceptional-set structure:** 37 of the 38 members of {n ≤ 10^5 :
  z(n) ≥ 4} are even; the only odd member is 1155 = 3·5·7·11. All residues
  mod 4 occur; no single congruence characterization.
- **z = 2 share rises with N:** 0.4174 (N=10^3) → 0.4398 (10^4) → 0.4705
  (10^5), residues mod 4 all present — consistent with the {1,2}-alphabet
  set having Hausdorff dimension 0.5312… > 1/2.
- **First witnesses pin the Cantor-set edge:** over n ≤ 10^7,
  min m_first/n = 0.170820393 (the infimum (3√5−5)/10 to 9 decimals),
  mean = 0.170834, max = 5/6 (n = 6). The smallest witness is, in almost
  every case, essentially the first rational above the left edge of the
  K ≤ 5 set. Mean tries from the infimum grows slowly: 16.1 (first 10^6)
  → 33.9 (tenth 10^6); worst single n in 10^7: 1308 tries (n = 8,604,750).

**Not claimed:** anything about n > 10^8 (S1′), n > 10^7 with certificates
(S1), or n > 10^5 (S2); finiteness of {n : z(n) ≥ 4}; that 6234 is its
largest member (only: largest below 10^5); any fitted growth law.

## Why it failed / what survived

Nothing failed in the refutation sense — but a census cannot decide the
conjecture, and its most suggestive output is exactly what a finite window
can never settle: the exceptional set {n : z(n) ≥ 4} stops at 6234 and stays
empty for the remaining ~94% of the census range. A gap that long is either
finiteness or a very slowly thinning sequence, and the two are
indistinguishable from inside any finite range. `SPECULATION`: the set is
finite (equivalently z(n) ≤ 3 for all n > 6234); the counting heuristic —
the number of {1..4}-witnesses per denominator grows like n^{2δ₄−1} with
2δ₄−1 ≈ 0.57 — says exceptions should die out, but nothing here approaches a
proof route.

The real obstruction found is quantitative and reusable: naive witness
scanning is dominated by the **dead zone below the Cantor-set infimum**, and
after removing it, by the gap structure of the K ≤ 5 set just above its left
edge. First witnesses concentrate at that edge (mean m/n within 1.4×10^−5 of
the infimum), so per-n scanning does ~30–40 cheap aborts each; that is fine
to 10^8 and hopeless well before 10^10. The scaling tool for lead 2 is the
other implementation: generation-and-mark is output-linear and embarrassingly
sliceable.

Reusable: the kernel (`zscan.c`), the reference CLI (`zaremba.py`), the
generation cross-checker (`crosscheck_tree.py`, also the right skeleton for
a bitmap mark-all scan), the witness verifier, and the exceptional-set data
file with certificates.

## Leads generated

1. **Exact census to 10^6** (`zscan exact 2 1000000`, ~100× the 10^5 cost;
   slice it per `harness/common/run_slices.sh` conventions). Definite either
   way: any member of {z ≥ 4} above 6234 kills the finiteness speculation as
   stated; none extends the empty gap by a decade.
2. **Mark-all scan to 10^9:** port `crosscheck_tree.py` to C with a bitmap
   (1 bit per n = 125 MB), count unreached n per decade. Decides z ≤ 5
   coverage at 10^9 in one output-linear pass, and yields the z ≤ 3
   hitting-density curve for free.
3. **Odd exceptional members:** 37/38 are even yet 1155 exists. Exact z on
   odd n only, to 10^6 (half-range, and the scan over m is unchanged).
   Falsifiable: a second odd member, or a decade of none.
4. **z = 2 residue structure vs the corrected local–global conjecture:**
   check the z = 2 set against moduli 8, 12, 16 (it missed nothing mod 4). A
   missed class would be an obstruction datum for the {1,2} alphabet; none
   is a datum for Bourgain–Kontorovich's corrected conjecture.
5. **Hard-n anatomy:** the per-slice worst n (858840, …, 9495570) and the
   exceptional set share heavy divisibility by small primes. Compute
   tries(n) vs ω(n)/smoothness on the 10^7 witness log and test whether
   "hard for the scan" and "large z" are the same phenomenon or two.

## References

- J. Bourgain, A. Kontorovich, *On Zaremba's conjecture*, Ann. of Math. 180
  (2014); arXiv:1107.3776.
- S. Huang, *An improvement to Zaremba's conjecture*, GAFA 25 (2015);
  arXiv:1310.3772. [T]
- A. Kontorovich, *From Apollonius to Zaremba*, Bull. AMS 50 (2013);
  arXiv:1208.5460. [T]
- H. Niederreiter, *Dyadic fractions with small partial quotients*,
  Monatsh. Math. 101 (1986). [T] — attribution via secondary sources.
- This repo: `problems/zaremba/PROBLEM.md` (contract);
  `harness/zaremba/zscan.c`, `harness/zaremba/zaremba.py`,
  `problems/zaremba/explore/crosscheck_tree.py`,
  `problems/zaremba/explore/verify_witnesses.py`,
  `problems/zaremba/data/z-exceptional-1e5.txt`.
