# Singmaster — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`. Full records: `attempts/`.

## Editorial view of the attack surface

- Extend the search for entries with multiplicity ≥ 8; search in (n, k) space
  via collision-finding among binomials with small k.
- Build the table of Diophantine equations C(n, j) = C(m, k) for small fixed
  j < k — each pair is a curve; which are resolved, which open?
- Boundary vs interior: map exactly what parameter region the MRST-T result
  leaves open.

## Attempts

### 001 — Binomial-coefficient collision search to 2.5×10^29 · `VERIFIED`

Complete census of multiplicities ≥ 5 for all values ≤ 2.5×10^29 (~4× beyond
C(104,39) ≈ 6.12×10^28):

| multiplicity | count | values |
|---|---|---|
| 8 | 1 | 3003 = C(78,2) = C(15,5) = C(14,6) |
| 7 | 0 | — |
| 6 | 7 | 120, 210, 1540, 7140, 11628, 24310, C(104,39) = C(103,40) |
| 5 | 0 | — |

**Completeness argument** (this is what makes it a census rather than a
search). For fixed k, C(n,k) is strictly increasing in n ≥ k, so V has at most
one canonical cell per k. One canonical cell gives multiplicity ≤ 4. Hence
mult(V) ≥ 5 requires ≥ 2 canonical cells with distinct k, at most one of which
has k = 2 — so every V ≤ B with mult ≥ 5 has a canonical cell with k ≥ 3, and
all of its canonical cells have value ≤ B. Enumerating canonical cells with
k ≥ 3 and value ≤ B therefore suffices.

Scanned 11,447,142,421 cells at k = 3, 49,492,314 at k = 4, 2,304,208 at
k ∈ [5,50]. 2,066 s on 4 workers, peak ~1.0 GB. Self-test plus independent
brute force to 10^7, double re-verification of every hit, and agreement with
de Weger's coincidence list.

**Multiplicity convention used:** mult(V) = #{(n,k) : 0 ≤ k ≤ n, C(n,k) = V},
counting symmetric cells separately and including the trivial cells C(V,1),
C(V,V−1). Under this convention 3003 has multiplicity 8.

## Open lines

- Diophantine curve table: all in-range coincidences come from known
  families, so searching deeper is now low value; the equation table is the
  better use of budget.
