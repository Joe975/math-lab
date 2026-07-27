"""Independent re-derivations of every number this route claims.

Each block below recomputes a claim by a *different formulation* from the one
that produced it, so that agreement is evidence rather than a repeated run.
"""

import os
import sys
from itertools import combinations, product

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)

import srg

FAILED = []


def claim(label, a, b):
    ok = (a == b)
    if not ok:
        FAILED.append(label)
    print(f"  {'OK  ' if ok else 'FAIL'} {label}: {a!r} vs {b!r}")


print("=== 1. multiplicities of (99,14,1,2) ===")
# route A: solve the 2x2 linear system (harness)
sA = srg.spectrum(99, 14, 1, 2)
# route B: from the closed form f,g = 1/2[(n-1) -+ (2k+(n-1)(lam-mu))/sqrt(D)]
n, k, lam, mu = 99, 14, 1, 2
D = (lam - mu) ** 2 + 4 * (k - mu)
t = srg.exact_isqrt(D)
num = 2 * k + (n - 1) * (lam - mu)
assert num % t == 0
fB = ((n - 1) - num // t) // 2
gB = ((n - 1) + num // t) // 2
claim("multiplicity of eigenvalue 3", sA["f"], fB)
claim("multiplicity of eigenvalue -4", sA["g"], gB)
# route C: the two trace identities, solved by brute force
solC = [(f, n - 1 - f) for f in range(n)
        if k + 3 * f - 4 * (n - 1 - f) == 0
        and k * k + 9 * f + 16 * (n - 1 - f) == n * k]
claim("trace-identity brute force", solC, [(54, 44)])
print(f"  NOTE: the tier-0 statement gives 56 and 42; those fail tr(A)=0 "
      f"({k + 3*56 - 4*42}, need 0) and tr(A^2)=nk "
      f"({k*k + 9*56 + 16*42}, need {n*k}).")

print()
print("=== 2. the lam=1, mu=2 family ===")
# route A: harness scan over k (enumerate_feasible)
famA = [(x[0], x[1]) for x in srg.enumerate_feasible(1, 2, 2_000_000)]
# route B: parameterise by the eigenvalue gap t, 4k-7 = t^2, and test
# multiplicity integrality from the trace equation directly
famB = []
tt = 1
while True:
    kk = (tt * tt + 7) // 4
    if (tt * tt + 7) % 4 == 0 and kk >= 2:
        nn = (kk * kk + 2) // 2
        if nn > 2_000_000:
            break
        r, s = (tt - 1) // 2, (-1 - tt) // 2
        # k + f r + g s = 0, f + g = n-1
        den = r - s
        numr = -kk - (nn - 1) * s
        if den and numr % den == 0:
            f = numr // den
            g = nn - 1 - f
            if f >= 0 and g >= 0 and nn > kk + 1:
                famB.append((nn, kk))
    tt += 2
    if tt > 5000:
        break
claim("feasible (n,k) with lam=1,mu=2, n<=2e6", famA, famB)

print()
print("=== 3. the forced pair model's arithmetic for k=14 ===")
# route A: counts produced by the propagator (see partner_regular.py output)
# route B: direct combinatorial arithmetic, written independently here
npos = 7
nblocks = npos * (npos - 1) // 2               # 21
nd = nblocks * 4                               # 84
claim("|D|", nd, 14 * 13 // 2 - 7)
total_pairs = nd * (nd - 1) // 2
claim("total D-pair variables", total_pairs, 3486)
in_block = nblocks * (4 * 3 // 2)              # 6 pairs per block
# blocks meeting a given block (sharing exactly one position)
meet_per_block = 2 * (npos - 2)                # 10
meeting_pairs = nd * (meet_per_block * 4) // 2
disjoint_blocks = (npos - 2) * (npos - 3) // 2  # C(5,2) = 10
disjoint_pairs = nd * (disjoint_blocks * 4) // 2
claim("in-block + meeting + disjoint == total",
      in_block + meeting_pairs + disjoint_pairs, total_pairs)
claim("pairs decided by partner-regularity (in-block 126 + meeting 1680)",
      in_block + meeting_pairs, 1806)
claim("pairs left open", disjoint_pairs, 1680)
claim("edges forced (2 per D-vertex)", nblocks * 4, 84)
claim("neighbours still needed per D-vertex", 12 - 2, 10)
claim("candidates per D-vertex", disjoint_blocks * 4, 40)

print()
print("=== 4. the n_B multiplicity profiles, second algorithm ===")
# route A (probe): brute force over 3^10 assignments
R = list(range(5))
blocks = list(combinations(R, 2))
solA = []
for cur in product((0, 1, 2), repeat=10):
    if sum(cur) != 10:
        continue
    degs = [0] * 5
    for m, (a, b) in zip(cur, blocks):
        degs[a] += m
        degs[b] += m
    if all(d == 4 for d in degs):
        solA.append(cur)
# route B: choose the multiplicity-2 support first, then solve for the
# multiplicity-1 support as a subgraph with prescribed remaining degrees
solB = []
for r2 in range(0, 6):
    for supp2 in combinations(range(10), r2):
        deg2 = [0] * 5
        for i in supp2:
            a, b = blocks[i]
            deg2[a] += 2
            deg2[b] += 2
        if any(d > 4 for d in deg2):
            continue
        rem = [4 - d for d in deg2]
        need1 = sum(rem) // 2
        if sum(rem) % 2 or 10 - 2 * r2 != need1:
            continue
        others = [i for i in range(10) if i not in supp2]
        for supp1 in combinations(others, need1):
            deg1 = list(deg2)
            for i in supp1:
                a, b = blocks[i]
                deg1[a] += 1
                deg1[b] += 1
            if all(d == 4 for d in deg1):
                v = [0] * 10
                for i in supp2:
                    v[i] = 2
                for i in supp1:
                    v[i] = 1
                solB.append(tuple(v))
import hashlib
claim("labelled n_B profiles (sha256 of the sorted list)",
      hashlib.sha256(repr(sorted(solA)).encode()).hexdigest()[:16],
      hashlib.sha256(repr(sorted(solB)).encode()).hexdigest()[:16])
claim("number of labelled n_B profiles", len(solA), 73)

print()
print("=== 5. automorphism arithmetic ===")
claim("99 mod 7 (fixed points of an order-7 automorphism)", 99 % 7, 1)
claim("99 mod 5 (fixed points of an order-5 automorphism, mod 5)", 99 % 5, 4)
claim("order of the fixed subgraph for p=5, from k=4,lam=1,mu=2",
      [f for f in range(1, 200) if srg.basic_identity(f, 4, 1, 2)], [9])
claim("9 mod 5 matches 99 mod 5", 9 % 5, 99 % 5)
claim("99 - 9 divisible by 5", (99 - 9) % 5, 0)
claim("99 = 1 + 14*7", 1 + 14 * 7, 99)
claim("primes dividing 99", sorted({3, 11}), sorted({p for p in (2,3,5,7,11,13)
                                                     if 99 % p == 0}))

print()
print("ALL CROSS-CHECKS AGREE" if not FAILED else f"DISAGREEMENTS: {FAILED}")
sys.exit(1 if FAILED else 0)
