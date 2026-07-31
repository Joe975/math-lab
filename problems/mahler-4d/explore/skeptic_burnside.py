#!/usr/bin/env python3
"""SKEPTIC re-implementation (written from scratch for the adversarial review,
sharing no code with census.c, burnside_check.py or hanner_match.py).

1. Burnside orbit counts of k-subsets of the 40 antipodal pairs of nonzero
   {0,+-1}^4 points under B4 (order 384), for k = 1..K.
2. Kernel of the action on pairs (expected {+-I}, order 2).
3. Canonicality audit of level files: a mask is canonical iff it equals the
   minimum of its 40-bit images over the group.  Full audit for small levels,
   seeded random sample for large ones.  Plus a uniqueness check.

Usage:
  skeptic_burnside.py counts K
  skeptic_burnside.py audit LEVELFILE [SAMPLE] [SEED]
"""
import itertools
import random
import sys

# --- the 40 pair representatives, in the documented fixed order:
# iterate (c0,c1,c2,c3), each in -1,0,1 ascending, c0 outermost; skip 0;
# keep iff first nonzero coordinate is +1.
PAIRS = []
for c0 in (-1, 0, 1):
    for c1 in (-1, 0, 1):
        for c2 in (-1, 0, 1):
            for c3 in (-1, 0, 1):
                p = (c0, c1, c2, c3)
                nz = [x for x in p if x]
                if nz and nz[0] == 1:
                    PAIRS.append(p)
assert len(PAIRS) == 40
PIDX = {p: i for i, p in enumerate(PAIRS)}


def norm_pair(p):
    for x in p:
        if x:
            return p if x > 0 else tuple(-y for y in p)
    raise ValueError("zero point")


def group_pair_perms():
    """All 384 elements of B4 as permutations of the 40 pairs.
    Element = (sigma, eps): x -> (eps_i * x_sigma(i))_i."""
    out = []
    for sigma in itertools.permutations(range(4)):
        for eps in itertools.product((1, -1), repeat=4):
            perm = tuple(
                PIDX[norm_pair(tuple(eps[i] * p[sigma[i]] for i in range(4)))]
                for p in PAIRS)
            out.append(perm)
    return out


def cycle_lengths(perm):
    seen = [False] * len(perm)
    out = []
    for i in range(len(perm)):
        if seen[i]:
            continue
        l, j = 0, i
        while not seen[j]:
            seen[j] = True
            j = perm[j]
            l += 1
        out.append(l)
    return out


def burnside_counts(K):
    perms = group_pair_perms()
    ident = tuple(range(40))
    kernel = sum(1 for g in perms if g == ident)
    tot = [0] * (K + 1)
    for g in perms:
        # generating function prod over cycles (1 + x^len), truncated at K
        poly = [0] * (K + 1)
        poly[0] = 1
        for l in cycle_lengths(g):
            if l <= K:
                for k in range(K, l - 1, -1):
                    poly[k] += poly[k - l]
        for k in range(K + 1):
            tot[k] += poly[k]
    assert all(t % len(perms) == 0 for t in tot[1:]), "Burnside not integral!"
    return kernel, [t // len(perms) for t in tot]


# --- canonicality audit ------------------------------------------------

def make_tables():
    """Per distinct pair-permutation, 5 byte-lookup tables for fast mask
    image computation."""
    uniq = sorted(set(group_pair_perms()))
    tabs = []
    for g in uniq:
        t = []
        for byte in range(5):
            row = [0] * 256
            for v in range(256):
                m = 0
                for b in range(8):
                    if v >> b & 1:
                        m |= 1 << g[8 * byte + b]
                row[v] = m
            t.append(row)
        tabs.append(t)
    return tabs


def is_canonical(mask, tabs):
    b = [(mask >> (8 * i)) & 255 for i in range(5)]
    for t in tabs:
        img = t[0][b[0]] | t[1][b[1]] | t[2][b[2]] | t[3][b[3]] | t[4][b[4]]
        if img < mask:
            return False
    return True


def audit(path, sample, seed):
    tabs = make_tables()
    with open(path) as f:
        masks = [int(line.split()[0], 16) for line in f if line.strip()]
    n = len(masks)
    dup = n - len(set(masks))
    if sample and sample < n:
        rng = random.Random(seed)
        idx = rng.sample(range(n), sample)
        checked = [(i, masks[i]) for i in idx]
        scope = f"sample {sample}/{n} seed={seed}"
    else:
        checked = list(enumerate(masks))
        scope = f"ALL {n}"
    bad = [(i, m) for i, m in checked if not is_canonical(m, tabs)]
    print(f"{path}: {n} masks, {dup} duplicates, canonicality {scope}: "
          f"{len(bad)} NON-CANONICAL")
    for i, m in bad[:10]:
        print(f"  line {i + 1}: {m:x} NOT canonical")
    return not bad and dup == 0


if __name__ == "__main__":
    if sys.argv[1] == "counts":
        K = int(sys.argv[2])
        kernel, c = burnside_counts(K)
        print(f"kernel size of action on pairs: {kernel} (expect 2)")
        for k in range(1, K + 1):
            print(f"k={k:2d}  orbits={c[k]}")
    elif sys.argv[1] == "audit":
        sample = int(sys.argv[2 + 1]) if len(sys.argv) > 3 else 0
        seed = int(sys.argv[4]) if len(sys.argv) > 4 else 12345
        ok = audit(sys.argv[2], sample, seed)
        sys.exit(0 if ok else 1)
