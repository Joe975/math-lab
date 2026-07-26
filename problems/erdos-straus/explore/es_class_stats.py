#!/usr/bin/env python3
"""Class-enrichment statistics for Erdos-Straus representation counts f(p).

Reads the checkpointed chunk files produced by es_fcount_run.sh
(attempts/erdos-straus/data/f1mod24_*.txt: lines "p f(p)", primes
p ≡ 1 mod 24, p < 10^6) and answers attempt 002's questions:

  * For each residue class mod 840 / 168 / 120 (restricted to p ≡ 1 mod 24),
    how over-represented is the class among the bottom 0.5% / 2% of f-values,
    vs its share of the population?  Exact counts, enrichment ratios and
    one-sided binomial tail p-values.
  * Two selection modes for "bottom": raw f (population share is still a fair
    null because the classes are equidistributed in p), and size-normalized
    rank (rank of f within a dyadic band of p, so small primes cannot
    dominate the bottom set).
  * Class 601 mod 840 is the anomaly under test (6 of bottom 50 at 10^5).

Usage: python3 es_class_stats.py [--data DIR] [--json OUT]
"""
import argparse
import glob
import json
import math
import os
import sys

QR840 = sorted({u * u % 840 for u in range(840) if math.gcd(u, 840) == 1})


def load(data_dir):
    data = []
    for fn in sorted(glob.glob(os.path.join(data_dir, "f1mod24_*.txt"))):
        for line in open(fn):
            p, f = map(int, line.split())
            data.append((p, f))
    data.sort()
    assert all(p % 24 == 1 for p, _ in data)
    assert len({p for p, _ in data}) == len(data)
    return data


def binom_sf(k, n, q):
    """P[Bin(n, q) >= k], exact-ish via lgamma logs."""
    if k <= 0:
        return 1.0
    lq, l1q = math.log(q), math.log1p(-q)
    tot = 0.0
    for i in range(k, n + 1):
        lp = (math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
              + i * lq + (n - i) * l1q)
        if lp < -745:
            if i > k + 10 and tot > 0:
                break
            continue
        tot += math.exp(lp)
    return min(tot, 1.0)


def bottom_sets(data, fracs, norm):
    """Return {frac: set of primes in the bottom `frac` of f-values}.

    norm=False: order by raw f (ties by p).
    norm=True:  order by percentile-rank of f within dyadic bands of p
                (band = floor(log2 p)), so every size range contributes
                proportionally.
    """
    if not norm:
        order = sorted(data, key=lambda t: (t[1], t[0]))
        keyed = [(p, None) for p, _ in order]
    else:
        bands = {}
        for p, f in data:
            bands.setdefault(p.bit_length(), []).append((f, p))
        rank = {}
        for b, v in bands.items():
            v.sort()
            n = len(v)
            for i, (f, p) in enumerate(v):
                rank[p] = (i + 0.5) / n
        keyed = sorted(((p, rank[p]) for p, _ in data), key=lambda t: t[1])
    out = {}
    for fr in fracs:
        k = max(1, round(fr * len(data)))
        out[fr] = set(p for p, _ in keyed[:k])
    return out


def enrichment_table(data, bottom, mod, label, min_class_n=1):
    N = len(data)
    classes = sorted({p % mod for p, _ in data})
    pop = {c: sum(1 for p, _ in data if p % mod == c) for c in classes}
    k = len(bottom)
    rows = []
    for c in classes:
        if pop[c] < min_class_n:
            continue
        q = pop[c] / N
        obs = sum(1 for p in bottom if p % mod == c)
        exp = k * q
        enr = obs / exp if exp else float("nan")
        pval = binom_sf(obs, k, q)
        rows.append((c, pop[c], obs, exp, enr, pval))
    rows.sort(key=lambda r: -r[4])
    print(f"\n--- mod {mod}, {label} (bottom set size {k} of {N}) ---")
    print(f"{'class':>6} {'QR840?':>6} {'pop':>6} {'obs':>4} {'exp':>7} "
          f"{'enrich':>7} {'binom p':>10}")
    for c, n_c, obs, exp, enr, pval in rows:
        qr = ""
        if mod == 840:
            qr = "QR" if c in QR840 else "-"
        elif mod == 120:
            qr = "QR" if c in {u * u % 120 for u in range(120)
                              if math.gcd(u, 120) == 1} else "-"
        elif mod == 168:
            qr = "QR" if c in {u * u % 168 for u in range(168)
                              if math.gcd(u, 168) == 1} else "-"
        print(f"{c:>6} {qr:>6} {n_c:>6} {obs:>4} {exp:>7.2f} "
              f"{enr:>7.2f} {pval:>10.3g}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--data", default=os.path.join(
        here, "..", "attempts", "erdos-straus", "data"))
    ap.add_argument("--json", default=None,
                    help="dump machine-readable results here")
    args = ap.parse_args()

    data = load(args.data)
    print(f"loaded f(p) for {len(data)} primes ≡ 1 mod 24, "
          f"p in [{data[0][0]}, {data[-1][0]}]")

    results = {}
    for norm in (False, True):
        tag = "normalized-rank" if norm else "raw-f"
        bots = bottom_sets(data, (0.005, 0.02), norm)
        for fr, bset in sorted(bots.items()):
            label = f"bottom {fr*100:g}% by {tag}"
            for mod in (840, 168, 120):
                rows = enrichment_table(data, bset, mod, label)
                results[f"{tag}|{fr}|{mod}"] = rows
    if args.json:
        with open(args.json, "w") as f:
            json.dump({k: [list(r) for r in v] for k, v in results.items()},
                      f, indent=1)
        print(f"\njson dumped to {args.json}")


if __name__ == "__main__":
    main()
