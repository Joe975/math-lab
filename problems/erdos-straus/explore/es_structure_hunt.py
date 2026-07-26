#!/usr/bin/env python3
"""Structure hunt on low-f Erdos-Straus primes (attempt 002, step 3).

Loads the f(p) dataset for p ≡ 1 mod 24, p < 10^6 (chunk files from
es_fcount_run.sh), takes the bottom-100 primes (size-normalized rank and raw),
and compares them against a size-matched random control sample of the same
population on:

  * factorization shape of p-1 and p+1 (d(.), largest prime factor);
  * quadratic-residue class memberships (mod 5, 7, 840; classes mod
    24/40/56/120/168);
  * N_B(p): the number of small Elsholtz-Tao Type I parameter triples
    (a <= b <= B, c | a+b) whose identity applies to p, i.e.
    4ab | c*p + a + b.  This generalizes "how many identity families cover
    p's class" to a per-prime count (001, lead 1).

Also computes, per residue class mod 840 (the 24 classes ≡ 1 mod 24), the
number of full-class Type I covering families (ab | 210 sweep from 001) and
the class's mean/min f on [10^5, 10^6), to test the class-level version of
the same hypothesis.

Usage: python3 es_structure_hunt.py [--data DIR] [--nbot 100] [--abound 30]
"""
import argparse
import glob
import math
import os
import random
import statistics as st


def load(data_dir):
    data = []
    for fn in sorted(glob.glob(os.path.join(data_dir, "f1mod24_*.txt"))):
        for line in open(fn):
            p, f = map(int, line.split())
            data.append((p, f))
    data.sort()
    return data


def factorize(n):
    fac = {}
    t = n
    for q in range(2, int(n ** 0.5) + 1):
        while t % q == 0:
            fac[q] = fac.get(q, 0) + 1
            t //= q
        if t == 1:
            break
    if t > 1:
        fac[t] = fac.get(t, 0) + 1
    return fac


def ndiv(fac):
    r = 1
    for e in fac.values():
        r *= e + 1
    return r


def fac_str(fac):
    return "*".join(f"{q}^{e}" if e > 1 else str(q)
                    for q, e in sorted(fac.items()))


def n_typeI(p, B):
    """# of (a,b,c): a<=b<=B, c | a+b, 4ab | c*p + a + b."""
    cnt = 0
    for a in range(1, B + 1):
        for b in range(a, B + 1):
            m = 4 * a * b
            s = a + b
            for c in range(1, s + 1):
                if s % c == 0 and (c * p + s) % m == 0:
                    cnt += 1
    return cnt


def class_cover_count():
    """Per class r mod 840: # of full-class Type I families (001 sweep)."""
    cover = {r: 0 for r in range(840)}
    for a in range(1, 211):
        for b in range(a, 211):
            if 210 % (a * b):
                continue
            m4 = 4 * a * b
            s = a + b
            for c in range(1, s + 1):
                if s % c:
                    continue
                r0 = (-(s // c)) % m4
                for r in range(r0, 840, m4):
                    cover[r] += 1
    return cover


def normalized_bottom(data, k):
    bands = {}
    for p, f in data:
        bands.setdefault(p.bit_length(), []).append((f, p))
    rank = {}
    for b, v in bands.items():
        v.sort()
        for i, (f, p) in enumerate(v):
            rank[p] = (i + 0.5) / len(v)
    return [p for p, _ in sorted(((p, rank[p]) for p, _ in data),
                                 key=lambda t: t[1])[:k]], rank


def summarize(name, primes, fmap, B):
    qr5 = {1, 4}
    qr7 = {1, 2, 4}
    rows = dict(
        d_pm1=[], d_pp1=[], lpf_pm1=[], lpf_pp1=[], ntyp=[],
        qr5=0, qr7=0, qr840=0)
    for p in primes:
        f1 = factorize(p - 1)
        f2 = factorize(p + 1)
        rows["d_pm1"].append(ndiv(f1))
        rows["d_pp1"].append(ndiv(f2))
        rows["lpf_pm1"].append(max(f1))
        rows["lpf_pp1"].append(max(f2))
        rows["ntyp"].append(n_typeI(p, B))
        rows["qr5"] += p % 5 in qr5
        rows["qr7"] += p % 7 in qr7
        rows["qr840"] += (p % 5 in qr5) and (p % 7 in qr7)
    n = len(primes)
    print(f"\n=== {name} (n={n}) ===")
    for k in ("d_pm1", "d_pp1", "lpf_pm1", "lpf_pp1", "ntyp"):
        v = rows[k]
        print(f"  {k:8s}: mean={st.mean(v):9.1f} median={st.median(v):8.1f} "
              f"min={min(v):7d} max={max(v):8d}")
    print(f"  QR mod5: {rows['qr5']}/{n}  QR mod7: {rows['qr7']}/{n}  "
          f"QR840: {rows['qr840']}/{n}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--data", default=os.path.join(
        here, "..", "attempts", "erdos-straus", "data"))
    ap.add_argument("--nbot", type=int, default=100)
    ap.add_argument("--abound", type=int, default=30)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    data = load(args.data)
    fmap = dict(data)
    print(f"loaded {len(data)} primes ≡ 1 mod 24")

    bot_norm, rank = normalized_bottom(data, args.nbot)
    order_raw = sorted(data, key=lambda t: (t[1], t[0]))
    bot_raw = [p for p, _ in order_raw[:args.nbot]]

    # size-matched control: same dyadic band multiset as bot_norm
    rng = random.Random(args.seed)
    byband = {}
    for p, _ in data:
        byband.setdefault(p.bit_length(), []).append(p)
    control = []
    botset = set(bot_norm)
    for p in bot_norm:
        band = byband[p.bit_length()]
        while True:
            c = rng.choice(band)
            if c not in botset:
                control.append(c)
                break

    print("\nBottom-30 primes by size-normalized rank "
          "(p, f, p%840, p-1, p+1):")
    for p in bot_norm[:30]:
        print(f"  p={p:7d} f={fmap[p]:5d} p%840={p % 840:3d} "
              f"p%120={p % 120:3d}  p-1={fac_str(factorize(p - 1))}  "
              f"p+1={fac_str(factorize(p + 1))}")

    b = summarize(f"bottom-{args.nbot} (normalized rank)", bot_norm, fmap,
                  args.abound)
    summarize(f"bottom-{args.nbot} (raw f)", bot_raw, fmap, args.abound)
    c = summarize("size-matched control", control, fmap, args.abound)

    # N_typeI separation diagnostics
    for thresh in range(0, 8):
        tp = sum(1 for v in b["ntyp"] if v <= thresh)
        fp = sum(1 for v in c["ntyp"] if v <= thresh)
        print(f"  N_typeI(B={args.abound}) <= {thresh}: "
              f"bottom {tp}/{len(b['ntyp'])}, control {fp}/{len(c['ntyp'])}")

    # Spearman-ish check on a size window: correlation of f with N_typeI
    win = [(p, f) for p, f in data if 5 * 10**5 <= p < 10**6]
    sub = rng.sample(win, min(800, len(win)))
    xs = [n_typeI(p, args.abound) for p, _ in sub]
    ys = [f for _, f in sub]
    def ranks(v):
        s = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(s):
            j = i
            while j + 1 < len(s) and v[s[j + 1]] == v[s[i]]:
                j += 1
            for t in range(i, j + 1):
                r[s[t]] = (i + j) / 2
            i = j + 1
        return r
    rx, ry = ranks(xs), ranks(ys)
    mx, my = st.mean(rx), st.mean(ry)
    num = sum((a - mx) * (b2 - my) for a, b2 in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx)
                    * sum((b2 - my) ** 2 for b2 in ry))
    print(f"\nSpearman corr(f, N_typeI(B={args.abound})) on {len(sub)} primes "
          f"in [5e5,1e6): {num / den:.3f}")

    # class-level: covering families vs class f stats
    cover = class_cover_count()
    big = [(p, f) for p, f in data if p >= 10**5]
    bycls = {}
    for p, f in big:
        bycls.setdefault(p % 840, []).append(f)
    bset = set(bot_norm)
    print("\nclass mod840 | #TypeI-covering-families | mean f | min f | "
          "#in bottom-100 (norm)")
    rows = []
    for cls in sorted(bycls):
        v = bycls[cls]
        nb = sum(1 for p in bot_norm if p % 840 == cls)
        rows.append((cover[cls], cls, st.mean(v), min(v), nb))
    for cov, cls, mf, mnf, nb in sorted(rows):
        print(f"  {cls:4d} | {cov:3d} | {mf:7.1f} | {mnf:4d} | {nb}")
    # correlation of class mean f with covering count
    xs = [r[0] for r in rows]
    ys = [r[2] for r in rows]
    mx, my = st.mean(xs), st.mean(ys)
    num = sum((a - mx) * (b2 - my) for a, b2 in zip(xs, ys))
    den = math.sqrt(sum((a - mx) ** 2 for a in xs)
                    * sum((b2 - my) ** 2 for b2 in ys))
    print(f"Pearson corr(class mean f, class covering count): {num/den:.3f}")


if __name__ == "__main__":
    main()
