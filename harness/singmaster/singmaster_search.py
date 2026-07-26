#!/usr/bin/env python3
"""Collision search for Singmaster's conjecture: find binomial coefficients of
high multiplicity in Pascal's triangle.

Multiplicity convention
-----------------------
mult(V) = #{(n,k) : 0 <= k <= n, C(n,k) = V}, counting the symmetric cells
C(n,k) and C(n,n-k) separately, and counting the trivial cells C(V,1) and
C(V,V-1).  Equivalently, for V >= 3:

    mult(V) = 2 + sum over *canonical* cells (n,k) with 2 <= k <= n/2 of
              (1 if n == 2k else 2)

Under this convention 3003 has multiplicity 8 via the cells
(3003,1),(3003,3002),(78,2),(78,76),(15,5),(15,10),(14,6),(14,8).

Completeness claim
------------------
For fixed k, C(n,k) is strictly increasing in n (n >= k), so a value has at
most one canonical cell per k.  One canonical cell gives mult <= 4; hence
mult(V) >= 5 requires >= 2 canonical cells with *distinct* k, at most one of
which has k = 2.  Therefore every V <= B with mult >= 5 has a canonical cell
with 3 <= k, and all its canonical cells have value V <= B.  The search
enumerates ALL canonical cells with k >= 3 and value <= B:

  * k = 3 and k = 4 are streamed as two sorted sequences (incremental
    big-int recurrences) and two-pointer merged, catching C(a,3) = C(b,4);
  * k >= 5 values are precomputed into a hash set S5 (small: O(B^{1/5}));
    collisions among k >= 5 and triangular membership are found at build;
  * every streamed value is tested for membership in S5 and for being
    triangular (C(m,2) = V  <=>  8V+1 is a perfect square), via a two-level
    quadratic-residue prefilter + math.isqrt.

So the output list of values with mult >= 5 is complete up to the bound.
Values with exactly one canonical cell (mult 3 if central, else 4) form
generic infinite families and are only counted, not listed.

Re-runnable; per-worker checkpoints make long runs restartable
(--resume, default on).  The k>=5 table is regenerated deterministically in
seconds, so it is not persisted.

Usage:
  singmaster_search.py --bound 2.5e29 [--procs 4] [--outdir DIR]
  singmaster_search.py --selftest          # validate on known facts (B=1e8)
  singmaster_search.py --brute 10000000    # independent brute-force cross-check
"""

import argparse
import json
import os
import sys
import time
from collections import Counter
from math import comb, isqrt

# ----------------------------------------------------------------------------
# square-test prefilter: t is a perfect square only if t is a QR mod M1 and M2
M1 = 693            # 9 * 7 * 11
M2 = 65 * 64        # 4160


def _sqmask(m):
    ba = bytearray(m)
    for i in range(m):
        ba[(i * i) % m] = 1
    return bytes(ba)


T1 = _sqmask(M1)
T2 = _sqmask(M2)

# globals shared with forked workers
G = {}


# ----------------------------------------------------------------------------
def max_n_for_k(bound, k):
    """Largest n with C(n,k) <= bound, or 0 if none with n >= 2k."""
    if comb(2 * k, k) > bound:
        return 0
    hi = 2 * k
    while comb(hi, k) <= bound:
        hi *= 2
    lo = 2 * k
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if comb(mid, k) <= bound:
            lo = mid
        else:
            hi = mid - 1
    return lo


def first_n_at_least(value, k):
    """Smallest n >= 2k with C(n,k) >= value."""
    lo, hi = 2 * k, 2 * k
    while comb(hi, k) < value:
        hi *= 2
    while lo < hi:
        mid = (lo + hi) // 2
        if comb(mid, k) < value:
            lo = mid + 1
        else:
            hi = mid
    return lo


def solve_nk(V, k):
    """Return n >= 2k with C(n,k) == V (canonical cell), or None.
    Uses math.comb directly -- independent of the incremental scan arithmetic."""
    if comb(2 * k, k) > V:
        return None
    n = first_n_at_least(V, k)
    return n if comb(n, k) == V else None


def exact_cells(V):
    """All canonical cells (n,k), 2 <= k <= n/2, with C(n,k) == V.
    Exact big-int verification via math.comb."""
    cells = []
    k = 2
    while comb(2 * k, k) <= V:
        n = solve_nk(V, k)
        if n is not None:
            cells.append((n, k))
        k += 1
    return cells


def multiplicity(V):
    """Full multiplicity under the stated convention (V >= 3)."""
    return 2 + sum(1 if n == 2 * k else 2 for n, k in exact_cells(V))


# ----------------------------------------------------------------------------
def build_s5(bound):
    """Enumerate all C(n,k), k >= 5, 2k <= n, value <= bound.
    Returns (set_of_values, hits, per_k_counts).  hits are (value, tag)
    for collisions inside k>=5 and for triangular members found here."""
    hits = []
    counts = {}
    s5 = set()
    k = 5
    while comb(2 * k, k) <= bound:
        n = 2 * k
        c = comb(n, k)
        cnt = 0
        while c <= bound:
            t = 8 * c + 1
            if T1[t % M1] and T2[t % M2]:
                s = isqrt(t)
                if s * s == t:
                    hits.append((c, "tri@k=%d" % k))
            if c in s5:
                hits.append((c, "s5-collision@k=%d" % k))
            s5.add(c)
            cnt += 1
            # C(n+1,k) = C(n,k) * (n+1) // (n+1-k)
            n += 1
            c = c * n // (n - k)
        counts[k] = cnt
        k += 1
    return s5, hits, counts


# ----------------------------------------------------------------------------
def scan_chunk(widx, vlo, vhi, bound, outdir, resume):
    """Stream-merge the k=3 and k=4 sequences over values in [vlo, vhi).
    Tests every value for: triangular (k=2), membership in S5 (k>=5),
    and k3 == k4 equality.  Checkpoints (n3, n4, hits) atomically.
    Returns list of (value, tag) hits."""
    s5 = G["s5"]
    ckpath = os.path.join(outdir, "ck_b%s_w%d.json" % (G["btag"], widx))

    n3 = first_n_at_least(vlo, 3)
    n4 = first_n_at_least(vlo, 4)
    hits = []
    if resume and os.path.exists(ckpath):
        try:
            with open(ckpath) as f:
                ck = json.load(f)
            if ck.get("vlo") == str(vlo) and ck.get("vhi") == str(vhi):
                n3, n4 = ck["n3"], ck["n4"]
                hits = [(int(v), tag) for v, tag in ck["hits"]]
                if ck.get("done"):
                    return hits
        except (ValueError, KeyError):
            pass

    c3, d3 = comb(n3, 3), comb(n3, 2)          # c3=C(n3,3), d3=C(n3,2)
    c4, e4, f4 = comb(n4, 4), comb(n4, 3), comb(n4, 2)

    t1, t2 = T1, T2                             # local refs for speed
    last_ck = time.time()

    def save(done=False):
        tmp = ckpath + ".tmp"
        with open(tmp, "w") as f:
            json.dump({"vlo": str(vlo), "vhi": str(vhi), "n3": n3, "n4": n4,
                       "done": done,
                       "hits": [[str(v), tag] for v, tag in hits]}, f)
        os.replace(tmp, ckpath)

    while True:
        lim = c4 if c4 < vhi else vhi
        # ---- hot loop: k=3 lane -------------------------------------------
        while c3 < lim:
            t = 8 * c3 + 1
            if t1[t % M1] and t2[t % M2]:
                s = isqrt(t)
                if s * s == t:
                    hits.append((c3, "tri@k=3"))
            if c3 in s5:
                hits.append((c3, "s5@k=3"))
            c3 += d3
            d3 += n3
            n3 += 1
        # -------------------------------------------------------------------
        if c4 >= vhi:
            break                                # k=3 lane also reached vhi
        # checkpoint occasionally (before processing c4: state consistent)
        now = time.time()
        if now - last_ck > 30.0:
            save()
            last_ck = now
        if c3 == c4:
            hits.append((c4, "k3=k4"))
        t = 8 * c4 + 1
        if t1[t % M1] and t2[t % M2]:
            s = isqrt(t)
            if s * s == t:
                hits.append((c4, "tri@k=4"))
        if c4 in s5:
            hits.append((c4, "s5@k=4"))
        c4 += e4
        e4 += f4
        f4 += n4
        n4 += 1

    save(done=True)
    return hits


def _worker(args):
    return scan_chunk(*args)


# ----------------------------------------------------------------------------
def run_search(bound, procs, outdir, resume, quiet=False):
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    n3max = max_n_for_k(bound, 3)
    if n3max == 0:
        raise SystemExit("bound too small")

    if not quiet:
        print("bound B = %.4g (%d); k=3 stream: n up to %d" % (bound, bound, n3max))
    s5, s5_hits, s5_counts = build_s5(bound)
    n4max = max_n_for_k(bound, 4)
    if not quiet:
        print("k=4 stream: n up to %d (%d values, streamed)" % (n4max, max(0, n4max - 7)))
        print("k>=5 table: %d values in set (%s)  [build: %.1fs]"
              % (len(s5), {k: c for k, c in s5_counts.items() if c}, time.time() - t0))
    G["s5"] = s5
    G["btag"] = "%x" % bound          # short tag for checkpoint filenames

    # partition k=3 n-range into equal chunks; chunk boundaries define value
    # cutoffs so that k=4 work and equality detection stay chunk-local.
    total3 = n3max + 1 - 6
    nchunks = procs if total3 > 4_000_000 else 1
    bnds = [6 + (total3 * i) // nchunks for i in range(nchunks)]
    specs = []
    for i in range(nchunks):
        vlo = comb(bnds[i], 3)
        vhi = comb(bnds[i + 1], 3) if i + 1 < nchunks else bound + 1
        specs.append((i, vlo, vhi, bound, outdir, resume))

    if nchunks == 1:
        all_hits = list(scan_chunk(*specs[0]))
    else:
        import multiprocessing as mp
        ctx = mp.get_context("fork")     # workers inherit G["s5"] copy-on-write
        with ctx.Pool(nchunks) as pool:
            results = pool.map_async(_worker, specs)
            while not results.ready():
                results.wait(60)
                if not quiet and not results.ready():
                    done = 0
                    for i in range(nchunks):
                        p = os.path.join(outdir, "ck_b%s_w%d.json" % (G["btag"], i))
                        try:
                            with open(p) as f:
                                ck = json.load(f)
                            frac = (ck["n3"] - bnds[i]) / max(1, (bnds[i + 1] if i + 1 < nchunks else n3max + 1) - bnds[i])
                            done += min(1.0, frac)
                        except Exception:
                            pass
                    el = time.time() - t0
                    pct = 100.0 * done / nchunks
                    eta = el / max(pct, 1e-9) * (100 - pct)
                    print("  ... %5.1f%% done, %.0fs elapsed, ETA %.0fs" % (pct, el, eta))
            all_hits = [h for sub in results.get() for h in sub]
    all_hits.extend(s5_hits)
    elapsed = time.time() - t0

    # ---- post-process: exact verification of every hit --------------------
    tags = {}
    for v, tag in all_hits:
        tags.setdefault(v, set()).add(tag)
    records = []
    for v in sorted(tags):
        cells = exact_cells(v)           # independent math.comb verification
        mult = 2 + sum(1 if n == 2 * k else 2 for n, k in cells)
        assert len(cells) >= 2, "hit %d has <2 canonical cells?!" % v
        records.append({"value": v, "cells": cells, "multiplicity": mult,
                        "tags": sorted(tags[v])})
    census = Counter(r["multiplicity"] for r in records)

    result = {
        "bound": bound,
        "n3max": n3max,
        "n4max": n4max,
        "s5_size": len(s5),
        "s5_counts": s5_counts,
        "procs": nchunks,
        "elapsed_sec": round(elapsed, 1),
        "hit_count": len(records),
        "census_mult_ge5": dict(sorted(census.items())),
        "mult_ge6": [r for r in records if r["multiplicity"] >= 6],
        "all_hits": records,
    }
    return result


# ----------------------------------------------------------------------------
def brute_force(limit):
    """Independent brute force: full multiplicity of every value <= limit by
    direct row-by-row enumeration with math.comb.  Returns {V: mult} for all
    V <= limit with mult >= 5."""
    canon = Counter()      # V -> sum of (1 if central else 2)
    n = 4
    while comb(n, 2) <= limit:
        for k in range(2, n // 2 + 1):
            v = comb(n, k)
            if v > limit:
                break
            canon[v] += 1 if n == 2 * k else 2
        n += 1
    return {v: 2 + c for v, c in canon.items() if 2 + c >= 5}


def selftest():
    """Validate against known facts at bound 1e8 plus brute-force cross-check."""
    print("=== selftest: search at B=1e8 ===")
    res = run_search(10 ** 8, procs=1, outdir="/tmp/singmaster_selftest",
                     resume=False, quiet=True)
    found = {r["value"]: r for r in res["all_hits"]}

    v = found.get(3003)
    assert v, "3003 not found"
    assert v["multiplicity"] == 8, v
    assert set(map(tuple, v["cells"])) == {(78, 2), (15, 5), (14, 6)}, v
    print("3003: multiplicity 8, cells (78,2),(15,5),(14,6)  OK")

    for known in (120, 210, 1540, 7140, 11628, 24310):
        v = found.get(known)
        assert v and v["multiplicity"] == 6, (known, v)
        print("%d: multiplicity 6 %s  OK" % (known, v["cells"]))

    # brute-force cross-check at 1e7: same mult>=5 sets, same multiplicities
    bf = brute_force(10 ** 7)
    pipe = {r["value"]: r["multiplicity"] for r in res["all_hits"]
            if r["value"] <= 10 ** 7}
    assert bf == pipe, (set(bf) ^ set(pipe),
                        {v: (bf.get(v), pipe.get(v)) for v in set(bf) | set(pipe)
                         if bf.get(v) != pipe.get(v)})
    print("brute force <= 1e7: %d values with mult>=5, identical to pipeline  OK"
          % len(bf))
    print("census at 1e8 (mult>=5): %s" % res["census_mult_ge5"])
    print("selftest PASSED")


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bound", type=float, default=None,
                    help="search bound B on the value C(n,k), e.g. 2.5e29")
    ap.add_argument("--procs", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--outdir", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "attempts", "singmaster", "data"))
    ap.add_argument("--no-resume", action="store_true",
                    help="ignore existing checkpoints")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--brute", type=float, default=None,
                    help="run independent brute force up to this limit and print")
    args = ap.parse_args()

    if args.selftest:
        selftest()
        return
    if args.brute:
        bf = brute_force(int(args.brute))
        for v in sorted(bf):
            print(v, bf[v], exact_cells(v))
        return
    if not args.bound:
        ap.error("--bound required (or --selftest / --brute)")

    bound = int(args.bound)
    outdir = os.path.abspath(args.outdir)
    res = run_search(bound, args.procs, outdir, resume=not args.no_resume)

    respath = os.path.join(outdir, "results_b%s.json" % ("%x" % bound))
    with open(respath, "w") as f:
        json.dump(res, f, indent=1, default=str)
    print("\n=== results (%.0fs) ===" % res["elapsed_sec"])
    print("values <= %.4g with multiplicity >= 5: %d" % (bound, res["hit_count"]))
    print("census: %s" % res["census_mult_ge5"])
    print("multiplicity >= 6:")
    for r in res["mult_ge6"]:
        print("  mult %d: %d  cells %s" % (r["multiplicity"], r["value"], r["cells"]))
    best = max((r["multiplicity"] for r in res["all_hits"]), default=0)
    if best > 8:
        print("*** MULTIPLICITY > 8 FOUND -- UNVERIFIED, needs adversarial review ***")
    print("written: %s" % respath)


if __name__ == "__main__":
    main()
