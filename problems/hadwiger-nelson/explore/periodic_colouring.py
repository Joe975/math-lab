"""Look for a periodic 4-colouring of the whole group, which would cap a family.

The unit-distance graph on G_m is a Cayley graph on the coefficient group Z^(2m)
with the group's unit vectors as generators.  A proper colouring of the finite
quotient Cayley graph on Z^(2m)/L pulls back to a proper colouring of the whole
infinite graph, provided no generator lies in L (which would be a self-loop).

So one 4-colourable loopless quotient is a *proof* that every finite subgraph
of G_m is 4-colourable -- that is, that no construction inside G_m can ever
reach five colours.  A failure to find one is only evidence.

Two quotient families are tried:

* diagonal, L = N.Z^(2m), i.e. all coefficients read modulo N;
* per-power diagonal, L = diag(N_0, N_0, N_1, N_1, ...), which lets the
  lattice and its rotated copies wrap at different rates.

Run:  python problems/hadwiger-nelson/explore/periodic_colouring.py --m 2 --max-n 8
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def _harness_dir():
    d = HERE
    while True:
        cand = os.path.join(d, "harness", "hadwiger-nelson")
        if os.path.isdir(cand):
            return cand
        parent = os.path.dirname(d)
        if parent == d:
            raise RuntimeError("could not locate harness/hadwiger-nelson")
        d = parent


sys.path.insert(0, _harness_dir())
sys.path.insert(0, HERE)

import colouring_fast as F  # noqa: E402
from rho_group import RhoGroup  # noqa: E402

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))


def flatten(coeffs):
    out = []
    for a, b in coeffs:
        out.extend((a, b))
    return tuple(out)


def quotient_graph(unit_vecs, moduli):
    """Cayley graph on prod Z_{moduli[i]}, with the given generators.

    Returns (n, edges, has_loop).  Vertices are mixed-radix encoded.
    """
    gens = [tuple(x % m for x, m in zip(flatten(u), moduli)) for u in unit_vecs]
    has_loop = any(all(g == 0 for g in gen) for gen in gens)
    size = 1
    for m in moduli:
        size *= m

    def encode(v):
        idx = 0
        for x, m in zip(v, moduli):
            idx = idx * m + x
        return idx

    edges = set()
    for combo in itertools.product(*[range(m) for m in moduli]):
        i = encode(combo)
        for gen in gens:
            d = tuple((x + g) % m for x, g, m in zip(combo, gen, moduli))
            j = encode(d)
            if i != j:
                edges.add((min(i, j), max(i, j)))
    return size, sorted(edges), has_loop


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, default=2)
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--max-n", type=int, default=8, help="largest modulus")
    ap.add_argument("--k", type=int, default=4, help="number of colours to try")
    ap.add_argument("--search-patch", type=int, default=12)
    ap.add_argument("--out", default=os.path.join(OUT, "periodic_colouring.json"))
    args = ap.parse_args()

    group = RhoGroup(args.m, n=args.n)
    uv = group.unit_vectors([args.search_patch] * args.m)
    print(f"G_{args.m} (rho from n={args.n}): {len(uv)} unit vectors "
          f"in the coefficient box of norm <= {args.search_patch}")

    rows = []
    hdr = f"{'moduli':>18} {'|V|':>7} {'|E|':>8} {'loop':>5} {'k-col':>6} {'s':>7}"
    print(hdr)
    print("-" * len(hdr))

    # per-power moduli, non-decreasing to avoid testing both (2,3) and (3,2)
    candidates = []
    for combo in itertools.product(range(2, args.max_n + 1), repeat=args.m):
        if list(combo) != sorted(combo):
            continue
        size = 1
        for c in combo:
            size *= c * c
        if size > 60000:
            continue
        candidates.append(combo)
    candidates.sort(key=lambda c: [x * x for x in c].__len__() and
                    __import__("math").prod(x * x for x in c))

    found = None
    for combo in candidates:
        moduli = []
        for c in combo:
            moduli.extend((c, c))
        t0 = time.time()
        n, edges, loop = quotient_graph(uv, moduli)
        if loop:
            rows.append({"moduli": moduli, "vertices": n, "self_loop": True})
            print(f"{str(moduli):>18} {n:>7} {'-':>8} {'yes':>5} {'-':>6} "
                  f"{time.time() - t0:>7.1f}")
            continue
        col = F.k_colourable(n, edges, args.k)
        ok = col is not None
        verified = F.verify(n, edges, col, args.k) if ok else None
        rows.append({
            "moduli": moduli, "vertices": n, "edges": len(edges),
            "self_loop": False, f"{args.k}_colourable": ok,
            "colouring_verified": verified,
            "colouring": col if ok else None,
            "seconds": round(time.time() - t0, 2),
        })
        print(f"{str(moduli):>18} {n:>7} {len(edges):>8} {'no':>5} "
              f"{str(ok):>6} {time.time() - t0:>7.1f}")
        if ok and found is None:
            found = rows[-1]

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump({"m": args.m, "n": args.n, "unit_vectors": len(uv),
                   "k": args.k, "rows": rows}, fh, indent=2)
    print(f"\nwrote {args.out}")
    if found:
        print(f"FOUND a periodic {args.k}-colouring at moduli {found['moduli']}: "
              f"every finite subgraph of G_{args.m} is {args.k}-colourable.")
    else:
        print(f"No periodic {args.k}-colouring found in the quotients tried. "
              f"That is evidence, not proof.")


if __name__ == "__main__":
    main()
