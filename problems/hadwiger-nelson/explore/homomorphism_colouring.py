"""Cap a whole level of the ring with one periodic colouring.

H_k = 3^-k (E + rho E) is, as a unit-distance graph, the Cayley graph on the
rank-4 group E + rho E whose generators are the vectors of squared length 9^k.
Any group homomorphism h : E + rho E -> A into a finite-or-not group, together
with a proper colouring of the Cayley graph on A with the *images* of those
generators, pulls back to a proper colouring of the whole infinite graph -- so
long as no generator lands on the identity.

The homomorphisms tried here are h(z, w) = alpha z + beta w into E itself.  A
proper 4-colouring of the resulting Cayley graph on E is then found as a
periodic colouring of a quotient E / L, which is a finite graph.

One success is a *proof* that every finite unit-distance graph on H_k is
4-colourable, i.e. that no construction at that level can reach five colours.
Failure is only evidence about the homomorphisms tried.

Run:
  python problems/hadwiger-nelson/explore/homomorphism_colouring.py --levels 0,1,2
"""

from __future__ import annotations

import argparse
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
from rho_group import RhoGroup, ball, edges_from_unit_vectors  # noqa: E402
from unit_distance import eisenstein_norm  # noqa: E402

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))


def eis_mul(p, q):
    """(a + b zeta)(c + d zeta) with zeta^2 = zeta - 1."""
    a, b = p
    c, d = q
    # (a + b z)(c + d z) = ac + (ad + bc) z + bd z^2 = (ac - bd) + (ad + bc + bd) z
    return (a * c - b * d, a * d + b * c + b * d)


def image_generators(gens, alpha, beta):
    """h(z, w) = alpha z + beta w applied to each generator."""
    out = []
    for z, w in gens:
        out.append(tuple(x + y for x, y in zip(eis_mul(alpha, z), eis_mul(beta, w))))
    return out


def quotient_colouring(image_gens, k, P, Q, shear):
    """Proper k-colouring of the Cayley graph on E / L, L = <(P,0), (shear,Q)>.

    Returns (n, edges, colouring) or None. A generator congruent to 0 is a
    self-loop and rules the quotient out.
    """
    def red(a, b):
        bb = b % Q
        aa = (a - ((b - bb) // Q) * shear) % P
        return (aa, bb)

    n = P * Q
    edges = set()
    for a in range(P):
        for b in range(Q):
            i = red(a, b)[0] * Q + red(a, b)[1]
            for da, db in image_gens:
                r = red(a + da, b + db)
                j = r[0] * Q + r[1]
                if i == j:
                    return None
                edges.add((min(i, j), max(i, j)))
    edges = sorted(edges)
    col = F.k_colourable(n, edges, k)
    if col is None or not F.verify(n, edges, col, k):
        return None
    return n, edges, col


def verify_on_ball(group, level, gens, alpha, beta, P, Q, shear, colouring,
                   depth=3):
    """Independent check: colour an actual ball by the rule and test every edge.

    This does not reuse the quotient graph at all -- it walks the real
    unit-distance graph, applies h and the periodic colouring to each vertex,
    and checks the colouring edge by edge.
    """
    V = ball(group, gens, depth)
    E = edges_from_unit_vectors(V, gens)

    def red(a, b):
        bb = b % Q
        aa = (a - ((b - bb) // Q) * shear) % P
        return (aa, bb)

    cols = []
    for c in V:
        z, w = c[0], c[1]
        h = tuple(x + y for x, y in zip(eis_mul(alpha, z), eis_mul(beta, w)))
        r = red(*h)
        cols.append(colouring[r[0] * Q + r[1]])
    bad = [(a, b) for a, b in E if cols[a] == cols[b]]
    return {
        "depth": depth, "vertices": len(V), "edges": len(E),
        "monochromatic_edges": len(bad),
        "colours_used": len(set(cols)),
        "proper": not bad,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", default="0,1,2")
    ap.add_argument("--k", type=int, default=4)
    ap.add_argument("--max-period", type=int, default=14)
    ap.add_argument("--alphas", default="1,0")
    ap.add_argument("--betas", default="-1,0")
    ap.add_argument("--verify-depth", type=int, default=3)
    ap.add_argument("--out", default=os.path.join(OUT, "homomorphism_colouring.json"))
    args = ap.parse_args()

    group = RhoGroup(2)
    alpha = tuple(int(x) for x in args.alphas.split(","))
    beta = tuple(int(x) for x in args.betas.split(","))

    results = []
    for level in [int(x) for x in args.levels.split(",")]:
        gens = group.vectors_of_norm(9 ** level)
        img = image_generators(gens, alpha, beta)
        zero_hits = [g for g in img if g == (0, 0)]
        norms = sorted({eisenstein_norm(*g) for g in img})
        print(f"level {level}: {len(gens)} unit vectors -> image generators with "
              f"norms {norms}; {len(zero_hits)} map to zero")
        rec = {"level": level, "alpha": list(alpha), "beta": list(beta),
               "generators": len(gens), "image_norms": norms,
               "generators_mapping_to_zero": len(zero_hits)}
        if zero_hits:
            rec["capped"] = False
            rec["reason"] = "homomorphism kills a generator"
            results.append(rec)
            print("  -> this homomorphism cannot work at this level\n")
            continue

        found = None
        t0 = time.time()
        for P in range(2, args.max_period + 1):
            for Q in range(2, args.max_period + 1):
                for shear in range(Q):
                    r = quotient_colouring(img, args.k, P, Q, shear)
                    if r is not None:
                        found = (P, Q, shear, r)
                        break
                if found:
                    break
            if found:
                break
        if not found:
            rec["capped"] = None
            rec["reason"] = (f"no periodic {args.k}-colouring with period up to "
                             f"{args.max_period}")
            print(f"  -> none found up to period {args.max_period} "
                  f"({time.time() - t0:.1f}s)\n")
            results.append(rec)
            continue

        P, Q, shear, (n, edges, col) = found
        chk = verify_on_ball(group, level, gens, alpha, beta, P, Q, shear, col,
                             depth=args.verify_depth)
        rec.update({
            "capped": True,
            "period": [P, Q, shear], "quotient_vertices": n,
            "quotient_edges": len(edges),
            "colouring": col,
            "ball_check": chk,
            "seconds": round(time.time() - t0, 1),
        })
        results.append(rec)
        print(f"  -> PERIODIC {args.k}-COLOURING with period ({P},{Q},shear={shear}), "
              f"quotient {n} vertices")
        print(f"     independent ball check at depth {chk['depth']}: "
              f"{chk['vertices']} vertices, {chk['edges']} edges, "
              f"{chk['monochromatic_edges']} monochromatic edges, "
              f"proper={chk['proper']}\n")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
