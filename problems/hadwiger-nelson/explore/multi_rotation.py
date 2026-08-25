"""Unions of many rotated copies of one lattice patch.

Two extra degrees of freedom beyond ``rotation_census``:

  * the rotation centre.  Rotating the patch about a lattice point c gives
    c + rho(P - c), a translate of rho(P) by c - rho(c); since rho(P) is not
    lattice-periodic, different centres give genuinely different interfaces.
  * several different rotations rho_n at once, and their inverses.

Each copy is specified as "n:cx,cy:e" meaning rotate about the lattice point
(cx, cy) by rho_n^e.  The base patch is always copy 0.

Run:
  python explore/multi_rotation.py --patch 13 --copies "3:0,0:1 3:1,0:1 3:0,1:1"
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
    """Locate harness/hadwiger-nelson by walking up from this file."""
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

import colouring as C  # noqa: E402
from exact_field import field_for  # noqa: E402
from lattice import (  # noqa: E402
    eisenstein_patch,
    rotation_for_norm,
)
from unit_distance import (  # noqa: E402
    Point,
    build_graph,
    certify,
    eisenstein_norm,
    eisenstein_point,
)

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))


def parse_copy(spec):
    n, centre, power = spec.split(":")
    cx, cy = centre.split(",")
    return int(n), (int(cx), int(cy)), int(power)


def build(patch_norm, copy_specs):
    ns = sorted({n for n, _, _ in copy_specs})
    field = field_for(3, *[4 * n - 1 for n in ns])
    rots = {n: rotation_for_norm(n, field) for n in ns}
    patch = eisenstein_patch(patch_norm)
    base = [eisenstein_point(field, a, b) for a, b in patch]

    points, tags, seen = [], [], {}
    for p, ab in zip(base, patch):
        seen[p] = len(points)
        points.append(p)
        tags.append((0, ab))
    for idx, (n, c, e) in enumerate(copy_specs, start=1):
        centre = eisenstein_point(field, *c)
        rot = rots[n]
        if e < 0:  # inverse rotation: same cosine, opposite sine
            rot = type(rot)(field, rot.c, -rot.s, label=f"n={n} inverse")
        for p, ab in zip(base, patch):
            q = Point(p.x - centre.x, p.y - centre.y)
            for _ in range(abs(e)):
                q = rot(q)
            q = Point(q.x + centre.x, q.y + centre.y)
            if q in seen:
                continue
            seen[q] = len(points)
            points.append(q)
            tags.append((idx, ab))
    return field, points, tags


def analyse(field, points, tags, label, certify_full=True, core_limit=22):
    t0 = time.time()
    edges, stats = build_graph(points)
    if certify_full:
        rep = certify(points, edges)
        if not rep["ok"]:
            raise AssertionError(f"certificate failed for {label}: {rep}")
    cross = sum(1 for a, b in edges if tags[a][0] != tags[b][0])
    chi = C.chromatic_number(len(points), edges, lo=1, hi=6)
    col = C.k_colouring(len(points), edges, chi)
    row = {
        "label": label,
        "field": repr(field),
        "vertices": len(points),
        "edges": len(edges),
        "cross_copy_edges": cross,
        "chi": chi,
        "colouring_verified": C.verify_colouring(len(points), edges, col),
        "seconds": round(time.time() - t0, 2),
    }
    if chi >= 4:
        kept, sub = C.vertex_critical_subgraph(len(points), edges, chi - 1)
        row["core_vertices"] = len(kept)
        row["core_edges"] = len(sub)
        row["core_tags"] = [list(tags[i]) for i in kept]
        if len(kept) <= core_limit:
            row["core_ie_not_k_minus_1_colourable"] = not C.is_k_colourable_ie(
                len(kept), sub, chi - 1)
    row["seconds"] = round(time.time() - t0, 2)
    return row


def sweep_centres(patch_norm, n, max_centre_norm, out_path, powers=(1,), ns=None):
    """Base patch plus rho_n about every lattice centre of small norm.

    Copies are added one at a time so the chi-versus-size curve is visible.
    """
    centres = [(a, b) for a in range(-4, 5) for b in range(-4, 5)
               if eisenstein_norm(a, b) <= max_centre_norm]
    centres.sort(key=lambda ab: (eisenstein_norm(*ab), ab))
    centres = [(nn, c, e) for c in centres for nn in (ns or [n]) for e in powers]
    rows = []
    specs = []
    hdr = f"{'copies':>7} {'|V|':>5} {'|E|':>6} {'cross':>6} {'chi':>4} {'core':>5} {'s':>7}"
    print(hdr)
    print("-" * len(hdr))
    for spec in centres:
        specs.append(spec)
        field, points, tags = build(patch_norm, specs)
        row = analyse(field, points, tags,
                      f"n={n} patch<={patch_norm} copies={len(specs)}")
        row.update({"n": n, "patch_norm": patch_norm,
                    "copy_specs": [[s[0], list(s[1]), s[2]] for s in specs]})
        rows.append(row)
        print(f"{len(specs):>7} {row['vertices']:>5} {row['edges']:>6} "
              f"{row['cross_copy_edges']:>6} {row['chi']:>4} "
              f"{row.get('core_vertices', '-'):>5} {row['seconds']:>7}")
        with open(out_path, "w") as fh:
            json.dump(rows, fh, indent=2)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patch", type=int, default=7)
    ap.add_argument("--copies", default="")
    ap.add_argument("--sweep-centres", type=int, default=0,
                    help="if > 0, sweep rotation centres of norm <= this")
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--ns", default="", help="comma list of rotations to mix in the sweep")
    ap.add_argument("--powers", default="1", help="comma list, may include -1")
    ap.add_argument("--out", default=os.path.join(OUT, "multi_rotation.json"))
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    if args.sweep_centres:
        rows = sweep_centres(
            args.patch, args.n, args.sweep_centres, args.out,
            powers=tuple(int(x) for x in args.powers.split(",")),
            ns=[int(x) for x in args.ns.split(",")] if args.ns else None,
        )
        print(f"\nmax chi seen: {max(r['chi'] for r in rows)}")
        print(f"wrote {args.out}")
        return

    specs = [parse_copy(s) for s in args.copies.split()]
    field, points, tags = build(args.patch, specs)
    row = analyse(field, points, tags, args.copies)
    print(json.dumps({k: v for k, v in row.items() if k != "core_tags"}, indent=2))
    with open(args.out, "w") as fh:
        json.dump([row], fh, indent=2)


if __name__ == "__main__":
    main()
