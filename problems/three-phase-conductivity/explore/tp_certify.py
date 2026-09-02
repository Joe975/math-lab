#!/usr/bin/env python3
"""Exact certification of a screened three-phase laminate tree.

`tp_search.py` (float, Nelder-Mead) writes candidate laminate trees to
`problems/three-phase-conductivity/data/screen/*.json` whose node fractions
are floating-point strings that only *approximately* hit the target phase
fractions and only *approximately* isotropic tensors.  This tool turns one
screened tree into an EXACT record: exact-in-Q node fractions that hit the
target phase fractions *exactly*, the resulting exact effective tensor, and
a certified rational enclosure of the isotropy residual.

Method (exact fraction solve)
------------------------------
`fractions_of(tree)` (see harness/laminate.py) is, as a function of the
internal-node fractions m_1..m_r, MULTILINEAR: each m_i occurs to power 0 or
1 in every monomial (a leaf's contribution is the product of the fraction
weights along its root path, and a path visits each node at most once).  With
3 phases there are only 2 independent constraints (the third fraction is
`1 - f1 - f2`), so for rank r we have r-2 "extra" degrees of freedom.  We fix
r-2 of the node fractions to the rationalization of the screened float
(`Fraction(x).limit_denominator(N)`), and solve the remaining 2 EXACTLY.

The 2 solved-for nodes are chosen as a (parent P, cherry C) pair, where a
"cherry" is an internal node both of whose children are leaves (every binary
tree with >=1 internal node has one -- the deepest node must be a cherry).
Writing a = m_P, b = m_C, and fixing every other node's fraction to a
rational constant, `fractions_of(tree)` becomes, in each phase component, an
affine-in-each-variable polynomial of the special form

    A + B*a + C*b + D*a*b

Because C's two leaves cover only 2 of the 3 phases, the THIRD phase Z (not
a leaf of C) never depends on `b` -- its subtree contribution beneath C is
constant in b, so its coefficients C and D vanish identically for phase Z:
the equation for phase Z is affine in `a` ALONE.  Solve that for `a`.  Then
substitute `a` into the equation for one of C's own leaf phases X, which
(now that `a` is a number) is affine in `b` alone.  Solve that for `b`.  The
third phase's equation is then satisfied automatically, because the three
phase fractions sum to 1 identically in a, b (a polynomial identity, true
for any m's, not just the solution) and the target fractions were required
to sum to 1 too; we still check it as a consistency assertion.

If no (parent, cherry) pair yields a valid solve (both fractions in (0, 1),
the structural degeneracy checks pass) -- e.g. a cherry whose two leaves are
the same phase, in which case that node's own fraction is meaningless (see
`laminate.py`'s comment "laminating a material with itself changes nothing")
-- the tool falls back to rationalizing every node's fraction directly and
reports the exact residual error between the resulting phase fractions and
the target, without claiming an exact fraction match.

Isotropy residual
------------------
The residual lambda_max - lambda_min of tensor [[a, b], [b, c]] is
sqrt((a - c)^2 + 4 b^2), generally irrational.  We certify a rational
enclosure [lo, hi] with hi - lo < tol (default 1e-12) by bisecting the
"is mid^2 < D" PSD-style test on D = (a - c)^2 + 4 b^2 (equivalent to
`laminate.eig_range_within` narrowed around the tensor's own mean), which
gives an exact, checkable upper bound `hi` on the residual.

Records and sidecars
---------------------
Each certified record is written in exactly the format `laminate.py --json`
produces (tree, phases, tensor, fractions, wiener, hs, all as exact strings)
to `problems/three-phase-conductivity/data/certified/<name>.json`, verified
against `harness/three-phase-conductivity/verify_laminate.py`, plus a
`<name>.meta.json` sidecar with source, denominators, eigenvalue enclosure,
isotropy residual bound, HS bounds and the attainment gap.

Usage:
    python tp_certify.py --all [--denom 1000000] [--tol 1e-12]
    python tp_certify.py SCREEN.json [--index 0] [--out NAME]
    python tp_certify.py --selftest

Standard library only.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import subprocess
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
PROBLEM_DIR = os.path.dirname(HERE)
REPO_ROOT = os.path.dirname(os.path.dirname(PROBLEM_DIR))
HARNESS_DIR = os.path.join(REPO_ROOT, "harness", "three-phase-conductivity")
SCREEN_DIR = os.path.join(PROBLEM_DIR, "data", "screen")
CERTIFIED_DIR = os.path.join(PROBLEM_DIR, "data", "certified")

sys.path.insert(0, HARNESS_DIR)
import laminate  # noqa: E402

DEFAULT_PHASES = {"p1": Fr(1), "p2": Fr(2), "p3": Fr(5)}
ALL_PHASES = {"p1", "p2", "p3"}


# ---------------------------------------------------------------------------
# tiny 2-variable multilinear polynomial: c00 + c10*a + c01*b + c11*a*b
# ---------------------------------------------------------------------------


class Bi:
    __slots__ = ("c00", "c10", "c01", "c11")

    def __init__(self, c00=Fr(0), c10=Fr(0), c01=Fr(0), c11=Fr(0)):
        self.c00, self.c10, self.c01, self.c11 = c00, c10, c01, c11

    def __repr__(self):
        return f"Bi({self.c00},{self.c10},{self.c01},{self.c11})"


def bi_add(x: Bi, y: Bi) -> Bi:
    return Bi(x.c00 + y.c00, x.c10 + y.c10, x.c01 + y.c01, x.c11 + y.c11)


def bi_scale(x: Bi, s: Fr) -> Bi:
    return Bi(x.c00 * s, x.c10 * s, x.c01 * s, x.c11 * s)


def bi_mul_var(x: Bi, var: str) -> Bi:
    """Multiply x (which must not already depend on `var`) by the symbol var."""
    if var == "a":
        assert x.c10 == 0 and x.c11 == 0, "polynomial already depends on a"
        return Bi(Fr(0), x.c00, Fr(0), x.c01)
    assert var == "b"
    assert x.c01 == 0 and x.c11 == 0, "polynomial already depends on b"
    return Bi(Fr(0), Fr(0), x.c00, x.c10)


def combine(left: dict[str, Bi], right: dict[str, Bi], weight) -> dict[str, Bi]:
    """weight*left + (1-weight)*right, weight a Fraction or ('var','a'/'b')."""
    out: dict[str, Bi] = {}
    for p in ALL_PHASES:
        L = left.get(p, Bi())
        R = right.get(p, Bi())
        if isinstance(weight, Fr):
            out[p] = bi_add(bi_scale(L, weight), bi_scale(R, 1 - weight))
        else:
            _, var = weight
            diff = Bi(L.c00 - R.c00, L.c10 - R.c10, L.c01 - R.c01, L.c11 - R.c11)
            out[p] = bi_add(R, bi_mul_var(diff, var))
    return out


def build_bi(tree: dict, var_of: dict[int, str], fixed_of: dict[int, Fr]) -> dict[str, Bi]:
    if "phase" in tree:
        return {p: (Bi(Fr(1)) if p == tree["phase"] else Bi()) for p in ALL_PHASES}
    L = build_bi(tree["layers"][0], var_of, fixed_of)
    R = build_bi(tree["layers"][1], var_of, fixed_of)
    nid = id(tree)
    if nid in var_of:
        weight = ("var", var_of[nid])
    else:
        weight = fixed_of[nid]
    return combine(L, R, weight)


# ---------------------------------------------------------------------------
# tree structure helpers
# ---------------------------------------------------------------------------


def internal_nodes(tree: dict) -> list[dict]:
    if "phase" in tree:
        return []
    return [tree] + internal_nodes(tree["layers"][0]) + internal_nodes(tree["layers"][1])


def parent_map(tree: dict, parent=None, out=None) -> dict[int, dict | None]:
    if out is None:
        out = {}
    if "phase" not in tree:
        out[id(tree)] = parent
        parent_map(tree["layers"][0], tree, out)
        parent_map(tree["layers"][1], tree, out)
    return out


def is_cherry(node: dict) -> bool:
    return "phase" in node["layers"][0] and "phase" in node["layers"][1]


# ---------------------------------------------------------------------------
# exact fraction solve
# ---------------------------------------------------------------------------


def rationalize(x, denom: int) -> Fr:
    return Fr(str(x)).limit_denominator(denom)


def solve_exact_fractions(tree: dict, target: dict[str, Fr], denom: int) -> tuple[dict[int, Fr], str] | None:
    """Return (fraction-by-node-id, method-note) with fractions_of(tree) == target exactly, or None."""
    nodes = internal_nodes(tree)
    if not nodes:
        return None
    if len(nodes) == 1:
        # rank 1: single node, two leaves, one free constraint (target of one leaf phase).
        n = nodes[0]
        leaves = (n["layers"][0]["phase"], n["layers"][1]["phase"])
        if leaves[0] == leaves[1]:
            return None
        m = target[leaves[0]]
        if not (0 < m < 1):
            return None
        return {id(n): m}, "rank-1 direct solve"

    pmap = parent_map(tree)

    def ancestors(n: dict) -> list[dict]:
        out = []
        p = pmap[id(n)]
        while p is not None:
            out.append(p)
            p = pmap[id(p)]
        return out

    cherries = [n for n in nodes if is_cherry(n)]
    for C in cherries:
        leaves_C = (C["layers"][0]["phase"], C["layers"][1]["phase"])
        if leaves_C[0] == leaves_C[1]:
            continue  # degenerate: this node's own fraction is meaningless
        for P in ancestors(C):  # nearest ancestor first, then further out
            others = [n for n in nodes if n is not C and n is not P]
            var_of = {id(P): "a", id(C): "b"}
            fixed_of = {id(n): rationalize(n["fraction"], denom) for n in others}
            bidict = build_bi(tree, var_of, fixed_of)

            Z = next(iter(ALL_PHASES - set(leaves_C)))
            bz = bidict[Z]
            if bz.c01 != 0 or bz.c11 != 0 or bz.c10 == 0:
                continue  # Z doesn't (yet) depend on 'a' through this ancestor -- try a farther one
            a_val = (target[Z] - bz.c00) / bz.c10
            if not (0 < a_val < 1):
                continue

            X, W = leaves_C
            bx = bidict[X]
            coeff_b = bx.c01 + bx.c11 * a_val
            const_ = bx.c00 + bx.c10 * a_val
            if coeff_b == 0:
                continue
            b_val = (target[X] - const_) / coeff_b
            if not (0 < b_val < 1):
                continue

            bw = bidict[W]
            w_val = (bw.c00 + bw.c10 * a_val) + (bw.c01 + bw.c11 * a_val) * b_val
            assert w_val == target[W], "sum-to-1 identity violated -- bug in the Bi solve"

            frac_by_id = dict(fixed_of)
            frac_by_id[id(P)] = a_val
            frac_by_id[id(C)] = b_val
            note = (
                f"solved cherry {leaves_C} via ancestor: phase {Z} (ancestor's fraction) then "
                f"phase {X} (cherry's fraction); {len(others)} other node(s) rationalized to "
                f"denominator<={denom}"
            )
            return frac_by_id, note
    return None


def apply_fractions(tree: dict, frac_by_id: dict[int, Fr]) -> None:
    for n in internal_nodes(tree):
        n["fraction"] = str(frac_by_id[id(n)])


def fallback_rationalize_all(tree: dict, denom: int) -> None:
    for n in internal_nodes(tree):
        n["fraction"] = str(rationalize(n["fraction"], denom))


# ---------------------------------------------------------------------------
# exact isotropy-residual enclosure
# ---------------------------------------------------------------------------


def residual_enclosure(t: laminate.Tensor, tol: Fr) -> tuple[Fr, Fr]:
    """Rational [lo, hi] enclosing lambda_max - lambda_min = sqrt((a-c)^2+4b^2)."""
    a, b, c = t
    D = (a - c) * (a - c) + 4 * b * b
    if D == 0:
        return Fr(0), Fr(0)
    lo, hi = Fr(0), max(Fr(1), D)
    while hi * hi < D:
        hi *= 2
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if mid * mid < D:
            lo = mid
        else:
            hi = mid
    return lo, hi


def eig_bounds(t: laminate.Tensor, tol: Fr) -> tuple[Fr, Fr, Fr, Fr]:
    """(lambda_min_lo, lambda_min_hi, lambda_max_lo, lambda_max_hi)."""
    a, b, c = t
    res_lo, res_hi = residual_enclosure(t, tol)
    mean = (a + c) / 2
    return mean - res_hi / 2, mean - res_lo / 2, mean + res_lo / 2, mean + res_hi / 2


# ---------------------------------------------------------------------------
# certify one screen file
# ---------------------------------------------------------------------------


def target_from_name(name: str) -> tuple[list[Fr], str]:
    """Fallback: parse target fractions and side from a rN_a-b_c-d_e-f_side name."""
    parts = name.split("_")
    side = parts[-1]
    fracs = [Fr(int(p.split("-")[0]), int(p.split("-")[1])) for p in parts[1:-1] if "-" in p]
    return fracs, side


def target_from_log(log_path: str) -> tuple[list[Fr], list[Fr], str] | None:
    if not os.path.exists(log_path):
        return None
    with open(log_path, encoding="utf-8") as fh:
        first = fh.readline()
    # f=1/8,1/8,6/8 sigma=1,2,5 rank=3 side=lower topologies=11520
    fields = dict(tok.split("=", 1) for tok in first.split())
    f = [Fr(x) for x in fields["f"].split(",")]
    sigma = [Fr(x) for x in fields["sigma"].split(",")]
    side = fields["side"]
    return f, sigma, side


def certify_file(screen_path: str, index: int, denom: int, tol: Fr, out_name: str | None) -> dict:
    with open(screen_path, encoding="utf-8") as fh:
        candidates = json.load(fh)
    entry = candidates[index]
    tree = entry["tree"]

    base = os.path.splitext(os.path.basename(screen_path))[0]
    log_path = os.path.join(os.path.dirname(screen_path), base + ".log")
    parsed = target_from_log(log_path)
    if parsed is not None:
        f_list, sigma_list, side = parsed
    else:
        f_list, side = target_from_name(base)
        sigma_list = [DEFAULT_PHASES["p1"], DEFAULT_PHASES["p2"], DEFAULT_PHASES["p3"]]

    phases = {"p1": sigma_list[0], "p2": sigma_list[1], "p3": sigma_list[2]}
    target = {"p1": f_list[0], "p2": f_list[1], "p3": f_list[2]}
    assert sum(target.values()) == 1, f"target fractions must sum to 1: {target}"

    solved = solve_exact_fractions(tree, target, denom)
    if solved is not None:
        frac_by_id, method_note = solved
        apply_fractions(tree, frac_by_id)
        exact_fraction_match = True
    else:
        fallback_rationalize_all(tree, denom)
        method_note = f"fallback: rationalized every node fraction to denominator<={denom} (no exact solve found)"
        exact_fraction_match = False

    t = laminate.effective(tree, phases)
    fr = laminate.fractions_of(tree)
    names = sorted(fr)
    fracs = [fr[n] for n in names]
    sigs = [phases[n] for n in names]
    wiener = laminate.wiener_bounds(fracs, sigs)
    hs = laminate.hs_bounds(fracs, sigs)
    duality = laminate.keller_check(tree, phases)
    assert duality, "Keller-Dykhne identity failed on the exact tree"

    fraction_error = {p: fr[p] - target[p] for p in ALL_PHASES}
    if exact_fraction_match:
        assert all(v == 0 for v in fraction_error.values()), fraction_error

    lam_min_lo, lam_min_hi, lam_max_lo, lam_max_hi = eig_bounds(t, tol)
    res_lo, res_hi = residual_enclosure(t, tol)

    if side == "lower":
        gap_lo = lam_max_lo - hs[0]
        gap_hi = lam_max_hi - hs[0]
    else:
        gap_lo = hs[1] - lam_min_hi
        gap_hi = hs[1] - lam_min_lo

    name = out_name or base
    record = {
        "tree": tree,
        "phases": {k: str(v) for k, v in phases.items()},
        "tensor": [str(x) for x in t],
        "fractions": {k: str(v) for k, v in fr.items()},
        "wiener": [str(x) for x in wiener],
        "hs": [str(x) for x in hs],
    }
    meta = {
        "source_screen_file": os.path.relpath(screen_path, PROBLEM_DIR),
        "candidate_index": index,
        "target_fractions": {k: str(v) for k, v in target.items()},
        "side": side,
        "denominator_limit": denom,
        "exact_fraction_match": exact_fraction_match,
        "fraction_error": {k: str(v) for k, v in fraction_error.items()},
        "solve_method": method_note,
        "eigenvalue_enclosure": {
            "lambda_min": [str(lam_min_lo), str(lam_min_hi)],
            "lambda_max": [str(lam_max_lo), str(lam_max_hi)],
            "tol": str(tol),
        },
        "isotropy_residual_bound": {
            "lambda_max_minus_lambda_min": [str(res_lo), str(res_hi)],
        },
        "hs_bounds": {"lo": str(hs[0]), "hi": str(hs[1])},
        "gap_to_hs": {"side": side, "enclosure": [str(gap_lo), str(gap_hi)]},
        "command_line": " ".join(sys.argv),
    }
    return {"name": name, "record": record, "meta": meta}


def write_and_verify(result: dict) -> str:
    os.makedirs(CERTIFIED_DIR, exist_ok=True)
    name = result["name"]
    record_path = os.path.join(CERTIFIED_DIR, name + ".json")
    meta_path = os.path.join(CERTIFIED_DIR, name + ".meta.json")
    with open(record_path, "w", encoding="utf-8") as fh:
        json.dump(result["record"], fh, indent=1)
    with open(meta_path, "w", encoding="utf-8") as fh:
        json.dump(result["meta"], fh, indent=1)

    verify_script = os.path.join(HARNESS_DIR, "verify_laminate.py")
    proc = subprocess.run(
        [sys.executable, verify_script, record_path], capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(f"verify_laminate.py FAILED for {record_path}:\n{proc.stdout}\n{proc.stderr}")
    return proc.stdout.strip()


# ---------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------


def selftest() -> None:
    F = Fr
    s1, s2 = F(1), F(4)
    tree = {
        "fraction": "1/2",
        "normal": [1, 1],
        "layers": [
            {
                "fraction": "1/3",
                "normal": [1, 0],
                "layers": [{"phase": "p1"}, {"phase": "p2"}],
            },
            {
                "fraction": "1/3",
                "normal": [0, 1],
                "layers": [{"phase": "p1"}, {"phase": "p2"}],
            },
        ],
    }
    phases = {"p1": s1, "p2": s2}
    t = laminate.effective(tree, phases)
    assert t == (F(49, 20), F(1, 20), F(49, 20)), t
    assert t[0] == t[2], "hand-built construction must give a == c exactly"

    lo, hi = residual_enclosure(t, F(1, 10**12))
    assert hi - lo < F(1, 10**12)
    exact_residual = F(1, 10)  # sqrt((a-c)^2 + 4b^2) = 2|b| since a == c
    assert lo <= exact_residual <= hi, (lo, exact_residual, hi)

    fr = laminate.fractions_of(tree)
    assert fr == {"p1": F(1, 3), "p2": F(2, 3)}
    assert laminate.keller_check(tree, {"p1": s1, "p2": s2})
    print("selftest: exact hand-built rank-2 tree certified (a == c exactly; "
          f"isotropy residual enclosed in [{lo}, {hi}], contains {exact_residual})")

    # also exercise the general solver on a synthetic 3-phase, rank-3 screened tree
    F3 = {"p1": Fr(1), "p2": Fr(2), "p3": Fr(5)}
    target = {"p1": Fr(1, 8), "p2": Fr(1, 8), "p3": Fr(6, 8)}
    screened_tree = {
        "fraction": "0.18310642055900347",
        "normal": [1, 0],
        "layers": [
            {"fraction": "0.3173368819160453", "normal": [1, 0],
             "layers": [{"phase": "p1"}, {"phase": "p2"}]},
            {"fraction": "0.0818877527312502", "normal": [0, 1],
             "layers": [{"phase": "p1"}, {"phase": "p3"}]},
        ],
    }
    solved = solve_exact_fractions(screened_tree, target, 10**6)
    assert solved is not None
    frac_by_id, note = solved
    apply_fractions(screened_tree, frac_by_id)
    fr3 = laminate.fractions_of(screened_tree)
    assert fr3 == target, fr3
    assert laminate.keller_check(screened_tree, F3)
    print(f"selftest: general cherry-solve reproduced target fractions exactly ({note})")


# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(description="Exact certification of a screened three-phase laminate tree")
    ap.add_argument("screen", nargs="?", help="path to a data/screen/*.json file")
    ap.add_argument("--index", type=int, default=0, help="which candidate in the screen file (default: best, 0)")
    ap.add_argument("--denom", type=int, default=10**6, help="limit_denominator bound for rationalized fractions")
    ap.add_argument("--tol", type=float, default=1e-12, help="eigenvalue-enclosure width target")
    ap.add_argument("--out", help="output record name (default: derived from the screen filename)")
    ap.add_argument("--all", action="store_true", help="certify every file in data/screen/")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        selftest()
        return

    tol = Fr(args.tol).limit_denominator(10**15)

    if args.all:
        files = sorted(glob.glob(os.path.join(SCREEN_DIR, "*.json")))
        if not files:
            print("no screen files found", file=sys.stderr)
            sys.exit(1)
        rows = []
        for path in files:
            result = certify_file(path, args.index, args.denom, tol, None)
            verify_out = write_and_verify(result)
            rows.append((result["name"], result["meta"]))
            print(f"{result['name']}: certified and verified ({verify_out})")
            print(f"  method: {result['meta']['solve_method']}")
            print(f"  isotropy residual bound: {result['meta']['isotropy_residual_bound']}")
            print(f"  HS bounds: {result['meta']['hs_bounds']}  gap: {result['meta']['gap_to_hs']}")
        return

    if not args.screen:
        ap.error("need SCREEN.json (or --all / --selftest)")
    result = certify_file(args.screen, args.index, args.denom, tol, args.out)
    verify_out = write_and_verify(result)
    print(f"{result['name']}: certified and verified ({verify_out})")
    print(json.dumps(result["meta"], indent=1))


if __name__ == "__main__":
    main()
