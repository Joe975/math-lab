#!/usr/bin/env python3
"""Rank-5 (and rank-6) AXIS-NORMAL-ONLY laminate screen for the three-phase
isotropic lower bound, extending attempt 001's rank-3/rank-4 census
(problems/three-phase-conductivity/attempts/001-coated-laminate-census.md).

Why axis-only is cheap
-----------------------
With normals restricted to {e1, e2}, every tensor in the tree stays diagonal
(harness/three-phase-conductivity/laminate.py's general 2x2-matrix lamination
formula reduces, for diagonal A, B and normal e1, to a weighted HARMONIC mean
in the e1 (lamination-axis) component and a weighted ARITHMETIC mean in the
e2 (transverse) component -- symmetric for normal e2).  So a whole
axis-normal tree can be evaluated with a pair of floats (t1, t2) per node
instead of the general (a, b, c) tensor + rotation machinery tp_search.py
uses for the four-normal set.  This is the same rule tp_shapes.py already
uses for its closed forms (see its harm()/arith(); this file is the
brute-force sibling of that closed-form analysis, one level up in rank).

Global-flip symmetry
----------------------
Swapping EVERY normal in the tree (e1<->e2) transposes the final tensor
(t1, t2) -> (t2, t1); this changes nothing about an isotropic point (t1==t2)
or the objective (symmetric in t1, t2).  So the root node's normal is fixed
to e1 and only the remaining (rank-1) normals are enumerated, an exact 2x
reduction in the normal-assignment count, not an approximation.

Search space at rank R (Catalan(R) shapes x 3^(R+1) labellings, filtered to
require all three phases present, x 2^(R-1) normals after the root-fix):
rank 5 -> 42 x 540 x 16 = 362,880 topologies; rank 6 -> 132 x 1806 x 32 =
7,629,984 (see --calibrate for whether this is affordable this session).

Nothing here is exact.  Nelder-Mead screens for isotropic-ish, fraction-
correct trees; problems/three-phase-conductivity/explore/tp_certify.py
(unmodified, imported as a module) turns a screened tree into a certified
exact-in-Q record, re-verified by harness/three-phase-conductivity/
verify_laminate.py, written under data/rank5/ instead of data/certified/ so
this attempt's records are distinguishable from attempt 001's rank<=4 ones.

Usage:
    python tp_rank5.py --selftest
    python tp_rank5.py --calibrate --rank 5 --sample 3000
    python tp_rank5.py --f 1/8,1/8,6/8 --sigma 1,2,5 --rank 5 --side lower \
        --checkpoint ../data/rank5/ckpt_1-1-6.json --out ../data/rank5/screen_1-1-6.json
    python tp_rank5.py --certify ../data/rank5/screen_1-1-6.json --name r5_1-8_1-8_6-8_lower
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import random
import sys
import time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
PROBLEM_DIR = os.path.dirname(HERE)
REPO_ROOT = os.path.dirname(os.path.dirname(PROBLEM_DIR))
HARNESS_DIR = os.path.join(REPO_ROOT, "harness", "three-phase-conductivity")
RANK5_DIR = os.path.join(PROBLEM_DIR, "data", "rank5")

sys.path.insert(0, HARNESS_DIR)
sys.path.insert(0, HERE)
import laminate as L  # noqa: E402
import tp_shapes as SH  # noqa: E402
import tp_certify as CERT  # noqa: E402
from tp_search import shapes, internal_nodes  # noqa: E402


# ---------------------------------------------------------------------------
# axis-only scalar evaluation
# ---------------------------------------------------------------------------

def harmf(m, p, q):
    return 1.0 / (m / p + (1 - m) / q)


def arithf(m, p, q):
    return m * p + (1 - m) * q


class AxisTopology:
    """A binary-tree topology restricted to axis normals, evaluated as a pair
    of floats (t1, t2) -- the diagonal entries of the (always-diagonal)
    effective tensor -- instead of the general (a, b, c) + rotation form."""

    def __init__(self, shape, labels, normals):
        self.shape, self.labels, self.normals = shape, labels, normals
        self.nodes = internal_nodes(shape)          # root first (preorder)
        self.node_index = {n: i for i, n in enumerate(self.nodes)}
        self.n_params = len(self.nodes)

    def evaluate(self, params, sig):
        it = iter(params)

        def rec(t):
            if isinstance(t, int):
                s = sig[self.labels[t]]
                fr = [0.0, 0.0, 0.0]
                fr[self.labels[t]] = 1.0
                return (s, s), fr
            m = next(it)
            idx = self.node_index[t]
            (a1, a2), fa = rec(t[0])
            (b1, b2), fb = rec(t[1])
            fr = [m * x + (1 - m) * y for x, y in zip(fa, fb)]
            if self.normals[idx] == 0:      # e1: axis component harmonic
                out = (harmf(m, a1, b1), arithf(m, a2, b2))
            else:                            # e2: axis component harmonic in slot 2
                out = (arithf(m, a1, b1), harmf(m, a2, b2))
            return out, fr
        return rec(self.shape)

    def to_tree(self, params_frac):
        it = iter(params_frac)

        def rec(t):
            if isinstance(t, int):
                return {"phase": f"p{self.labels[t] + 1}"}
            m = next(it)
            idx = self.node_index[t]
            normal = [1, 0] if self.normals[idx] == 0 else [0, 1]
            return {"fraction": str(m), "normal": normal,
                    "layers": [rec(t[0]), rec(t[1])]}
        return rec(self.shape)


def has_degenerate_cherry(shape, labels):
    """True if some internal node's two children are both leaves with the
    SAME phase label (that node's fraction is then meaningless -- laminating
    a phase with itself -- see laminate.py's own comment on this case).
    Pruned to cut wasted search volume, not for correctness."""
    if isinstance(shape, int):
        return False
    left, right = shape
    if isinstance(left, int) and isinstance(right, int) and labels[left] == labels[right]:
        return True
    return has_degenerate_cherry(left, labels) or has_degenerate_cherry(right, labels)


# ---------------------------------------------------------------------------
# optimisation (float Nelder-Mead, mirrors tp_search.py)
# ---------------------------------------------------------------------------

def sigmoid(x):
    return 1 / (1 + math.exp(-x)) if x > -700 else 0.0


def objective(top, x, sig, f_target, sign, W=200.0):
    params = [sigmoid(v) for v in x]
    (t1, t2), fr = top.evaluate(params, sig)
    tr = t1 + t2
    aniso = abs(t1 - t2) / tr
    ferr = sum(abs(p - q) for p, q in zip(fr, f_target))
    return sign * tr / 2 + W * (aniso + ferr), tr / 2, aniso, ferr


def nelder_mead(fn, x0, step=1.0, iters=200, tol=1e-12):
    n = len(x0)
    pts = [list(x0)]
    for i in range(n):
        p = list(x0)
        p[i] += step
        pts.append(p)
    vals = [fn(p)[0] for p in pts]
    for _ in range(iters):
        order = sorted(range(n + 1), key=lambda i: vals[i])
        pts = [pts[i] for i in order]
        vals = [vals[i] for i in order]
        if abs(vals[-1] - vals[0]) < tol:
            break
        cen = [sum(p[j] for p in pts[:-1]) / n for j in range(n)]
        worst = pts[-1]
        refl = [cen[j] + (cen[j] - worst[j]) for j in range(n)]
        fr_ = fn(refl)[0]
        if fr_ < vals[0]:
            exp = [cen[j] + 2 * (cen[j] - worst[j]) for j in range(n)]
            fe = fn(exp)[0]
            if fe < fr_:
                pts[-1], vals[-1] = exp, fe
            else:
                pts[-1], vals[-1] = refl, fr_
        elif fr_ < vals[-2]:
            pts[-1], vals[-1] = refl, fr_
        else:
            con = [cen[j] + 0.5 * (worst[j] - cen[j]) for j in range(n)]
            fc = fn(con)[0]
            if fc < vals[-1]:
                pts[-1], vals[-1] = con, fc
            else:
                best = pts[0]
                for i in range(1, n + 1):
                    pts[i] = [best[j] + 0.5 * (pts[i][j] - best[j]) for j in range(n)]
                    vals[i] = fn(pts[i])[0]
    i = min(range(n + 1), key=lambda i: vals[i])
    return pts[i], vals[i]


# ---------------------------------------------------------------------------
# topology enumeration
# ---------------------------------------------------------------------------

def all_topologies(rank, f_target, prune_degenerate=True):
    n_leaf = rank + 1
    required = set(i for i, f in enumerate(f_target) if f > 0)
    for shape in shapes(rank):
        for labels in itertools.product(range(3), repeat=n_leaf):
            if required - set(labels):
                continue
            if prune_degenerate and has_degenerate_cherry(shape, labels):
                continue
            for rest in itertools.product((0, 1), repeat=rank - 1):
                normals = (0,) + rest   # root normal fixed to e1 (global-flip symmetry)
                yield shape, labels, normals


def count_topologies(rank, f_target, prune_degenerate=True):
    return sum(1 for _ in all_topologies(rank, f_target, prune_degenerate))


# ---------------------------------------------------------------------------
# search with checkpointing
# ---------------------------------------------------------------------------

def load_checkpoint(path):
    if path and os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    return None


def save_checkpoint(path, done, total, best, keep=25):
    if not path:
        return
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".tmp"
    payload = {
        "done": done,
        "total": total,
        "best": [
            {"mean": mean, "aniso": aniso, "ferr": ferr,
             "shape": top.shape, "labels": top.labels, "normals": top.normals,
             "params": params}
            for mean, aniso, ferr, top, params in best[:keep]
        ],
    }
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh)
    os.replace(tmp, path)


def resume_best(payload, rank):
    out = []
    for e in payload.get("best", []):
        shape = tuple_ify(e["shape"])
        labels = tuple(e["labels"])
        normals = tuple(e["normals"])
        top = AxisTopology(shape, labels, normals)
        out.append((e["mean"], e["aniso"], e["ferr"], top, e["params"]))
    return out


def tuple_ify(x):
    if isinstance(x, list):
        return tuple(tuple_ify(e) for e in x)
    return x


def search(f_target, sig, rank, side, restarts, seed, tol_report,
           checkpoint=None, checkpoint_every=5000, time_budget=None,
           resume=None):
    rng = random.Random(seed)
    sign = 1.0 if side == "lower" else -1.0
    best = list(resume) if resume else []
    start_index = resume_done = 0
    if resume is not None:
        pass
    t_start = time.time()
    done = 0
    total = count_topologies(rank, f_target)
    stopped_early = False
    for shape, labels, normals in all_topologies(rank, f_target):
        done += 1
        top = AxisTopology(shape, labels, normals)

        def fn(x, top=top):
            return objective(top, x, sig, f_target, sign)

        for _ in range(restarts):
            x0 = [rng.uniform(-2, 2) for _ in range(top.n_params)]
            x, val = nelder_mead(fn, x0, iters=200)
            x, val = nelder_mead(fn, x, step=0.05, iters=200)
            _, mean, aniso, ferr = fn(x)
            if aniso < tol_report and ferr < tol_report:
                best.append((mean, aniso, ferr, top, [sigmoid(v) for v in x]))
        if done % checkpoint_every == 0:
            best.sort(key=lambda r: sign * r[0])
            save_checkpoint(checkpoint, done, total, best)
            elapsed = time.time() - t_start
            print(f"  ... {done}/{total} ({100*done/total:.1f}%) elapsed={elapsed:.0f}s "
                  f"best={best[0][0] if best else None}", file=sys.stderr)
        if time_budget and (time.time() - t_start) > time_budget:
            stopped_early = True
            break
    best.sort(key=lambda r: sign * r[0])
    save_checkpoint(checkpoint, done, total, best)
    return best, done, total, stopped_early


# ---------------------------------------------------------------------------
# writing screen-format output (compatible with tp_certify.py) + certifying
# ---------------------------------------------------------------------------

def write_screen(out_path, f_target, sigma, rank, side, total_topologies, best, top_n=10):
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    seen = set()
    shown = []
    for mean, aniso, ferr, top, params in best:
        key = (top.shape, top.labels, top.normals)
        if key in seen:
            continue
        seen.add(key)
        shown.append((mean, aniso, ferr, top, params))
        if len(shown) >= top_n:
            break
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump([{"mean": m, "aniso": an, "ferr": fe,
                     "tree": t.to_tree([repr(p) for p in pr])}
                    for m, an, fe, t, pr in shown], fh, indent=1)
    log_path = os.path.splitext(out_path)[0] + ".log"
    with open(log_path, "w", encoding="utf-8") as fh:
        f_str = ",".join(str(Fr(x).limit_denominator(10**6)) for x in f_target)
        s_str = ",".join(str(Fr(x).limit_denominator(10**6)) for x in sigma)
        fh.write(f"f={f_str} sigma={s_str} rank={rank} side={side} "
                  f"topologies={total_topologies}\n")
        fh.write("# axis-normal-only rank5 screen (tp_rank5.py)\n")
    return out_path


def certify_screen(screen_path, out_dir, name, denom=10**6, tol=1e-12, index=0):
    os.makedirs(out_dir, exist_ok=True)
    result = CERT.certify_file(screen_path, index, denom, Fr(tol).limit_denominator(10**15), name)
    record_path = os.path.join(out_dir, result["name"] + ".json")
    meta_path = os.path.join(out_dir, result["name"] + ".meta.json")
    with open(record_path, "w", encoding="utf-8") as fh:
        json.dump(result["record"], fh, indent=1)
    with open(meta_path, "w", encoding="utf-8") as fh:
        json.dump(result["meta"], fh, indent=1)
    import subprocess
    verify_script = os.path.join(HARNESS_DIR, "verify_laminate.py")
    proc = subprocess.run([sys.executable, verify_script, record_path],
                           capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"verify_laminate.py FAILED for {record_path}:\n{proc.stdout}\n{proc.stderr}")
    return record_path, meta_path, proc.stdout.strip(), result["meta"]


# ---------------------------------------------------------------------------
# selftest: reproduces attempt 001's rank-3 and rank-4 numbers, and checks the
# axis-only scalar evaluator against harness.effective() on concrete trees
# ---------------------------------------------------------------------------

def selftest():
    F = Fr
    phases = {"p1": F(1), "p2": F(2), "p3": F(5)}
    sig = [1.0, 2.0, 5.0]

    # --- axis-only scalar evaluator matches harness.effective() exactly on a
    # hand-built rank-3 tree (reuse tp_shapes' own rank-3 shape) ---
    f1, f2, f3 = F(1, 8), F(1, 8), F(6, 8)
    m_top = F(3, 20)
    m_L, m_R = SH.rank3_from_fractions(m_top, f2, f3)
    tree = SH.rank3_tree(m_top, m_L, m_R)
    t_exact = L.effective(tree, phases)
    assert t_exact[1] == 0

    # same tree as an AxisTopology: shape ((0,1),2), labels leaf0=p1,leaf1=p2,leaf2=p3? map
    # L12 = lam(p1,p2,e1) fraction m_L; L13 = lam(p1,p3,e2) fraction m_R; top = lam(L12,L13,e1) m_top
    # shape: node0=(node1,node2), node1=(leaf0,leaf1), node2=(leaf2,leaf3)
    shape = ((0, 1), (2, 3))
    labels = (0, 1, 0, 2)   # p1,p2,p1,p3  (0-indexed phase)
    normals = (0, 0, 1)     # root e1, node1 e1, node2 e2
    top = AxisTopology(shape, labels, normals)
    (t1, t2), fr = top.evaluate([float(m_top), float(m_L), float(m_R)], sig)
    assert abs(t1 - float(t_exact[0])) < 1e-9 and abs(t2 - float(t_exact[2])) < 1e-9, (t1, t2, t_exact)
    assert abs(fr[0] - float(1 - f2 - f3)) < 1e-9 and abs(fr[1] - float(f2)) < 1e-9 and abs(fr[2] - float(f3)) < 1e-9

    # --- known rank-3 / rank-4 numbers from attempt 001, sigma=(1,2,5), f=(1,1,6)/8 ---
    lo_bound = SH.hs_lo((f1, f2, f3))
    assert abs(float(lo_bound) - 37 / 11) < 1e-12

    a_star, d_star, v_star = SH.rank4_optimize(f1, f2)
    gap4 = float(v_star) - float(lo_bound)
    assert abs(gap4 - 5.456e-4) < 2e-6, gap4   # attempt 001's exact rank-4 optimum

    r3_meta_path = os.path.join(PROBLEM_DIR, "data", "certified", "r3_1-8_1-8_6-8_lower.meta.json")
    with open(r3_meta_path, encoding="utf-8") as fh:
        r3_meta = json.load(fh)
    lo3, hi3 = (Fr(x) for x in r3_meta["gap_to_hs"]["enclosure"])
    gap3 = float((lo3 + hi3) / 2)
    assert abs(gap3 - 2.116e-2) < 2e-4, gap3   # attempt 001's certified rank-3 gap

    # --- degenerate-cherry pruning is conservative (never removes a topology
    # whose cherry leaves differ) ---
    assert has_degenerate_cherry((0, (1, 2)), (0, 0, 1)) is False  # leaves 1,2 differ (labels[1]=0,labels[2]=1)
    assert has_degenerate_cherry((0, (1, 2)), (0, 1, 1)) is True   # leaves 1,2 both label 1

    # --- topology counts match the header's arithmetic ---
    f_full = (1.0, 1.0, 1.0)  # all three phases required (dummy, only used for the `required` filter)
    n5 = count_topologies(5, f_full, prune_degenerate=False)
    assert n5 == 42 * 540 * 16, n5

    print("selftest: all checks passed "
          f"(rank-4 gap {gap4:.4e} ~ 5.456e-4, rank-3 gap {gap3:.4e} ~ 2.116e-2, "
          f"rank-5 raw topology count {n5})")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--f", default="1/8,1/8,6/8")
    ap.add_argument("--sigma", default="1,2,5")
    ap.add_argument("--rank", type=int, default=5)
    ap.add_argument("--side", choices=["lower", "upper"], default="lower")
    ap.add_argument("--restarts", type=int, default=1)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--out", help="write best screened trees as JSON (tp_search-compatible)")
    ap.add_argument("--checkpoint", help="checkpoint file (resumable, killable)")
    ap.add_argument("--checkpoint-every", type=int, default=5000)
    ap.add_argument("--time-budget", type=float, help="seconds; stop (and checkpoint) after this")
    ap.add_argument("--resume", action="store_true", help="resume from --checkpoint's best list")
    ap.add_argument("--calibrate", action="store_true", help="time a random sample, estimate full run")
    ap.add_argument("--sample", type=int, default=3000)
    ap.add_argument("--certify", help="certify a screen JSON file (path)")
    ap.add_argument("--name", help="output record name for --certify")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()

    if a.selftest:
        selftest()
        return

    if a.certify:
        record_path, meta_path, verify_out, meta = certify_screen(a.certify, RANK5_DIR, a.name)
        print(f"certified -> {record_path}")
        print(f"verify: {verify_out}")
        print(json.dumps(meta, indent=1))
        return

    f = [float(Fr(x)) for x in a.f.split(",")]
    sig = [float(Fr(x)) for x in a.sigma.split(",")]

    if a.calibrate:
        total = count_topologies(a.rank, f)
        gen = all_topologies(a.rank, f)
        sample = list(itertools.islice(gen, a.sample))
        rng = random.Random(a.seed)
        t0 = time.time()
        for shape, labels, normals in sample:
            top = AxisTopology(shape, labels, normals)

            def fn(x, top=top):
                return objective(top, x, sig, f, 1.0)
            x0 = [rng.uniform(-2, 2) for _ in range(top.n_params)]
            x, val = nelder_mead(fn, x0, iters=200)
            x, val = nelder_mead(fn, x, step=0.05, iters=200)
        dt = time.time() - t0
        per = dt / len(sample)
        print(f"rank={a.rank} total_topologies={total} sample={len(sample)} "
              f"time={dt:.2f}s ({per*1000:.3f} ms/topology, restarts=1)")
        print(f"estimated full run (restarts={a.restarts}): {per*total*a.restarts:.0f}s "
              f"= {per*total*a.restarts/60:.1f} min")
        return

    resume = None
    if a.resume and a.checkpoint:
        payload = load_checkpoint(a.checkpoint)
        if payload:
            resume = resume_best(payload, a.rank)
            print(f"resumed {len(resume)} candidates from {a.checkpoint} "
                  f"(previously done {payload['done']}/{payload['total']})", file=sys.stderr)

    best, done, total, stopped_early = search(
        f, sig, a.rank, a.side, a.restarts, a.seed, a.tol,
        checkpoint=a.checkpoint, checkpoint_every=a.checkpoint_every,
        time_budget=a.time_budget, resume=resume,
    )
    F_ = [Fr(x) for x in a.f.split(",")]
    S_ = [Fr(x) for x in a.sigma.split(",")]
    lo_c, hi_c = min(S_), max(S_)
    hs_lo = 1 / sum(p / (s + lo_c) for p, s in zip(F_, S_)) - lo_c
    hs_hi = 1 / sum(p / (s + hi_c) for p, s in zip(F_, S_)) - hi_c
    print(f"f={a.f} sigma={a.sigma} rank={a.rank} side={a.side} "
          f"topologies={total} done={done} stopped_early={stopped_early}")
    print(f"HS_lo={float(hs_lo):.10f} HS_hi={float(hs_hi):.10f}")
    for mean, aniso, ferr, top, params in best[:5]:
        print(f"  {mean:.9f}  aniso={aniso:.1e} ferr={ferr:.1e} gap={mean-float(hs_lo):.6e} "
              f"labels={top.labels} normals={top.normals} shape={top.shape}")
    if a.out:
        write_screen(a.out, f, sig, a.rank, a.side, total, best)
        print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
