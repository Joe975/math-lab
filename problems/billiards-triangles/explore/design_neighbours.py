#!/usr/bin/env python3
"""Attempt 006: design hunt for stable-orbit word structures past 135 deg.

Queue item 12 (exploratory track). The W(a,b) = (0(12)^a(02)^b)^2 family's
death law gamma_d = 180 - 90(a+b)/(a(b+1)) is now a per-member theorem
(003/004). This tool INVERTS the machinery: it enumerates the structure
class that the glide-reflection reduction applies to -- ALL doubled odd
half-words w = u^2 (any odd word u composed of reflections is
orientation-reversing; u^2 always has trivial rotation part, so every such
w is automatically a translation word) -- and hunts for members whose alive
windows reach angle regions the W family does not cover.

Frontier semantics (pinned; see attempt 006 record): the W family itself is
alive far past 135 (certified members to gamma ~ 159.25 in 001). The actual
135-deg stall is: (i) no word of length <= 26 is known alive at any
gamma >= 135 (001's census, sample-bounded); (ii) the pinch gap
(135.000, 135.0486) between death(W(3,3)) = death(W(4,2)) = 135 and
birth(W(4,3)) ~ 135.0486 has NOTHING known alive at any length; (iii) above
135 the known alive set is the union of the W-family windows
[birth(a,b), gamma_d(a,b)] -- thin isolated windows, not area. A new
structure advances the frontier iff it is alive (a) inside a pinch gap,
(b) at gamma >= 135 with length <= 26, or (c) at some gamma with a shorter
word than the W family provides.

Float screening only -- nothing here is a certificate; exact claims go
through the deathlaw_exact.py pattern. Corridor evaluation is by the
numeric glide reduction (half gates + axis offset m, the 003 reduction
confirmed by 004), cross-checked in --selftest against skeptic_family's
independent full-corridor width (they agree up to the exact factor |tau|).

Subcommands:
  selftest                       reduction vs full corridor; translation
                                 facts; W-family windows reproduced
  screen --max-half L [--klo K --khi K2] [--arcs ...] --out F
                                 enumerate + screen a class:
                                 mode all:  every odd half word |u| <= L
                                 mode s1:   2-alternating class
                                            u = y0 y1 2 y2 2 ... yk 2
  window --word W --out F        adaptive birth/death bisection of one word
  triple --amax A --out F        death table of the three-block family
                                 u = 0(12)^a(02)^b(12)^c
  gapscan --words F --gamma G    dense scan of one arc for listed words

Deterministic (no randomness in enumeration/screen paths). Stdlib only.
"""

import argparse
import cmath
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from skeptic_family import width as sk_width  # noqa: E402

ENDS = ((1, 2), (2, 0), (0, 1))


# ---------------------------------------------------------------- geometry

def apex_from_angles(al_deg, be_deg):
    ta = math.tan(math.radians(al_deg))
    tb = math.tan(math.radians(be_deg))
    return tb / (ta + tb), ta * tb / (ta + tb)


def half_glide_numeric(half, ax, ay):
    """Unfold the HALF word; return (gates, mu, t) with the composed
    orientation-reversing map z -> mu*conj(z) + t, or None if half even."""
    A, B, C = 0j, 1 + 0j, complex(ax, ay)
    tri = [A, B, C]
    gates = []
    for s in half:
        i, j = ENDS[s]
        u, v = tri[i], tri[j]
        gates.append((u, v))
        d = v - u
        e = d * d / (d.real * d.real + d.imag * d.imag)
        tri = [u + e * (z - u).conjugate() for z in tri]
    # composed map from vertex images: z -> mu conj(z) + t, A = 0 -> tri[0]
    t = tri[0]
    mu = (tri[1] - tri[0])  # / conj(B - A) = /1
    return gates, mu, t


def reduced_width(word, ax, ay):
    """Corridor width of the doubled word u^2 via the glide reduction,
    normalized (unit axis direction). Returns (width, |tau|).
    width = min(hi, 2m-lo) - max(lo, 2m-hi) over half-gate projections.
    Requires word = u*2 with |u| odd (checked)."""
    n = len(word)
    half = word[:n // 2]
    if n % 2 or (n // 2) % 2 == 0 or word != half * 2:
        raise ValueError("reduced_width needs a doubled odd half word")
    gates, mu, t = half_glide_numeric(half, ax, ay)
    tau = t + mu * t.conjugate()
    at = abs(tau)
    if at < 1e-14:
        return -math.inf, at
    delta = tau / at
    dc = delta.conjugate()
    m = (dc * t).imag / 2.0
    lo, hi = -math.inf, math.inf
    for u, v in gates:
        pu = (dc * u).imag
        pv = (dc * v).imag
        if pu > pv:
            pu, pv = pv, pu
        if pu > lo:
            lo = pu
        if pv < hi:
            hi = pv
    a = lo if lo > 2 * m - hi else 2 * m - hi
    b = hi if hi < 2 * m - lo else 2 * m - lo
    return b - a, at


# ---------------------------------------------------------------- sampling

def arc_alphas(g, uniform=60, jmax=12, depth=10):
    """Sample alphas on the arc gamma = g: uniform grid over the WHOLE arc
    (both halves) plus geometric accumulation at every candidate window edge
    alpha = 90/j and beta = 90/j (the W family's window edges sit at such
    values; a new family's may too -- the uniform grid is the fallback)."""
    th = 180.0 - g
    als = [th * (k + 0.5) / uniform for k in range(uniform)]
    for j in range(1, jmax + 1):
        e = 90.0 / j
        if e >= th:
            continue
        for d in range(1, depth + 1):
            eps = 0.5 * 2.0 ** -d
            for al in (e - eps, e + eps, th - e + eps, th - e - eps):
                if 0 < al < th:
                    als.append(al)
    return als


def best_alive(word, g, uniform=60, jmax=12, depth=10):
    """(best reduced width, (alpha, beta) argmax) on the arc gamma = g."""
    th = 180.0 - g
    best, arg = -math.inf, None
    for al in arc_alphas(g, uniform, jmax, depth):
        be = th - al
        x, y = apex_from_angles(al, be)
        v, _ = reduced_width(word, x, y)
        if v > best:
            best, arg = v, (al, be)
    return best, arg


ALIVE_TOL = 1e-12


# ------------------------------------------------------------- enumeration

def canonical(word):
    """Canonical form of the doubled word under cyclic rotation, reversal,
    and the 0<->1 letter swap (mirror; gamma-coverage-equivalent)."""
    n = len(word)
    best = None
    for w in (word, word[::-1]):
        for sw in (False, True):
            ww = tuple((1 - c if c < 2 else 2) for c in w) if sw else w
            for r in range(n):
                cand = ww[r:] + ww[:r]
                if best is None or cand < best:
                    best = cand
    return best


def is_power(u):
    """True if u = v^k for some k >= 2 (then u^2 duplicates v^2's corridor
    when |v| odd, or is a repeat of an even translation word)."""
    n = len(u)
    for p in range(1, n):
        if n % p == 0 and u == u[:p] * (n // p):
            return True
    return False


def gen_all_odd(L):
    """All odd half words u, |u| = L, over {0,1,2}, no equal adjacent
    letters in u^2 (i.e. within u and u[-1] != u[0]), u not a proper power;
    yields deduped doubled words via canonical()."""
    seen = set()
    stack = [(c,) for c in (0, 1, 2)]
    while stack:
        u = stack.pop()
        if len(u) == L:
            if u[-1] == u[0] or is_power(u):
                continue
            w = u * 2
            key = canonical(w)
            if key not in seen:
                seen.add(key)
                yield w
            continue
        for c in (0, 1, 2):
            if c != u[-1]:
                stack.append(u + (c,))


def gen_s1(k):
    """The 2-alternating class: u = y0 y1 2 y2 2 ... yk 2, y in {0,1}^{k+1},
    y0 != y1 (one double cell, side 2 struck every other bounce elsewhere).
    W(a,b) is y = 0 1^a 0^b. |u| = 2k+1, |w| = 4k+2."""
    seen = set()
    for bits in range(2 ** (k + 1)):
        y = [(bits >> i) & 1 for i in range(k + 1)]
        if y[0] == y[1]:
            continue
        u = (y[0], y[1], 2)
        for i in range(2, k + 1):
            u += (y[i], 2)
        if is_power(u):
            continue
        w = u * 2
        key = canonical(w)
        if key not in seen:
            seen.add(key)
            yield w


# ------------------------------------------------------------- subcommands

DEFAULT_ARCS = [135.02, 135.5, 136.5, 138.0, 140.0, 142.5, 145.0,
                148.0, 151.0, 155.0, 159.0]


def cmd_screen(args):
    arcs = [float(g) for g in args.arcs.split(",")] if args.arcs \
        else DEFAULT_ARCS
    words = []
    if args.mode == "all":
        for L in range(3, args.max_half + 1, 2):
            words.extend(gen_all_odd(L))
    else:
        for k in range(args.klo, args.khi + 1):
            words.extend(gen_s1(k))
    print(f"# {len(words)} canonical doubled words, arcs {arcs}")
    hits, near = [], []
    for idx, w in enumerate(words):
        best_g, best_v, best_arg = None, -math.inf, None
        for g in arcs:
            v, arg = best_alive(w, g, args.uniform, depth=args.depth)
            if v > best_v:
                best_v, best_g, best_arg = v, g, arg
            if v > ALIVE_TOL:
                hits.append({"word": "".join(map(str, w)), "len": len(w),
                             "gamma": g, "width": v,
                             "alpha": arg[0], "beta": arg[1]})
        if ALIVE_TOL >= best_v > -args.near_tol:
            near.append({"word": "".join(map(str, w)), "len": len(w),
                         "gamma": best_g, "width": best_v,
                         "alpha": best_arg[0], "beta": best_arg[1]})
        if (idx + 1) % 500 == 0:
            print(f"# {idx + 1}/{len(words)} screened, "
                  f"{len(hits)} arc-hits so far", flush=True)
    # cross-check every hit against the independent full corridor
    for h in hits:
        w = tuple(int(c) for c in h["word"])
        x, y = apex_from_angles(h["alpha"], h["beta"])
        h["sk_width"] = sk_width(w, x, y)
        h["sk_agrees"] = bool(h["sk_width"] > 0)
    res = {"mode": args.mode, "max_half": args.max_half,
           "klo": args.klo, "khi": args.khi,
           "n_words": len(words), "arcs": arcs,
           "uniform": args.uniform, "depth": args.depth,
           "hits": hits, "near_misses": near}
    json.dump(res, open(args.out, "w"), indent=1)
    print(f"# wrote {args.out}: {len(hits)} arc-hits "
          f"({len(set(h['word'] for h in hits))} words), "
          f"{len(near)} near misses")


def measure_window(w, g_seed=None, uniform=400, depth=16, iters=44):
    """Adaptive (birth, death) bisection for word w. Coarse walk first."""
    alive = []
    g = 90.25
    while g < 179.6:
        if best_alive(w, g, uniform=120, depth=12)[0] > ALIVE_TOL:
            alive.append(g)
        g += 0.25
    if not alive and g_seed is not None:
        if best_alive(w, g_seed, uniform, depth=depth)[0] > ALIVE_TOL:
            alive = [g_seed]
    if not alive:
        return None
    # death: bisect above the last alive arc
    lo, hi = alive[-1], alive[-1] + 0.25
    arg_d = None
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        v, ar = best_alive(w, mid, uniform, depth=depth)
        if v > ALIVE_TOL:
            lo, arg_d = mid, ar
        else:
            hi = mid
    death = (lo, hi)
    # birth: bisect below the first alive arc
    lo2, hi2 = alive[0] - 0.25, alive[0]
    arg_b = None
    for _ in range(iters):
        mid = 0.5 * (lo2 + hi2)
        v, ar = best_alive(w, mid, uniform, depth=depth)
        if v > ALIVE_TOL:
            hi2, arg_b = mid, ar
        else:
            lo2 = mid
    birth = (lo2, hi2)
    return {"birth_between": birth, "death_between": death,
            "alive_coarse": [alive[0], alive[-1]],
            "death_argmax": arg_d, "birth_argmin": arg_b}


def cmd_window(args):
    w = tuple(int(c) for c in args.word)
    r = measure_window(w, uniform=args.uniform, iters=args.iters)
    res = {"word": args.word, "len": len(w), "window": r}
    if args.out:
        json.dump(res, open(args.out, "w"), indent=1)
    print(json.dumps(res))


def cmd_triple(args):
    rows = []
    for a in range(1, args.amax + 1):
        for b in range(1, args.amax + 1):
            for c in range(1, args.amax + 1):
                u = (0,) + (1, 2) * a + (0, 2) * b + (1, 2) * c
                if is_power(u):
                    continue
                w = u * 2
                r = measure_window(w, uniform=args.uniform, iters=30)
                d = None if r is None else 0.5 * (r["death_between"][0]
                                                 + r["death_between"][1])
                rows.append({"a": a, "b": b, "c": c, "len": len(w),
                             "window": r, "death_mid": d})
                print(f"({a},{b},{c}) len {len(w)}: "
                      f"{'DEAD everywhere sampled' if r is None else r}",
                      flush=True)
    json.dump({"family": "u=0(12)^a(02)^b(12)^c doubled",
               "amax": args.amax, "rows": rows},
              open(args.out, "w"), indent=1)
    print(f"# wrote {args.out}")


def cmd_gapscan(args):
    words = [tuple(int(c) for c in s) for s in
             json.load(open(args.words))] if args.words.endswith(".json") \
        else [tuple(int(c) for c in args.words)]
    out = []
    for w in words:
        v, arg = best_alive(w, args.gamma, uniform=args.uniform,
                            depth=args.depth)
        out.append({"word": "".join(map(str, w)), "gamma": args.gamma,
                    "best_width": v,
                    "alpha": None if arg is None else arg[0]})
        if v > ALIVE_TOL:
            print(f"ALIVE {out[-1]}")
    print(json.dumps({"n": len(words),
                      "n_alive": sum(1 for o in out
                                     if o["best_width"] > ALIVE_TOL),
                      "max_width": max(o["best_width"] for o in out)}))
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)


def cmd_selftest(_args):
    ok = True
    # 1. doubled odd words are translation words: full-word composed map has
    #    linear part 1 (checked numerically at a generic apex)
    test_words = [((0,) + (1, 2) * 2 + (0, 2) * 1) * 2,
                  (0, 1, 2) * 2,
                  (0, 2, 1, 2, 0, 2, 1) * 2]
    for w in test_words:
        A, B = 0j, 1 + 0j
        tri = [A, B, complex(0.41, 0.093)]
        lin, orient = 1 + 0j, 1
        for s in w:
            i, j = ENDS[s]
            u, v = tri[i], tri[j]
            d = v - u
            e = d * d / abs(d) ** 2
            tri = [u + e * (z - u).conjugate() for z in tri]
            lin = e * lin.conjugate()
            orient = -orient
        good = orient == 1 and abs(lin - 1) < 1e-9
        ok &= good
        print(f"[translation {''.join(map(str, w))}]",
              "OK" if good else "FAIL")
    # 2. reduced width * |tau| == full corridor width (skeptic_family), on
    #    the W family and random-ish other doubled words at many apexes
    words = [((0,) + (1, 2) * a + (0, 2) * b) * 2
             for a, b in ((2, 1), (3, 3), (4, 2))]
    words += [(0, 1, 2, 1, 2, 1, 2, 0, 2, 1, 2) * 2,
              (0, 2, 1, 2, 1, 2, 0, 2, 0, 2, 1) * 2]
    import random
    rng = random.Random(20260731)
    worst = 0.0
    n_alive_agree = 0
    for w in words:
        for _ in range(120):
            g = rng.uniform(91.0, 178.0)
            th = 180.0 - g
            al = rng.uniform(0.03 * th, 0.97 * th)
            x, y = apex_from_angles(al, th - al)
            rw, at = reduced_width(w, x, y)
            fw = sk_width(w, x, y)
            if math.isinf(rw) or math.isinf(fw):
                continue
            err = abs(rw * at - fw) / max(1.0, abs(fw))
            worst = max(worst, err)
            if (rw > 1e-12) == (fw > 1e-12):
                n_alive_agree += 1
    good = worst < 1e-8
    ok &= good
    print(f"[reduced vs full corridor: rel err {worst:.2e}, "
          f"verdict agreement {n_alive_agree}/600]", "OK" if good else "FAIL")
    # 3. reproduce two W-family windows (001/003 values)
    w33 = ((0,) + (1, 2) * 3 + (0, 2) * 3) * 2
    r = measure_window(w33, uniform=200, iters=30)
    good = r is not None and abs(r["death_between"][0] - 135.0) < 2e-3 \
        and abs(r["birth_between"][1] - 127.5064) < 2e-2
    ok &= good
    print(f"[W(3,3) window {r['birth_between'][1]:.4f}.."
          f"{r['death_between'][0]:.4f} vs 127.5064..135.0]",
          "OK" if good else "FAIL")
    w42 = ((0,) + (1, 2) * 4 + (0, 2) * 2) * 2
    r = measure_window(w42, uniform=200, iters=30)
    good = r is not None and abs(r["death_between"][0] - 135.0) < 2e-3
    ok &= good
    print(f"[W(4,2) death {r['death_between'][0]:.4f} vs 135.0 "
          f"(mirror-half window)]", "OK" if good else "FAIL")
    print("SELFTEST", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest").set_defaults(fn=cmd_selftest)
    c = sub.add_parser("screen")
    c.add_argument("--mode", choices=["all", "s1"], default="all")
    c.add_argument("--max-half", type=int, default=13)
    c.add_argument("--klo", type=int, default=1)
    c.add_argument("--khi", type=int, default=7)
    c.add_argument("--arcs")
    c.add_argument("--uniform", type=int, default=60)
    c.add_argument("--depth", type=int, default=10)
    c.add_argument("--near-tol", type=float, default=1e-3)
    c.add_argument("--out", required=True)
    c.set_defaults(fn=cmd_screen)
    c = sub.add_parser("window")
    c.add_argument("--word", required=True)
    c.add_argument("--uniform", type=int, default=400)
    c.add_argument("--iters", type=int, default=44)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_window)
    c = sub.add_parser("triple")
    c.add_argument("--amax", type=int, default=3)
    c.add_argument("--uniform", type=int, default=200)
    c.add_argument("--out", required=True)
    c.set_defaults(fn=cmd_triple)
    c = sub.add_parser("gapscan")
    c.add_argument("--words", required=True)
    c.add_argument("--gamma", type=float, required=True)
    c.add_argument("--uniform", type=int, default=2000)
    c.add_argument("--depth", type=int, default=20)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_gapscan)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
