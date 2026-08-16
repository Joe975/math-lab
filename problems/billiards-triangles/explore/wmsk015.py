#!/usr/bin/env python3
"""Skeptic engine for attempt 015 (window-law margin system) — attempt 016.

Independent verification layers (default stance: refute), stdlib only,
deterministic (seeded RNG):

  selftest   (a) margin closed forms re-implemented from scratch and checked
                 against the *geometry*: sign of min signed margin vs the
                 corridor width of skeptic_family.py (002's blind-era
                 implementation, independent of the whole 003-015 stack),
                 near-window AND far-from-window triangles;
             (b) identity checks on the closed-form system:
                 AC0[death] = -N1  (015 SS1 says -N1-N2: WRONG, see record),
                 BC0[death] = N4 (I4, all-positive segment form),
                 AB_a = -PB (sine addition),
                 tau: Re(conj(delta) tau) = -4 sin(th) sin(a al) sin((b+1)be)
                 (my closed form; checked against deathlaw_symbolic's exact
                 tau, 40 members, random points);
             (c) ring bijection of every closed form vs the symbolic engine
                 on members 015 never checked (all a <= 10 plus large spot
                 members incl. a > 2b+3 shapes).
  lemmas     adversarial scans of the seven S-lemmas AND of every named
             intermediate step in the 016 proofs (the steps are what the
             proofs rest on; a false step must show up as a negative
             margin here). Grid: all (a,b), a <= 40, dense t; plus random
             large-a spot checks to a = 3000.
  t2         adversarial search for alive() points with a.alpha >= 90 or
             (b+1).beta >= 90 (T2 says none exist), targeted at branch
             windows the record's uniform scan under-samples.

Angles in degrees throughout.  T = a(b+1); death range t in (0, 22.5/T]
(k = 1/2), birth range the same.
"""

import argparse
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

D2R = math.pi / 180.0


def sind(x):
    return math.sin(x * D2R)


def cosd(x):
    return math.cos(x * D2R)


# ---------------------------------------------------------------- margins
# My own implementation of the 015 closed forms (written from the 015
# derivation, section 1; independent code, same mathematics).

def margins(a, b, al, be):
    th = al + be
    N2 = cosd(a * al) * sind((b + 1) * be) * sind(th)
    PB = cosd((b + 1) * be) * sind(a * al) * sind(th)
    out = [("N2", N2, +1), ("PB", PB, +1)]
    for s in range(a):
        out.append((f"AC{s}",
                    N2 + sind(be) * sind((a - 2 * s - 1) * al - (b + 1) * be),
                    -1))
    for s in range(1, a + 1):
        out.append((f"AB{s}",
                    N2 + sind(th) * sind((a - 2 * s) * al - (b + 1) * be),
                    -1))
    for j in range(b):
        out.append((f"BC{j}",
                    -PB + sind(al) * sind(a * al + (b - 2 * j) * be), +1))
    for j in range(1, b + 1):
        out.append((f"BA{j}",
                    -PB + sind(th) * sind(a * al + (b + 1 - 2 * j) * be), +1))
    return out


def min_signed(a, b, al, be):
    return min(v * s for _, v, s in margins(a, b, al, be))


def alive(a, b, al, be):
    return min_signed(a, b, al, be) > 0


# ------------------------------------------------------------- segments

def death_point(a, b, t):
    A = 90.0 / a
    c = 90.0 * (a - 1) / (a * (b + 1))
    return A - t, c + 2 * t


def birth_point(a, b, t):
    A = 90.0 / a
    H = 90.0 / (b + 1)
    return A - t, H - t


# --------------------------------------------------------------- selftest

def check_geometry(n_samples=1200, amax=12, seed=20260817):
    """Sign of min signed margin vs corridor width of skeptic_family (002)."""
    from skeptic_family import width as sk_width, family_word
    rng = random.Random(seed)
    agree = tried = 0
    worst = None
    while tried < n_samples:
        a = rng.randint(2, amax)
        b = rng.randint(1, a)
        # mix: near the death corner, near the birth corner, random obtuse
        mode = rng.random()
        A = 90.0 / a
        if mode < 0.4:
            ac, bc = A, 90.0 * (a - 1) / (a * (b + 1))
        elif mode < 0.8:
            ac, bc = A, 90.0 / (b + 1)
        else:
            ac, bc = rng.uniform(1, 60), rng.uniform(1, 60)
        al = ac + rng.uniform(-0.8, 0.8)
        be = bc + rng.uniform(-0.8, 0.8)
        if not (0 < al and 0 < be and al + be < 90):
            continue
        ta, tb = math.tan(al * D2R), math.tan(be * D2R)
        x, y = tb / (ta + tb), ta * tb / (ta + tb)
        w = sk_width(family_word(a, b), x, y)
        mine = min_signed(a, b, al, be)
        if abs(w) < 1e-9 or abs(mine) < 1e-9:
            continue
        tried += 1
        ok = (w > 0) == (mine > 0)
        agree += ok
        if not ok and worst is None:
            worst = (a, b, al, be, w, mine)
    return agree, tried, worst


def check_identities(seed=20260817):
    """Float identity checks on the closed-form system (my formulas)."""
    rng = random.Random(seed)
    max_ac0 = max_bc0 = max_aba = 0.0
    for _ in range(4000):
        a = rng.randint(2, 40)
        b = rng.randint(1, a)
        T = a * (b + 1)
        t = rng.uniform(1e-6, 22.5 / T)
        al, be = death_point(a, b, t)
        th = al + be
        m = dict((lbl, v) for lbl, v, _ in margins(a, b, al, be))
        # N1 (brief/plaw_suffice form) on the death segment
        N1 = (cosd(a * t) * sind(be) * sind((2 * b + 1) * t)
              - sind(a * t) * sind(al) * sind(b * be))
        max_ac0 = max(max_ac0, abs(m["AC0"] + N1))
        # N4 (I4 all-positive segment form) = BC0
        N4 = (sind(a * t) * sind(al) * sind(b * be)
              + cosd(a * t) * sind(be) * sind((2 * b + 1) * t))
        max_bc0 = max(max_bc0, abs(m["BC0"] - N4))
        # AB_a = -PB
        max_aba = max(max_aba, abs(m[f"AB{a}"] + m["PB"]))
        # the 015 SS1 claim AC0 = -(N1+N2) would need |AC0 + N1 + N2| ~ 0:
        # tracked implicitly: if AC0 = -N1 holds to 1e-12 and N2 is O(t),
        # the -(N1+N2) version is falsified whenever N2 is not tiny.
    return max_ac0, max_bc0, max_aba


def check_tau(seed=20260818):
    """My closed form R = -4 sin(th) sin(a al) sin((b+1) be) against the
    exact symbolic tau of deathlaw_symbolic (2 random points x 40 members)."""
    import deathlaw_symbolic as ds
    rng = random.Random(seed)
    worst = 0.0
    members = [(a, b) for a in range(2, 9) for b in range(1, a + 1)]
    members += [(10, 3), (11, 1), (12, 12), (13, 5), (15, 7), (16, 2)]
    for (a, b) in members:
        gd = ds.glide_data(a, b)
        pr = ds.pre(ds.pmul(ds.pconj(gd["delta"]), gd["tau"]))
        for _ in range(2):
            al = rng.uniform(0.5, 80.0 / a)
            be = rng.uniform(0.5, 80.0 / (b + 1))
            if al + be >= 90:
                continue
            got = ds.peval(pr, al * D2R, be * D2R).real
            mine = -4 * sind(al + be) * sind(a * al) * sind((b + 1) * be)
            worst = max(worst, abs(got - mine))
    return worst, len(members)


def check_ring(members=None):
    """Ring bijection of the closed forms on members 015 did not check."""
    from wmargins import _ring_check
    if members is None:
        members = [(a, b) for a in range(2, 11) for b in range(1, a + 1)]
        members += [(11, 2), (12, 7), (13, 13), (14, 1), (16, 3),
                    (17, 16), (19, 1), (20, 7)]
    bad = [m for m in members if not _ring_check(*m)]
    return bad, len(members)


def cmd_selftest(args):
    ok = True
    agree, tried, worst = check_geometry()
    print(f"[geometry] closed-form alive vs skeptic_family width: "
          f"{agree}/{tried} agree")
    if worst:
        print("  FIRST DISAGREEMENT:", worst)
    ok &= agree == tried

    e1, e2, e3 = check_identities()
    print(f"[identity] AC0 = -N1 (death):   max |err| {e1:.3e}")
    print(f"[identity] BC0 = N4 (death):    max |err| {e2:.3e}")
    print(f"[identity] AB_a = -PB:          max |err| {e3:.3e}")
    ok &= max(e1, e2, e3) < 1e-12

    et, nm = check_tau()
    print(f"[tau] R = -4 sin(th) sin(a al) sin((b+1) be) vs symbolic tau: "
          f"max |err| {et:.3e} over {nm} members")
    ok &= et < 1e-9

    bad, n = check_ring()
    print(f"[ring] closed-form bijection on {n} members "
          f"(fresh set): {'ALL EXACT' if not bad else f'FAIL {bad}'}")
    ok &= not bad

    print("SELFTEST", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


# ----------------------------------------------------------- lemma scans
# Each entry: (name, callable(a, b, t) -> min margin over the family index).
# All must be strictly positive on 0 < t <= 22.5/T.  These are BOTH the
# S-lemma statements themselves and the intermediate steps my proofs use.

def scan_items_death(a, b, t):
    A = 90.0 / a
    T = a * (b + 1)
    c = 90.0 * (a - 1) / T
    al, be = death_point(a, b, t)
    th = al + be
    thd = A + c
    sbb = sind((b + 1) * be)
    out = []

    # --- [D-AB]  sin(L_s) > sin(at) sin((b+1)be), L_s = (2s-1)A+(a+2b+2-2s)t
    m_full = min(sind((2 * s - 1) * A + (a + 2 * b + 2 - 2 * s) * t)
                 - sind(a * t) * sbb for s in range(1, a))
    out.append(("D-AB", m_full))
    # steps: L_s in (at, 180-at)
    m_lo = min((2 * s - 1) * A + (a + 2 * b + 2 - 2 * s) * t - a * t
               for s in range(1, a))
    m_hi = min(180 - a * t - ((2 * s - 1) * A + (a + 2 * b + 2 - 2 * s) * t)
               for s in range(1, a))
    out.append(("D-AB/step-window", min(m_lo, m_hi)))

    # --- [D-BA]  cos(at - r be) > sin(A-2(b+1)t) cos(at), r = b+1-2j
    m_full = min(cosd(a * t - (b + 1 - 2 * j) * be)
                 - sind(A - 2 * (b + 1) * t) * cosd(a * t)
                 for j in range(1, b + 1))
    out.append(("D-BA", m_full))
    out.append(("D-BA/step-angle", (90 - A) - (a * t + (b - 1) * be)))

    # --- [D-BC]  sin(al) cos(at - r be) > PB, r = b-2j, j=1..b-1
    PB = sind(A - 2 * (b + 1) * t) * cosd(a * t) * sind(th)
    if b >= 2:
        m_full = min(sind(al) * cosd(a * t - (b - 2 * j) * be) - PB
                     for j in range(1, b))
        out.append(("D-BC", m_full))
        out.append(("D-BC/step-angle",
                    b * c - t - (a * t + (b - 2) * be)))

    # --- [D-AC]  sin(be) sin(M_s) > N2, M_s = 2sA + (a+2b+1-2s)t
    N2 = sind(a * t) * sbb * sind(th)
    m_full = min(sind(be) * sind(2 * s * A + (a + 2 * b + 1 - 2 * s) * t)
                 - N2 for s in range(1, a))
    out.append(("D-AC", m_full))
    m0 = 2 * A + (a - 2 * b - 3) * t
    # step: all M_s within [m0, 180-m0] (CLOSED interval: s = a-1 attains
    # the upper end exactly), and m0 in (0, 90]
    m_win = min(min(2 * s * A + (a + 2 * b + 1 - 2 * s) * t - m0,
                    180 - m0 - (2 * s * A + (a + 2 * b + 1 - 2 * s) * t))
                for s in range(1, a))
    out.append(("D-AC/step-window", min(m_win + 1e-9, m0, 90 - m0 + 1e-12)))
    # step: core  sin(be) sin(m0) > sin(at) sin(th)
    out.append(("D-AC/core", sind(be) * sind(m0) - sind(a * t) * sind(th)))
    # step: pairing pieces of the core:
    #   pair2 = sin(3A/2 + c + (a-b)t) * sin(A/2 - (b+1)t)  must be > 0
    #   phi1  = c - A/2 + (b+3-a)t  claimed >= 0 for all (a,b) != (2,2)
    pair2 = (sind(1.5 * A + c + (a - b) * t) * sind(A / 2 - (b + 1) * t))
    out.append(("D-AC/step-pair2", pair2))
    if (a, b) != (2, 2):
        out.append(("D-AC/step-phi1", c - A / 2 + (b + 3 - a) * t))
    # Delta itself (= 2 * core, via the pairing) as a cross-check
    delta = (2 * sind(1.5 * A - (b + 2) * t)
             * sind(c - A / 2 + (b + 3 - a) * t) + 2 * pair2)
    out.append(("D-AC/pairing-identity",
                1e-9 - abs(delta - 2 * (sind(be) * sind(m0)
                                        - sind(a * t) * sind(th)))))

    # --- [D-N1]  N1 > 0 and my proof chain
    N1 = (cosd(a * t) * sind(be) * sind((2 * b + 1) * t)
          - sind(a * t) * sind(al) * sind(b * be))
    out.append(("D-N1", N1))
    # step: sin(lx) >= l sin(x) cos(lx) with x = at, l = (2b+1)/a
    lam = (2 * b + 1) / a
    out.append(("D-N1/step-lambda",
                sind((2 * b + 1) * t)
                - lam * sind(a * t) * cosd((2 * b + 1) * t)))
    # step: sin(b be) <= b sin(be)  (equality iff b = 1: skip that case)
    if b >= 2:
        out.append(("D-N1/step-nsin", b * sind(be) - sind(b * be)))
    # step: (2b+1) cos(at) cos((2b+1)t) > a b sin(A)
    out.append(("D-N1/step-final",
                (2 * b + 1) * cosd(a * t) * cosd((2 * b + 1) * t)
                - a * b * sind(A)))

    # --- N4 death row: positivity of the all-positive form + range checks
    N4 = (sind(a * t) * sind(al) * sind(b * be)
          + cosd(a * t) * sind(be) * sind((2 * b + 1) * t))
    out.append(("D-N4", N4))
    out.append(("D-N4/step-range", min(180 - b * be, 180 - (2 * b + 1) * t)))
    return out


def scan_items_birth(a, b, t):
    A = 90.0 / a
    H = 90.0 / (b + 1)
    al, be = birth_point(a, b, t)
    th = al + be
    N2 = sind(a * t) * cosd((b + 1) * t) * sind(th)
    PB = sind((b + 1) * t) * cosd(a * t) * sind(th)
    out = []

    # --- [B-AB] sin(G_s) > sin(at) cos((b+1)t)
    m_full = min(sind(2 * s * A + (a - b - 1 - 2 * s) * t)
                 - sind(a * t) * cosd((b + 1) * t) for s in range(1, a))
    out.append(("B-AB", m_full))
    m_lo = min(2 * s * A + (a - b - 1 - 2 * s) * t - a * t
               for s in range(1, a))
    m_hi = min(180 - a * t - (2 * s * A + (a - b - 1 - 2 * s) * t)
               for s in range(1, a))
    out.append(("B-AB/step-window", min(m_lo, m_hi)))

    # --- [B-AC] sin(be) sin(J_s) > N2, J_s = (2s+1)A + (a-b-2-2s)t
    m_full = min(sind(be) * sind((2 * s + 1) * A + (a - b - 2 - 2 * s) * t)
                 - N2 for s in range(a))
    out.append(("B-AC", m_full))
    # step: sin(J_s) >= sin(J_0) for s >= 1 (s = 0 is equality by definition)
    sJ0 = sind(A + (a - b - 2) * t)
    if a >= 2:
        out.append(("B-AC/step-J0min",
                    min(sind((2 * s + 1) * A + (a - b - 2 - 2 * s) * t) - sJ0
                        for s in range(1, a))))
    # step: E-identity  sin(be) sin(J_0) - N2 ==
    #        sin(al) sin(H-(a-b)t) - sin(th) cos(at) sin((b+1)t)
    E1 = sind(be) * sJ0 - N2
    core = (sind(al) * sind(H - (a - b) * t)
            - sind(th) * cosd(a * t) * sind((b + 1) * t))
    out.append(("B-AC/E-identity", 1e-9 - abs(E1 - core)))
    # step: the core inequality itself
    out.append(("B-AC/core", core))
    # step: chord-domain  H-(a-b)t <= th  and  H-(a-b)t >= (3/4)H
    out.append(("B-AC/step-chord",
                min(th - (H - (a - b) * t),
                    (H - (a - b) * t) - 0.75 * H)))
    # step: final constants (b < a uses 3H/4; b = a uses H)
    if b < a:
        out.append(("B-AC/step-const", 2.3625 * H - th))
    else:
        out.append(("B-AC/step-const-ba", 3.15 * H - th))

    # --- [B-BC] sin(al) cos(rH - (a+r)t) > PB, r = b-2j, j=0..b-1
    m_full = min(sind(al) * cosd((b - 2 * j) * H - (a + b - 2 * j) * t) - PB
                 for j in range(b))
    out.append(("B-BC", m_full))
    # step: |X_r| <= 90 - H
    out.append(("B-BC/step-Xr",
                min(90 - H - abs((b - 2 * j) * H - (a + b - 2 * j) * t)
                    for j in range(b))))

    # --- [B-BA] cos(at - r'(H-t)) > sin((b+1)t) cos(at)
    m_full = min(cosd(a * t - (b + 1 - 2 * j) * (H - t))
                 - sind((b + 1) * t) * cosd(a * t) for j in range(1, b + 1))
    out.append(("B-BA", m_full))
    out.append(("B-BA/step-angle",
                (2 * H - a * t + (b - 1) * t) - (b + 1) * t))

    # --- AB_a (P row): N2 > 0 via addition formula, trivial positivity
    out.append(("B-ABa", sind(th) * cosd(a * t) * sind((b + 1) * t)))
    return out


def cmd_lemmas(args):
    worst = {}

    def feed(items, a, b, t):
        for name, v in items:
            if name not in worst or v < worst[name][0]:
                worst[name] = (v, a, b, t)

    # dense grid
    for a in range(2, args.amax + 1):
        for b in range(1, a + 1):
            T = a * (b + 1)
            top = 22.5 / T
            for i in range(1, args.steps + 1):
                t = top * i / args.steps
                feed(scan_items_death(a, b, t), a, b, t)
                feed(scan_items_birth(a, b, t), a, b, t)
    # random large-a spot checks
    rng = random.Random(20260819)
    for _ in range(args.rand):
        a = int(math.exp(rng.uniform(math.log(2), math.log(3000))))
        a = max(2, a)
        b = rng.randint(1, a)
        T = a * (b + 1)
        t = rng.uniform(1e-9, 22.5 / T)
        feed(scan_items_death(a, b, t), a, b, t)
        feed(scan_items_birth(a, b, t), a, b, t)

    bad = 0
    for name in sorted(worst):
        v, a, b, t = worst[name]
        flag = "ok" if v > 0 else "VIOLATED"
        if v <= 0:
            bad += 1
        print(f"{name:24s} min {v:+.3e}  at W({a},{b}) t={t:.6g}  {flag}")
    print(f"{len(worst)} tracked inequalities, {bad} violated")
    sys.exit(1 if bad else 0)


# ------------------------------------------------------------------- t2

def cmd_t2(args):
    """Adversarial search for alive points beyond the birth walls."""
    rng = random.Random(args.seed)
    hits = []
    tried = 0
    for _ in range(args.samples):
        a = rng.randint(2, 14)
        b = rng.randint(1, a)
        branch = rng.random()
        if branch < 0.5:
            # a al slightly above 90 (the cos(a al) < 0 branch)
            al = (90.0 + rng.uniform(0, 30)) / a
            be = rng.uniform(0.1, 90.0 / (b + 1) * rng.choice([0.9, 1.4]))
        else:
            # (b+1) be above 90, including higher branches (360m + x)
            m = rng.choice([0, 0, 0, 1])
            be = (90.0 * rng.uniform(1.0, 1.8) + 360.0 * m) / (b + 1)
            al = rng.uniform(0.1, 90.0 / a)
        if not (al > 0 and be > 0 and al + be < 90):
            continue
        if a * al < 90 and (b + 1) * be < 90:
            continue  # not a violation candidate
        tried += 1
        if alive(a, b, al, be):
            hits.append((a, b, al, be))
    print(f"[t2] {tried} candidate points beyond the walls, "
          f"{len(hits)} alive (T2 predicts 0)")
    for h in hits[:10]:
        print("  ALIVE BEYOND WALL:", h)
    sys.exit(1 if hits else 0)


# ------------------------------------------------------------------ n4leg
# The missing certificate leg of wmargins_certify.py: the death-segment
# BC0 row (= N4).  Closed exactly, per member, in two exact steps:
#   (1) ring identity  BC0_form == N4_form  (pure trig identity in
#       Q[e^{i al}, e^{i be}]; this is 005's I4 restated in the margin
#       language; also proven for ALL (a,b) by plaw_general/Lean),
#   (2) rational range checks making every factor of N4's death-segment
#       form  sin(at) sin(al) sin(b be) + cos(at) sin(be) sin((2b+1)t)
#       strictly positive on t in (0, 22.5/T].
# Together: BC0 > 0 on the whole death segment, exact, no floats.

def n4_member_exact(a, b):
    from fractions import Fraction as F
    import deathlaw_symbolic as ds
    cos_of = lambda mm, nn: ds.pscale(
        ds.padd(ds.mono(mm, nn), ds.mono(-mm, -nn)), (F(1, 2), F(0)))
    sin_of = ds.sin_of
    # BC0 closed form: sin(al) sin(a al + b be) - PB,
    # PB = cos((b+1) be) sin(a al) sin(al+be)     [angle unit: arg/2]
    PB = ds.pmul(ds.pmul(cos_of(0, 2 * (b + 1)), sin_of(2 * a, 0)),
                 sin_of(2, 2))
    BC0 = ds.psub(ds.pmul(sin_of(2, 0), sin_of(2 * a, 2 * b)), PB)
    # N4 (I4): cos(a al) sin(al) sin(b be) - sin(a al) sin(be) cos(al+(b+1)be)
    N4 = ds.psub(
        ds.pmul(ds.pmul(cos_of(2 * a, 0), sin_of(2, 0)), sin_of(0, 2 * b)),
        ds.pmul(ds.pmul(sin_of(2 * a, 0), sin_of(0, 2)),
                cos_of(2, 2 * (b + 1))))
    ring_ok = ds.pclean(ds.psub(BC0, N4)) == {}
    # rational ranges on the death segment, t in (0, T], T = 45/(2 a (b+1)):
    # at in (0,90); al = A - t in (0,90); b be = b c + 2bt in (0,180);
    # (2b+1) t in (0,180); be in (0,180).  All endpoints rational.
    T = F(45, 2 * a * (b + 1))
    A = F(90, a)
    c = F(90 * (a - 1), a * (b + 1))
    rng_ok = (a * T < 90 and 0 < A - T and A < 90
              and b * c + 2 * b * T < 180 and (2 * b + 1) * T < 180
              and c + 2 * T < 180)
    return ring_ok, rng_ok


def cmd_n4leg(args):
    ok = True
    n = 0
    for a in range(2, args.amax + 1):
        for b in range(1, a + 1):
            ring_ok, rng_ok = n4_member_exact(a, b)
            n += 1
            if not (ring_ok and rng_ok):
                ok = False
                print(f"W({a},{b}): ring={ring_ok} ranges={rng_ok}  FAIL")
    print(f"[n4leg] death BC0 = N4 identity (exact ring) + positive-factor "
          f"ranges (exact rational): {n} members a <= {args.amax}: "
          f"{'ALL CERTIFIED' if ok else 'FAILURES'}")
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest").set_defaults(fn=cmd_selftest)
    c = sub.add_parser("n4leg")
    c.add_argument("--amax", type=int, default=14)
    c.set_defaults(fn=cmd_n4leg)
    c = sub.add_parser("lemmas")
    c.add_argument("--amax", type=int, default=40)
    c.add_argument("--steps", type=int, default=96)
    c.add_argument("--rand", type=int, default=20000)
    c.set_defaults(fn=cmd_lemmas)
    c = sub.add_parser("t2")
    c.add_argument("--samples", type=int, default=500000)
    c.add_argument("--seed", type=int, default=20260820)
    c.set_defaults(fn=cmd_t2)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
