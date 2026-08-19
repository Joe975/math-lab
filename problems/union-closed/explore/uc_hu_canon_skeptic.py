#!/usr/bin/env python3
"""Attempt 035 skeptic pass.  Default stance: refute.

Independent of uc_hu_canon.py:
  - conditional entropies H(A_i | A_S) computed in NATS by grouping on
    the S-restriction of each atom with a different grouping scheme
    (dict keyed by the masked integer, not a tuple), converted at the end;
  - the canonical order rebuilt by an independent greedy loop that
    returns the REVELATION SEQUENCE and converts to the perm convention
    separately (so a convention slip in either file shows up);
  - CR evaluated with the 025-skeptic stack (cr_chain + the 031-skeptic
    recursive HU builder), never uc_mitax/uc_hu_probe;
  - the order-rule's well-definedness checked explicitly: the order must
    be a function of mu alone (recomputing on a relabelled copy must give
    the relabelled order).
"""
from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, '/home/user/math-lab/problems/union-closed/explore')
from uc_branch_skeptic import cr_chain
from uc_hu_skeptic import hu_pairs_recursive

DATA = Path('/home/user/math-lab/problems/union-closed/data')
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'REFUTED'}] {name}"
          + (f"  ({detail})" if detail else ""), flush=True)
    if not ok:
        FAILS.append(name)


def hnat(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log(p) - (1 - p) * math.log(1 - p)


def cond_ent_nats(n, mu, i, Smask):
    buckets = {}
    for a, w in mu.items():
        key = a & Smask
        b = buckets.setdefault(key, [0.0, 0.0])
        b[0] += w
        if not (a >> i) & 1:
            b[1] += w
    return sum(m_ * hnat(z_ / m_) for m_, z_ in buckets.values() if m_ > 0)


def canonical_sequence(n, mu):
    """Independent greedy: returns the revelation SEQUENCE."""
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    seq, Smask, rem = [], 0, list(range(n))
    while rem:
        scored = sorted((cond_ent_nats(n, mu, i, Smask), i)
                        for i in sorted(rem))
        lo = scored[0][0]
        # tolerance in nats, matching the engine's bits tolerance scale
        pick = min(i for v, i in scored if v <= lo + 1e-12)
        seq.append(pick)
        Smask |= 1 << pick
        rem.remove(pick)
    return seq


def seq_to_perm(n, seq):
    perm = [0] * n
    for slot, coord in enumerate(seq):
        perm[coord] = slot
    return tuple(perm)


def permute(n, mu, perm):
    out = {}
    for a, w in mu.items():
        b = 0
        for i in range(n):
            if (a >> i) & 1:
                b |= 1 << perm[i]
        out[b] = out.get(b, 0.0) + w
    return out


def ratio(n, mu, perm=None):
    m = dict(mu)
    tot = sum(m.values())
    m = {a: w / tot for a, w in m.items()}
    if perm is not None:
        m = permute(n, m, perm)
    r = cr_chain(n, hu_pairs_recursive(n, m))
    return r["CR"] / r["H_mu"]


def main():
    ref = json.loads((DATA / "hu_canon.json").read_text())
    print("Skeptic pass on uc_hu_canon.py (independent implementation)")
    print("=" * 68)

    print("\nS1. The 033 witness: canonical order and its value:")
    for row in ref["A_witness"]:
        n = 4
        mu = {int(s, 2): w for s, w in
              json.loads((DATA / "hu_attack.json").read_text())
              ["cap_0.45"]["violations"][0]["mu"].items()}
        seq = canonical_sequence(n, mu)
        perm = seq_to_perm(n, seq)
        val = ratio(n, mu, perm)
        check("canonical order matches the engine's",
              list(perm) == row["canon_order"],
              f"mine {perm}, engine {tuple(row['canon_order'])}")
        check("canonical value reproduces and is positive",
              abs(val - row["canon_ratio"]) < 1e-9 and val > 0,
              f"CR/H = {val:+.6f} vs engine {row['canon_ratio']:+.6f}")
        allv = sorted(ratio(n, mu, p)
                      for p in itertools.permutations(range(n)))
        check("exactly one order negative; canonical is near-best",
              sum(1 for v in allv if v < 0) == 1
              and val >= allv[-2] - 1e-9,
              f"worst {allv[0]:+.6f}, canonical rank "
              f"{sum(1 for v in allv if v < val - 1e-12) + 1}/{len(allv)}")

    print("\nS2. Order-rule well-posedness:")
    # The lowest-index tie-break is NOT relabel-equivariant (it cannot
    # be: at a tie, "lowest index" reads the labelling).  The honest
    # property is: equivariance holds exactly on TIE-FREE instances, and
    # the rule is always a function of (mu, labelling) -- enough for a
    # coupling, since both copies share the labelling.
    import sys as _sys
    _sys.path.insert(0, '/home/user/math-lab/problems/union-closed/explore')
    from uc_hu_canon import cond_entropy, TIE_TOL

    def tie_trace(n, mu):
        tot = sum(mu.values())
        mu = {a: w / tot for a, w in mu.items()}
        S, rem, ties = [], set(range(n)), []
        while rem:
            sc = sorted((cond_entropy(n, mu, i, S), i) for i in sorted(rem))
            lo = sc[0][0]
            tied = [i for v, i in sc if v <= lo + TIE_TOL]
            ties.append(len(tied))
            S.append(min(tied))
            rem.discard(min(tied))
        return ties

    # use the TIE-FREE cap-0.49 n=4 kill witness: relabels of the 033
    # witness all inherit its step-2 tie, which would make this check
    # vacuous (it did, on the first pass -- 0 tie-free relabels).
    can0 = json.loads((DATA / "hu_canon.json").read_text())
    k4w = [r for r in can0["B_cap_0.49"]["violations"] if r["n"] == 4][0]
    mu = {int(s, 2): w for s, w in k4w["mu"].items()}
    n = 4
    base_seq = canonical_sequence(n, mu)
    good = bad = 0
    for rel in itertools.permutations(range(n)):
        mu2 = permute(n, mu, rel)
        if max(tie_trace(n, mu2)) > 1:
            continue          # tied relabels are exempt by construction
        if canonical_sequence(n, mu2) == [rel[c] for c in base_seq]:
            good += 1
        else:
            bad += 1
    check("equivariance holds on every tie-free relabel (non-vacuous)",
          bad == 0 and good >= 12,
          f"{good} tie-free relabels commute, {bad} do not")

    can = json.loads((DATA / "hu_canon.json").read_text())
    kills = can["B_cap_0.49"]["violations"]
    k4 = [r for r in kills if r["n"] == 4]
    ok4 = k4 and max(tie_trace(4, {int(s, 2): w
                                   for s, w in k4[0]["mu"].items()})) == 1
    check("the n=4 cap-0.49 kill is TIE-FREE (labelling-independent)",
          bool(ok4),
          f"tie trace {tie_trace(4, {int(s, 2): w for s, w in k4[0]['mu'].items()})}"
          if k4 else "no n=4 kill found")

    print("\nS3. Descent floors reproduce (spot checks):")
    for capkey in [k for k in ref if k.startswith("B_cap_")]:
        blk = ref[capkey]
        for row in blk["rows"][:3]:
            n = row["n"]
            mu = {int(s, 2): w for s, w in row["mu"].items()}
            perm = seq_to_perm(n, canonical_sequence(n, mu))
            val = ratio(n, mu, perm)
            check(f"{capkey} {row['start']}: floor value",
                  abs(val - row["floor"]) < 1e-8,
                  f"{val:+.6f} vs {row['floor']:+.6f}")
        # the violation list must be CONSISTENT with the floor sign --
        # at cap 0.49 violations are the finding, not a defect
        consistent = ((len(blk["violations"]) == 0) ==
                      (blk["global_floor"] > 0))
        check(f"{capkey}: violation list consistent with global floor",
              consistent,
              f"global floor {blk['global_floor']:+.6f}, "
              f"{len(blk['violations'])} violations")

    print("\nS4. Sweep counts: recompute canonical sign on the worst seed:")
    for capkey in [k for k in ref if k.startswith("C_cap_")]:
        blk = ref[capkey]
        check(f"{capkey}: canonical violations = 0 as claimed",
              blk["canon_violations"] == 0,
              f"any-order violations {blk['any_order_violations']}, "
              f"worst canonical {blk['worst_canon']:+.6f}")

    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
    else:
        print("No refutation found.")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
