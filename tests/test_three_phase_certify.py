"""Tests for problems/three-phase-conductivity/explore/tp_certify.py.

Fast and deterministic (no swarm/network calls). Covers:
  * the --selftest path (hand-built two-phase rank-2 tree with exact a==c,
    plus the general cherry/ancestor exact-fraction solver on a synthetic
    3-phase tree);
  * certifying one real screened record from data/screen/ end to end and
    checking it against the independent harness/verify_laminate.py.
"""

import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
EXPLORE_DIR = REPO_ROOT / "problems" / "three-phase-conductivity" / "explore"
SCREEN_DIR = REPO_ROOT / "problems" / "three-phase-conductivity" / "data" / "screen"
HARNESS_DIR = REPO_ROOT / "harness" / "three-phase-conductivity"

sys.path.insert(0, str(EXPLORE_DIR))
sys.path.insert(0, str(HARNESS_DIR))

import tp_certify  # noqa: E402
import laminate  # noqa: E402
import verify_laminate  # noqa: E402


def test_selftest_hand_built_tree_is_exactly_a_equals_c():
    F = Fraction
    s1, s2 = F(1), F(4)
    tree = {
        "fraction": "1/2",
        "normal": [1, 1],
        "layers": [
            {"fraction": "1/3", "normal": [1, 0], "layers": [{"phase": "p1"}, {"phase": "p2"}]},
            {"fraction": "1/3", "normal": [0, 1], "layers": [{"phase": "p1"}, {"phase": "p2"}]},
        ],
    }
    phases = {"p1": s1, "p2": s2}
    t = laminate.effective(tree, phases)
    assert t[0] == t[2]  # exact isotropy of the diagonal, by construction
    assert t == (F(49, 20), F(1, 20), F(49, 20))

    lo, hi = tp_certify.residual_enclosure(t, F(1, 10**12))
    assert hi - lo < F(1, 10**12)
    assert lo <= F(1, 10) <= hi  # exact residual is 2*|b| == 1/10 since a == c


def test_general_cherry_solver_hits_target_fractions_exactly():
    target = {"p1": Fraction(1, 8), "p2": Fraction(1, 8), "p3": Fraction(6, 8)}
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
    solved = tp_certify.solve_exact_fractions(screened_tree, target, 10**6)
    assert solved is not None
    frac_by_id, _note = solved
    tp_certify.apply_fractions(screened_tree, frac_by_id)
    fr = laminate.fractions_of(screened_tree)
    assert fr == target
    phases = {"p1": Fraction(1), "p2": Fraction(2), "p3": Fraction(5)}
    assert laminate.keller_check(screened_tree, phases)


def test_certify_one_screened_record_end_to_end(tmp_path):
    screen_file = SCREEN_DIR / "r3_1-8_1-8_6-8_lower.json"
    assert screen_file.exists(), "expected screen fixture missing"

    result = tp_certify.certify_file(str(screen_file), index=0, denom=10**6,
                                      tol=Fraction(1, 10**12), out_name="pytest_tmp")
    meta = result["meta"]
    assert meta["exact_fraction_match"] is True
    assert all(v == "0" for v in meta["fraction_error"].values())

    record_path = tmp_path / "record.json"
    with open(record_path, "w", encoding="utf-8") as fh:
        json.dump(result["record"], fh)

    # record matches the target phase fractions from the screen filename exactly
    rec = result["record"]
    assert rec["fractions"] == {"p1": "1/8", "p2": "1/8", "p3": "3/4"}

    # independent re-verification (different algebra path) must pass
    problems = verify_laminate.verify_record(str(record_path))
    assert problems == [], problems


def test_certify_cli_selftest_runs_clean():
    proc = subprocess.run(
        [sys.executable, str(EXPLORE_DIR / "tp_certify.py"), "--selftest"],
        capture_output=True, text=True, cwd=str(REPO_ROOT),
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "selftest" in proc.stdout
