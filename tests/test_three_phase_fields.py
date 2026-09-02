"""Fast deterministic tests for the exact laminate field tool.

These fail if the field propagation, the Hashin-Shtrikman attainment-field
formula, or the Keller-Dykhne duality reduction is broken.
"""

from __future__ import annotations

import sys
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "harness" / "three-phase-conductivity"))
sys.path.insert(0, str(ROOT / "problems" / "three-phase-conductivity" / "explore"))

import laminate as L  # noqa: E402
import tp_coated as TC  # noqa: E402
import tp_fields as TF  # noqa: E402

SIGMA = [Fr(1), Fr(2), Fr(5)]
PHASES = {"p1": SIGMA[0], "p2": SIGMA[1], "p3": SIGMA[2]}
E0 = (Fr(1), Fr(1))

TREE = {
    "fraction": "1/3", "normal": [1, 0],
    "layers": [{"phase": "p1"},
               {"fraction": "3/7", "normal": [1, 2],
                "layers": [{"phase": "p2"}, {"phase": "p3"}]}],
}


def test_field_identities_exact():
    """Average field, average flux and energy all reproduce, in Q."""
    for e in [(Fr(1), Fr(0)), (Fr(0), Fr(1)), (Fr(1), Fr(1)), (Fr(3), Fr(-2))]:
        TF.check_identities(TREE, PHASES, e)


def test_milton_assemblage_realises_attainment_fields():
    """Phases 2 and 3 carry exactly E_i = (s_HS+s1)/(s_i+s1) E0."""
    f = [Fr(3, 8), Fr(2, 8), Fr(3, 8)]
    tree, ph, _ = TC.milton(f, SIGMA, 0)
    assert tree is not None
    hs_lo, targets = TF.hs_target_fields(f, SIGMA, E0)
    got = TF.per_phase(tree, ph, E0)
    for i, name in [(1, "p2"), (2, "p3")]:
        for _, E in got[name]:
            assert E == targets[i]
    # phase 1 is the comparison medium: free pointwise, but its volume
    # average is still pinned to the target.
    tot = sum(v for v, _ in got["p1"])
    avg = tuple(sum(v * E[k] for v, E in got["p1"]) / tot for k in (0, 1))
    assert avg == targets[0]


def test_middle_phase_field_condition_is_miltons_threshold():
    """E_2 <= E0 exactly when HS_lo <= sigma2, i.e. m1 >= 2*Theta*(1-m2)."""
    s1, s2, s3 = SIGMA
    theta = s1 * (s3 - s2) / ((s2 + s1) * (s3 - s1))
    for i in range(1, 8):
        for j in range(1, 8 - i):
            f = [Fr(i, 8), Fr(j, 8), Fr(8 - i - j, 8)]
            hs_lo, targets = TF.hs_target_fields(f, SIGMA, E0)
            ratio = targets[1][0] / E0[0]  # |E_2| / |E0|
            assert (ratio <= 1) == (hs_lo <= s2)
            assert (hs_lo <= s2) == (f[0] >= 2 * theta * (1 - f[1]))


def test_duality_reduces_upper_side_to_lower_side():
    """HS_hi(f, sigma) inverts to HS_lo(f, 1/sigma); isotropic trees follow."""
    f = [Fr(1, 8), Fr(1, 8), Fr(6, 8)]
    lo, hi = L.hs_bounds(f, SIGMA)
    dlo, dhi = L.hs_bounds(f, [1 / s for s in SIGMA])
    assert dhi == 1 / lo and dlo == 1 / hi
    tree, ph, _ = TC.milton([Fr(3, 8), Fr(2, 8), Fr(3, 8)], SIGMA, 0)
    sig = L.effective(tree, ph)
    dual = L.effective(tree, {k: 1 / v for k, v in ph.items()})
    assert sig[1] == 0 and sig[0] == sig[2]
    assert dual == (1 / sig[0], Fr(0), 1 / sig[0])
