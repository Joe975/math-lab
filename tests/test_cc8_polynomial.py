"""Adversarial tests of the independent lifted-polynomial verifier."""
import copy
import importlib.util
from pathlib import Path

import pytest

PATH = Path(__file__).resolve().parents[1] / 'harness/central-configurations-8/verify_polynomial.py'
SPEC = importlib.util.spec_from_file_location('cc8_polynomial', PATH)
checker = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(checker)


@pytest.fixture(scope='module')
def cert():
    pytest.importorskip('numpy')
    # Two equal masses have 2a=(2)^(1/3) at lambda=1.
    a = 2**(-2/3)
    return checker.certify([[a, 0], [-a, 0]])


def test_valid_and_roundtrip(cert):
    import json
    assert checker.verify(json.loads(json.dumps(cert)))


@pytest.mark.parametrize('attack', ['sign', 'anchor', 'center', 'inverse', 'collision',
                                   'mass', 'radius', 'dimension'])
def test_tampering(cert, attack):
    bad = copy.deepcopy(cert)
    if attack == 'sign':
        bad['center'][-1] = '-' + bad['center'][-1]
    elif attack == 'anchor':
        bad['center'][0] = '0'
    elif attack == 'center':
        bad['center'][1] = '-1'
    elif attack == 'inverse':
        bad['inverse'] = [['0']*4 for _ in range(4)]
    elif attack == 'collision':
        bad['center'][1] = bad['center'][0]
    elif attack == 'mass':
        bad['masses'][0] = '2'
    elif attack == 'radius':
        bad['radius'] = '0'
    else:
        bad['center'].pop()
    assert not checker.verify(bad)


def test_reject_malformed():
    assert not checker.verify({})
    assert not checker.verify(None)


def test_reject_proposal_collision():
    pytest.importorskip('numpy')
    with pytest.raises(ValueError, match='Collision'):
        checker.certify([[1, 0], [1, 0]])


def test_eight_body_polygon_and_pure_stdlib_verification(monkeypatch):
    import builtins
    import math
    pytest.importorskip('numpy')
    radius = (sum(1/math.sin(math.pi*k/8) for k in range(1, 8))/4)**(1/3)
    points = [[radius*math.cos(math.pi*k/4), radius*math.sin(math.pi*k/4)]
              for k in range(8)]
    certificate = checker.certify(points)
    original_import = builtins.__import__

    def restricted_import(name, *args, **kwargs):
        if name.split('.')[0] in {'numpy', 'scipy', 'mpmath'}:
            raise AssertionError('Numerical package used during verification')
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', restricted_import)
    assert checker.verify(certificate)


def test_interval_jacobian_contains_difference_quotients():
    # Exact rational secants, separate from any floating numerical derivative.
    F = checker.F
    x = [F(3, 5), F(-3, 5), F(1, 20), F(4, 5)]
    h = F(1, 10000)
    _, jac = checker.system([(t-h, t+h) for t in x], 2)
    for j in range(4):
        left, right = x.copy(), x.copy()
        left[j] -= h
        right[j] += h
        fl, _ = checker.system([(t, t) for t in left], 2)
        fr, _ = checker.system([(t, t) for t in right], 2)
        for i in range(4):
            slope = (fr[i][0]-fl[i][0])/(2*h)
            lo, hi = jac[i].get(j, checker.ZERO)
            assert lo <= slope <= hi


def test_independent_distance_audit_gauge_bounds():
    path = PATH.parents[1].parent / 'problems/central-configurations-8/explore/audit_certificates.py'
    spec = importlib.util.spec_from_file_location('audit', path)
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    # Synthetic coordinate box only, deliberately not a root certificate.
    # Seven pairs touch the fixed y0=0 anchor, so max distance squared is 5;
    # all other pairs can differ by 2 in both coordinates, attaining 8.
    lo, hi = audit.distance_orders({'center': ['0']*15, 'radius': '1'})
    assert lo == [0]*28
    assert hi == [5]*7 + [8]*21
