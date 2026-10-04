import importlib.util
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('cc8_direct', ROOT/'harness/central-configurations-8/central.py')
cc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cc)


def square():
    r = ((1+2*math.sqrt(2))/4)**(1/3)
    return [[r, 0], [0, r], [-r, 0], [0, -r]]


def test_exact_square_root_enclosure():
    for value in ['0', '2', '1/10', '123456789/73']:
        a = cc.Interval(value).sqrt()
        assert a.lo*a.lo <= cc.Q(value) <= a.hi*a.hi
    assert cc.Interval(-1, 2).square().lo == 0


def test_known_square_and_adversaries():
    pytest.importorskip('numpy')
    import copy
    cert = cc.certify(square())
    assert cc.verify(cert)
    for mutation in ('shift', 'inverse', 'gauge', 'collision', 'negative-radius', 'normalization'):
        bad = copy.deepcopy(cert)
        if mutation == 'shift':
            bad['center'][0] = str(cc.Q(bad['center'][0])+cc.Q(1, 100))
        elif mutation == 'inverse':
            bad['inverse'] = [['0']*7 for _ in range(7)]
        elif mutation == 'gauge':
            bad['center'][0] = '0'
        elif mutation == 'collision':
            bad['center'][1:3] = [bad['center'][0], '0']
        elif mutation == 'negative-radius':
            bad['radius'] = '-1'
        else:
            bad['normalization'] = 'I=1'
        assert not cc.verify(bad), mutation


def test_derivative_against_finite_difference():
    np = pytest.importorskip('numpy')
    x = np.array([1.1, .1, .9, -.8, .2, -.1, -.85])
    _, j = cc.interval_system([cc.Interval(str(v)) for v in x])
    for k in range(len(x)):
        eps = np.zeros(len(x)); eps[k] = 1e-6
        fp, _ = cc.interval_system([cc.Interval(str(v)) for v in x+eps])
        fm, _ = cc.interval_system([cc.Interval(str(v)) for v in x-eps])
        for i in range(len(x)):
            fd = (float(fp[i].lo)-float(fm[i].lo))/2e-6
            assert abs(fd-float(j[i][k].lo)) < 1e-6
