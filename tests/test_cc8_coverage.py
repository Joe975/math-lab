import copy
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('cc8_coverage', ROOT/'problems/central-configurations-8/explore/coverage.py')
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


def test_coverage_and_tampering():
    cert = c.build(3)
    assert c.verify(cert)
    bad = copy.deepcopy(cert); del bad['tree']['left']
    assert not c.verify(bad)
    bad = copy.deepcopy(cert); bad['tree']['cut'] = '999'
    assert not c.verify(bad)
    bad = copy.deepcopy(cert); bad['tree'] = {'leaf': 'center-of-mass'}
    bad['leaf_counts'] = {'center-of-mass': 1}
    assert not c.verify(bad)
    bad = copy.deepcopy(cert); bad['complete_classification'] = True
    assert not c.verify(bad)


def test_exclusion_rules():
    box = [(c.Q(1), c.Q(2))]+[(c.Q(1), c.Q(2))]*14
    assert c.reason(box) == 'center-of-mass'
    box = [(c.Q(1), c.Q(2))]+[(c.Q(-14), c.Q(14))]*14
    box[1] = (c.Q(3), c.Q(4))
    assert c.reason(box) == 'farthest-anchor'
