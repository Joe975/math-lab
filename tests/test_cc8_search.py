import importlib.util
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'problems/central-configurations-8/explore'))
np = pytest.importorskip('numpy')
pytest.importorskip('scipy')
spec = importlib.util.spec_from_file_location('cc8_search', ROOT/'problems/central-configurations-8/explore/search.py')
s = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)
import symmetric
import continuation


def test_symmetry_canonicalization_and_homometric_adversary():
    q = s.ring(8)
    assert s.equivalent(q, q[::-1] @ np.array([[0., 1.], [1., 0.]])+[3, 2])
    a = np.array([[v, 0.] for v in [0, 1, 4, 10, 12, 17]])
    b = np.array([[v, 0.] for v in [0, 1, 8, 11, 13, 17]])
    assert np.array_equal(np.sort(s.distances(a).ravel()), np.sort(s.distances(b).ravel()))
    assert not s.equivalent(a, b)


def test_solver_and_numeric_derivative():
    q = s.solve(s.gauge(s.ring(8)))
    z = np.delete(q.ravel(), 1)
    f, j = s.system(z)
    assert np.max(np.abs(f)) < 1e-12
    assert np.linalg.norm(q.mean(axis=0)) < 1e-12
    z = z+np.linspace(.001, .02, 15)
    for masses in [None, [2., 1., 1., 1., 1., 1., 1., 1.]]:
        _, j = s.system(z, masses)
        for k in range(15):
            delta = np.zeros(15); delta[k] = 1e-6
            diff = (s.system(z+delta, masses)[0]-s.system(z-delta, masses)[0])/2e-6
            assert np.max(np.abs(diff-j[:, k])) < 1e-6


def test_aligned_scalar_matches_direct_force_and_bracket():
    for t in [1.1, 1.5, 2., 3., 10.]:
        p = np.vstack([s.ring(4), s.ring(4, t)])
        d = p[:, None]-p[None, :]
        r = np.linalg.norm(d, axis=2); np.fill_diagonal(r, np.inf)
        force = np.sum(d/r[:, :, None]**3, axis=1)
        actual = force[0, 0]-force[4, 0]/t
        c = (1+2*np.sqrt(2))/4
        expected = c*(1-t**-3)-(1+1/t)/(t-1)**2+(1-1/t)/(t+1)**2
        assert abs(actual-expected) < 1e-10
    cert = symmetric.bracket()
    assert symmetric.verify_bracket(cert)
    assert not symmetric.verify_bracket({'lo': '2', 'hi': '21/10'})
    assert not symmetric.verify_bracket({'lo': '3', 'hi': '2'})


def test_mass_derivative_and_arclength_step():
    q = s.solve(s.gauge(s.ring(8)))
    w = np.r_[np.delete(q.ravel(), 1), 1.]
    f, j = continuation.augmented(w)
    delta = np.zeros(16); delta[-1] = 1e-6
    fd = (continuation.augmented(w+delta)[0]-continuation.augmented(w-delta)[0])/2e-6
    assert np.max(np.abs(fd-j[:, -1])) < 1e-7
    v = continuation.tangent(w)
    assert np.linalg.norm(j @ v) < 1e-12
    nxt = continuation.step(w, v, .01)
    assert nxt is not None
    assert np.max(np.abs(continuation.augmented(nxt)[0])) < 1e-10
