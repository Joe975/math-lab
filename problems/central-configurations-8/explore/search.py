"""Reproducible, explicitly non-exhaustive central-configuration searches.

NumPy/SciPy discover candidates. Independent exact checkers certify them.
Run with --help; all output goes to the explicit --output directory.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import root, brentq

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'harness/central-configurations-8'))
import central
import verify_polynomial


def system(z, masses=None):
    q = np.insert(np.asarray(z), 1, 0.).reshape(-1, 2)
    n = len(q)
    m = np.ones(n) if masses is None else np.asarray(masses)
    diff = q[:, None, :]-q[None, :, :]
    r2 = np.sum(diff*diff, axis=2)
    np.fill_diagonal(r2, np.inf)
    inv3 = r2**(-1.5)
    f = q-np.sum(diff*inv3[:, :, None]*m[None, :, None], axis=1)
    h = inv3[:, :, None, None]*(np.eye(2)-3*diff[:, :, :, None]*diff[:, :, None, :]/r2[:, :, None, None])
    j = np.zeros((n, 2, n, 2))
    for i in range(n):
        j[i, :, i, :] = np.eye(2)-np.sum(h[i]*m[:, None, None], axis=0)
        for k in range(n):
            if k != i:
                j[i, :, k, :] = m[k]*h[i, k]
    keep = [i for i in range(2*n) if i != 1]
    return f.ravel()[keep], j.reshape(2*n, 2*n)[np.ix_(keep, keep)]


def gauge(q, masses=None, rescale=True):
    q = np.array(q, dtype=float)
    m = np.ones(len(q)) if masses is None else np.asarray(masses)
    q -= np.average(q, axis=0, weights=m)
    anchor = int(np.argmax(np.sum(q*q, axis=1)))
    q[[0, anchor]] = q[[anchor, 0]]
    # Reordering unequal masses is not allowed; continuation never uses this.
    if masses is not None and not np.all(m == m[0]) and anchor != 0:
        raise ValueError('cannot reorder unequal masses')
    theta = math.atan2(q[0, 1], q[0, 0])
    q = q @ np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
    q[0, 1] = 0
    if rescale:
        r = np.linalg.norm(q[:, None]-q[None, :], axis=2)
        u = sum(m[i]*m[j]/r[i, j] for i in range(len(q)) for j in range(i))
        inertia = np.sum(m[:, None]*q*q)
        q *= (u/inertia)**(1/3)
    return q


def solve(q, masses=None):
    z = np.delete(np.asarray(q).ravel(), 1)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        sol = root(lambda v: system(v, masses)[0], z,
                   jac=lambda v: system(v, masses)[1], tol=1e-10)
    if not np.all(np.isfinite(sol.x)):
        return None
    f, j = system(sol.x, masses)
    if np.max(np.abs(f)) > 2e-9 or abs(sol.x[0]) < 1e-5:
        return None
    # Newton polishing is still a proposal, not certification.
    for _ in range(3):
        try:
            sol.x -= np.linalg.solve(j, f)
            f, j = system(sol.x, masses)
        except np.linalg.LinAlgError:
            return None
    q = np.insert(sol.x, 1, 0.).reshape(-1, 2)
    return q if np.max(np.abs(f)) < 2e-10 else None


def distances(q):
    return np.linalg.norm(q[:, None]-q[None, :], axis=2)


def equivalent(a, b, tol=1e-6):
    """Numerical graph isomorphism with consistent vertex permutation.

    This is only discovery deduplication, never an exact equivalence certificate.
    """
    a, b = distances(a), distances(b)
    if len(a) != len(b):
        return False
    n = len(a)
    candidates = [[j for j in range(n) if np.max(np.abs(np.sort(a[i])-np.sort(b[j]))) < tol]
                  for i in range(n)]
    order = sorted(range(n), key=lambda i: len(candidates[i]))
    mapping = {}

    def extend(k):
        if k == n:
            return True
        i = order[k]
        for j in candidates[i]:
            if j not in mapping.values() and all(abs(a[i, ii]-b[j, jj]) < tol for ii, jj in mapping.items()):
                mapping[i] = j
                if extend(k+1):
                    return True
                del mapping[i]
        return False
    return extend(0)


def ring(count, radius=1., phase=0.):
    t = np.arange(count)*2*np.pi/count+phase
    return radius*np.stack([np.cos(t), np.sin(t)], axis=1)


def seeds(rng, count):
    yield 'regular-octagon', ring(8)
    yield 'heptagon-center', np.vstack([ring(7), [0., 0.]])
    yield 'collinear', np.stack([np.linspace(-2, 2, 8), np.zeros(8)], axis=1)
    for ratio in [.3, .5, .7]:
        for phase in [0., np.pi/4]:
            yield 'two-squares', np.vstack([ring(4), ring(4, ratio, phase)])
    for k in [2, 3, 4, 5, 6]:
        for ratio in [.25, .55, .8]:
            yield 'two-rings', np.vstack([ring(k), ring(8-k, ratio, .23)])
    for i in range(count):
        if i % 3 == 0:
            q = rng.normal(size=(8, 2))
        elif i % 3 == 1:
            q = ring(8)+rng.normal(scale=.25, size=(8, 2))
        else:
            q = rng.normal(size=(8, 2)); q[:, 1] *= .15
        yield 'random', q


def run(output, count=300, seed=20260906, continuation=True):
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    known = []
    rows = []
    attempts = []

    def accept(q, origin):
        q = gauge(q, rescale=False)
        q = solve(q)
        if q is None:
            return None
        for idx, old in enumerate(known):
            if equivalent(q, old):
                return idx
        ident = f'cc-{len(known)+1:03d}'
        try:
            direct = central.certify(q)
            lifted = verify_polynomial.certify(q)
        except ValueError:
            attempts.append({'origin': origin, 'status': 'UNRESOLVED-certification', 'points': q.tolist()})
            return None
        for suffix, cert in [('direct', direct), ('polynomial', lifted)]:
            (output/f'{ident}-{suffix}.json').write_text(json.dumps(cert, indent=2)+'\n')
        known.append(q)
        rows.append({'id': ident, 'origin': origin, 'points': q.tolist(),
                     'direct_certificate': f'{ident}-direct.json',
                     'polynomial_certificate': f'{ident}-polynomial.json'})
        (output/'configurations.json').write_text(json.dumps(rows, indent=2)+'\n')
        print(f'certified {ident} from {origin}', flush=True)
        return len(known)-1

    for k, (origin, start) in enumerate(seeds(rng, count)):
        q = solve(gauge(start))
        idx = accept(q, origin) if q is not None else None
        attempts.append({'round': 'seeds', 'index': k, 'origin': origin,
                         'matched': idx, 'converged': q is not None})
        if k % 50 == 0:
            print(f'seeds {k}: {len(known)} certified classes', flush=True)
    # Continue each discovered branch away from, then back to equal masses.
    # Perturb the anchor mass only. Fix its label and gauge during continuation.
    if continuation:
        initial = list(known)
        for idx, q0 in enumerate(initial):
            for target in [.5, 2.]:
                q = q0.copy()
                path = np.r_[np.linspace(1, target, 11)[1:], np.linspace(target, 1, 11)[1:]]
                min_sv = float('inf')
                steps = 0
                for mass in path:
                    m = np.ones(8); m[0] = mass
                    nxt = solve(q, m)
                    if nxt is None:
                        break
                    q = nxt; steps += 1
                    min_sv = min(min_sv, float(np.linalg.svd(system(np.delete(q.ravel(), 1), m)[1], compute_uv=False)[-1]))
                result = accept(q, 'mass-continuation') if steps == len(path) else None
                attempts.append({'round': 'continuation', 'start': idx, 'mass_target': target,
                                 'steps': steps, 'min_sampled_singular_value': min_sv if steps else None, 'matched': result})
                # Independent branch-breaking proposals at the equal-mass endpoint.
                if steps == len(path):
                    for _ in range(3):
                        perturbed = gauge(q+rng.normal(scale=.08, size=q.shape))
                        qq = solve(perturbed)
                        result = accept(qq, 'continuation-perturbation') if qq is not None else None
                        attempts.append({'round': 'continuation-perturbation', 'start': idx, 'matched': result})
    certs = [json.loads((output/row['direct_certificate']).read_text()) for row in rows]
    separated = all(central.distinct_by_distances(a, b) for a, b in itertools.combinations(certs, 2))
    summary = {'seed': seed, 'random_starts': count, 'certified_classes': len(rows),
               'pairwise_inequivalence_certified': separated, 'complete': False,
               'normalization': 'unit masses, lambda=1, y0=0', 'attempts': attempts,
               'numpy': np.__version__}
    (output/'search-summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'attempts'}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--random-starts', type=int, default=300)
    p.add_argument('--seed', type=int, default=20260906)
    p.add_argument('--no-continuation', action='store_true')
    a = p.parse_args()
    run(a.output, a.random_starts, a.seed, not a.no_continuation)
