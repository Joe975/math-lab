"""Pseudo-arclength push through folds missed by fixed-mass continuation.

Evidence only: a sampled curve is not an interval-certified continuation.
Only equal-mass endpoints that pass both root checkers count as objects.
"""
import argparse
import json
from pathlib import Path

import search as s


def augmented(w):
    z, mass = w[:-1], w[-1]
    m = s.np.ones(8); m[0] = mass
    f, j = s.system(z, m)
    p = s.np.insert(z, 1, 0.).reshape(8, 2)
    col = s.np.zeros((8, 2))
    for i in range(1, 8):
        d = p[i]-p[0]
        col[i] = -d/s.np.linalg.norm(d)**3
    return f, s.np.column_stack([j, s.np.delete(col.ravel(), 1)])


def tangent(w, previous=None, direction=1):
    _, j = augmented(w)
    v = s.np.linalg.svd(j, full_matrices=True)[2][-1]
    if previous is not None:
        if v @ previous < 0:
            v = -v
    elif v[-1]*direction < 0:
        v = -v
    return v


def step(w, v, ds):
    prediction = w+ds*v

    def fun(x):
        f, j = augmented(x)
        return s.np.r_[f, (x-prediction) @ v], s.np.vstack([j, v])

    with s.np.errstate(invalid='ignore', divide='ignore', over='ignore'):
        sol = s.root(lambda x: fun(x)[0], prediction, jac=lambda x: fun(x)[1], tol=1e-10)
    if s.np.all(s.np.isfinite(sol.x)) and s.np.max(s.np.abs(fun(sol.x)[0])) < 1e-9:
        return sol.x
    return None


def run(source, output, steps=140):
    output.mkdir(parents=True, exist_ok=True)
    known = json.loads((source/'configurations.json').read_text())
    originals = list(known)
    paths = []
    for row in originals:
        for direction in [-1, 1]:
            q = s.np.array(row['points'])
            w = s.np.r_[s.np.delete(q.ravel(), 1), 1.]
            v = tangent(w, direction=direction)
            trace = []
            ds = .025
            crossings = []
            stop = 'step-budget'
            for k in range(steps):
                new = step(w, v, ds)
                if new is None:
                    ds /= 2
                    if ds < .0001:
                        stop = 'corrector-failure'; break
                    continue
                if new[-1] <= .05 or new[-1] > 4 or abs(new[0]) < .01:
                    stop = 'mass-or-gauge-boundary'; break
                if k > 1 and (w[-1]-1)*(new[-1]-1) < 0:
                    proposal = s.np.insert(new[:-1], 1, 0.).reshape(8, 2)
                    endpoint = s.solve(proposal)
                    if endpoint is not None:
                        endpoint = s.solve(s.gauge(endpoint, rescale=False))
                    if endpoint is not None:
                        matches = [old['id'] for old in known if s.equivalent(endpoint, s.np.array(old['points']))]
                        if not matches:
                            ident = f'arc-{len(known)-len(originals)+1:03d}'
                            try:
                                c = s.central.certify(endpoint)
                                p = s.verify_polynomial.certify(endpoint)
                            except ValueError:
                                crossings.append({'step': k, 'status': 'UNRESOLVED-certificate'})
                            else:
                                for label, cert in [('direct', c), ('polynomial', p)]:
                                    (output/f'{ident}-{label}.json').write_text(json.dumps(cert, indent=2)+'\n')
                                known.append({'id': ident, 'origin': 'pseudo-arclength', 'points': endpoint.tolist()})
                                matches = [ident]
                                print(f'new certified endpoint {ident}', flush=True)
                        crossings.append({'step': k, 'matches': matches})
                w = new
                v = tangent(w, previous=v)
                trace.append({'step': k, 'mass': float(w[-1]), 'dm_ds': float(v[-1]), 'step_size': ds})
            paths.append({'source': row['id'], 'direction': direction, 'stop': stop, 'trace': trace, 'crossings': crossings})
            print(f"{row['id']} direction {direction}: {len(trace)} steps, {stop}, {len(crossings)} crossings", flush=True)
            (output/'summary.json').write_text(json.dumps({'steps_requested': steps, 'new_endpoints': known[len(originals):],
                                                          'paths': paths, 'complete': False}, indent=2)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--steps', type=int, default=140)
    a = p.parse_args()
    run(a.source, a.output, a.steps)
