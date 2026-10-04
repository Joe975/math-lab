import importlib.util
import json
import math
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXPLORE = ROOT/'problems/central-configurations-8/explore'
spec = importlib.util.spec_from_file_location('cc8_animation', EXPLORE/'build_animation.py')
animation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(animation)


def test_animation_data_and_render_reproducibility():
    rows = animation.configuration_data()
    assert len(rows) == 20
    assert all(len(row['points']) == 8 for row in rows)
    for row in rows:
        assert max(abs(sum(p[k] for p in row['points'])) for k in range(2)) < 1e-8
    template = (EXPLORE/'rotation-template.html').read_text(encoding='utf-8')
    expected = template.replace('__CONFIGURATION_DATA__', json.dumps(rows)).replace(
        '__ROTATION_MATH__', (EXPLORE/'rotation.js').read_text(encoding='utf-8'))
    actual = (EXPLORE.parent/'data/rotating-configurations.html').read_text(encoding='utf-8')
    assert actual == expected
    assert len(actual.encode()) < 1_000_000


def test_real_animation_math_preserves_newtonian_relative_equilibrium():
    node = shutil.which('node')
    if not node:
        pytest.skip('Node required to execute actual animation JavaScript')
    script = r"""
const fs = require('node:fs');
const {rotateConfiguration, advanceRotation} = require(process.argv[1]);
const rows = JSON.parse(fs.readFileSync(0, 'utf8'));
const out = rows.map(r => ({before:r.points, after:rotateConfiguration(r.points, 1.234)}));
console.log(JSON.stringify({out, paused:advanceRotation(1,3,2,false),
  turn:advanceRotation(0,10,1,true), half:advanceRotation(0,5,1,true)}));
"""
    result = subprocess.run([node, '-e', script, str(EXPLORE/'rotation.js')],
                            input=json.dumps(animation.configuration_data()), capture_output=True,
                            text=True, check=True)
    data = json.loads(result.stdout)
    assert data['paused'] == 1
    assert abs(data['turn']) < 1e-12
    assert abs(data['half']-math.pi) < 1e-12
    for case in data['out']:
        p, q = case['before'], case['after']
        for i in range(8):
            for j in range(i):
                assert abs(math.dist(p[i], p[j])-math.dist(q[i], q[j])) < 1e-12
            acceleration = [sum((q[j][k]-q[i][k])/math.dist(q[i],q[j])**3
                                for j in range(8) if i != j) for k in range(2)]
            assert max(abs(acceleration[k]+q[i][k]) for k in range(2)) < 1e-8
