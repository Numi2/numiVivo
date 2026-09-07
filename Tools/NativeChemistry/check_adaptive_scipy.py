#!/usr/bin/env python3
"""Independent SciPy matrix-exponential/Frechet oracle for native test output.

SciPy is an optional external test dependency, never a NumiVivo runtime.
"""
import json
import pathlib
import sys
import numpy as np
import scipy
from scipy.linalg import expm, expm_frechet

path = pathlib.Path(sys.argv[1])
data = json.loads(path.read_text())
# These values are independently declared by the fixture, not fitted to native
# probabilities, derivatives, eigenvectors, residuals, or convergence flags.
q = np.array([[-0.8, 0.5], [0.2, -0.6]], dtype=np.float64)
b = np.array([[-0.5, 0.5], [0.0, 0.0]], dtype=np.float64)
p = np.array([0.6, 0.4], dtype=np.float64)
checks = []
for point in data['result']['observations']:
    time = point['timeSeconds']
    reference = p @ expm(q * time)
    derivative = p @ expm_frechet(q * time, b * time, compute_expm=False)
    actual = np.array(point['probabilityByState'])
    actual_d = np.array(point['derivatives'][0]['probabilityByState'])
    err = float(np.max(np.abs(actual-reference)))
    derr = float(np.max(np.abs(actual_d-derivative)))
    assert err < 3e-11, (time, err)
    assert derr < 3e-10, (time, derr)
    checks.append(dict(timeSeconds=time, maximumProbabilityError=err,
                       maximumDerivativeError=derr))
assert len(checks) == 4
report = dict(schema='numivivo.org/kinetic-scipy-conformance/v1',
              numpyVersion=np.__version__, scipyVersion=scipy.__version__,
              source='independently declared two-state generator and rate derivative',
              checks=checks, status='passed',
              scope='finite-state numerical oracle, not molecular or experimental validation')
path.with_name('scipy-kinetic-comparison.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
