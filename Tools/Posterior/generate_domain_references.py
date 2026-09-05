"""Independent finite-pool ODE reference, not a replacement simulation runtime."""
import json
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp

rates = dict(kon=1e6, koff=0.2, kinact=0.1, turnover=0.003,
             clearance=0.01, metabolism=0.02, metaboliteClearance=0.005, baseline=2e-6)
initial = [1e-6, 2e-6, 0, 0, 0, 0, 0, 0]

def rhs(t, s):
    d, f, r, c, m, eliminated, synthesized, removed = s
    association=rates['kon']*d*f
    dissociation=rates['koff']*r
    conversion=rates['kinact']*r
    turnover=rates['turnover']
    return [-association+dissociation-(rates['clearance']+rates['metabolism'])*d,
            -association+dissociation-turnover*f+turnover*rates['baseline'],
            association-dissociation-conversion-turnover*r,
            conversion-turnover*c,
            rates['metabolism']*d+turnover*(r+c)-rates['metaboliteClearance']*m,
            rates['clearance']*d+rates['metaboliteClearance']*m,
            turnover*rates['baseline'], turnover*(f+r+c)]

values={}
for method in ['DOP853', 'Radau']:
    result=solve_ivp(rhs,[0,100],initial,method=method,rtol=1e-11,atol=1e-19)
    assert result.success
    values[method]=result.y[:,-1].tolist()
assert np.max(np.abs(np.array(values['DOP853'])-np.array(values['Radau']))) < 1e-16
(root:=Path(__file__).resolve().parent).joinpath('domain-reference.json').write_text(json.dumps({
    'generator':'SciPy solve_ivp, independent ODE equations in generate_domain_references.py',
    'scipy':scipy.__version__,'numpy':np.__version__,'rates':rates,'initial':initial,'durationSeconds':100,
    'expectedStates':values},indent=2)+'\n')
