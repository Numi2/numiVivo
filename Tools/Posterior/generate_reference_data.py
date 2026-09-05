"""Independent numerical fixtures, not biological measurements.
Requires mpmath and NumPy; production NumiVivo does not depend on either.
"""
import json, math
from pathlib import Path
import mpmath as mp
import numpy as np
mp.mp.dps = 100
cdf = lambda z: mp.erfc(-mp.mpf(z)/mp.sqrt(2))/2
normal = [(x,float(mp.log(cdf(x)))) for x in [-100,-40,-12,-10,-9,-1,0,1,9]]
intervals=[]
for a,b in [(-12,-10),(10,11),(-1,1),(0,1e-14),(-1e-8,1e-8),(30,30.0000000001)]:
    if a>0: v=mp.log(cdf(-a)-cdf(-b))
    else: v=mp.log(cdf(b)-cdf(a))
    intervals.append([a,b,float(v)])
times=np.array([0.,1.,1.]); values=np.array([.4,.2,.7]); means=np.array([.3,.1,.5]); sds=np.array([.1,.2,.3])
sd=np.hypot(sds*1.2,.03)
r=.7*np.exp(-abs(times[:,None]-times[None,:])/2.)
np.fill_diagonal(r,1)
cov=sd[:,None]*r*sd[None,:]
z=values-means-.02
ll=float(-.5*z@np.linalg.solve(cov,z)-.5*np.linalg.slogdet(cov)[1]-1.5*np.log(2*np.pi))
Path('Tools/Posterior/reference-data.json').write_text(json.dumps({'normalLogCDF':normal,'intervalLogProbability':intervals,'correlatedLogLikelihood':ll,'generators':{'mpmath':mp.__version__,'numpy':np.__version__,'precisionDecimalDigits':100}},indent=2)+'\n')
# A separate, explicitly synthetic end-to-end problem with measured precision
# left unknown and an estimable residual SD. Observations from analytic binding,
# not from the native forward implementation under test.
evidence={'source':'synthetic analytic binding fixture','locator':'Tools/Posterior/generate_reference_data.py'}
context={'compound':'synthetic-ligand','target':'synthetic-target','targetVariant':'reference','site':'one-site','chemicalState':'synthetic-state','hostContext':'well-mixed-fixture','temperatureK':300,'pH':7,'ionicStrengthM':.1}
param=lambda x,u:dict(value=x,unit=u,origin='assumed',uncertainty={'unknown':{}},evidence=evidence)
kinetics=dict(schemaVersion=1,identifier='synthetic-assay',context=context,association=param(1e5,'M^-1 s^-1'),dissociation=param(.2,'1/s'),inactivation=param(0,'1/s'),targetTurnover=param(0,'1/s'),baselineTarget=param(1e-7,'M'),maximumUnboundDrugM=1e-3,exposureTreatment='externallyMaintainedUnbound',turnoverModel='equalStateLossAndConstantSynthesis')
cases=[]
# Fixed deviations are a synthetic fixture, not draws used for SBC.
errors=[.002,-.007,.004,-.003,.009]
for idx,(name,role,L) in enumerate([('fit-low','calibration',1e-6),('fit-high','calibration',2e-6),('heldout','test',1.5e-6)]):
    times=[1,2,5,10,20]; a=1e5*L
    f=lambda t: a/(a+.2)*(-math.expm1(-(a+.2)*min(t,10)))*math.exp(-.2*max(0,t-10))
    trace=dict(schemaVersion=1,context=context,interpolation='rightContinuousPiecewiseConstant',knots=[dict(timeSeconds=0,unboundDrugM=L),dict(timeSeconds=10,unboundDrugM=0),dict(timeSeconds=20,unboundDrugM=0)],origin='assumed',evidence=evidence)
    exp=dict(schemaVersion=1,kinetics=kinetics,exposure=trace,initial=dict(free=1,reversible=0,covalent=0,competitor=0),sampleTimesSeconds=[0]+times)
    obs=[dict(identifier=f'o{t}',timeSeconds=t,observable='drugOccupancy',value=f(t)+errors[(j+idx)%5],origin='assumed',evidence=evidence) for j,t in enumerate(times)]
    cases.append(dict(identifier=name,leakageGroup=name,partition=role,experiment=exp,observations=obs))
ids=[c['identifier'] for c in cases]
bind=lambda name,field,unit,lo,hi:dict(parameter=dict(identifier=name,unit=unit,lower=lo,upper=hi,prior='uniformLogPhysical'),field=field,caseIdentifiers=ids,priorEvidence=evidence)
problem=dict(schemaVersion=1,study=dict(schemaVersion=1,identifier='synthetic-unknown-assay-sd',cases=cases),bindings=[bind('association-rate','association','M^-1 s^-1',3e4,3e5),bind('dissociation-rate','dissociation','1/s',.03,1),bind('assay-sd','assayNoiseFloor','fraction',.002,.04)],sampler=dict(particleCount=128,maximumStages=128,mutationSweeps=6,targetESSFraction=.8,proposalScale=1,covarianceRegularization=1e-6,independentPriorProposalProbability=.1,minimumTemperatureIncrement=1e-12,maximumLikelihoodEvaluations=50000,seed=20260904),numerics=dict(maximumMatrixSquarings=48,maximumPropagations=1048576,fractionConservationTolerance=1e-8),parallelEvaluations=4)
Path('Examples/target-engagement/synthetic-assay-inference.json').write_text(json.dumps(problem,separators=(',',':'))+'\n')
