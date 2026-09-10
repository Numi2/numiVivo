#!/usr/bin/env python3
"""Verify every frozen arm and aggregate the moderation qualification evidence."""
import argparse, gzip, hashlib, json, re
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root',type=Path,required=True)
args=parser.parse_args(); root=args.root
sha=lambda raw:hashlib.sha256(raw).hexdigest()
load=lambda p:json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
runs=load(root/'family/complete.json'); inputs=load(root/'inputs/complete.json')
assert len(runs)==58 and {(r['case'],r['method']) for r in runs}=={(r['case'],r['method']) for r in inputs}
maxima={}; features=0; boundaries={'lower':0,'upper':0,'interior':0}; profiles=0; robust=0; work=0; peak=0
default_error=0; default_df_error=0; diagnostics=[]
for run in runs:
    directory=root/'family'/run['case']/run['method']
    assert run['status']=='passed' and run['referenceReturnCode']==0
    for name,key in [('input.json.gz','inputSHA256'),('native.json.gz','nativeSHA256'),('reference.json.gz','referenceSHA256')]:
        assert sha((directory/name).read_bytes())==run[key]
    inp=load(directory/'input.json.gz'); native=load(directory/'native.json.gz'); ref=load(directory/'reference.json.gz')
    fit=native['fit']; assert inp['featureIndices']==native['featureIndices']
    n=len(inp['featureIndices']); assert n==run['features']
    for name in ['posteriorVariances','priorScales','priorDegreesOfFreedom']:
        assert len(fit[name])==n
    features+=n; profiles+=len(fit['profiles']); robust+=int(len(fit['profiles'])==2); work+=fit['featureEvaluations']
    assert ref['status']=='passed' and all(c['passed'] and c['error']<=c['tolerance'] for c in ref['checks'])
    for c in ref['checks']:
        name=re.sub(r'-[12]$','',c['name']); maxima[name]=max(maxima.get(name,0),c['error'])
    for p in fit['profiles']:
        boundaries['lower' if p['lowerBoundary'] else 'upper' if p['upperBoundary'] else 'interior']+=1
    default_error=max(default_error,ref['defaultPosteriorRelativeError'])
    default_df_error=max(default_df_error,abs(fit['commonPriorDegreesOfFreedom']-ref['defaultReference']['df2'])/max(1,abs(ref['defaultReference']['df2'])))
    memory=re.search(r'(\d+)\s+maximum resident set size',(directory/'native.log').read_text()); assert memory
    peak=max(peak,int(memory[1]))
    diagnostics.append(dict(case=run['case'],method=run['method'],features=n,seconds=run['seconds'],
        defaultPosteriorRelativeError=ref['defaultPosteriorRelativeError'],
        nativeCommonPriorDF=fit['commonPriorDegreesOfFreedom'],defaultCommonPriorDF=ref['defaultReference']['df2']))
numerics=load(root/'numerics-attempt-2/checks.json')
assert numerics['status']=='passed' and len(numerics['checks'])==265
numerical_maxima={}
for kind in ['smoother','shape','tail']:
    rows=[c for c in numerics['checks'] if c['kind']==kind]; assert all(c['status']=='passed' for c in rows)
    numerical_maxima[kind]=dict(cases=len(rows),maximumError=max(c.get('maximumRelativeError',c.get('maximumAbsoluteLogError')) for c in rows))
tests=(root/'swift-tests-family.log').read_text(); assert 'Test run with 43 tests in 8 suites passed' in tests
result=dict(status='passed',baseCommit='6cb43383fce3de26b9a48361370c3cfea29c73bf',arms=58,features=features,
    profiles=profiles,robustRefits=robust,boundaries=boundaries,maximumErrors=maxima,numericalChecks=numerical_maxima,
    featureEvaluations=work,maximumResidentBytes=peak,maximumDefaultPosteriorRelativeError=default_error,
    maximumDefaultCommonDFRelativeError=default_df_error,swiftTests=43,swiftSuites=8,
    armsDiagnostics=diagnostics,
    scope='Native unequal-DF robust QL prior and posterior variance, with native abundance integration. No cohort QL hypothesis, Poisson bound, varying-support borrowing, FDR calibration, biological truth, GPU or million-cell claim.')
(root/'summary.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='armsDiagnostics'},sort_keys=True,indent=2))
