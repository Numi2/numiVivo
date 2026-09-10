#!/usr/bin/env python3
"""Combine bounded objective checks, native numerics and reference sensitivity."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--prior-sensitivity',type=Path);a=p.parse_args()
prior_root=a.prior_sensitivity or (a.root/'initial-complete-sensitivity' if (a.root/'initial-complete-sensitivity').exists() else None)
native=json.loads((a.root/'native-checks.json').read_text())
assert native['status']=='passed-selected-valid-profile-numerics'
references=json.loads((a.root/'runs.json').read_text());sensitivities=json.loads((a.root/'sensitivity/runs.json').read_text())
assert len(references)==len(sensitivities)==20
assert all(x.get('status')=='checked-every-eligible-reference-profile' for x in references)
assert all(x['result']['status']=='completed' for x in sensitivities)
def bh(p):
 p=np.asarray(p);order=np.argsort(p);result=np.empty_like(p)
 result[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return result
rows=[]
for record in references:
 study=record['study'];seed=record['seed'];d=a.root/study/str(seed)
 checked=pd.read_csv(d/'checked.tsv.gz',sep='\t');upper=json.loads((d/'status.json').read_text())['gridUpper']
 valid=checked.gridDispersion.between(1e-8,upper);improved=valid&(checked.gridObjectiveGain>1e-4);boundary=improved&(checked.originalDispersion<1e-6)
 sensitivity=next(s for s in sensitivities if s['study']==study and s['seed']==seed)['result']
 assert sensitivity['results']['changedGenes']==int(improved.sum())
 for method in ['baseline','candidate']:
  f=a.root/'sensitivity'/study/str(seed)/(method+'.tsv.gz')
  table=pd.read_csv(f,sep='\t',float_precision='round_trip')
  if prior_root is not None:
   prior=pd.read_csv(prior_root/study/str(seed)/(method+'.tsv.gz'),sep='\t',float_precision='round_trip')
   assert table.equals(prior),'Convergence diagnostic collection changed fitted outputs'
  assert list(table.featureID)==list(checked.featureID)
  tested=table.pValue.notna();assert np.allclose(bh(table.loc[tested,'pValue']),table.loc[tested,'adjustedPValue'],rtol=1e-10,atol=1e-12)
  assert sensitivity['results'][method]['testedGenes']==int(tested.sum())
 rows.append(dict(study=study,seed=seed,eligibleGenes=len(checked),boundedImproved=int(improved.sum()),boundedBoundaryImproved=int(boundary.sum()),
  boundedInitialFloorImproved=int((boundary&(checked.initialDispersion==1e-8)).sum()),
  largestBoundedObjectiveGain=float(checked.loc[improved,'gridObjectiveGain'].max()),outsideBounds=int((~valid).sum()),sensitivity=sensitivity))
studies=[]
protocol_failures=[]
for study in ['kang','hagai']:
 selected=[r for r in rows if r['study']==study]
 value=dict(study=study,splits=10,eligibleGeneSplitTests=sum(r['eligibleGenes'] for r in selected),
  boundedImprovedRange=[min(r['boundedImproved'] for r in selected),max(r['boundedImproved'] for r in selected)],
  boundedBoundaryImprovedRange=[min(r['boundedBoundaryImproved'] for r in selected),max(r['boundedBoundaryImproved'] for r in selected)],
  largestBoundedObjectiveGain=max(r['largestBoundedObjectiveGain'] for r in selected))
 for method in ['baseline','candidate']:
  cases=[r['sensitivity']['results'][method] for r in selected]
  value[method]=dict(totalBH005=sum(c['bhBelow005'] for c in cases),splitsWithAnyBH005=sum(c['bhBelow005']>0 for c in cases),
   nonconvergedCoefficients=sum(c['nonconvergedCoefficients'] for c in cases),trendMethods=sorted(set(c['trendMethod'] for c in cases)),
   parametricProtocolSplits=sum(c['trendMethod']=='parametric' for c in cases))
  for row in selected:
   if row['sensitivity']['results'][method]['trendMethod']!='parametric':
    protocol_failures.append(dict(study=study,seed=row['seed'],method=method,reason='DESeq2 automatically substituted a local trend; the declared parametric sensitivity did not complete for this case',messages=row['sensitivity']['messages']))
 studies.append(value)
result=dict(status='completed-with-parametric-sensitivity-failures' if protocol_failures else 'completed-profile-and-optimization-sensitivity',studies=studies,perSplit=rows,protocolFailures=protocol_failures,
 nativeProfiles=len(native['checks']),nativeValidPoints=sum(c['points'] for c in native['checks']),
 nativeMaximumErrors={key:max(c['maxErrors'][key] for c in native['checks']) for key in native['limits']},
 nativeMaximumGridGap=max(c['independentGridOverNativeProfile'] for c in native['checks']),expectedDomainRejections=native['expectedDomainRejections'],
 nativeNumericalFailures=native['failures'],referenceWarnings=sum(len(r['sensitivity']['warnings']) for r in rows),
 baselineMaximumReproductionError=max(v for r in rows for v in r['sensitivity']['results']['baselineMaximumErrors'].values()),
 convergenceDiagnosticCollectionPreservedAllTables=True if prior_root is not None else None,
 qualification='Default reference optimization is contradicted by better in-bounds objectives; finite-grid improvements do not prove global optimality, calibrated inference or power. Native inference unchanged.')
(a.root/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='perSplit'},indent=2))
