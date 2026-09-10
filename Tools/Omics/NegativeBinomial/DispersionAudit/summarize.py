#!/usr/bin/env python3
"""Summarize the completed stage audit without changing fitted inference."""
import argparse,gzip,json
from pathlib import Path
import numpy as np
import pandas as pd
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args()
runs=json.loads((a.root/'runs.json').read_text());curves=json.loads((a.root/'trend-decomposition/trends.json').read_text())
assert len(runs)==20 and all(r.get('status')=='reproduced-and-audited' for r in runs)
assert len(curves)==80
summary=[]
for study in ['kang','hagai']:
 selected=[r for r in runs if r['study']==study]
 def span(values):return [float(min(values)),float(max(values))]
 value=dict(study=study,splits=len(selected))
 for name,values in [
   ('nativePriorRange',[r['nativeTrend']['priorLogVariance'] for r in selected]),
   ('referencePriorRange',[r['reference']['priorLogVariance'] for r in selected]),
   ('nativeResidualVarianceRange',[r['nativeTrend']['robustLogResidualVariance'] for r in selected]),
   ('referenceResidualVarianceRange',[r['reference']['residualLogVariance'] for r in selected]),
   ('medianReferenceOverNativeTrendRange',[r['dispersionRatios']['trendDispersion']['median'] for r in selected]),
   ('medianReferenceOverNativeFinalRange',[r['dispersionRatios']['finalDispersion']['median'] for r in selected])]:value[name]=span(values)
 value['stageCoverage']=[];value['trendDecomposition']=[]
 for seed in range(1,11):
  directory=a.root/study/str(seed);n=json.loads(gzip.decompress((directory/'native-stages.json.gz').read_bytes()))
  nt=pd.DataFrame(n['features']).set_index('featureID')
  rt=pd.read_csv(directory/'stages.tsv.gz',sep='\t',index_col='featureID')
  joint=nt.loc[rt.index];available=joint.geneWiseDispersion.notna()
  value['stageCoverage'].append(dict(seed=seed,eligibleGenes=len(rt),nativeGeneWiseFits=int(available.sum()),
   nativeReferenceGenes=len(n['trend']['referenceFeatureIndices']),referenceAboveFloor=int((rt.geneWiseDispersion>=1e-6).sum()),
   nativeAboveFloor=int((joint.geneWiseDispersion>=1e-6).sum()),
   referenceBelowFloorButNativeAbove=int(((rt.geneWiseDispersion<1e-6)&(joint.geneWiseDispersion>=1e-6)).sum())))
  for c in [c for c in curves if c['study']==study and c['seed']==seed]:
   row=dict(seed=seed,case=c['case'],status=c['status'],genes=c['genes'])
   if c['status']=='completed':
    coefficients=c['coefficients'];f=coefficients['asymptDisp']+coefficients['extraPois']/rt.meanNormalizedCount
    native=n['trend']['intercept']+n['trend']['inverseMeanCoefficient']/rt.meanNormalizedCount
    row.update(medianCurveOverNative=float(np.median(f/native)),medianCurveOverReference=float(np.median(f/rt.trendDispersion)),maximumCurveRelativeDifference=c['maximumCurveRelativeDifference'])
   value['trendDecomposition'].append(row)
 summary.append(value)
calls=json.loads((a.root/'null-calls.json').read_text())
result=dict(status='completed-stage-diagnosis',studies=summary,referenceReproductionMaximumAbsoluteError=max(v for r in runs for v in r['referenceReproductionErrors'].values()),
 referenceWarnings=sum(len(r['reference']['warnings']) for r in runs),trendWarnings=sum(len(c['warnings']) for c in curves),
 trendFailures=[c for c in curves if c['status']!='completed'],calls=len(calls),nativeCallOutliers=sum(c['nativeOutlier'] for c in calls),referenceCallOutliers=sum(c['referenceOutlier'] for c in calls),
 qualification='Post-score diagnosis of frozen overlapping null splits. No revised inference, general FDR or power claim.')
(a.root/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='studies'},indent=2))
