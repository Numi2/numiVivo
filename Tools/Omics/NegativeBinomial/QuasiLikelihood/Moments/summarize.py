#!/usr/bin/env python3
"""Keep numerical correctness, full-family coverage and work failures separate."""
import argparse,json,re
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root;final=root/'final'
runs=json.loads((final/'complete.json').read_text());independent=json.loads((final/'independent-summary.json').read_text());grid=json.loads((final/'grid-check.json').read_text());sensitivity=json.loads((root/'work-sensitivity/summary.json').read_text());equality=json.loads((final/'requalification-equality.json').read_text())
assert len(runs)==58 and independent['arms']==58 and len(equality)==58 and all(r['exactlyEqual'] for r in equality)
assert grid['status']=='passed'
numerical=[];unavailable=[];memory=[]
for r in runs:
 d=final/r['case']/r['method'];check=json.loads((d/'independent-check.json').read_text());assert check['outputSHA256']==r['outputSHA256']
 for f in check['failures']:
  if 'error' in f and 'exhausted' in f['error']:unavailable.append(dict(case=r['case'],method=r['method'],**f))
  else:numerical.append(dict(case=r['case'],method=r['method'],**f))
 for f in r['failures']:
  if not ('error' in f and 'exhausted' in f['error']):numerical.append(dict(case=r['case'],method=r['method'],**f))
 match=re.search(r'(\d+)\s+maximum resident set size',(d/'stderr.log').read_text())
 if match:memory.append(int(match[1]))
result=dict(status='completed-with-resource-limit-failures' if unavailable else 'completed',attemptedArms=58,fullyAvailableArms=sum(not r['failures'] for r in runs),defaultUnavailableGeneArms=len(unavailable),defaultUnavailable=unavailable,defaultMoments=independent['moments'],defaultSupportEvaluations=independent['supportTerms'],maximumRelativeMomentError=independent['maximumRelativeMomentError'],maximumReferencePMFMassError=independent['maximumPMFMassError'],maximumLeverageError=max(r['maximumLeverageError'] for r in runs),maximumRelativeTruncationBound=max(r['maximumRelativeTruncationBound'] for r in runs),maximumRelativeReferenceDevianceDifference=max(r['relativeDevianceDifferenceQuantiles'][-1] for r in runs),maximumRelativeReferenceDFDifference=max(r['relativeDFDifferenceQuantiles'][-1] for r in runs),independentNumericalFailures=numerical,gridCases=grid['cases'],gridMaximumRelativeError=grid['maximumRelativeError'],higherWorkSensitivity={k:v for k,v in sensitivity.items() if k!='checks'},nativePeakResidentBytesRange=[min(memory),max(memory)] if memory else None,requalifiedByteIdenticalArms=58,reusedExactOutputChecks=sum(r['independentCheckReused'] for r in equality),nativeBinarySHA256=runs[0]['binarySHA256'],qualification='Native conditional moment and residual adjustment only. Twenty primary work-limit failures retained; supplemental higher-work results do not change default coverage. One higher-work fit remains unavailable. Global QL scale, robust prior and cohort test remain unimplemented; no calibration or performance superiority claim.')
if numerical:result['status']='completed-with-numerical-and-resource-failures'
(root/'summary.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['defaultUnavailable','higherWorkSensitivity','independentNumericalFailures']},indent=2))
