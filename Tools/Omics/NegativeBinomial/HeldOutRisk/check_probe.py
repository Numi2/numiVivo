#!/usr/bin/env python3
"""Exact measurement-harness agreement with the qualified full-cohort product."""
import argparse,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root
def load(p):return json.loads(gzip.decompress(p.read_bytes()))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
protocol=json.loads((r/'protocol.json').read_text())
for name in ['productReport','sourceReport']:assert sha(Path(protocol[name]))==protocol[name+'SHA256']
model=load(r/'probe-model.json.gz');empirical=load(Path(protocol['productReport']))['contrasts'][0];fixed=load(Path(protocol['sourceReport']))['contrasts'][0]
assert model['design']==empirical['design'] and model['trend']==empirical['negativeBinomial']['trend'] and model['prior']==empirical['negativeBinomial']['effectPriorEstimate']
count=0
for gene,feature,diag,old in zip(model['genes'],empirical['features'],empirical['negativeBinomial']['features'],fixed['negativeBinomial']['features'],strict=True):
 assert gene['status']==feature['status'] and gene['mean']==feature['meanNormalizedCount']
 if gene['status']!='tested':continue
 for method,reference in [('mle',diag['finalFit']),('empirical',diag['effectShrinkageFit']),('fixed',old['effectShrinkageFit'])]:
  assert gene[method]['coefficients']==reference['coefficients'] and gene[method]['scaledScore']==reference['maximumScaledScore']
  assert gene[method]['standardDeviation']==reference['standardError' if method=='mle' else 'posteriorStandardDeviation']
  assert gene[method]['converged']==reference['converged'];count+=1
result=dict(status='exact-production-owner-agreement',testedGenes=count//3,coefficientFits=count,designTrendAndPriorExactlyEqual=True,
 sourceProductReports=[protocol['sourceReportSHA256'],protocol['productReportSHA256']],probeSHA256=sha(r/'probe-model.json.gz'),checkerSHA256=sha(Path(__file__)))
(r/'probe-check.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps(result))
