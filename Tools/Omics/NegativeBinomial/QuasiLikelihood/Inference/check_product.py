#!/usr/bin/env python3
"""Bind the real pseudobulk owner result to original counts and independent R tests."""
import argparse,gzip,hashlib,json,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',type=Path,required=True);p.add_argument('--frozen-input',type=Path,required=True)
p.add_argument('--source-sha256',required=True)
a=p.parse_args();load=lambda p:json.loads(gzip.decompress(p.read_bytes()))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(a.root/'source.json.gz')==a.source_sha256
source=load(a.root/'source.json.gz');out=load(a.root/'output.json.gz');frozen=load(a.frozen_input)
assert all(out[k] for k in ['sourceDesignExactlyEqual','sourceTrendExactlyEqual','sourceEvidenceExactlyEqual'])
result=out['expression'];q=result['negativeBinomial']['quasiLikelihood'];indices=q['featureIndices']
assert q['inference']['completed'] and not q['failures']
assert indices==frozen['featureIndices'] and result['design']['rows']==frozen['design']
assert result['design']['contrast']==frozen['contrast']
matrix=source['pseudobulk']['matrix'];counts=[[0]*len(result['design']['rows']) for _ in indices];position={g:j for j,g in enumerate(indices)}
for j,row in enumerate(result['design']['sourcePseudobulkIndices']):
 for k in range(matrix['rowOffsets'][row],matrix['rowOffsets'][row+1]):
  g=matrix['featureIndices'][k]
  if g in position:counts[position[g]][j]=matrix['counts'][k]
assert counts==frozen['counts']
assert result['testedFeatures']==len(indices)
for j,g in enumerate(indices):
 feature=result['features'][g];test=q['inference']['tests'][j]
 assert feature['status']=='tested' and feature['pValue']==test['pValue'] and feature['adjustedPValue']==test['adjustedPValue']
 assert feature['fStatistic']==test['fStatistic']
 assert all(feature.get(k) is None for k in ['zStatistic','tStatistic','standardError','intervalLower','intervalUpper'])
comparison=dict(sourceSHA256=a.source_sha256,outputSHA256=sha(a.root/'output.json.gz'),frozenInputSHA256=sha(a.frozen_input),
 countsExactlyEqual=True,designExactlyEqual=True,featureIndicesExactlyEqual=True,
 features=len(result['features']),testedFeatures=result['testedFeatures'],
 calls=sum(f.get('adjustedPValue') is not None and f['adjustedPValue']<=.05 for f in result['features']),seconds=out['seconds'])
with (a.root/'reference.log').open('wb') as log:
 run=subprocess.run(['/opt/homebrew/bin/Rscript',str(Path(__file__).with_name('product_reference.R')),
  str(a.root/'output.json.gz'),str(a.frozen_input),str(a.root/'reference.json')],stdout=log,stderr=log,
  env={**os.environ,'R_LIBS_USER':'/Users/home/numivivo-r-library-20260909'})
comparison['referenceReturnCode']=run.returncode
reference=json.loads((a.root/'reference.json').read_text());comparison['referenceStatus']=reference['status']
comparison['status']='passed' if run.returncode==0 and reference['status']=='passed' else 'failed'
(a.root/'checks.json').write_text(json.dumps(comparison,sort_keys=True,indent=2)+'\n')
print(json.dumps(comparison,sort_keys=True,indent=2))
raise SystemExit(0 if comparison['status']=='passed' else 1)
