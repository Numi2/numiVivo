#!/usr/bin/env python3
"""Verify all native QL families and retain descriptive original-reference comparisons."""
import argparse,gzip,hashlib,json,re
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--ql-root',type=Path,required=True)
a=p.parse_args();load=lambda p:json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
runs=load(a.root/'family-final/complete.json');inputs=load(a.root/'inputs/complete.json')
assert len(runs)==58 and {(r['case'],r['method']) for r in runs}=={(r['case'],r['method']) for r in inputs}
maxima={};features=0;maximum_rss=0;comparisons=[];shams={};old_reference={};source_inputs={}
for run in runs:
 assert run['status']=='passed' and run['referenceReturnCode']==0
 directory=a.root/'family-final'/run['case']/run['method']
 for name,key in [('input.json.gz','inputSHA256'),('native.json.gz','nativeSHA256'),('reference.json.gz','referenceSHA256')]:assert sha(directory/name)==run[key]
 native=load(directory/'native.json.gz');ref=load(directory/'reference.json.gz');inp=load(directory/'input.json.gz')
 assert native['completed'] and all(native[k] for k in ['refitsExactlyEqual','residualDFExactlyEqual','abundanceExactlyEqual','moderationExactlyEqual'])
 assert native['featureIndices']==inp['featureIndices'] and len(native['inference']['testedIndices'])==run['features']
 assert ref['status']=='passed' and all(c['passed'] and c['error']<=c['tolerance'] for c in ref['checks'])
 features+=run['features']
 for c in ref['checks']:maxima[c['name']]=max(maxima.get(c['name'],0),c['error'])
 memory=re.search(r'(\d+)\s+maximum resident set size',(directory/'native.log').read_text());assert memory
 maximum_rss=max(maximum_rss,int(memory[1]))
 if run['case'] not in old_reference:
  # Keep only this case's tables, not the much larger historical means/unit matrices.
  old=load(a.ql_root/run['case']/'reference.json.gz')
  old_reference[run['case']]={k:v['table'] for k,v in old['results'].items()}
  source_inputs[run['case']]=load(a.ql_root/run['case']/'input.json.gz')
 original=old_reference[run['case']][run['method']];data=source_inputs[run['case']]
 assert [r['featureIndex'] for r in original]==native['featureIndices']
 tests=native['inference']['tests'];calls=[g for g,t in zip(native['featureIndices'],tests) if t['adjustedPValue']<=.05]
 old_calls=[r['featureIndex'] for r in original if r['BH']<=.05]
 comp=dict(case=run['case'],method=run['method'],features=run['features'],seconds=run['seconds'],
  nativeCalls=len(calls),originalDefaultReferenceCalls=len(old_calls),originalWaldCalls=data['nativeWaldCalls'],
  originalLRTCalls=data['nativeLRTCalls'],nativeCalledIndices=calls,originalReferenceCalledIndices=old_calls,
  maximumOriginalReferencePValueDifference=max(abs(t['pValue']-r['pValue']) for t,r in zip(tests,original)),
  originalReferenceSHA256=sha(a.ql_root/run['case']/'reference.json.gz'))
 comparisons.append(comp)
 if '-null-' in run['case']:
  key=run['case'].split('-')[0]+'/'+run['method'];s=shams.setdefault(key,dict(cases=0,nativeCalls=0,originalDefaultReferenceCalls=0,originalWaldCalls=0))
  s['cases']+=1
  for k in ['nativeCalls','originalDefaultReferenceCalls','originalWaldCalls']:s[k]+=comp[k]
product=load(a.root/'product-kang/checks.json');assert product['status']=='passed'
assert sha(a.root/'product-kang/output.json.gz')==product['outputSHA256']
assert 'Test run with 50 tests in 9 suites passed' in (a.root/'swift-tests-final.log').read_text()
result=dict(status='passed',baseCommit='9da4ee0d50aea97d6b8312c3acdd4796225c124a',arms=58,features=features,
 maximumErrors=maxima,maximumResidentBytes=maximum_rss,swiftTests=50,swiftSuites=9,
 shamCallInventory=shams,product=product,
 scope='Complete native adjusted QL chain and one-constraint tests on 58 frozen full-support families, plus Kang through the existing internal pseudobulk expression owner and controlled public dataset-route tests. No legacy Poisson-bound, varying-support, effect-interval, FDR calibration, biological truth, GPU or million-cell claim.')
(a.root/'original-reference-comparisons.json').write_text(json.dumps(comparisons,sort_keys=True,indent=2)+'\n')
(a.root/'summary.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
print(json.dumps(result,sort_keys=True,indent=2))
