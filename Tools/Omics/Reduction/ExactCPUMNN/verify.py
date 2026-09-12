"""Verify retained evidence and complete cohort accounting, not fresh fitting."""
from pathlib import Path
import hashlib,json,tarfile
r=Path(__file__).parent;m=json.loads((r/'manifest.json').read_text())
assert hashlib.sha256((r/'evidence.tar.gz').read_bytes()).hexdigest()==m['archiveSHA256']
with tarfile.open(r/'evidence.tar.gz') as t:
 names=[x.name for x in t.getmembers() if x.isfile()]
 assert len(names)==len(set(names)) and set(names)==set(m['members'])
 data={n:t.extractfile(n).read() for n in names}
for n,b in data.items():assert hashlib.sha256(b).hexdigest()==m['members'][n],n
read=lambda n:json.loads(data[n])
assert read('terminal.json')['returnCode']==read('checks-terminal.json')['returnCode']==0
s=read('summary.json');assert s==json.loads((r/'summary.json').read_text())
f=read('candidate/freeze.json');out=read('candidate/outputs-frozen.json')
assert out['allOriginalCells'] and out['allReplayArraysExact'] and out['metricsRead'] is False
for n,h in f['sourceSHA256'].items():assert hashlib.sha256(data['fit-source/'+n]).hexdigest()==h
cpp=data['fit-source/LocalMNN/LocalNeighbors.cpp'].decode()
assert cpp.count('index.setEf(512)')==2 and 'index.setEf(128)' not in cpp
assert f['librarySHA256']==m['members']['libExact.dylib']
for n,h in out['files'].items():
 path='candidate/'+n
 assert h==(m['members'][path] if path in m['members'] else m['externalArrays'][path]['SHA256'])
expected={'hagai':13863,'kang':24673,'ding':44031}
assert {x['cohort']:x['cells'] for x in s['cohorts']}==expected
for x in s['cohorts']:
 e=read('evaluation/'+x['cohort']+'/checks.json')
 assert x['gates']==e['gates'] and x['cells']==e['cells']
 n=next(v for v in read('exact-check.json')['results'] if v['cohort']==x['cohort'])
 assert n['allAnchorsExact'] and n['orderExact'] and n['maximumCoordinateError']<1e-10
 assert n['matchingDistanceEvaluations']==x['matchingDistanceEvaluations']
 assert all(value for key,value in e['gates'].items() if key!='completeClassifierStrata')
assert read('evaluation/kang/checks.json')['gates']['completeClassifierStrata'] is False
assert len(read('evaluation/kang/checks.json')['metrics']['missingClassifierStrata'])==4
assert read('exact-check.json')['checkerSHA256']==m['members']['check_exact.py']
assert read('exact-check.json')['outputsFreezeSHA256']==m['members']['candidate/outputs-frozen.json']
print('PASS all source/output identities and complete 82,567-cell accounting; inspect individual biological gates')
