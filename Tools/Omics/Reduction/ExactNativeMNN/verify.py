from pathlib import Path
import hashlib,json,tarfile
r=Path(__file__).parent;repo=r.parents[3];m=json.loads((r/'manifest.json').read_text())
assert hashlib.sha256((r/'evidence.tar.gz').read_bytes()).hexdigest()==m['archiveSHA256']
with tarfile.open(r/'evidence.tar.gz') as t:
 names=[x.name for x in t.getmembers() if x.isfile()]
 assert len(names)==len(set(names)) and set(names)==set(m['members'])
 data={n:t.extractfile(n).read() for n in names}
for n,b in data.items():assert hashlib.sha256(b).hexdigest()==m['members'][n]
v=json.loads(data['verification.json']);assert v==json.loads((r/'verification.json').read_text())
assert v['sanitizerTestsPassed'] and v['swiftCancellationPassed']
assert {x['cohort']:x['cells'] for x in v['results']}=={'hagai':13863,'kang':24673,'ding':44031}
for x in v['results']:
 assert x['scoresBitwiseExact'] and x['reportExact'] and x['maximumCoordinateError']==0
 for name,h in x['outputs'].items():assert h==x['inputs']['mnn/'+name]
for n,h in v['sourceSHA256'].items():
 assert m['members']['source/'+n]==h
 assert hashlib.sha256((repo/n).read_bytes()).hexdigest()==h,n
assert v['cliSHA256']==m['externalArtifacts']['build/numivivo-omics']['SHA256']
assert v['ownerDriverSHA256']==m['externalArtifacts']['owner-check']['SHA256']
assert hashlib.sha256((r/'Test.cpp').read_bytes()).hexdigest()==v['testSources']['Test.cpp']
assert 'PCA score read tile' in data['driver-attempt-1/owner.log'].decode()
print('PASS current owner source, archive, CLI identity, complete byte-equivalence and retained test/failure accounting')
