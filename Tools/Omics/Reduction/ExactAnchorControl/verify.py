from pathlib import Path
import hashlib,json,tarfile
r=Path(__file__).parent;m=json.loads((r/'manifest.json').read_text())
assert hashlib.sha256((r/'evidence.tar.gz').read_bytes()).hexdigest()==m['archiveSHA256']
with tarfile.open(r/'evidence.tar.gz') as t:
 names=[p.name for p in t.getmembers() if p.isfile()]
 assert len(names)==len(set(names)) and set(names)==set(m['members'])
 data={n:t.extractfile(n).read() for n in names}
for n,b in data.items():assert hashlib.sha256(b).hexdigest()==m['members'][n]
read=lambda n:json.loads(data[n])
for n in ['results.json','changed-anchors.json']:assert read(n)==json.loads((r/n).read_text())
assert read('terminal.json')['returnCode']==0
f=read('freeze.json')
for n,h in f['sourceSHA256'].items():assert m['members'][n]==h
assert f['protocolSHA256']==m['members']['protocol.md']
s=read('results.json');assert s['allCohortsChecked'] and s['allCoordinatesPassed'] and s['freezeSHA256']==m['members']['freeze.json']
assert {x['cohort']:x['cells'] for x in s['results']}=={'hagai':13863,'kang':24673,'ding':44031}
assert all(x['passed'] and x['orderExact'] and x['maximumCoordinateError']<1e-10 for x in s['results'])
for c in read('changed-anchors.json')['results']:
 for kind in ['allExact','missing','added']:
  d=c[kind];assert sum(d['endpointOccurrences'].values())==2*d['anchors']
  assert sum(d['uniqueEndpointLabels'].values())==d['uniqueEndpointCells']
print('PASS all retained hashes, full-coordinate comparison results and endpoint accounting')
