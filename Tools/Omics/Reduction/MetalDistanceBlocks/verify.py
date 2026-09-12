from pathlib import Path
import tarfile,json,hashlib
r=Path(__file__).parent;p=r/'evidence.tar.gz';m=json.loads((r/'manifest.json').read_text())
assert p.stat().st_size==m['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==m['sha256']
with tarfile.open(p) as t:
 def raw(n):return t.extractfile(n).read()
 def obj(n):return json.loads(raw(n))
 entries={(e['label']+'/'+e['path']):e for e in obj('external-manifest.json')}
 for member in t.getmembers():
  if member.name=='external-manifest.json':continue
  e=entries[member.name];b=raw(member.name);assert len(b)==e['bytes'] and hashlib.sha256(b).hexdigest()==e['sha256']
 for name in ['Main.swift','run.py','check.py']:assert (r/name).read_bytes()==raw('hoisted/'+name)
 for label in ['initial','hoisted']:
  q=obj(label+'/verification.json');assert q['allCells']==13863 and q['cpuMetalIndicesExact'] and len(q['timings'])==6
  assert all(z['missingNeighbors']==0 and z['changedRows']==1 for z in q['comparison'].values())
  assert all(z['maximumScaledSquaredDistanceError']<1e-5 for z in q['comparison'].values())
 assert obj('hoisted/verification.json')['medianProcessSeconds']['metal']<obj('hoisted/verification.json')['medianProcessSeconds']['cpu']
print('PASS: retained full-cohort research evidence; product and biological qualification remain open')
