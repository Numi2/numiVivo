from pathlib import Path
import tarfile,json,hashlib,statistics
r=Path(__file__).parent;p=r/'evidence.tar.gz';m=json.loads((r/'manifest.json').read_text())
assert p.stat().st_size==m['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==m['sha256']
with tarfile.open(p) as t:
 def raw(n):return t.extractfile(n).read()
 def obj(n):return json.loads(raw(n))
 files={e['path']:e for e in obj('external-manifest.json')['files']}
 for member in t.getmembers():
  if member.name=='external-manifest.json':continue
  e=files[member.name];b=raw(member.name);assert len(b)==e['bytes'] and hashlib.sha256(b).hexdigest()==e['sha256']
 for n in ['run.py','check.py','benchmark.py']:assert (r/n).read_bytes()==raw(n)
 v=obj('verification.json');assert hashlib.sha256(raw('protocol.json')).hexdigest()==v['protocolSHA256']
 assert v['status']=='PASS' and v['cells']==24673 and v['missingNeighbors']==v['rowsWithOrderingOrMembershipChange']==0 and not v['biologicalGeneralizationQualified']
 assert v['annotationConcordance']['cpu']==v['annotationConcordance']['metal']
 for seed in [7,19,42]:
  a=obj('cpu-'+str(seed)+'/result.json');b=obj('metal-'+str(seed)+'/result.json');assert a['cells']==b['cells'] and len(a['labels'])==len(b['labels'])==24673
  forward={};reverse={}
  for x,y in zip(a['labels'],b['labels']):assert forward.setdefault(x,y)==y and reverse.setdefault(y,x)==x
 assert len(obj('terminal.json')['results'])==15 and all(x['exitCode']==0 for x in obj('terminal.json')['results'])
 bench=obj('benchmark/terminal.json');assert len(bench['results'])==6 and all(x['exitCode']==0 and x['allEightFilesMatch'] for x in bench['results'])
 for mode in ['cpu','metal']:assert statistics.median(x['processWallSeconds'] for x in bench['results'] if x['backend']==mode)==bench['medianSeconds'][mode]
print('PASS: complete Kang preservation and repeated native CLI evidence; generalization unqualified')
