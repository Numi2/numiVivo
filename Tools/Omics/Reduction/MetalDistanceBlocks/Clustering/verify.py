from pathlib import Path
import tarfile,json,hashlib
r=Path(__file__).parent;m=json.loads((r/'manifest.json').read_text());p=r/'evidence.tar.gz'
assert p.stat().st_size==m['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==m['sha256']
with tarfile.open(p) as t:
 def raw(n):return t.extractfile(n).read()
 def obj(n):return json.loads(raw(n))
 files={e['path']:e for e in obj('external-manifest.json')['files']}
 for member in t.getmembers():
  if member.name=='external-manifest.json':continue
  e=files[member.name];b=raw(member.name);assert len(b)==e['bytes'] and hashlib.sha256(b).hexdigest()==e['sha256']
 for n in ['run.py','check.py']:assert raw(n)==(r/n).read_bytes()
 protocol=obj('protocol.json');v=obj('verification.json');assert hashlib.sha256(raw('protocol.json')).hexdigest()==v['protocolSHA256']
 assert len(obj('terminal.json')['results'])==12 and all(x['exitCode']==0 for x in obj('terminal.json')['results'])
 assert protocol['seeds']==[7,19,42] and v['status']=='PASS' and not v['biologicalAccuracyQualified']
 for seed in protocol['seeds']:
  a=obj('cpu-'+str(seed)+'/result.json');b=obj('metal-'+str(seed)+'/result.json');assert a['cells']==b['cells'] and len(a['labels'])==len(b['labels'])==13863
  forward={};reverse={}
  for x,y in zip(a['labels'],b['labels']):
   assert forward.setdefault(x,y)==y and reverse.setdefault(y,x)==x
 assert all(x['ARI']==1 and x['matchedAssignmentAgreement']==1 for x in v['comparisons'])
 assert all(x['modularityOracleError']<1e-10 for x in v['nativeChecks'])
print('PASS: all three complete CPU/Metal partitions equivalent; biological accuracy remains unqualified')
