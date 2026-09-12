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
 for n in ['Graph.swift','run.py','scanpy_graph.py','check.py']:assert raw(n)==(r/n).read_bytes()
 v=obj('verification.json');assert hashlib.sha256(raw('protocol.json')).hexdigest()==v['protocolSHA256']
 assert v['cells']==24673 and v['repetitionsExact'] and not v['fullCLISpeedupClaim'] and not v['biologicalQualification']
 for c in v['comparisons'].values():assert c['missingScanpyNeighbors']==0 and c['scanpyNonSelfNeighborEntries']==468787 and c['nativeEdges']==c['scanpyEdges']==c['commonEdges']==706016
 runs=obj('terminal.json')['results'];assert len(runs)==9 and all(x['exitCode']==0 for x in runs)
 for backend in ['cpu','metal','scanpy']:
  assert statistics.median(x['graphStageSeconds'] for x in runs if x['backend']==backend)==v['timings'][backend]['medianGraphStageSeconds']
print('PASS: retained complete exact-Scanpy graph comparison; stage timing only')
