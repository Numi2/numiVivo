from pathlib import Path
import hashlib,json,tarfile
r=Path(__file__).parent;repo=r.parents[4];p=r/'evidence.tar.gz';m=json.loads((r/'manifest.json').read_text())
assert p.stat().st_size==m['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==m['sha256']
with tarfile.open(p) as t:
 def raw(n):return t.extractfile(n).read()
 def obj(n):return json.loads(raw(n))
 manifest=obj('manifest.json');files={e['path']:e for e in manifest['files']}
 for member in t.getmembers():
  if member.name=='manifest.json':continue
  e=files[member.name];b=raw(member.name);assert len(b)==e['bytes'] and hashlib.sha256(b).hexdigest()==e['sha256']
 for n,h in manifest['sourceHashes'].items():assert hashlib.sha256((repo/n).read_bytes()).hexdigest()==h,n
 for n in ['Check.swift','Controls.swift','cli.py','check_graph.py']:assert (r/n).read_bytes()==raw(n)
 c=obj('cli-terminal.json');assert len(c['results'])==10 and all(x['exitCode']==(65 if x['step'].startswith('reject') else 0) for x in c['results'])
 g=obj('graph-check.json');assert g['cells']==13863 and g['allNeighborMembershipExact'] and g['allEdgeCoordinatesExact'] and g['metalBinaryJSONNeighborsExact'] and not g['biologicalQualification']
 assert g['maximumBinaryJSONWeightDifference']==0 and g['maximumFuzzyWeightDifferenceVsCPU']<1e-5
 assert obj('owner-output/reference-check.json')['maximumDistanceError']==0
 assert 'cancellation PASS' in raw('owner.log').decode()
 assert 'duplicate ties and partial query/candidate tiles PASS' in raw('controls.log').decode()
print('PASS: current native Metal sources and retained full-cohort owner/CLI qualification')
