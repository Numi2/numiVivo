"""Verify retained execution evidence; does not refit or rerun biological scoring."""
from pathlib import Path
import json,gzip,hashlib,tarfile,statistics
r=Path(__file__).parent
selection=json.loads(gzip.decompress((r/'selection.json.gz').read_bytes()))
for filename,key in [('protocol.json','protocolSHA256'),('select_regularization.py','selectorSHA256')]:
 assert hashlib.sha256((r/filename).read_bytes()).hexdigest()==selection[key]
assert len(selection['folds'])==150
for f in selection['folds']:
 assert len(f['losses'])==5
 for loss in f['losses']:
  assert len(loss['innerTargets'])==29 and f['target'] not in loss['innerTargets']
  assert len(loss['innerRMSE'])==29
  assert abs(statistics.mean(loss['innerRMSE'])-loss['meanRMSE'])<1e-12
 assert min(f['losses'],key=lambda x:(x['meanRMSE'],-x['regularization']))['regularization']==f['selectedRegularization']==1
m=json.loads((r/'manifest.json').read_text());p=r/'evidence.tar.gz'
assert p.stat().st_size==m['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==m['sha256']
with tarfile.open(p) as t:
 def raw(n):return t.extractfile(n).read()
 def obj(n):return json.loads(raw(n))
 external={x['path']:x for x in obj('external-manifest.json')['files']}
 for member in t.getmembers():
  if member.name=='external-manifest.json':continue
  e=external[member.name];b=raw(member.name);assert len(b)==e['bytes'] and hashlib.sha256(b).hexdigest()==e['sha256']
 result=obj('result.json');assert result['byteIdenticalReports']==150 and result['nativeCommands']==601 and result['comparedVectors']==750
 assert obj('completion.json')['exitCode']==0
 freeze=obj('native/prediction-freeze.json');assert freeze['allNativeReplaysPassed'] and freeze['firstRepeatedPredictionExact'] and not freeze['scoringStarted']
 assert len(result['reports'])==150 and len(result['primary'])==5
 assert all(x['beatsAllThreeComparators'] for x in result['primary'])
 assert sum(len(x['kernelWorseThanMean']) for x in result['primary'])==37
 assert sum(len(x['kernelWorseThanNoChange']) for x in result['primary'])==34
print('PASS: complete retained nested-selection and native execution evidence; unchanged predictor')
