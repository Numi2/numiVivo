from pathlib import Path
import json,hashlib,tarfile
p=Path(__file__).parent;m=json.loads((p/'manifest.json').read_text());a=p/'evidence.tar.gz'
assert hashlib.sha256(a.read_bytes()).hexdigest()==m['archiveSHA256']
with tarfile.open(a,'r:gz') as t:
 files={x.name:t.extractfile(x).read() for x in t.getmembers() if x.isfile()}
 assert set(files)==set(m['members'])
 for n,b in files.items():assert hashlib.sha256(b).hexdigest()==m['members'][n],n
 def read(n):return json.loads(files[n])
 s=read('summary.json');assert s==json.loads((p/'summary.json').read_text())
 c=read('evaluation/rigid.json')['comparisons'];e=[x for x in c if x['sufficientSupport'] and x['controlSensitive']]
 assert len(c)==s['comparisons']==471 and len(e)==s['eligible']==146
 assert sum(x['marginFailure'] for x in e)==s['failures']==31
 assert sum(x['rare'] and x['marginFailure'] for x in e)==s['rareFailures']==9
 assert read('complete.json')['allNumericalChecksPassed']
 for name in ['response','mixed-programs','within-programs']:
  q=read('evaluation/'+name+'.json')['comparisons'];v=[x for x in q if x.get('controlSensitive',True)]
  assert sum(all(x['gates'].values()) for x in v)==s[name]['passedAmongSensitive']
 assert read('output-freeze.json')['files']['rigid.bin']==m['externalCoordinate']['SHA256']
print('PASS complete archive hashes and full-cohort result accounting; preservation FAIL retained')
