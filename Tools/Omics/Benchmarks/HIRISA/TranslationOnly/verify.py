"""Verify retained full-cohort evidence, not a fresh execution."""
from pathlib import Path
import hashlib,json,tarfile
root=Path(__file__).parent
m=json.loads((root/'manifest.json').read_text())
assert hashlib.sha256((root/'evidence.tar.gz').read_bytes()).hexdigest()==m['archiveSHA256']
with tarfile.open(root/'evidence.tar.gz') as t:
    names=[x.name for x in t.getmembers() if x.isfile()]
    assert len(names)==len(set(names)) and set(names)==set(m['members'])
    data={n:t.extractfile(n).read() for n in names}
for name,payload in data.items():assert hashlib.sha256(payload).hexdigest()==m['members'][name],name
def read(n):return json.loads(data[n])
s=read('summary.json');assert s==json.loads((root/'summary.json').read_text())
assert read('terminal.json')['returnCode']==0 and read('complete.json')['allNumericalChecksPassed']
a=read('evaluation/translation.json');e=[x for x in a['comparisons'] if x['sufficientSupport'] and x['controlSensitive']]
assert a['cells']==s['cells']==1612594 and a['foldCount']==s['folds']==114
assert len(a['comparisons'])==s['comparisons']==471 and len(e)==s['eligible']==146
assert sum(x['marginFailure'] for x in e)==s['failures']
assert sum(x['rare'] and x['marginFailure'] for x in e)==s['rareFailures']
for name in ['response','mixed-programs','within-programs']:
    c=read('evaluation/'+name+'.json')['comparisons'];v=[x for x in c if x.get('controlSensitive',True)]
    assert len(c)==s[name]['comparisons'] and len(v)==s[name]['sensitive']
    assert sum(all(x['gates'].values()) for x in v)==s[name]['passedAmongSensitive']
assert read('output-freeze.json')['files']['translation.bin']==m['externalCoordinate']['SHA256']
print('PASS archive, complete-cohort accounting and output binding; annotation failures:',s['failures'])
