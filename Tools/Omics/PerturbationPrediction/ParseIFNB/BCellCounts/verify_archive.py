from pathlib import Path
import json,hashlib,tarfile
p=Path(__file__).parent;m=json.loads((p/'manifest.json').read_text());a=p/'preparation.tar.gz'
assert hashlib.sha256(a.read_bytes()).hexdigest()==m['archiveSHA256']
with tarfile.open(a,'r:gz') as t:
 files={x.name:t.extractfile(x).read() for x in t.getmembers() if x.isfile()}
 assert set(files)==set(m['members'])
 for n,b in files.items():assert hashlib.sha256(b).hexdigest()==m['members'][n],n
 f=json.loads(files['donor-plans/manifest.json'])
 assert sum(x['cells'] for x in f['parts'])==72446 and sum(x['records'] for x in f['parts'])==124909573
 for part in f['parts']:
  root='donor-plans/'+part['donor']+'/'
  for name,key in [('plan.json','planSHA256'),('runs.json','runsSHA256'),('parent-selected-rows.npy','rowMapSHA256')]:assert hashlib.sha256(files[root+name]).hexdigest()==part[key]
 assert json.loads(files['preparation-verification.json'])['allSelectedRowsExactlyOnce']
 assert json.loads(files['tests.json'])['status']=='pass'
print('PASS frozen preparation archive; full native execution/replay completion remains separate')
