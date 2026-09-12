from pathlib import Path
import gzip,hashlib,json,tarfile
p=Path(__file__).parent;m=json.loads((p/'manifest.json').read_text());archive=p/'evidence.tar.gz'
assert hashlib.sha256(archive.read_bytes()).hexdigest()==m['archiveSHA256']
with tarfile.open(archive,'r:gz') as t:
 files={x.name:t.extractfile(x).read() for x in t.getmembers() if x.isfile()}
 assert set(files)==set(m['members'])
 for n,b in files.items():assert hashlib.sha256(b).hexdigest()==m['members'][n],n
 report=json.loads(files['checks.json'])
 for rec in report['runs']:
  b=files['checks/'+rec['name']+'-output.json.gz'];assert hashlib.sha256(b).hexdigest()==rec['compressedSHA256'];assert hashlib.sha256(gzip.decompress(b)).hexdigest()==rec['outputSHA256']
 for method in ['ridge','mnn']:
  old=gzip.decompress(files['checks/'+method+'-old-output.json.gz']);new=gzip.decompress(files['checks/'+method+'-new-output.json.gz']);assert old==new
  a=json.loads(new);b=json.loads(gzip.decompress(files['checks/'+method+'-protected-output.json.gz']));assert b['options'].pop('protectedSampleGroups')==json.loads(files['checks/'+method+'-protected.json'])['protectedSampleGroups'];assert a==b
 assert report['status']=='pass' and len(report['rejections'])==9
 cli=json.loads(files['cli-checks.json']);assert cli['status']=='pass' and len(cli['runs'])==6
print('PASS archive hashes, retained exact native regression and report scope')
