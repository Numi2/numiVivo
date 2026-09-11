#!/usr/bin/env python3
"""Run full legacy-derived count/pseudobulk owners on a second native host."""
import argparse,hashlib,json,os,shutil,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);a=p.parse_args();root=a.root
reference=json.loads((root/'native-reference.json').read_text())
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
assert sha(a.binary)==reference['binarySHA256'] and sha(root/'annotated.h5ad')==reference['annotatedSHA256']
records=[]
for name,publish,verify,plan in [('count-store','count-store','verify-count-store','mapping.json'),('aggregate','aggregate','verify-aggregate','aggregate-plan.json')]:
 dest=root/name;assert not dest.exists()
 free=shutil.disk_usage(root).free;required=2*(sum(x['bytes'] for x in reference['files'][name].values())+(root/'annotated.h5ad').stat().st_size)+64*2**20
 assert free>=required,(name,'insufficient native disk admission',free,required)
 for label,args in [('publish',[publish,str(root/'annotated.h5ad'),str(root/plan),str(dest)]),('verify',[verify,str(dest)])]:
  start=time.monotonic()
  with (root/(name+'-'+label+'.stdout')).open('wb') as out,(root/(name+'-'+label+'.log')).open('wb') as log:
   r=subprocess.run(['/usr/bin/time','-l',str(a.binary),*args],stdout=out,stderr=log)
  records.append(dict(stage=name,operation=label,returncode=r.returncode,seconds=time.monotonic()-start));assert r.returncode==0,(name,label,r.returncode)
 for filename,expected in reference['files'][name].items():
  q=dest/filename;assert q.stat().st_size==expected['bytes'] and sha(q)==expected['SHA256'],(name,filename)
 print(name,'all payloads exact; native reconstruction PASS',flush=True)
(root/'native-checks.json').write_text(json.dumps(dict(status='passed',binarySHA256=sha(a.binary),checkerSHA256=sha(Path(__file__)),referenceSHA256=sha(root/'native-reference.json'),commands=records,allCountStoreAndAggregatePayloadsExact=True,currentAvailableBytes=shutil.disk_usage(root).free),indent=2)+'\n')
