#!/usr/bin/env python3
"""Replay hash-bound legacy fixtures and complete cohorts on a native host."""
import argparse,hashlib,json,os,shutil,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
reference=json.loads((a.reference/'manifest.json').read_text());a.out.mkdir(parents=True,exist_ok=False);records=[]
for case in reference['cases']:
 source=a.reference/case['source'];plan=a.reference/case['plan'];name=case['name'];out=a.out/name
 assert sha(source)==case['files']['original.h5ad'] and sha(plan)==case['files']['plan.json']
 available=shutil.disk_usage(a.out).free;required=2*(source.stat().st_size+case['outputBytes'])+64*2**20
 assert available>=required,(name,'insufficient disk before publication and reconstruction',available,required)
 with (a.out/(name+'.log')).open('w') as log:
  started=time.monotonic();subprocess.run([str(a.binary),'project',str(source),str(plan),str(out)],check=True,stdout=log,stderr=subprocess.STDOUT)
  published=time.monotonic()-started
  subprocess.run([str(a.binary),'verify-project',str(out)],check=True,stdout=log,stderr=subprocess.STDOUT)
 checks={k:sha(out/k) for k in case['files']};assert checks==case['files'],name
 records.append(dict(name=name,files=checks,allFourPayloadsExact=True,publicationSeconds=published,initialAvailableBytes=available,admittedAdditionalBytes=required))
 print(name,'exact cross-host replay PASS',flush=True)
(a.out/'summary.json').write_text(json.dumps(dict(status='passed',binarySHA256=sha(a.binary),checkerSHA256=sha(Path(__file__)),referenceSHA256=sha(a.reference/'manifest.json'),cases=records,finalAvailableBytes=shutil.disk_usage(a.out).free),indent=2)+'\n')
