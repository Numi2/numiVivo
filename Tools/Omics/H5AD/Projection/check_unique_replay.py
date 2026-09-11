#!/usr/bin/env python3
"""Republish every prior unique-selection bundle and require exact original bytes."""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--previous',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();previous=json.loads((a.previous/'summary.json').read_text());cases=[]
for record in previous['cases']:
 name=record['name'];old=a.previous/name;output=a.out/name
 assert sha(old/'original.h5ad')==record['sourceSHA256'] and sha(old/'projected.h5ad')==record['outputSHA256']
 with (a.out/(name+'.log')).open('w') as log:
  subprocess.run([str(a.binary),'project',str(old/'original.h5ad'),str(old/'plan.json'),str(output)],stdout=log,stderr=subprocess.STDOUT,check=True)
  subprocess.run([str(a.binary),'verify-project',str(output)],stdout=log,stderr=subprocess.STDOUT,check=True)
 files={}
 for f in ('original.h5ad','projected.h5ad','plan.json','report.json'):
  assert (old/f).read_bytes()==(output/f).read_bytes(),(name,f);files[f]=sha(output/f)
 cases.append(dict(name=name,allFourPayloadsExact=True,files=files));print(name,'original bytes exact',flush=True)
(a.out/'summary.json').write_text(json.dumps(dict(status='passed',cases=cases,binarySHA256=sha(a.binary),previousSummarySHA256=sha(a.previous/'summary.json'),checkerSHA256=sha(Path(__file__)),scope='Every prior declared unique-selection fixture and real output; actual receipt implementation stays separately bound to the current scoped binary.'),indent=2)+'\n')
