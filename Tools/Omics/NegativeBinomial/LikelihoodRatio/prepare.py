#!/usr/bin/env python3
"""Bind existing measured model archives before LRT evaluation; no count rewrite."""
import argparse,hashlib,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--out',type=Path,required=True);p.add_argument('--host',default='macmini')
p.add_argument('--null-root',default='/Users/n/numivivo-null-benchmark-20260910')
p.add_argument('--treatment-root',default='/Users/n/numivivo-empirical-prior-20260910')
a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
def read(path):return subprocess.check_output(['ssh',a.host,'cat',str(path)])
cases=[]
for study in ['kang','hagai']:
 for seed in range(1,11):
  for policy in ['default','active']:
   root=Path(a.null_root)/study/str(seed)/policy
   archive=json.loads(read(root/'archive.json'));e=next(e for e in archive['entries'] if e['logicalPath']=='report.json')
   cases.append(dict(id=f'{study}-null-{seed:02}-{policy}',family='null',study=study,source=str(root/'report.json.gz'),sourceSHA256=e['sha256'],logicalSHA256=e['logicalSHA256']))
for rec in json.loads(read(Path(a.treatment_root)/'protocol.json'))['records']:
 root=Path(a.treatment_root)/rec['id']/'native';archive=json.loads(read(root/'archive.json'));e=next(e for e in archive['entries'] if e['logicalPath']=='report.json')
 cases.append(dict(id=rec['id'],family='treatment',study=rec['id'].split('-')[0],source=str(root/'report.json.gz'),sourceSHA256=e['sha256'],logicalSHA256=e['logicalSHA256']))
protocol=dict(host=a.host,baselineCommit='e5c5d4ee15f51c5354ade241e491c7e621e9fabf',protocolSHA256=hashlib.sha256(Path(__file__).with_name('PROTOCOL.md').read_bytes()).hexdigest(),cases=cases)
(a.out/'protocol.json').write_text(json.dumps(protocol,sort_keys=True,indent=2)+'\n')
print(json.dumps(dict(cases=len(cases),protocolSHA256=protocol['protocolSHA256'])))
