#!/usr/bin/env python3
"""Native full Baron donor-held-out projection and legacy-reference regression; stdlib only."""
import argparse,hashlib,json,re,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);commands=[]
def run(label,args):
    command=['/usr/bin/time','-l',str(a.binary),*map(str,args)]
    r=subprocess.run(command,capture_output=True,text=True);(a.out/(label+'.log')).write_text(r.stdout+r.stderr)
    rss=re.search(r'(\d+)\s+maximum resident set size',r.stderr);elapsed=re.search(r'([\d.]+) real',r.stderr)
    commands.append(dict(label=label,command=command,exitCode=r.returncode,seconds=float(elapsed[1]) if elapsed else None,maximumResidentBytes=int(rss[1]) if rss else None))
    (a.out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n');assert r.returncode==0,r.stderr
for donor in ['human1','human2','human3','human4']:
    inputs=a.inputs/donor;root=a.out/donor;root.mkdir()
    run(donor+'-fit',['singlecell-h5ad-pca',inputs/'train.h5ad','--plan',inputs/'fit.json','--output',root/'reference'])
    run(donor+'-query',['singlecell-h5ad-pca-query',inputs/'query.h5ad','--plan',inputs/'query.json','--reference',root/'reference','--output',root/'query'])
    run(donor+'-verify',['singlecell-h5ad-pca-query-verify',root/'query'])
    run(donor+'-legacy-fit',['singlecell-reference-fit',inputs/'train.h5ad','--plan',inputs/'legacy-fit.json','--output',root/'legacy-reference'])
    run(donor+'-legacy-map',['singlecell-reference-map',inputs/'query.h5ad','--plan',inputs/'query.json','--reference',root/'legacy-reference','--output',root/'legacy-query'])
    assert json.loads((root/'query/report.json').read_text())['overlappingDonorIDs']==[]
    print(donor+' native lifecycle passed',flush=True)
(a.out/'checks.json').write_text(json.dumps(dict(status='passed',folds=4,commands=len(commands)),indent=2)+'\n')
