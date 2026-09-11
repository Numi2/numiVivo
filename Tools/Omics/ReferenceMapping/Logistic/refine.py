#!/usr/bin/env python3
"""Freeze stricter numerical convergence while reusing every original source byte."""
import argparse,copy,json,os
from pathlib import Path
from prepare import sha,write
p=argparse.ArgumentParser();p.add_argument('--original',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();old=json.loads((a.original/'input-freeze.json').read_text())
for n,h in old['files'].items():assert sha(a.original/n)==h,n
a.out.mkdir(parents=True,exist_ok=False);sources=copy.deepcopy(old['sources']);files={}
for donor in ['human1','human2','human3','human4']:
 root=a.out/donor;root.mkdir();plan=json.loads((a.original/donor/'fit.json').read_text());assert plan['logistic']['gradientTolerance']==1e-7;plan['logistic']['gradientTolerance']=1e-8;write(root/'fit.json',plan)
 (root/'query.json').write_bytes((a.original/donor/'query.json').read_bytes())
 for group in ['train','query']:
  spec=sources[donor+'/'+group];spec['gzipPath']=os.path.relpath(a.original/spec['gzipPath'],a.out);files[spec['gzipPath']]=spec['gzipSHA256']
for f in a.out.rglob('*.json'):files[str(f.relative_to(a.out))]=sha(f)
write(a.out/'input-freeze.json',dict(schemaVersion=1,originalInputFreezeSHA256=sha(a.original/'input-freeze.json'),originalProtocolSHA256=old['protocolSHA256'],refinementSHA256=sha(Path(__file__).with_name('NUMERICAL_REFINEMENT.md')),refinerSHA256=sha(__file__),sources=sources,files=files,allFourDonors=True,fitStarted=False,queryLabelsExcludedByMapping=True,numericalChange='gradientTolerance 1e-7 -> 1e-8 only'))
