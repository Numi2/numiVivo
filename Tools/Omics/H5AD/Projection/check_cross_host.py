#!/usr/bin/env python3
"""Replay selected completed projection bundles on another native host."""
import argparse,hashlib,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);rows=[]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for src in sorted(a.reference.iterdir()):
 if not src.is_dir():continue
 out=a.out/src.name
 with (a.out/(src.name+'.log')).open('w') as log:
  subprocess.run([str(a.binary),'project',str(src/'original.h5ad'),str(src/'plan.json'),str(out)],check=True,stdout=log,stderr=subprocess.STDOUT)
  subprocess.run([str(a.binary),'verify-project',str(out)],check=True,stdout=log,stderr=subprocess.STDOUT)
 checks={n:sha(out/n) for n in ('original.h5ad','projected.h5ad','plan.json','report.json')}
 for n,h in checks.items():assert h==sha(src/n),(src.name,n)
 rows.append(dict(case=src.name,allFourPayloadsExact=True,files=checks))
assert len(rows)>0
record=dict(status='passed',cases=rows,binarySHA256=sha(a.binary),checkerSHA256=sha(Path(__file__)),scope='Native cross-host projection/reconstruction of complete stored format bundles, including repeated axes and transfer boundaries. No biological claim.')
(a.out/'summary.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))
