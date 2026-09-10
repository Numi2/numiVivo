#!/usr/bin/env python3
"""Check archive restoration rejects damaged sources/payloads without publication."""
import argparse,json,shutil,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archive',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);checks=[]
for case in ['source','payload','schema']:
 root=a.out/case;archive=root/'study/1/default';archive.parent.mkdir(parents=True)
 shutil.copytree(a.archive,archive)
 shared=archive.parent.parent/'annotated.h5ad';shutil.copyfile(a.archive/'../../annotated.h5ad',shared)
 if case=='source':
  with shared.open('ab') as f:f.write(b'changed')
 elif case=='payload':
  with (archive/'report.json.gz').open('ab') as f:f.write(b'changed')
 else:
  m=json.loads((archive/'archive.json').read_text());m['entries'][0]['logicalPath']='outside.json';(archive/'archive.json').write_text(json.dumps(m))
 dest=root/'output';r=subprocess.run([sys.executable,str(Path(__file__).with_name('restore_bundle.py')),'--archive',str(archive),'--out',str(dest)],capture_output=True,text=True)
 (root/'rejection.log').write_text(r.stdout+r.stderr)
 assert r.returncode==1 and not dest.exists() and not list(root.glob('.numivivo-restore-*')),(case,r.returncode,r.stderr)
 checks.append(case)
(a.out/'checks.json').write_text(json.dumps(dict(status='passed',rejectedBeforePublication=checks),indent=2)+'\n')
print(json.dumps(dict(status='passed',checks=len(checks))))
