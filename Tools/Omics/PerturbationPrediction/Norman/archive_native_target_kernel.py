#!/usr/bin/env python3
"""Archive native target-kernel evidence with explicit exclusions and hash verification."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--work',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--verify',action='store_true')
a=p.parse_args();archiver=Path(__file__).parents[2]/'Reduction/archive_adaptive_integration.py'
if not a.verify:
    if a.work is None:p.error('--work required when packing')
    assert json.loads((a.work/'native/receipt.json').read_text())['status']=='passed'
    assert json.loads((a.work/'checks/checks.json').read_text())['status']=='passed'
    assert not (a.work/'native/scratch').exists()
    with tempfile.TemporaryDirectory(prefix='numivivo-native-kernel-evidence-') as directory:
        stage=Path(directory);root=stage/'run';root.mkdir()
        files=[]
        for label in ['native','integrity','checks','checks-repeat']:
            files += [path for path in (a.work/label).rglob('*') if path.is_file() and path != a.work/'native/training.json']
        files += [path for path in (a.work/'portable-replay').iterdir() if path.is_file()]
        files += [path for path in a.work.iterdir() if path.is_file() and path.suffix in ('.json','.log') and not path.name.startswith(('packaging','publication'))]
        files += [a.work/'check_clipping.py']
        files += [a.work/'source-bundle'/name for name in ['plan.json','receipt.json']]
        files += [a.work/'scoped/sources.sha256']
        for source in files:
            assert not source.is_symlink()
            destination=root/source.relative_to(a.work);destination.parent.mkdir(parents=True,exist_ok=True)
            os.link(source,destination)
        subprocess.run([sys.executable,str(archiver),'--out',str(a.out),'--external',str(a.work/'external-inputs.json'),'--source','run='+str(root)],check=True)
    manifest=json.loads((a.out/'manifest.json').read_text());readable=[]
    for label,source in [('native-checks.json',a.work/'checks/checks.json'),('summary.json',a.work/'checks/summary.json'),
                         ('coverage.json',a.work/'checks/coverage.json'),('integrity-checks.json',a.work/'integrity/checks.json'),
                         ('replay-checks.json',a.work/'replay-checks.json'),('platform-replay.json',a.work/'portable-replay/checks.json')]:
        raw=source.read_bytes();(a.out/label).write_bytes(raw);readable.append(dict(path=label,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
    manifest['humanReadableFiles']=readable
    manifest['retention']='Full source H5AD/report and native executable remain external. native/training.json.gz retains exact prepared training bytes. All folds retain prediction report JSON as gzip plus input plans, queries, receipts and hashes. Only first-model retains complete fitted model JSON. See NATIVE_TARGET_KERNEL.md for reconstruction of other models.'
    (a.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
subprocess.run([sys.executable,str(archiver),'--verify','--out',str(a.out)],check=True)
manifest=json.loads((a.out/'manifest.json').read_text())
assert max(entry['storedBytes'] for entry in manifest['files'])<100*1024*1024
