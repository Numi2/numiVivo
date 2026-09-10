#!/usr/bin/env python3
"""Exercise descriptor evidence corruption, including rehashed false annotations."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('capture','identities','out'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    pinned=sha(a.capture/'receipt.json');env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    results=[]
    for case in ['changed-annotation','rehashed-annotation','missing-capture','unlisted-capture','wrong-receipt','wrong-identities']:
        d=a.out/case;shutil.copytree(a.capture,d)
        receipt=json.loads((d/'receipt.json').read_text());expected=pinned;identities=a.identities
        if case in ['changed-annotation','rehashed-annotation']:
            path=d/'annotations.json';value=json.loads(path.read_text());value['MANF']['terms'].append('GO:9999999')
            path.write_text(json.dumps(value))
            if case=='rehashed-annotation':
                entry=next(x for x in receipt['files'] if x['path']=='annotations.json')
                entry.update(sha256=sha(path),bytes=path.stat().st_size)
                (d/'receipt.json').write_text(json.dumps(receipt));expected=sha(d/'receipt.json')
        elif case=='missing-capture': (d/'responses/ENSG00000145050.json').unlink()
        elif case=='unlisted-capture': (d/'unlisted.json').write_text('{}')
        elif case=='wrong-receipt': expected='0'*64
        else:
            identities=a.out/'changed-identities.json';value=json.loads(a.identities.read_text())
            value['features'][0]='ENSG00000000000';identities.write_text(json.dumps(value))
        destination=a.out/(case+'-unexpected.json')
        command=[sys.executable,str(Path(__file__).with_name('check_descriptors.py')),'--capture',str(d),
                 '--identities',str(identities),'--receipt-sha256',expected,'--out',str(destination)]
        result=subprocess.run(command,env=env,capture_output=True,text=True)
        (a.out/(case+'.log')).write_text(result.stdout+result.stderr)
        assert result.returncode!=0 and not destination.exists(),case
        results.append(dict(case=case,returnCode=result.returncode,unexpectedOutputAbsent=True))
        shutil.rmtree(d) # Generated mutation copy only; the captured source remains immutable.
    report=dict(status='passed',cases=results)
    (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
