#!/usr/bin/env python3
"""Run all frozen requests sequentially with immutable per-case checkpoints."""
import argparse, gzip, hashlib, json, os, subprocess, time
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
p.add_argument('--binary',type=Path,required=True);a=p.parse_args();root=a.root
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
check=json.loads((root/'source-check.json').read_text());assert check['status']=='passed'
assert check['compressedReportSHA256']==sha(root/'source-report.json.gz')
assert check['protocolSHA256']==sha(Path(__file__).with_name('PROTOCOL.md'))
binary=a.binary.resolve();binary_hash=sha(binary)
manifest=json.loads((root/'requests/manifest.json').read_text());assert len(manifest)==27
out=root/'native';out.mkdir(exist_ok=True)
lock=out/'driver.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
os.write(fd,str(os.getpid()).encode());os.close(fd)
try:
    for case in manifest:
        request=root/'requests'/case['path'];assert sha(request)==case['requestSHA256']
        name=case['cohort']+'-'+case['method'];directory=out/name
        identity=dict(binarySHA256=binary_hash,sourceReportSHA256=check['compressedReportSHA256'],
            protocolSHA256=check['protocolSHA256'],requestSHA256=case['requestSHA256'])
        if (directory/'receipt.json').exists():
            saved=json.loads((directory/'receipt.json').read_text())
            assert all(saved[k]==v for k,v in identity.items())
            assert saved['outputSHA256']==sha(directory/'output.json.gz')
            print(name,'retained checkpoint',flush=True);continue
        directory.mkdir(exist_ok=False)
        command=['/usr/bin/time','-l',str(binary),str(root/'source-report.json.gz'),str(request)]
        start=time.time()
        with (directory/'output.json').open('wb') as stdout,(directory/'stderr.log').open('wb') as stderr:
            completed=subprocess.run(command,stdout=stdout,stderr=stderr)
        payload=(directory/'output.json').read_bytes()
        (directory/'output.json.gz').write_bytes(gzip.compress(payload,mtime=0))
        error=None
        try:
            result=json.loads(payload)
            assert result['sourceReportSHA256']==identity['sourceReportSHA256']
            error=result.get('error')
        except Exception as e:error='Invalid or missing harness output: '+str(e)
        record=dict(**case,**identity,command=command,returnCode=completed.returncode,seconds=time.time()-start,
            error=error,outputSHA256=sha(directory/'output.json.gz'))
        (directory/'receipt.json').write_text(json.dumps(record,sort_keys=True,indent=2)+'\n')
        (directory/'output.json').unlink() # identical content retained losslessly in gzip
        print(name,'returnCode',completed.returncode,'seconds',round(record['seconds'],2),'error',error,flush=True)
finally:
    lock.unlink()
