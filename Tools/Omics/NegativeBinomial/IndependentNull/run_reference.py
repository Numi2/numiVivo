#!/usr/bin/env python3
"""Run every frozen Bioconductor cohort, retaining warnings and failures."""
import argparse, hashlib, json, os, subprocess, time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
p.add_argument('--rscript',type=Path,required=True);p.add_argument('--r-library',type=Path,required=True)
a=p.parse_args();root=a.root;sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
script=Path(__file__).with_name('reference.R');out=root/'reference';out.mkdir(exist_ok=True)
lock=out/'driver.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,str(os.getpid()).encode());os.close(fd)
try:
    for number in range(1,10):
        case=f'{number:02d}';inp=root/'reference-inputs'/case;directory=out/case
        meta=json.loads((inp/'input.json').read_text())
        for name,digest in meta['files'].items():assert sha(inp/name)==digest
        identity=dict(inputSHA256=sha(inp/'input.json'),scriptSHA256=sha(script))
        receipt=out/(case+'-receipt.json')
        if receipt.exists():
            saved=json.loads(receipt.read_text());assert all(saved[k]==v for k,v in identity.items())
            for name,digest in saved['files'].items():assert sha(directory/name)==digest
            print(case,'retained checkpoint',flush=True);continue
        assert not directory.exists()
        command=[str(a.rscript),str(script),str(inp),str(directory)];start=time.time()
        with (out/(case+'.log')).open('wb') as log:
            r=subprocess.run(command,stdout=log,stderr=log,env={**os.environ,'R_LIBS_USER':str(a.r_library),
                'VECLIB_MAXIMUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
        record=dict(**identity,command=command,returnCode=r.returncode,seconds=time.time()-start,
            files={p.name:sha(p) for p in sorted(directory.iterdir()) if p.is_file()} if directory.exists() else {})
        receipt.write_text(json.dumps(record,sort_keys=True,indent=2)+'\n')
        print(case,'returnCode',r.returncode,'seconds',round(record['seconds'],2),flush=True)
finally:lock.unlink()
