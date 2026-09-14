#!/usr/bin/env python3
"""Execute every fixed site; freeze outputs before qPCR scoring."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import time

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',type=Path,required=True)
p.add_argument('--binary',type=Path,required=True)
a=p.parse_args();r=a.root;out=r/'native';out.mkdir(exist_ok=False)
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
records=[];binary=h(a.binary);protocol=h(r/'PROTOCOL.md')
for entry in json.loads((r/'native-inputs/manifest.json').read_text()):
    site=entry['site'];source=r/'native-inputs'/(site+'.json.gz');request=r/'native-inputs'/(site+'-request.json')
    assert h(source)==entry['inputSHA256'] and h(request)==entry['requestSHA256']
    dest=out/(site+'.json');start=time.time()
    with dest.open('xb') as f,(out/(site+'.log')).open('xb') as log:
        result=subprocess.run([str(a.binary),str(source),str(request)],stdout=f,stderr=log)
    with dest.open('rb') as f,gzip.open(out/(site+'.json.gz'),'wb') as g:
        while block:=f.read(1024*1024):g.write(block)
    record=dict(site=site,inputSHA256=h(source),requestSHA256=h(request),outputSHA256=h(dest),gzipSHA256=h(out/(site+'.json.gz')),
                returnCode=result.returncode,seconds=time.time()-start)
    records.append(record);(out/(site+'-receipt.json')).write_text(json.dumps(record,indent=2))
    print(site,result.returncode,record['seconds'],flush=True)
(r/'native-freeze.json').write_text(json.dumps(dict(protocolSHA256=protocol,binarySHA256=binary,scored=False,records=records),indent=2))
