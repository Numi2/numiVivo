#!/usr/bin/env python3
"""Independent byte-for-byte record check across release and ASAN owners."""
import argparse, hashlib, json, struct
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args();rows=[]
for count in [0,1,65535,65536,65537,131075]:
 expected=b''.join(struct.pack('<IIQ',i,2**32-1-i,2**64-1-i*17) for i in range(count))
 for mode in ['writer-checks','writer-checks-asan']:
  observed=(a.root/mode/(str(count)+'.bin')).read_bytes()
  assert observed==expected,(mode,count)
 rows.append(dict(entries=count,bytes=len(expected),SHA256=hashlib.sha256(expected).hexdigest()))
for mode in ['writer-checks','writer-checks-asan']:
 checks=json.loads((a.root/mode/'checks.json').read_text())
 assert checks['status']=='passed' and checks['failedFlushAppendRejected']
 assert (a.root/mode/'closed.bin').stat().st_size==0
result=dict(status='passed',cases=rows,releaseAndASANExact=True,independentFormat='<IIQ',failedFlushAppendRejected=True)
(a.root/'writer-independent.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
