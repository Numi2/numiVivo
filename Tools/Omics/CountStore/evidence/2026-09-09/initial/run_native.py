#!/usr/bin/env python3
"""Run the complete Norman count-store qualification on the native build host."""
import argparse,hashlib,json,os,re,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',type=Path,required=True)
a=p.parse_args();root=a.root
binary=root/'numivivo';source=root/'original.h5ad';plan=root/'plan.json'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
assert sha(source)=='efde6f5301fe256725dce1d980f37bd96a13481a9a16135515897368e631affc'
commands=[]
def run(label,args):
 start=time.monotonic()
 r=subprocess.run(['/usr/bin/time','-l',str(binary),*map(str,args)],capture_output=True,text=True)
 (root/(label+'.log')).write_text(r.stdout+r.stderr)
 peak=re.search(r'(\d+)\s+maximum resident set size',r.stderr)
 commands.append(dict(label=label,arguments=list(map(str,args)),exitCode=r.returncode,seconds=time.monotonic()-start,maximumResidentBytes=int(peak[1]) if peak else None))
 (root/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
 assert r.returncode==0,(label,r.stderr[-3000:])
run('import',['singlecell-h5ad-store',source,'--plan',plan,'--output',root/'store'])
# Both operations reconstruct raw counts from the retained H5AD before mapping.
run('normalize',['singlecell-count-store-normalize',root/'store','--target','10000','--output',root/'normalized'])
run('normalize-verify',['singlecell-count-store-normalize-verify',root/'normalized','--store',root/'store'])
receipt=json.loads((root/'store/receipt.json').read_text())
assert receipt['entries']==361582621
(root/'native-checks.json').write_text(json.dumps(dict(status='passed-native-reconstruction',binarySHA256=sha(binary),sourceSHA256=sha(source),planSHA256=sha(plan),
 entries=receipt['entries'],countBytes=(root/'store/counts.bin').stat().st_size,normalizedBytes=(root/'normalized/values.bin').stat().st_size,
 windowBytes=16777216,writeBufferBytes=1048576,maximumResidentBytes=max(v['maximumResidentBytes'] or 0 for v in commands),commands=commands,
 checkerSHA256=sha(Path(__file__))),indent=2)+'\n')
