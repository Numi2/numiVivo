#!/usr/bin/env python3
"""Compare the native local smoother with pinned R stats LOWESS."""
import argparse,json,math,os,subprocess,hashlib
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--binary',required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
cases=[]
for n in [2,7,100,1000]:
 x=[i/(n-1)*5 for i in range(n)]
 for pattern in ['constant','linear','curve','outlier','ties']:
  xx=[float(i//3) for i in range(n)] if pattern=='ties' else x
  y=[3.0 if pattern=='constant' else 2+4*v if pattern=='linear' else math.sin(v)+.1*math.cos(37*v) for v in xx]
  if pattern=='outlier':y[n//2]+=20
  cases.append(dict(id=f'{pattern}-{n}',x=xx,y=y))
for d in [0,.1,1]:cases.append(dict(id=f'delta-{d}',x=[i/10 for i in range(100)],y=[math.cos(i/10) for i in range(100)],deltaFraction=d))
cases.append(dict(id='constant-x',x=[0.0]*15,y=[float(i%4) for i in range(15)]))
cases.append(dict(id='reversed',x=list(reversed(cases[12]['x'])),y=list(reversed(cases[12]['y']))))
raw=json.dumps(cases,separators=(',',':')).encode();(a.out/'input.json').write_bytes(raw)
r=subprocess.run(['ssh','macmini',a.binary],input=raw,stdout=subprocess.PIPE,stderr=subprocess.PIPE);(a.out/'native.log').write_bytes(r.stderr);assert r.returncode==0;(a.out/'native.json').write_bytes(r.stdout)
script=Path(__file__).with_name('lowess_reference.R')
subprocess.run(['/opt/homebrew/bin/Rscript',str(script),str(a.out/'input.json'),str(a.out/'reference.json')],check=True,env={**os.environ,'R_LIBS_USER':'/Users/home/numivivo-r-library-20260909'})
native=json.loads(r.stdout);reference=json.loads((a.out/'reference.json').read_text());checks=[]
for i,n,ref in zip(cases,native,reference['results'],strict=True):
 errors=[abs(v-w)/max(1,abs(w)) for v,w in zip(n['fit']['fitted'],ref['fitted'],strict=True)] if 'fit' in n else [float('inf')]
 checks.append(dict(id=i['id'],maximumRelativeError=max(errors),status='passed' if max(errors)<=2e-7 else 'failed',native=n,reference=ref))
summary=dict(status='passed' if all(c['status']=='passed' for c in checks) else 'failed',checks=checks,maximumRelativeError=max(c['maximumRelativeError'] for c in checks),binarySHA256=subprocess.check_output(['ssh','macmini','shasum','-a','256',a.binary]).decode().split()[0],protocolSHA256=hashlib.sha256(Path(__file__).with_name('PROTOCOL.md').read_bytes()).hexdigest())
(a.out/'checks.json').write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='checks'},indent=2));print([(c['id'],c['maximumRelativeError']) for c in checks if c['status']!='passed'])
