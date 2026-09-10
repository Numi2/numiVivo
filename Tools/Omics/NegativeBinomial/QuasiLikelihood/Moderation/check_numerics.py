#!/usr/bin/env python3
"""Check QL numerical owners against pinned R and high-precision shape functions."""
import argparse,hashlib,json,math,os,subprocess
from pathlib import Path
import mpmath as mp
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--binary',required=True);args=parser.parse_args()
assert not args.out.exists();args.out.mkdir(parents=True)
sha=lambda b:hashlib.sha256(b).hexdigest();source=Path(__file__).parent
smoothers=[];labels=[]
for n in [3,7,15,100,1000]:
    for pattern in ['equal','graded','clipped','ties','dominant']:
        x=[float(k//3) if pattern=='ties' else 10*k/(n-1) for k in range(n)]
        y=[math.sin(v)+.2*math.cos(7*v) for v in x]
        w=[1.0 if pattern=='equal' else (0 if k%5==0 else 200 if k%7==0 else .1+k/n) if pattern=='clipped' else (50 if k==n//2 else .1) if pattern=='dominant' else .1+3*(k%7)/7 for k in range(n)]
        smoothers.append(dict(x=x,y=y,weights=w,span=.6));labels.append(f'{pattern}-{n}')
shapes=[1e-8,1e-6,.0001,.005,.01,.1,.5,1,2,5,15.99,16,100,5000,1e6]
tails=[dict(logStatistic=v,numeratorDF=a,denominatorDF=b) for a in [.01,1,5,20,1000] for b in [.01,1,10,100,9998] for v in [-700,-100,-20,-5,0,5,20,100,700]]
payload=json.dumps(dict(smoothers=smoothers,shapes=shapes,tails=tails),separators=(',',':')).encode();(args.out/'input.json').write_bytes(payload)
run=subprocess.run(['ssh','-4','macmini',args.binary],input=payload,capture_output=True);(args.out/'native.log').write_bytes(run.stderr);(args.out/'native.json').write_bytes(run.stdout);assert run.returncode==0
with (args.out/'reference.log').open('wb') as log:
    subprocess.run(['/opt/homebrew/bin/Rscript',str(source/'numerics_reference.R'),str(args.out/'input.json'),str(args.out/'reference.json')],check=True,stdout=log,stderr=log,env={**os.environ,'R_LIBS_USER':'/Users/home/numivivo-r-library-20260909'})
n=json.loads(run.stdout);r=json.loads((args.out/'reference.json').read_text());checks=[]
for label,row,native,ref in zip(labels,smoothers,n['smoothers'],r['smoothers'],strict=True):
    error=max(abs(a-b)/max(1,abs(b)) for a,b in zip(native['fit']['fitted'],ref,strict=True)) if 'fit' in native else None
    checks.append(dict(kind='smoother',id=label,maximumRelativeError=error,status='passed' if error is not None and error<=2e-7 else 'failed',native=native,reference=ref))
mp.mp.dps=80
for x,native in zip(shapes,n['shapes'],strict=True):
    xx=mp.mpf(x);h=mp.log(xx)-mp.digamma(xx);t=mp.polygamma(1,xx)
    error=max(abs(native['logMinusDigamma']/float(h)-1),abs(native['trigamma']/float(t)-1))
    checks.append(dict(kind='shape',x=x,maximumRelativeError=error,status='passed' if error<=2e-10 else 'failed',native=native,reference=dict(logMinusDigamma=mp.nstr(h,75),trigamma=mp.nstr(t,75))))
for row,native,ref in zip(tails,n['tails'],r['tails'],strict=True):
    error=max(abs(native['fit'][k]-ref[k]) for k in ['lower','upper']) if 'fit' in native and all(isinstance(ref[k],(int,float)) for k in ['lower','upper']) else None
    checks.append(dict(kind='tail',input=row,maximumAbsoluteLogError=error,status='passed' if error is not None and error<=2e-7 else 'failed',native=native,reference=ref))
summary=dict(status='passed' if all(c['status']=='passed' for c in checks) else 'failed',checks=checks,binarySHA256=subprocess.check_output(['ssh','-4','macmini','shasum','-a','256',args.binary]).decode().split()[0],protocolSHA256=sha((source/'PROTOCOL.md').read_bytes()),inputSHA256=sha(payload),nativeSHA256=sha(run.stdout))
(args.out/'checks.json').write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n')
for kind in ['smoother','shape','tail']:
    cases=[c for c in checks if c['kind']==kind];print(kind,'total',len(cases),'failed',sum(c['status']!='passed' for c in cases))
for c in checks:
    if c['status']!='passed':print({k:v for k,v in c.items() if k not in ['native','reference']})
raise SystemExit(0 if summary['status']=='passed' else 1)
