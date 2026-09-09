#!/usr/bin/env python3
"""Native standalone PCA lifecycle and tamper controls on numerical fixtures."""
import argparse,copy,hashlib,json,shutil,subprocess
from pathlib import Path
from pca_bundle_reference import decode
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--fixtures',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
old=json.loads((a.fixtures/'plan.json').read_text());plan={k:old[k] for k in ['schemaVersion','mapping','reduction']}
def write(path,value):path.write_text(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False))
planfile=a.out/'plan.json';write(planfile,plan);commands=[]
def run(args,reject=False,expected=None):
    r=subprocess.run([str(a.binary),*map(str,args)],capture_output=True,text=True)
    (a.out/(str(len(commands))+'.log')).write_text(r.stdout+r.stderr)
    commands.append(dict(arguments=list(map(str,args)),exitCode=r.returncode,expectedRejection=reject))
    write(a.out/'commands.json',commands)
    assert (r.returncode!=0) if reject else (r.returncode==0),r.stderr
    if expected:assert expected in r.stderr,r.stderr
results=[]
for encoding in ['csr','csc']:
    source=a.fixtures/(encoding+'.h5ad');out=a.out/encoding
    run(['singlecell-h5ad-pca',source,'--plan',planfile,'--output',out]);run(['singlecell-h5ad-pca-verify',out])
    run(['singlecell-h5ad-pca',source,'--plan',planfile,'--output',out],True)
    decoded,metadata,quality=decode(out);assert quality[2]['totalCounts']==0
    assert 'pseudobulk' not in decoded and not (out/'report.json').exists()
    results.append(decoded['reduction'])
assert results[0]==results[1]
repeat=a.out/'repeat';run(['singlecell-h5ad-pca',a.fixtures/'csr.h5ad','--plan',planfile,'--output',repeat])
assert (repeat/'receipt.json').read_bytes()==(a.out/'csr/receipt.json').read_bytes()
for name,edit in [('centers',{'pca':{'retainProjectionCenters':False}}),('cache',{'maximumCacheBytes':16}),('work',{'maximumEntryVisits':1})]:
    invalid=copy.deepcopy(plan);invalid['reduction'].update(edit);path=a.out/(name+'.json');write(path,invalid)
    run(['singlecell-h5ad-pca',a.fixtures/'csr.h5ad','--plan',path,'--output',a.out/('reject-'+name)],True)
for name in ['scores','model']:
    dst=a.out/('tamper-'+name);shutil.copytree(a.out/'csr',dst)
    receipt=json.loads((dst/'receipt.json').read_text())
    if name=='scores':
        path=dst/'scores.bin';value=bytearray(path.read_bytes());value[0]^=1;path.write_bytes(value)
    else:
        path=dst/'model.json';value=json.loads(path.read_text());value['projectionCenters'][0]+=1;write(path,value)
    receipt[name]={'bytes':list(hashlib.sha256(path.read_bytes()).digest())};write(dst/'receipt.json',receipt)
    run(['singlecell-h5ad-pca-verify',dst],True,'PCA source reconstruction differs')
assert not list(a.out.glob('.numivivo-pca-*'))
write(a.out/'checks.json',dict(status='passed',commands=len(commands),expectedRejections=sum(c['expectedRejection'] for c in commands),exactCSRCSCPCA=True,exactRepeat=True,rehashedTamperingRejected=True,binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest()))
