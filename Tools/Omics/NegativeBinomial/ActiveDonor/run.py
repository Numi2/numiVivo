#!/usr/bin/env python3
"""Run fresh default and active-donor products from an existing source bundle."""
import argparse,hashlib,json,os,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--source-bundle',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
plan=json.loads((a.source_bundle/'plan.json').read_text());assert len(plan['contrasts'])==1
assert 'zeroTotalDonorPolicy' not in plan['contrasts'][0]['negativeBinomialOptions']
summary=dict(binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),
    sourceSHA256=hashlib.sha256((a.source_bundle/'original.h5ad').read_bytes()).hexdigest(),runs=[])
def run(args,label):
    start=time.monotonic()
    with (a.out/(label+'.log')).open('w') as log:
        result=subprocess.run([str(a.binary.resolve())]+args,stdout=log,stderr=subprocess.STDOUT,
            env={**os.environ,'NUMIVIVO_HDF5_LIBRARY':'/opt/homebrew/opt/hdf5/lib/libhdf5.dylib'})
    summary['runs'].append(dict(command=args,seconds=time.monotonic()-start,exitCode=result.returncode))
    (a.out/'run.json').write_text(json.dumps(summary,indent=2)+'\n')
    assert result.returncode==0,(label,result.returncode)
for mode in ['default','active']:
    if mode=='active':plan['contrasts'][0]['negativeBinomialOptions']['zeroTotalDonorPolicy']='activeDonorProfile'
    path=a.out/(mode+'-plan.json');path.write_text(json.dumps(plan,sort_keys=True)+'\n')
    run(['singlecell-h5ad-pseudobulk',str(a.source_bundle/'original.h5ad'),'--plan',str(path),'--output',str(a.out/mode)],mode+'-publish')
    run(['singlecell-h5ad-pseudobulk-verify',str(a.out/mode)],mode+'-verify')
assert (a.out/'default/report.json').read_bytes()==(a.source_bundle/'report.json').read_bytes(),'default output drift'
summary['exactDefaultReport']=True
(a.out/'run.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
