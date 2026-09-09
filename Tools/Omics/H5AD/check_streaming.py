#!/usr/bin/env python3
"""Streaming correctness checks; synthetic fixtures are not biological evidence.

Run check_interop.py first to supply independent resident-path results.
"""
import argparse,json,shutil,subprocess
from pathlib import Path
import h5py
import numpy as np
import anndata as ad

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--interop',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
p.add_argument('--full-product',action='store_true')
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);checks=[]
def run(name,args,success=True):
    if a.full_product:
        args=(['singlecell-h5ad-pseudobulk',args[1],'--plan',args[2],'--output',args[3]] if args[0]=='aggregate' else ['singlecell-h5ad-pseudobulk-verify',args[1]])
    result=subprocess.run([str(a.binary.resolve()),*map(str,args)],capture_output=True,text=True)
    (a.out/(name+'.log')).write_text(result.stdout+result.stderr)
    assert result.returncode==(0 if success else 65),(name,result.returncode,result.stderr)
    checks.append(name)

mapping=json.loads((a.interop/'csr-plan.json').read_text())
plan=a.out/'plan.json';plan.write_text(json.dumps(dict(schemaVersion=1,mapping=mapping,contrasts=[])))
for layout in ['csr','csc','dense']:
    dest=a.out/layout
    run(layout,['aggregate',a.interop/(layout+'.h5ad'),plan,dest])
    report=json.loads((dest/'report.json').read_text())
    for field in ['quality','pseudobulk']:
        assert report[field]==json.loads((a.interop/layout/(field+'.json')).read_text()),(layout,field)
    assert report['canonicalNonzeros']==4
    run(layout+'-verify',['verify-aggregate',dest])
    run(layout+'-overwrite',['aggregate',a.interop/(layout+'.h5ad'),plan,dest],False)

# Duplicates crossing the HDF5 slice boundary must merge once per coordinate.
source=a.out/'boundary.h5ad';shutil.copyfile(a.interop/'csr.h5ad',source)
with h5py.File(source,'r+') as f:
    g=f['layers/counts']
    for key in ['data','indices','indptr']:del g[key]
    g.create_dataset('data',data=np.ones(65_539,dtype=np.uint64))
    g.create_dataset('indices',data=np.zeros(65_539,dtype=np.int32))
    g.create_dataset('indptr',data=np.array([0,65_539,65_539,65_539,65_539],dtype=np.int64))
run('boundary',['aggregate',source,plan,a.out/'boundary'])
r=json.loads((a.out/'boundary/report.json').read_text())
assert r['canonicalNonzeros']==1 and r['quality'][0]['totalCounts']==65_539 and r['quality'][0]['detectedFeatures']==1
assert r['pseudobulk']['matrix']['counts']==[65_539]
csc=a.out/'boundary-csc.h5ad';shutil.copyfile(source,csc)
with h5py.File(csc,'r+') as f:
    g=f['layers/counts'];g.attrs['encoding-type']='csc_matrix';del g['indptr']
    g.create_dataset('indptr',data=np.array([0,65_539,65_539,65_539],dtype=np.int64))
run('boundary-csc',['aggregate',csc,plan,a.out/'boundary-csc'])
rc=json.loads((a.out/'boundary-csc/report.json').read_text())
assert rc==r
empty=a.out/'empty.h5ad';ad.read_h5ad(a.interop/'csr.h5ad')[:0].copy().write_h5ad(empty)
run('empty',['aggregate',empty,plan,a.out/'empty'])
re=json.loads((a.out/'empty/report.json').read_text())
assert re['quality']==[] and re['pseudobulk']['groups']==[] and re['canonicalNonzeros']==0
oversized=json.loads(plan.read_text())
oversized['mapping']['samples']=[dict(mapping['samples'][0],id=str(i)+'x'*1020) for i in range(2500)]
huge=a.out/'oversized-plan.json';huge.write_text(json.dumps(oversized))
run('reject-oversized-plan',['aggregate',source,huge,a.out/'oversized-plan'],False)

for name in ['negative','nan','infinity','fractional','float-outside-exact-range']:
    run('reject-'+name,['aggregate',a.interop/(name+'.h5ad'),plan,a.out/('rejected-'+name)],False)
    assert not (a.out/('rejected-'+name)).exists()
for field in ['plan','report','source','implementation']:
    dest=a.out/('tamper-'+field);shutil.copytree(a.out/'csr',dest)
    if field=='source':
        with (dest/'original.h5ad').open('ab') as f:f.write(b'changed')
    elif field=='implementation':
        receipt=json.loads((dest/'receipt.json').read_text())
        receipt['implementation']=receipt['source']
        (dest/'receipt.json').write_text(json.dumps(receipt))
    else:
        with (dest/(field+'.json')).open('a') as f:f.write(' ')
    run('reject-tampered-'+field,['verify-aggregate',dest],False)

overflow=a.out/'overflow.h5ad';shutil.copyfile(source,overflow)
with h5py.File(overflow,'r+') as f:f['layers/counts/data'][0]=np.iinfo(np.uint64).max
run('reject-overflow',['aggregate',overflow,plan,a.out/'overflow'],False)
bad=a.out/'bad-index.h5ad';shutil.copyfile(source,bad)
with h5py.File(bad,'r+') as f:f['layers/counts/indices'][65_536]=3
run('reject-index-after-slice-boundary',['aggregate',bad,plan,a.out/'bad-index'],False)
assert not list(a.out.glob('.numivivo-stream-*'))
(a.out/'checks.json').write_text(json.dumps(dict(status='passed',fullProduct=a.full_product,checks=checks),indent=2)+'\n')
print(json.dumps(dict(status='passed',checks=len(checks))))
