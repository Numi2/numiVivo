#!/usr/bin/env python3
"""Persistent count-store CLI controls on AnnData CSR/CSC/dense inputs."""
import argparse,copy,hashlib,json,shutil,subprocess
from pathlib import Path
import anndata as ad
import numpy as np
from scipy import sparse

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--interop',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
checks=[]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_text(json.dumps(v,sort_keys=True,separators=(',',':'),ensure_ascii=False))
def run(label,args,success=True):
 r=subprocess.run([str(a.binary),*map(str,args)],capture_output=True,text=True)
 (a.out/(label+'.log')).write_text(r.stdout+r.stderr)
 checks.append(dict(label=label,arguments=list(map(str,args)),exitCode=r.returncode,expectedSuccess=success))
 write(a.out/'commands.json',checks)
 assert (r.returncode==0)==success,(label,r.stderr[-2500:])
record=np.dtype([('row','<u4'),('feature','<u4'),('count','<u8')])
values=np.dtype([('row','<u4'),('feature','<u4'),('value','<f8')])
for layout in ['csr','csc','dense']:
 source=a.interop/(layout+'.h5ad');plan=a.interop/(layout+'-plan.json');store=a.out/layout
 run(layout,['singlecell-h5ad-store',source,'--plan',plan,'--output',store])
 run(layout+'-verify',['singlecell-count-store-verify',store])
 run(layout+'-overwrite',['singlecell-h5ad-store',source,'--plan',plan,'--output',store],False)
 run(layout+'-normalize',['singlecell-count-store-normalize',store,'--target','10000','--output',a.out/(layout+'-normalized')])
 run(layout+'-normalized-verify',['singlecell-count-store-normalize-verify',a.out/(layout+'-normalized'),'--store',store])
 obj=ad.read_h5ad(source);matrix=sparse.csr_matrix(obj.layers['counts']);matrix.sum_duplicates();matrix.eliminate_zeros();matrix.sort_indices()
 rows,cols=matrix.nonzero();expected={(int(r),int(c)):int(matrix[r,c]) for r,c in zip(rows,cols)}
 raw=np.fromfile(store/'counts.bin',dtype=record)
 assert len(raw)==len(expected)
 assert {(int(v['row']),int(v['feature'])):int(v['count']) for v in raw}==expected
 qc=json.loads((store/'quality.json').read_text())
 assert qc['rowTotals']==np.asarray(matrix.sum(axis=1,dtype=np.uint64)).ravel().tolist()
 assert qc['featureTotals']==np.asarray(matrix.sum(axis=0,dtype=np.uint64)).ravel().tolist()
 assert qc['rowNonzeros']==matrix.getnnz(axis=1).tolist() and qc['featureNonzeros']==matrix.getnnz(axis=0).tolist()
 norm=np.fromfile(a.out/(layout+'-normalized')/'values.bin',dtype=values)
 np.testing.assert_array_equal(raw['row'],norm['row']);np.testing.assert_array_equal(raw['feature'],norm['feature'])
 reference=np.log1p(raw['count'].astype(np.float64)*(10000/np.asarray(qc['rowTotals'],dtype=np.float64)[raw['row']]))
 np.testing.assert_allclose(norm['value'],reference,rtol=1e-14,atol=1e-14)
 assert sha(source)==sha(store/'original.h5ad')
run('repeat',['singlecell-h5ad-store',a.interop/'csr.h5ad','--plan',a.interop/'csr-plan.json','--output',a.out/'repeat'])
for name in ['original.h5ad','counts.bin','metadata.json','quality.json','plan.json','receipt.json']:assert sha(a.out/'csr'/name)==sha(a.out/'repeat'/name)
for label in ['negative','fractional','nan','infinity','float-outside-exact-range']:
 run(label,['singlecell-h5ad-store',a.interop/(label+'.h5ad'),'--plan',a.interop/'csr-plan.json','--output',a.out/label],False)
 assert not (a.out/label).exists()
for label,target in [('zero-target','0'),('nan-target','nan'),('infinite-target','inf')]:
 run(label,['singlecell-count-store-normalize',a.out/'csr','--target',target,'--output',a.out/label],False)
 assert not (a.out/label).exists()
# Rehashed record mutation must fail source reconstruction, not only a digest check.
tamper=a.out/'rehashed-count';shutil.copytree(a.out/'csr',tamper)
x=np.fromfile(tamper/'counts.bin',dtype=record);x['count'][0]+=np.uint64(1);x.tofile(tamper/'counts.bin')
r=json.loads((tamper/'receipt.json').read_text());r['counts']={'bytes':list(bytes.fromhex(sha(tamper/'counts.bin')))};write(tamper/'receipt.json',r)
run('rehashed-count',['singlecell-count-store-verify',tamper],False)
tamper=a.out/'rehashed-normalized';shutil.copytree(a.out/'csr-normalized',tamper)
x=np.fromfile(tamper/'values.bin',dtype=values);x['value'][0]+=.1;x.tofile(tamper/'values.bin')
r=json.loads((tamper/'receipt.json').read_text());r['values']={'bytes':list(bytes.fromhex(sha(tamper/'values.bin')))};write(tamper/'receipt.json',r)
run('rehashed-normalized',['singlecell-count-store-normalize-verify',tamper,'--store',a.out/'csr'],False)
assert not list(a.out.glob('.numivivo-count-store-*'))
write(a.out/'checks.json',dict(status='passed',commands=len(checks),expectedRejections=sum(not v['expectedSuccess'] for v in checks),allCountsExact=True,normalizationReferenceTolerance=1e-14,repeatStoreBytesExact=True,binarySHA256=sha(a.binary),checkerSHA256=sha(Path(__file__))))
