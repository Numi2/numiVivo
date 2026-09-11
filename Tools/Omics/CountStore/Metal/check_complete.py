#!/usr/bin/env python3
"""Compare every normalized original-Kang record against FP64 Scanpy."""
import argparse,gzip,hashlib,importlib.metadata,json,platform,resource,shlex,subprocess,time
from pathlib import Path
import numpy as np
import scanpy as sc
from scipy import sparse
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
source=Path('/Users/home/numivivo-legacy-count-route-20260911/route/count-store/original.h5ad')
remote='/Users/n/numivivo-metal-normalization-20260911'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def fetch(name):
 b=subprocess.check_output(['ssh','-4','macmini','cat '+shlex.quote(remote+'/'+name)])
 out=a.out/name;out.parent.mkdir(parents=True,exist_ok=True);out.write_bytes(b);return json.loads(b)
assert sha(source)=='89c2c4d6aea0163676b45574aaafd3593525dce269ee315911680d8f31c014a0'
native=fetch('native-checks.json');assert native['status']=='passed'
t=time.monotonic();obj=sc.read_h5ad(source);load_seconds=time.monotonic()-t
assert obj.shape==(24673,15706)
obj.X=sparse.csc_matrix(obj.X,dtype=np.float64);obj.X.sum_duplicates();obj.X.sort_indices()
assert obj.X.nnz==14184532
before=obj.X.copy();t=time.monotonic();sc.pp.normalize_total(obj,target_sum=10000);sc.pp.log1p(obj);transform_seconds=time.monotonic()-t
reference=obj.X.tocsc();reference.sort_indices();assert reference.nnz==before.nnz
local_store=source.parent
metadata=json.loads((local_store/'metadata.json').read_text());quality=json.loads((local_store/'quality.json').read_text())
np.testing.assert_array_equal([c['barcode'] for c in metadata['cells']],obj.obs_names)
np.testing.assert_array_equal(quality['rowTotals'],np.asarray(before.sum(axis=1)).ravel())
np.testing.assert_array_equal(quality['featureTotals'],np.asarray(before.sum(axis=0)).ravel())
raw=np.memmap(local_store/'counts.bin',mode='r',dtype=[('row','<u4'),('feature','<u4'),('value','<u8')])
assert len(raw)==reference.nnz
for start in range(0,len(raw),262144):
 stop=min(start+262144,len(raw));columns=np.searchsorted(before.indptr[1:],np.arange(start,stop),side='right')
 np.testing.assert_array_equal(raw['row'][start:stop],before.indices[start:stop]);np.testing.assert_array_equal(raw['feature'][start:stop],columns)
 np.testing.assert_array_equal(raw['value'][start:stop],before.data[start:stop])
results=[]
for label in ['cpu-0','metal-0']:
 receipt=fetch(label+'/receipt.json');remote_metadata=fetch(label+'/metadata.json')
 assert remote_metadata==metadata
 proc=subprocess.Popen(['ssh','-4','macmini','gzip -1 -c '+shlex.quote(remote+'/'+label+'/values.bin')],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 digest=hashlib.sha256();position=0;maximum=0.;max_scaled=0.;violations=0
 with gzip.GzipFile(fileobj=proc.stdout) as f:
  while position<len(raw):
   n=min(262144,len(raw)-position);b=f.read(n*16);assert len(b)==n*16;digest.update(b)
   v=np.frombuffer(b,dtype=[('row','<u4'),('feature','<u4'),('value','<f8')]);end=position+n
   np.testing.assert_array_equal(v['row'],raw['row'][position:end]);np.testing.assert_array_equal(v['feature'],raw['feature'][position:end])
   expected=reference.data[position:end];error=np.abs(v['value']-expected)
   assert np.all(np.isfinite(v['value'])) and np.all(v['value']>0)
   tolerance=3e-6+3e-6*np.abs(expected) if label=='metal-0' else np.full(n,1e-14)
   violations+=int((error>tolerance).sum());maximum=max(maximum,float(error.max()));max_scaled=max(max_scaled,float((error/tolerance).max()))
   if label=='metal-0':np.testing.assert_array_equal(v['value'],v['value'].astype(np.float32).astype(np.float64))
   position=end
  assert f.read(1)==b''
 assert proc.wait()==0,proc.stderr.read().decode()
 assert digest.hexdigest()==bytes(receipt['values']['bytes']).hex() and violations==0
 results.append(dict(label=label,entries=position,SHA256=digest.hexdigest(),maximumAbsoluteError=maximum,maximumToleranceFraction=max_scaled,violations=violations))
result=dict(status='passed',sourceSHA256=sha(source),completeCells=24673,completeFeatures=15706,completeEntries=14184532,comparisons=results,allCoordinatesExact=True,metadataExact=True,reference=dict(platform=platform.platform(),loadH5ADSeconds=load_seconds,transformOnlySeconds=transform_seconds,peakProcessRSSBytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,packages={k:importlib.metadata.version(k) for k in ['numpy','scipy','anndata','scanpy','h5py']},scope='Different host and I/O boundary from native; descriptive reference timing, no speed comparison'),checkerSHA256=sha(Path(__file__)))
(a.out/'checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
