#!/usr/bin/env python3
"""Compare every remote count/normalized record against complete local Norman data.

The native host and reference host may differ: this is numerical interoperability
and native memory evidence, not a cross-host speed comparison. Large binary files
are streamed over SSH without a second local disk copy.
"""
import argparse,gzip,hashlib,importlib.metadata,json,platform,shlex,subprocess,time
from pathlib import Path
import h5py
import numpy as np
import scanpy as sc
from scipy import sparse

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source',type=Path,required=True)
p.add_argument('--remote-root',required=True)
p.add_argument('--host',default='macmini')
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
def sha(path):
 with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
assert sha(a.source)=='efde6f5301fe256725dce1d980f37bd96a13481a9a16135515897368e631affc'
def fetch(name):
 return subprocess.check_output(['ssh',a.host,'cat -- '+shlex.quote(a.remote_root+'/'+name)])
def get(name):
 raw=fetch(name);dest=a.out/name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(raw);return json.loads(raw)
status=get('native-checks.json');assert status['status']=='passed-native-reconstruction'
receipt=get('store/receipt.json');normalized=get('normalized/receipt.json');metadata=get('store/metadata.json');quality=get('store/quality.json');plan=get('store/plan.json')
get('commands.json')
for label in ['import','normalize','normalize-verify','run']:
 (a.out/(label+'.log')).write_bytes(fetch(label+'.log'))
started=time.monotonic();obj=sc.read_h5ad(a.source)
assert obj.shape==(111445,33694)
matrix=sparse.csc_matrix(obj.X);matrix.sum_duplicates();matrix.eliminate_zeros();matrix.sort_indices()
assert matrix.nnz==receipt['entries']==361582621
np.testing.assert_array_equal([v['barcode'] for v in metadata['cells']],obj.obs_names)
np.testing.assert_array_equal([v['sampleID'] for v in metadata['cells']],obj.obs[plan['sampleColumn']].astype(str))
np.testing.assert_array_equal([v['id'] for v in metadata['features']],obj.var[plan['featureIDColumn']] if plan.get('featureIDColumn') else obj.var_names)
np.testing.assert_array_equal([v['name'] for v in metadata['features']],obj.var[plan['featureNameColumn']] if plan.get('featureNameColumn') else obj.var_names)
np.testing.assert_array_equal(quality['rowTotals'],np.asarray(matrix.sum(axis=1,dtype=np.uint64)).ravel())
np.testing.assert_array_equal(quality['featureTotals'],np.asarray(matrix.sum(axis=0,dtype=np.uint64)).ravel())
np.testing.assert_array_equal(quality['rowNonzeros'],matrix.getnnz(axis=1))
np.testing.assert_array_equal(quality['featureNonzeros'],matrix.getnnz(axis=0))

def compare(path, expected, floating):
 dtype=np.dtype([('row','<u4'),('feature','<u4'),('value','<f8' if floating else '<u8')])
 proc=subprocess.Popen(['ssh',a.host,'gzip -c -- '+shlex.quote(a.remote_root+'/'+path)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 stream=gzip.GzipFile(fileobj=proc.stdout,mode='rb')
 h=hashlib.sha256();position=0;maximum=0.
 try:
  while position<expected.nnz:
   n=min(262144,expected.nnz-position);parts=[];remaining=n*16
   while remaining:
    block=stream.read(remaining)
    assert block,'truncated remote records'
    parts.append(block);remaining-=len(block)
   raw=b''.join(parts);h.update(raw);records=np.frombuffer(raw,dtype=dtype)
   np.testing.assert_array_equal(records['row'],expected.indices[position:position+n])
   columns=np.searchsorted(expected.indptr[1:],np.arange(position,position+n),side='right')
   np.testing.assert_array_equal(records['feature'],columns)
   if floating:
    ref=expected.data[position:position+n]
    maximum=max(maximum,float(np.max(np.abs(records['value']-ref))))
    np.testing.assert_allclose(records['value'],ref,rtol=1e-12,atol=1e-12)
   else:
    ref=expected.data[position:position+n]
    assert np.all(np.isfinite(ref)) and np.all(ref>=0) and np.all(ref==np.floor(ref))
    np.testing.assert_array_equal(records['value'],ref.astype(np.uint64))
   position+=n
  assert stream.read(1)==b'', 'extra remote records'
  assert proc.wait()==0,proc.stderr.read().decode()
 finally:
  if proc.poll() is None:proc.terminate();proc.wait()
 return dict(entriesCompared=position,sha256=h.hexdigest(),maximumAbsoluteError=maximum)
raw=compare('store/counts.bin',matrix,False)
assert raw['sha256']==bytes(receipt['counts']['bytes']).hex()
# Explicit FP64 reference. Scanpy's normalization and log1p define this comparison.
obj.X=matrix.astype(np.float64);del matrix
sc.pp.normalize_total(obj,target_sum=10000,inplace=True);sc.pp.log1p(obj)
reference=obj.X.tocsc();reference.sum_duplicates();reference.eliminate_zeros();reference.sort_indices()
assert reference.nnz==receipt['entries']
norm=compare('normalized/values.bin',reference,True)
assert norm['sha256']==bytes(normalized['values']['bytes']).hex()
(a.out/'checks.json').write_text(json.dumps(dict(status='passed',sourceSHA256=sha(a.source),native=status,raw=raw,normalized=norm,
 transport='gzip over SSH; digests and comparisons use every uncompressed record',allMetadataAxesExact=True,allCountMarginalsExact=True,scanpyFloat64NormalizationTolerance=1e-12,referenceSeconds=time.monotonic()-started,
 referencePlatform=platform.platform(),packages={n:importlib.metadata.version(n) for n in ['numpy','scipy','scanpy','anndata','h5py']},
 checkerSHA256=sha(Path(__file__)),scope='Full 111445-cell source; no million-cell, parallel-kernel or GPU qualification; hosts differ so no speedup claim'),indent=2)+'\n')
