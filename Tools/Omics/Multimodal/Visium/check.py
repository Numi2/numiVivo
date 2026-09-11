#!/usr/bin/env python3
"""Independently compare every native Visium count/spot/coordinate with the source."""
import argparse,gzip,hashlib,importlib.metadata,json,platform,tarfile,tempfile,shutil
from pathlib import Path
import h5py,numpy as np,pandas as pd,scanpy as sc
from scipy import sparse
from prepare import sha,write

def digest_stream(f):
 h=hashlib.sha256()
 for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--native',type=Path,required=True);p.add_argument('--prior',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 a.out.mkdir(parents=True,exist_ok=False);execution=json.loads((a.native/'execution.json').read_text());source=a.inputs/'outs/filtered_feature_bc_matrix.h5'
 with h5py.File(source) as f:
  m=f['matrix'];x=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).T.tocsr();x.sort_indices()
  ids=m['features/id'].asstr()[:].tolist();names=m['features/name'].asstr()[:].tolist();barcodes=m['barcodes'].asstr()[:].tolist()
 positions=pd.read_csv(a.inputs/'outs/spatial/tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','y','x'],index_col=0)
 assert positions.index.is_unique and set(positions.index[positions.in_tissue==1])==set(barcodes)
 xy=positions.loc[barcodes,['x','y']].to_numpy(dtype=np.float64)
 with tempfile.TemporaryDirectory(prefix='visium-reference-',dir=a.out) as d:
  d=Path(d);shutil.copy2(source,d/'filtered_feature_bc_matrix.h5');(d/'spatial').mkdir()
  with tarfile.open(a.prior/'source/spatial.tar.gz') as t:
   for name in ['spatial/tissue_positions_list.csv','spatial/scalefactors_json.json','spatial/tissue_hires_image.png','spatial/tissue_lowres_image.png']:
    m=t.getmember(name);assert m.isfile() and m.size<10000000;(d/name).write_bytes(t.extractfile(m).read())
  reference=sc.read_visium(d)
  assert reference.obs_names.tolist()==barcodes and reference.var.gene_ids.tolist()==ids and (reference.X!=x).nnz==0
  np.testing.assert_array_equal(reference.obsm['spatial'],xy)
 old={}
 for name in ['dataset.json','dataset.h5mu']:
  with gzip.open(a.prior/'bundle'/(name+'.gz'),'rb') as f:old[name]=digest_stream(f)
 results=[]
 for item in execution['archives']:
  path=a.native/item['path'];assert sha(path)==item['SHA256']
  with tarfile.open(path,'r:gz') as t:
   assert {m.name for m in t.getmembers()}==set(item['members'])
   for m in t.getmembers():
    assert m.isfile() and m.size==item['members'][m.name]['bytes']
    with t.extractfile(m) as f:assert digest_stream(f)==item['members'][m.name]['SHA256']
   native=json.load(t.extractfile('dataset.json'))
  assert len(native['observations'])==len(barcodes)
  assert all(o['kind']=='spot' for o in native['observations'])
  assert [o['identity']['barcode'] for o in native['observations']]==barcodes
  np.testing.assert_array_equal([o['position']['coordinates'] for o in native['observations']],xy)
  assay=native['assays'][0];mat=assay['matrix']
  assert [f['id'] for f in assay['features']]==ids and [f['name'] for f in assay['features']]==names
  assert assay['observationIndices']==list(range(len(barcodes)))
  np.testing.assert_array_equal(mat['rowOffsets'],x.indptr);np.testing.assert_array_equal(mat['featureIndices'],x.indices);np.testing.assert_array_equal(mat['counts'],x.data)
  assert native['samples']==[json.loads((a.inputs/'legacy-plan.json').read_text())['counts']['sample']]
  # Exact previous qualified dataset and H5MU bytes also bind all sample,
  # feature-space, frame and missing-value metadata, without ignored fields.
  for name,h in old.items():assert item['members'][name]['SHA256']==h,(item['label'],name)
  results.append(dict(format=item['label'],allCountsAndCoordinatesExact=True,priorQualifiedDatasetAndH5MUBytesExact=True,archiveSHA256=item['SHA256']))
 write(a.out/'checks.json',dict(status='passed',spots=x.shape[0],features=x.shape[1],nonzeros=x.nnz,totalUMIs=int(x.sum()),sourcePositionRows=len(positions),outsideTissueRows=int((positions.in_tissue==0).sum()),
  sourceSHA256=sha(source),inputFreezeSHA256=sha(a.inputs/'input-freeze.json'),executionSHA256=sha(a.native/'execution.json'),checkerSHA256=sha(__file__),results=results,
  independentScanpyCountsAndCoordinatesExact=True,priorQualifiedHashes=old,packages={n:importlib.metadata.version(n) for n in ['numpy','scipy','pandas','scanpy','h5py']},platform=platform.platform(),scope='Complete count and pixel-coordinate interchange; no spatial biological or predictive endpoint'))
 print((a.out/'checks.json').read_text())
if __name__=='__main__':main()
