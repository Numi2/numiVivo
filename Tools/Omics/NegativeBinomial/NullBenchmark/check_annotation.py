#!/usr/bin/env python3
"""Verify native sham annotations preserve every original HDF5 field and count."""
import argparse,hashlib,json
from pathlib import Path
import h5py,numpy as np
p=argparse.ArgumentParser(description=__doc__)
for name in ['source','annotated','plan','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();assert not a.out.exists();plan=json.loads(a.plan.read_text());checked=[]
def sha(path):
 with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def equal(x,y):
 x,y=np.asarray(x),np.asarray(y);assert x.shape==y.shape
 assert np.array_equal(x,y,equal_nan=True) if x.dtype.kind in 'fc' else np.array_equal(x,y)
assert sha(a.source)==bytes(plan['source']['bytes']).hex()
with h5py.File(a.source) as source,h5py.File(a.annotated) as annotated:
 def visit(name):
  x,y=source[name],annotated[name];assert isinstance(x,h5py.Dataset)==isinstance(y,h5py.Dataset)
  assert set(x.attrs)==set(y.attrs)
  for k,v in x.attrs.items():
   if name=='obs' and k=='column-order':equal(list(v)+[e['path'].split('/')[1] for e in plan['edits']],y.attrs[k])
   else:equal(v,y.attrs[k])
  if isinstance(x,h5py.Group):
   expected={e['path'].split('/')[1] for e in plan['edits']} if name=='obs' else ({'numivivo_edits'} if name=='uns' else set())
   assert set(y)==set(x)|expected
  else:
   assert x.shape==y.shape and x.dtype==y.dtype
   assert (x.chunks,x.compression,x.compression_opts,x.shuffle,x.fletcher32)==(y.chunks,y.compression,y.compression_opts,y.shuffle,y.fletcher32)
   if x.ndim==0:equal(x[()],y[()])
   else:
    for start in range(0,len(x),65536):equal(x[start:start+65536],y[start:start+65536])
   checked.append(name)
 visit('/');source.visit(visit)
 for edit in plan['edits']:
  g=annotated[edit['path']];value=edit['value']['categorical']
  equal(g['codes'][:],value['codes']);equal(g['categories'].asstr()[:],value['categories'])
  assert bool(g.attrs['ordered'])==value['ordered']
 journals=annotated['uns/numivivo_edits']
 assert any(json.loads(journals[k+'/plan_json'].asstr()[()])==plan for k in journals)
a.out.write_text(json.dumps(dict(status='passed-all-original-HDF5-fields-and-annotations',sourceSHA256=sha(a.source),annotatedSHA256=sha(a.annotated),originalDatasets=checked,addedColumns=len(plan['edits'])),indent=2)+'\n')
print(json.dumps(dict(status='passed',originalDatasets=len(checked),addedColumns=len(plan['edits']))))
