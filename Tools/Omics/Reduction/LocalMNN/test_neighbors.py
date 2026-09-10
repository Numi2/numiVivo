#!/usr/bin/env python3
import argparse,ctypes,json
from pathlib import Path
import numpy as np
from fit import neighbors,sha
p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();lib=ctypes.CDLL(str(a.library));rng=np.random.default_rng(20260911);cases=[]
for d in (1,20,64):
 source=rng.normal(size=(120,d));query=rng.normal(size=(19,d));levels=np.arange(120,dtype=np.uint32)%3
 ids,ds,r=neighbors(lib,source,query,64)
 expected=np.zeros((len(query),len(source)))
 for j in range(d):expected+=(query[:,j,None]-source[None,:,j])**2
 np.testing.assert_allclose(ds,np.take_along_axis(expected,ids,axis=1),rtol=1e-12,atol=1e-12)
 for i in range(len(query)):np.testing.assert_array_equal(ids[i],np.lexsort((np.arange(len(source)),expected[i]))[:64])
 ix,dd,rr=neighbors(lib,source,query,64);np.testing.assert_array_equal(ix,ids);np.testing.assert_array_equal(dd,ds)
 lowerupper,dist,r=neighbors(lib,source,source,20,levels)
 exact=np.zeros((len(source),len(source)))
 for j in range(d):exact+=np.abs(source[:,j,None]-source[None,:,j])
 for direction in range(2):
  for i in range(len(source)):
   candidates=np.flatnonzero(levels>levels[i] if direction else levels<levels[i]);got=lowerupper[direction,i]
   if not len(candidates):assert (got==-1).all();continue
   selected=candidates[np.lexsort((candidates,exact[i,candidates]))[:20]];np.testing.assert_array_equal(got,selected);np.testing.assert_allclose(dist[direction,i],exact[i,got],rtol=1e-12,atol=1e-12)
 cases.append(dict(dimensions=d,externalSquaredL2Exact=True,filteredL1Exact=True,replayExact=True))
for kwargs in ({'budget':1},):
 try:neighbors(lib,source,query,64,**kwargs)
 except AssertionError as e:assert '2' in str(e)
 else:raise AssertionError('budget accepted')
source[0,0]=np.nan
try:neighbors(lib,source,query,64)
except AssertionError as e:assert '4' in str(e)
else:raise AssertionError('nonfinite accepted')
a.out.write_text(json.dumps(dict(status='passed',cases=cases,budgetRejected=True,nonfiniteRejected=True,librarySHA256=sha(a.library),testSHA256=sha(Path(__file__)),scope='Small exhaustive numerical bridge tests; large-cohort HNSW remains approximate.'),indent=2)+'\n');print('passed: both metrics, filters, replay and admission')
