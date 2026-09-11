#!/usr/bin/env python3
"""Fixed eligible-level approximate matching with qualified all-anchor Gaussian."""
import argparse,ctypes as C,sys,time
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'LocalMNN'))
from fit import fit,sha,read,write,matrix
D=C.POINTER(C.c_double);U=C.c_uint32;L=C.c_uint64
class Gaussian:
 def __init__(self,path):
  self.lib=C.CDLL(str(path));self.f=self.lib.nvivo_omics_gaussian_bias
  self.f.restype=C.c_int32;self.f.argtypes=[D,D,U,D,U,U,C.c_double,L,L,L,D,L,D,L,D,L,C.POINTER(L),C.c_void_p,C.c_void_p]
 def __call__(self,source,bias,queries):
  start=time.perf_counter();source,bias,queries=[np.ascontiguousarray(x,dtype=np.float64) for x in (source,bias,queries)]
  na,d=source.shape;nq=len(queries);assert bias.shape==source.shape and queries.shape==(nq,d)
  delta=np.empty_like(queries);totals=np.empty(nq);workspace=np.empty(256*(d+32));terms=0;budget=2*na*nq*d
  ptr=lambda x:x.ctypes.data_as(D)
  for first in range(0,nq,32):
   q=queries[first:first+32];out=delta[first:first+32];total=totals[first:first+32];evaluated=L()
   status=self.f(ptr(source),ptr(bias),na,ptr(q),len(q),d,15.,source.size,q.size,budget-terms,ptr(out),out.size,ptr(total),len(total),ptr(workspace),workspace.size,C.byref(evaluated),None,None)
   assert status==0,('Gaussian status',status);terms+=evaluated.value
  assert na*nq*d<=terms<=budget and np.isfinite(delta).all() and np.isfinite(totals).all()
  return delta,totals,dict(distanceTerms=terms,maximumDistanceTerms=budget,workspaceBytes=workspace.nbytes,seconds=time.perf_counter()-start)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for n in ('source','manifest','library','gaussian-library','out','protocol'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);owner=Path(__file__).parent;manifest=read(a.manifest);lib=C.CDLL(str(a.library));kernel=Gaussian(a.gaussian_library);inputs={}
 for c in ('hagai','kang','ding'):
  entries=next(v['files'] for v in manifest if v['cohort']==c)
  for name in ('pca/scores.bin','pca/metadata.json','mnn/report.json','mnn/plan.json','mnn/scores.bin','mnn/anchors.bin'):
   entry=next(v for v in entries if v['path']==name);assert sha(a.source/c/name)==entry['sha256'];inputs[c+'/'+name]=entry['sha256']
 write(a.out/'freeze.json',dict(createdUnix=time.time(),protocolSHA256=sha(a.protocol),protocolName=a.protocol.name,sourceSHA256={str(p.relative_to(owner.parent)):sha(p) for p in list(owner.glob('*.py'))+list(owner.glob('*.md'))+[owner.parent/'LocalMNN/fit.py',owner.parent/'LocalMNN/LocalNeighbors.cpp']},librarySHA256=sha(a.library),gaussianLibrarySHA256=sha(a.gaussian_library),originalArtifactManifestSHA256=sha(a.manifest),inputs=inputs,metricsRead=False))
 for c in ('hagai','kang','ding'):
  out=a.out/c;out.mkdir();root=a.source/c;old=read(root/'mnn/report.json');meta=read(root/'pca/metadata.json');options=read(root/'mnn/plan.json')['mnn'];assert options['neighbors']==20 and options['sigma']==15 and options['minimumAlignment']==.1
  samples={s['id']:s for s in meta['samples']};field='batchID' if options['covariate']=='batch' else 'donorID';names=[samples[v['sampleID']][field] for v in meta['cells']];levels=sorted(set(names));batch=np.array([levels.index(v) for v in names],dtype=np.uint32);assert levels==old['levels'] and batch.tolist()==old['cellLevels']
  x=matrix(root/'pca/scores.bin',len(batch),20);r=fit(lib,x,batch,out,kernel=kernel);write(out/'report.json',r);print(c,'fit',r['fitSeconds'],flush=True)
  again=fit(lib,x,batch,out,r['arraySHA256'],kernel=kernel);write(out/'replay.json',again)
  for key in ('cells','components','scale','zeroDirectionCells','alignments','order','steps','panoramas','arraySHA256'):assert r[key]==again[key],key
  for first,second in zip(r['native'],again['native']):
   assert {k:v for k,v in first.items() if k!='seconds'}=={k:v for k,v in second.items() if k!='seconds'}
  print(c,'replay',again['fitSeconds'],flush=True)
 for n,h in inputs.items():assert sha(a.source/n)==h
 write(a.out/'outputs-frozen.json',dict(status='fit-and-replay-complete',allOriginalCells=True,allReplayArraysExact=True,librarySHA256=sha(a.library),files={str(p.relative_to(a.out)):sha(p) for p in a.out.rglob('*') if p.is_file()},metricsRead=False))
if __name__=='__main__':main()
