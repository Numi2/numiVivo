#!/usr/bin/env python3
"""Frozen full-cohort HNSW matching and local Gaussian development candidate."""
import argparse,ctypes as C,hashlib,json,time
from pathlib import Path
import numpy as np

class Report(C.Structure):
 _fields_=[(name,C.c_uint64) for name in ('construction','query','indexBytes')]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,x):p.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
def identity(arrays):
 h=hashlib.sha256()
 for key,x in sorted(arrays.items()):
  x=np.ascontiguousarray(x);h.update(json.dumps([key,x.dtype.str,x.shape],separators=(',',':')).encode());h.update(x.tobytes())
 return h.hexdigest()
def neighbors(lib,source,queries,k,levels=None,budget=500000000):
 source=np.ascontiguousarray(source,dtype=np.float64);queries=np.ascontiguousarray(queries,dtype=np.float64);n,d=source.shape;q=len(queries);kind=int(levels is None)
 levels=np.ascontiguousarray(levels,dtype=np.uint32) if levels is not None else np.zeros(1,dtype=np.uint32)
 shape=(2,q,k) if kind==0 else (q,k);ids=np.empty(shape,dtype=np.int64);ds=np.empty(shape,dtype=np.float64);r=Report();f=lib.local_neighbors
 f.argtypes=[C.POINTER(C.c_double),C.POINTER(C.c_double),C.POINTER(C.c_uint32),C.c_uint32,C.c_uint32,C.c_uint32,C.c_uint32,C.c_uint32,C.c_int,C.c_uint64,C.POINTER(C.c_int64),C.POINTER(C.c_double),C.c_uint64,C.POINTER(Report)];f.restype=C.c_int
 start=time.perf_counter();code=f(source.ctypes.data_as(C.POINTER(C.c_double)),queries.ctypes.data_as(C.POINTER(C.c_double)),levels.ctypes.data_as(C.POINTER(C.c_uint32)),n,q,d,k,int(levels.max())+1,kind,budget,ids.ctypes.data_as(C.POINTER(C.c_int64)),ds.ctypes.data_as(C.POINTER(C.c_double)),ids.size,C.byref(r));seconds=time.perf_counter()-start
 assert code==0,('native HNSW status',code)
 assert np.isfinite(ds).all() and (ds>=0).all() and (ids<n).all()
 if kind==1:assert (ids>=0).all()
 return ids,ds,dict(construction=int(r.construction),query=int(r.query),indexBytes=int(r.indexBytes),seconds=seconds)
def matrix(path,n,d):
 r=np.fromfile(path,dtype=[('row','<u4'),('column','<u4'),('value','<f8')]).reshape(n,d)
 assert np.array_equal(r['row'],np.broadcast_to(np.arange(n)[:,None],(n,d))) and np.array_equal(r['column'],np.broadcast_to(np.arange(d),(n,d)))
 x=r['value'].copy();assert np.isfinite(x).all();return x

def fit(lib,x,batch,out,expected=None):
 n,d=x.shape;B=int(batch.max())+1;k=20;start=time.perf_counter();witness_seconds=0.;hashes={};native=[]
 def emit(name,**arrays):
  nonlocal witness_seconds
  t=time.perf_counter();h=identity(arrays);hashes[name]=h
  if expected is None:np.savez_compressed(out/(name+'.npz'),**arrays)
  else:assert expected[name]==h,('replay array mismatch',name)
  witness_seconds+=time.perf_counter()-t
 norms=np.zeros(n)
 for j in range(d):norms+=x[:,j]*x[:,j]
 norms=np.sqrt(norms);ordered=np.sort(norms);scale=float(ordered[n//2] if n%2 else ordered[n//2-1]/2+ordered[n//2]/2);assert scale>0
 unit=np.divide(x,norms[:,None],out=np.zeros_like(x),where=norms[:,None]>0);corrected=x/scale
 ids,distances,r=neighbors(lib,unit,unit,k,batch);native.append(dict(stage='matching',**r));emit('matching',indices=ids,distances=distances)
 lower,upper=ids;pairs={};rows=[np.flatnonzero(batch==b) for b in range(B)]
 for first in range(0,n,4096):
  end=min(n,first+4096);i=np.arange(first,end);j=upper[first:end];mask=(j>=0)&np.any(lower[np.maximum(j,0)]==i[:,None,None],axis=2)
  aa=np.broadcast_to(i[:,None],j.shape)[mask];bb=j[mask]
  for bi in np.unique(batch[aa]):
   for bj in np.unique(batch[bb]):
    choose=(batch[aa]==bi)&(batch[bb]==bj)
    if choose.any():pairs.setdefault((int(bi),int(bj)),[]).append(np.stack([aa[choose],bb[choose]],axis=1))
 for key,pieces in list(pairs.items()):
  values=np.concatenate(pieces);values=values[np.lexsort((values[:,1],values[:,0]))];assert len(np.unique(values,axis=0))==len(values);pairs[key]=values
 emit('anchors',**{f'{a}-{b}':v for (a,b),v in sorted(pairs.items())})
 alignments=[]
 for (a,b),v in sorted(pairs.items()):alignments.append(dict(first=a,second=b,anchors=len(v),score=max(len(np.unique(v[:,0]))/len(rows[a]),len(np.unique(v[:,1]))/len(rows[b]))))
 order=sorted([i for i,a in enumerate(alignments) if a['score']>.1],key=lambda i:(-alignments[i]['score'],-alignments[i]['first'],-alignments[i]['second']))
 panoramas=[];visits=np.zeros(B,dtype=int);steps=[]
 def oriented(a,b):
  if a<b:return pairs.get((a,b),np.empty((0,2),dtype=np.int64))
  return pairs.get((b,a),np.empty((0,2),dtype=np.int64))[:,::-1]
 def apply(target,reference,match,alignment):
  selected=np.concatenate([rows[b] for b in target]);source=corrected[match[:,0]].copy();bias=corrected[match[:,1]]-source;query=corrected[selected].copy();count=min(64,len(match));assert count>0
  chosen,ds,r=neighbors(lib,source,query,count);native.append(dict(stage='kernel-'+str(len(steps)),**r))
  delta=np.empty_like(query);totals=np.empty(len(query))
  for first in range(0,len(query),512):
   end=min(len(query),first+512);ix=chosen[first:end];w=np.exp(-7.5*ds[first:end]);total=w.sum(axis=1);totals[first:end]=total
   dy=np.einsum('ij,ijd->id',w,bias[ix],optimize=False);delta[first:end]=np.divide(dy,total[:,None],out=np.zeros_like(dy),where=total[:,None]>0)
  assert np.isfinite(delta).all();emit('step-'+str(len(steps)),selected=selected,anchors=match,source=source,bias=bias,queries=query,indices=chosen,distances=ds,delta=delta,totalWeights=totals)
  corrected[selected]+=delta
  return dict(alignment=alignment,target=target,reference=reference,anchors=len(match),correctedCells=len(selected),zeroWeightCells=int(np.sum(totals==0)),minimumWeight=float(totals.min()),maximumWeight=float(totals.max()),correctionRMS=float(np.sqrt(np.sum(delta*delta)/len(selected))*scale),status='corrected')
 for alignment in order:
  i,j=alignments[alignment]['first'],alignments[alignment]['second'];visits[i]+=1;visits[j]+=1
  if visits[i]>3 and visits[j]>3:steps.append(dict(alignment=alignment,status='both-level-visit-counts-exceed-three'));continue
  pi=next((p for p,v in enumerate(panoramas) if i in v),None);pj=next((p for p,v in enumerate(panoramas) if j in v),None)
  if pi is None and pj is None:
   if len(rows[i])<len(rows[j]):i,j=j,i
   panoramas.append([i]);pi=len(panoramas)-1
  if pi is None:
   ref=panoramas[pj].copy();match=np.concatenate([oriented(i,b) for b in ref]);steps.append(apply([i],ref,match,alignment));panoramas[pj].append(i)
  elif pj is None:
   ref=panoramas[pi].copy();match=np.concatenate([oriented(j,b) for b in ref]);steps.append(apply([j],ref,match,alignment));panoramas[pi].append(j)
  else:
   target=panoramas[pi].copy();ref=panoramas[pj].copy();one=oriented(i,j);steps.append(apply(target,ref,np.concatenate([one,one]),alignment))
   if pi!=pj:panoramas[pi]+=ref;panoramas.pop(pj)
 for b in range(B):
  if not any(b in p for p in panoramas):panoramas.append([b])
 corrected*=scale;assert np.isfinite(corrected).all();emit('scores',scores=corrected)
 return dict(cells=n,components=d,scale=scale,zeroDirectionCells=int(np.sum(norms==0)),alignments=alignments,order=order,steps=steps,panoramas=panoramas,native=native,arraySHA256=hashes,fitSeconds=time.perf_counter()-start-witness_seconds,witnessSeconds=witness_seconds)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for n in ('source','manifest','library','out','protocol'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);owner=Path(__file__).parent;manifest=read(a.manifest);lib=C.CDLL(str(a.library));inputs={}
 for c in ('hagai','kang','ding'):
  entries=next(v['files'] for v in manifest if v['cohort']==c)
  for name in ('pca/scores.bin','pca/metadata.json','mnn/report.json','mnn/plan.json','mnn/scores.bin','mnn/anchors.bin'):
   entry=next(v for v in entries if v['path']==name);assert sha(a.source/c/name)==entry['sha256'];inputs[c+'/'+name]=entry['sha256']
 write(a.out/'freeze.json',dict(createdUnix=time.time(),protocolSHA256=sha(a.protocol),protocolName=a.protocol.name,sourceSHA256={p.name:sha(p) for p in owner.iterdir() if p.suffix in ('.py','.cpp','.md')},librarySHA256=sha(a.library),originalArtifactManifestSHA256=sha(a.manifest),inputs=inputs,metricsRead=False))
 for c in ('hagai','kang','ding'):
  out=a.out/c;out.mkdir();root=a.source/c;old=read(root/'mnn/report.json');meta=read(root/'pca/metadata.json');options=read(root/'mnn/plan.json')['mnn'];assert options['neighbors']==20 and options['sigma']==15 and options['minimumAlignment']==.1
  samples={s['id']:s for s in meta['samples']};field='batchID' if options['covariate']=='batch' else 'donorID';names=[samples[v['sampleID']][field] for v in meta['cells']];levels=sorted(set(names));batch=np.array([levels.index(v) for v in names],dtype=np.uint32);assert levels==old['levels'] and batch.tolist()==old['cellLevels']
  x=matrix(root/'pca/scores.bin',len(batch),20);r=fit(lib,x,batch,out);write(out/'report.json',r);print(c,'fit',r['fitSeconds'],flush=True)
  again=fit(lib,x,batch,out,r['arraySHA256']);write(out/'replay.json',again)
  for key in ('cells','components','scale','zeroDirectionCells','alignments','order','steps','panoramas','arraySHA256'):assert r[key]==again[key],key
  for first,second in zip(r['native'],again['native']):
   for key in ('stage','construction','query','indexBytes'):assert first[key]==second[key]
  print(c,'replay',again['fitSeconds'],flush=True)
 for n,h in inputs.items():assert sha(a.source/n)==h
 write(a.out/'outputs-frozen.json',dict(status='fit-and-replay-complete',allOriginalCells=True,allReplayArraysExact=True,librarySHA256=sha(a.library),files={str(p.relative_to(a.out)):sha(p) for p in a.out.rglob('*') if p.is_file()},metricsRead=False))
if __name__=='__main__':main()
