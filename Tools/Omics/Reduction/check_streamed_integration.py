#!/usr/bin/env python3
"""Independently reconstruct every ridge-integration output in bounded row batches."""
import argparse,contextlib,hashlib,json,platform,resource,time
from pathlib import Path
import ijson
import numpy as np
DT=np.dtype([('row','<u4'),('column','<u4'),('value','<f8')])

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def blocks(spec,n,tile=8192):
 with contextlib.ExitStack() as stack:
  handles={}
  for name,(path,width) in spec.items():
   assert path.is_file() and not path.is_symlink() and path.stat().st_size==n*width*16,(name,'size')
   handles[name]=stack.enter_context(path.open('rb'))
  for start in range(0,n,tile):
   end=min(n,start+tile);values={}
   for name,(path,width) in spec.items():
    raw=handles[name].read((end-start)*width*16);assert len(raw)==(end-start)*width*16,(name,'truncation')
    a=np.frombuffer(raw,dtype=DT).reshape(end-start,width)
    assert np.all(a['row']==np.arange(start,end,dtype=np.uint32)[:,None]),(name,'row coordinates')
    assert np.all(a['column']==np.arange(width,dtype=np.uint32)[None,:]),(name,'column coordinates')
    assert np.isfinite(a['value']).all(),(name,'nonfinite')
    values[name]=a['value']
   yield start,end,values
  assert all(f.read(1)==b'' for f in handles.values())

def check(bundle,out,expected_metadata=None):
 began=time.time();assert not out.exists();assert (bundle/'report.json').stat().st_size<=16777216
 report=json.loads((bundle/'report.json').read_text());plan=json.loads((bundle/'plan.json').read_text());assert plan.get('mnn') is None
 o=dict(clusters=88,covariate='donor',ridge=1.,temperature=.1,diversity=2.,maximumIterations=10,relativeTolerance=.01,ridgeScaling=None);o.update(plan['integration'])
 n,d,k=report['cells'],report['components'],report['clusters'];assert 1<=n<=2000000 and 1<=d<=64 and 2<=k<=100 and k==o['clusters']
 levels=report['levels'];assert levels==sorted(set(levels)) and 2<=len(levels)<=128;bcount=len(levels)
 assert len(report['cellLevels'])==n and all(type(v) is int for v in report['cellLevels'])
 batch=np.asarray(report['cellLevels'],dtype=np.int64);assert np.all((batch>=0)&(batch<bcount));sizes=np.bincount(batch,minlength=bcount).astype(np.float64);assert (sizes>0).all()
 metadata=bundle/'metadata.json';assert metadata.stat().st_size<=536870912
 meta_hash=sha(metadata);assert meta_hash==sha(bundle/'input/metadata.json')
 if expected_metadata:assert meta_hash==expected_metadata
 with metadata.open('rb') as f:samples=list(ijson.items(f,'samples.item'))
 lookup={s['id']:s for s in samples};assert len(lookup)==len(samples)
 cov='donorID' if o['covariate']=='donor' else 'batchID';assert o['covariate'] in ['donor','batch']
 count=0
 with metadata.open('rb') as f:
  for count,cell in enumerate(ijson.items(f,'cells.item'),1):
   assert count<=n and levels[batch[count-1]]==lookup[cell['sampleID']][cov],('metadata level',count)
 assert count==n
 paths={'x':(bundle/'input/scores.bin',d),'r':(bundle/'memberships.bin',k),'a':(bundle/'assignment-scores.bin',d),'y':(bundle/'scores.bin',d)}
 bindings={str(p.relative_to(bundle)):dict(bytes=p.stat().st_size,SHA256=sha(p)) for p,_ in paths.values()}
 for name in ['metadata.json','input/metadata.json','plan.json','report.json','receipt.json']:
  p=bundle/name;bindings[name]=dict(bytes=p.stat().st_size,SHA256=sha(p))
 centers=np.asarray(report['assignmentCenters'],dtype=np.float64);assert centers.shape==(k,d) and np.isfinite(centers).all()
 observed=np.zeros((k,bcount));weighted=np.zeros((bcount,k,d));objective=0.;row_error=0.
 for start,end,v in blocks({name:paths[name] for name in ['x','r','a']},n):
  x,r,a=v['x'],v['r'],v['a'];assert (r>=0).all() and (r<=1).all()
  row_error=max(row_error,float(np.max(np.abs(r.sum(axis=1)-1))));assert row_error<=1e-10
  distance=2*(1-a@centers.T);logr=np.zeros_like(r);np.log(r,out=logr,where=r>0)
  objective+=float(np.sum(r*distance+o['temperature']*r*logr,dtype=np.float64))
  for b in range(bcount):
   selection=batch[start:end]==b
   if not selection.any():continue
   rb=r[selection];observed[:,b]+=rb.sum(axis=0);weighted[b]+=rb.T@x[selection]
 masses=observed.sum(axis=1);expected=masses[:,None]*sizes[None,:]/n
 objective+=float(o['temperature']*o['diversity']*np.sum(observed*np.log((observed+expected+1)/(2*expected+1))))
 objective*=2000/n
 assert np.isfinite(objective) and np.isclose(objective,report['objectives'][-1],atol=1e-8,rtol=1e-9),('objective',objective,report['objectives'][-1])
 effects=np.zeros((bcount,k,d));solve_residual=0.
 penalties=np.full((k,bcount),o['ridge'])
 if o['ridgeScaling']=='expectedClusterBatchMass':
  penalties=o['ridge']*expected
  assert np.allclose(penalties,np.asarray(report['ridgePenalties']),atol=1e-8,rtol=1e-9)
 else:assert o['ridgeScaling'] is None and report.get('ridgePenalties') is None
 for c in range(k):
  active=np.flatnonzero(observed[c]/sizes>1e-5)
  if len(active)<2:continue
  m=observed[c,active];lam=penalties[c,active];s=weighted[active,c,:]
  # Independent dense normal-equation solve, not the native Schur implementation.
  matrix=np.zeros((len(active)+1,len(active)+1));matrix[0,0]=m.sum();matrix[0,1:]=m;matrix[1:,0]=m
  matrix[1:,1:]=np.diag(m+lam);rhs=np.vstack([s.sum(axis=0),s]);fit=np.linalg.solve(matrix,rhs)
  residual=float(np.max(np.abs(matrix@fit-rhs)/(1+np.abs(rhs))));solve_residual=max(solve_residual,residual);assert residual<1e-9
  effects[active,c,:]=fit[1:]
 max_error=0.;scale=0.
 for start,end,v in blocks({name:paths[name] for name in ['x','r','y']},n):
  predicted=np.array(v['x'],copy=True)
  for b in range(bcount):
   selection=batch[start:end]==b
   if selection.any():predicted[selection]-=v['r'][selection]@effects[b]
  error=np.abs(predicted-v['y']);assert np.all(error<=1e-8+1e-9*np.abs(v['y'])),('correction',float(error.max()))
  max_error=max(max_error,float(error.max()));scale=max(scale,float(np.abs(v['y']).max()))
 history=np.asarray(report['objectives']);improvements=-np.diff(history)/np.maximum(np.abs(history[:-1]),1e-12)
 assert len(history)>=2 and np.isfinite(history).all() and np.allclose(improvements,report['relativeImprovements'],atol=1e-12,rtol=1e-12)
 if report['stoppingReason']=='objective-increase':assert improvements[-1]<0
 elif report['stoppingReason']=='relative-objective-tolerance':assert 0<=improvements[-1]<o['relativeTolerance']
 else:assert report['stoppingReason']=='iteration-limit' and len(improvements)==o['maximumIterations']
 assert np.isfinite(report['maximumRidgeResidual']) and 0<=report['maximumRidgeResidual']<1e-10
 assert all((bundle/name).stat().st_size==v['bytes'] and sha(bundle/name)==v['SHA256'] for name,v in bindings.items()),'Input changed during check'
 result=dict(status='passed',cells=n,components=d,clusters=k,levels=levels,maximumRowBatch=8192,maximumMembershipRowSumError=row_error,independentObjective=objective,nativeObjective=report['objectives'][-1],maximumCorrectionError=max_error,maximumCorrectedMagnitude=scale,maximumIndependentSolveResidual=solve_residual,allMatrixCoordinatesAndValuesChecked=True,bindings=bindings,numpy=np.__version__,ijson=ijson.__version__,platform=platform.platform(),seconds=time.time()-began,maximumResidentBytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if platform.system()=='Darwin' else 1024),tolerances=dict(membershipSumAbsolute=1e-10,objectiveAbsolute=1e-8,objectiveRelative=1e-9,correctionAbsolute=1e-8,correctionRelative=1e-9),scope='Every published matrix record, final objective and categorical ridge correction checked independently in bounded row arrays. Does not reproduce the stochastic optimization trajectory, establish biological preservation or qualify prospective prediction.')
 out.parent.mkdir(parents=True,exist_ok=True)
 with out.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
 print(json.dumps({key:result[key] for key in ['status','cells','maximumCorrectionError','independentObjective','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bundle',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--expected-metadata-sha');a=p.parse_args();check(a.bundle,a.out,a.expected_metadata_sha)
