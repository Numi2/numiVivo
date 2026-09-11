#!/usr/bin/env python3
"""Independently reconstruct every full Gaussian update and all matching distances."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for n in ('source','root','previous','out'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();freeze=read(a.root/'outputs-frozen.json');assert freeze['allOriginalCells'] and freeze['allReplayArraysExact']
 for n,h in freeze['files'].items():assert sha(a.root/n)==h
 results=[]
 for c in ('hagai','kang','ding'):
  r=a.root/c;report=read(r/'report.json');n,d=report['cells'],report['components'];old=a.source/c
  raw=np.fromfile(old/'pca/scores.bin',dtype=[('row','<u4'),('column','<u4'),('value','<f8')]).reshape(n,d);assert (raw['row']==np.arange(n)[:,None]).all() and (raw['column']==np.arange(d)).all();x=raw['value'].copy();scale=read(old/'mnn/report.json')['medianRowNorm'];assert abs(scale-report['scale'])<1e-12
  current=x/report['scale'];batch=np.array(read(old/'mnn/report.json')['cellLevels']);rows=[np.flatnonzero(batch==b) for b in range(int(batch.max())+1)]
  squared=np.zeros(n)
  for j in range(d):squared+=x[:,j]*x[:,j]
  norm=np.sqrt(squared);unit=np.divide(x,norm[:,None],out=np.zeros_like(x),where=norm[:,None]>0)
  matching=np.load(r/'matching.npz');near=matching['indices'];dist=matching['distances'];max_match_error=0.;max_kernel_error=0.;max_delta_error=0.;query_rows=0
  for direction in range(2):
   for first in range(0,n,512):
    end=min(n,first+512);ix=near[direction,first:end];valid=ix>=0;expected_valid=np.broadcast_to(((batch[first:end]<batch.max()) if direction else (batch[first:end]>0))[:,None],ix.shape);assert np.array_equal(valid,expected_valid)
    safe=np.maximum(ix,0);assert ((batch[safe]>batch[first:end,None]) if direction else (batch[safe]<batch[first:end,None]))[valid].all()
    ordered=np.sort(ix,axis=1);assert (np.diff(ordered,axis=1)[valid[:,1:]]>0).all()
    dd=np.zeros(ix.shape)
    for j in range(d):dd+=np.abs(unit[first:end,j,None]-unit[safe,j])
    err=np.abs(dd[valid]-dist[direction,first:end][valid]);max_match_error=max(max_match_error,float(err.max(initial=0)));assert np.allclose(dd[valid],dist[direction,first:end][valid],rtol=1e-12,atol=1e-12)
  for number,step in enumerate(report['steps']):
   if step['status']!='corrected':continue
   z=np.load(r/('step-'+str(number)+'.npz'));selected=z['selected'];match=z['anchors'];source=z['source'];bias=z['bias'];query=z['queries'];delta=z['delta']
   assert np.array_equal(selected,np.concatenate([rows[b] for b in step['target']]))
   assert np.array_equal(source,current[match[:,0]]) and np.array_equal(bias,current[match[:,1]]-source) and np.array_equal(query,current[selected])
   for first in range(0,len(query),64):
    end=min(len(query),first+64);count=end-first;denom=np.zeros(count);numerator=np.zeros((count,d))
    for anchor in range(0,len(source),256):
     ss=source[anchor:anchor+256];dd=np.zeros((count,len(ss)))
     for j in range(d):
      difference=query[first:end,j,None]-ss[None,:,j];dd+=difference*difference
     weights=np.exp(-7.5*dd);denom+=weights.sum(axis=1);numerator+=weights@bias[anchor:anchor+256]
    expected=np.divide(numerator,denom[:,None],out=np.zeros_like(numerator),where=denom[:,None]>0)
    # Original cohorts must not silently test a different underflow convention.
    assert np.all((denom==0)|(denom>=1e-280)),'Subnormal full-cohort row requires explicit independent scalar oracle'
    np.testing.assert_allclose(denom,z['totalWeights'][first:end],rtol=1e-10,atol=1e-280)
    error=float(np.max(np.abs(expected-delta[first:end])));max_delta_error=max(max_delta_error,error)
    assert np.all(np.abs(expected-delta[first:end])<=1e-10*(1+np.abs(expected)))
   current[selected]+=delta;query_rows+=len(query)
  y=np.load(r/'scores.npz')['scores'];np.testing.assert_array_equal(current*report['scale'],y)
  exact=np.fromfile(old/'mnn/anchors.bin',dtype=[('first','<u4'),('second','<u4'),('distance','<f8')]);exact_keys=exact['first'].astype(np.int64)*n+exact['second'];candidate=np.load(r/'anchors.npz');candidate_pairs=np.concatenate([candidate[k] for k in candidate.files]);keys=candidate_pairs[:,0]*n+candidate_pairs[:,1];common=np.intersect1d(keys,exact_keys);precision=len(common)/len(keys);recall=len(common)/len(exact_keys)
  pairs=[]
  for i in range(len(rows)):
   for j in range(i+1,len(rows)):
    oldkeys=exact_keys[(batch[exact['first']]==i)&(batch[exact['second']]==j)];newkeys=keys[(batch[candidate_pairs[:,0]]==i)&(batch[candidate_pairs[:,1]]==j)];same=len(np.intersect1d(oldkeys,newkeys))
    pairs.append(dict(first=i,second=j,exactAnchors=len(oldkeys),candidateAnchors=len(newkeys),intersection=same,precision=same/len(newkeys) if len(newkeys) else None,recall=same/len(oldkeys) if len(oldkeys) else None))
  previous=read(a.previous/c/'report.json')
  for name in ('matching','anchors'):assert report['arraySHA256'][name]==previous['arraySHA256'][name],('Changed matching intervention',c,name)
  assemblyOrderOriginalExact=report['order']==read(old/'mnn/report.json')['assemblyOrder']
  base=np.fromfile(old/'mnn/scores.bin',dtype=raw.dtype)['value'].reshape(n,d);error=y-base
  results.append(dict(cohort=c,cells=n,kernelQueryRows=query_rows,maximumMatchingDistanceError=max_match_error,fullGaussianDistanceTerms=sum(s['distanceTerms'] for s in report['native'] if s['stage'].startswith('kernel-')),matchingAndAnchorsPreviousExact=True,assemblyOrderOriginalExact=assemblyOrderOriginalExact,maximumKernelDeltaError=max_delta_error,finalReconstructionExact=True,exactAnchors=len(exact_keys),candidateAnchors=len(keys),anchorPrecision=precision,anchorRecall=recall,anchorGatePassed=precision>=.95 and recall>=.95,anchorPairs=pairs,coordinateRMSE=float(np.sqrt(np.mean(error**2))),relativeFrobeniusError=float(np.linalg.norm(error)/np.linalg.norm(base)),maximumRowErrorInMedianNormUnits=float(np.max(np.linalg.norm(error,axis=1))/scale)))
  print(c,precision,recall,'all full Gaussian updates checked',flush=True)
 a.out.write_text(json.dumps(dict(status='checked',results=results,outputsFreezeSHA256=sha(a.root/'outputs-frozen.json'),checkerSHA256=sha(Path(__file__)),scope='Every selected matching distance and every full-anchor Gaussian update independently checked; all original complete-cohort anchors compared. Matching is approximate; no independent biological or million-cell qualification.'),indent=2)+'\n')
if __name__=='__main__':main()
