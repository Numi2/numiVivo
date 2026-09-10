#!/usr/bin/env python3
"""Test one protected-stratum regression with frozen native cluster memberships."""
import argparse,time
from pathlib import Path
import numpy as np
from check_integration_response import sha,RECORD
from prepare_annotation_retention import read,write


def records(path,n,d,rows=8192):
 assert path.stat().st_size==n*d*16 and 1<=d<=100
 with path.open('rb') as f:
  for first in range(0,n,rows):
   last=min(n,first+rows);raw=f.read((last-first)*d*16)
   assert len(raw)==(last-first)*d*16
   a=np.frombuffer(raw,dtype=RECORD).reshape(last-first,d)
   assert np.array_equal(a['row'],np.broadcast_to(np.arange(first,last)[:,None],a.shape))
   assert np.array_equal(a['column'],np.broadcast_to(np.arange(d),a.shape))
   x=a['value'].copy();assert np.isfinite(x).all();yield first,last,x
  assert not f.read(1)


def fit(mass,total):
 # Positive stratum intercepts are eliminated; absent strata never get invented.
 m=mass.sum(axis=1);valid=m>0;mass=mass[valid];total=total[valid];m=m[valid]
 group_mean=total.sum(axis=1)/m[:,None]
 gram=np.diag(mass.sum(axis=0))-np.einsum('sb,sc,s->bc',mass,mass,1/m)+np.eye(mass.shape[1])
 rhs=total.sum(axis=0)-mass.T@group_mean
 beta=np.linalg.solve(gram,rhs)
 intercept=group_mean-(mass@beta)/m[:,None]
 # Independent augmented least squares over weighted stratum/donor means.
 s,b=np.nonzero(mass>0);w=np.sqrt(mass[s,b]);design=np.zeros((len(s),len(m)+mass.shape[1]))
 design[np.arange(len(s)),s]=w;design[np.arange(len(s)),len(m)+b]=w
 response=total[s,b]/mass[s,b,None]*w[:,None]
 penalty=np.column_stack([np.zeros((mass.shape[1],len(m))),np.eye(mass.shape[1])])
 reference=np.linalg.lstsq(np.vstack([design,penalty]),np.vstack([response,np.zeros((mass.shape[1],total.shape[-1]))]),rcond=None)[0]
 err=float(np.max(np.abs(beta-reference[len(m):])))
 assert np.allclose(beta,reference[len(m):],rtol=1e-9,atol=1e-8)
 assert np.allclose(intercept,reference[:len(m)],rtol=1e-9,atol=1e-8)
 residual=float(np.max(np.abs(gram@beta-rhs)/(1+np.abs(rhs))))
 assert residual<1e-9
 return beta,reference[len(m):],dict(stratumIndices=np.flatnonzero(valid).tolist(),intercepts=intercept.tolist(),maximumSVDDiscrepancy=err,maximumNormalResidual=residual)


def append(stream,first,x):
 a=np.empty(x.shape,dtype=RECORD);a['row']=np.arange(first,first+len(x))[:,None];a['column']=np.arange(x.shape[1]);a['value']=x
 stream.write(a.tobytes())


def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ('study','out'):p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args();s=a.study;r=a.out;n=1612594;d=20;k=100
 assert not (r/'execution-freeze.json').exists()
 prior=s/'annotation-retention';freeze=read(prior/'freeze.json');ledger=read(prior/'ledger.json')
 for name,h in freeze['files'].items():assert sha(prior/name)==h
 rows=np.load(prior/'rows.npz');sc=rows['stratum'];dc=rows['donor'];S=len(ledger['strata']);B=len(ledger['donors'])
 assert len(sc)==len(dc)==n and S==23 and B==5
 x_path=s/'integration-full/pca/scores.bin';native=s/'integration-full/native';r_path=native/'memberships.bin'
 receipt=read(native/'receipt.json');plan=read(native/'plan.json')
 assert plan['integration']['ridge']==1 and not plan['integration'].get('ridgeScaling')
 assert read(native/'report.json')['levels']==ledger['donors']
 hashes={'original':freeze['matrices']['baseline']['SHA256'],'native':freeze['matrices']['native']['SHA256'],'memberships':bytes(receipt['memberships']['bytes']).hex()}
 inputs={'original':x_path,'native':native/'scores.bin','memberships':r_path}
 for name,path in inputs.items():assert sha(path)==hashes[name]
 write(r/'execution-freeze.json',dict(createdUnix=time.time(),protocolSHA256=sha(r/'protocol.md'),scriptSHA256=sha(Path(__file__)),annotationFreezeSHA256=sha(prior/'freeze.json'),inputs={name:dict(path=str(p.relative_to(s)),SHA256=hashes[name]) for name,p in inputs.items()},strata=ledger['strata'],donors=ledger['donors'],ridge=1,activeCriterion='cluster donor mass / original donor size > 1e-5',metricsRead=False))
 masses=np.zeros((k,S,B));totals=np.zeros((k,S,B,d));max_row_error=0.
 for (first,last,x),(i,j,weights) in zip(records(x_path,n,d),records(r_path,n,k)):
  assert (first,last)==(i,j) and weights.min()>=0
  max_row_error=max(max_row_error,float(np.max(np.abs(weights.sum(axis=1)-1))))
  codes=sc[first:last].astype(int)*B+dc[first:last]
  for code in np.unique(codes):
   ss,bb=divmod(int(code),B);mask=codes==code
   masses[:,ss,bb]+=weights[mask].sum(axis=0);totals[:,ss,bb]+=weights[mask].T@x[mask]
 assert max_row_error<1e-10
 sizes=np.bincount(dc,minlength=B);legacy=np.zeros((k,B,d));protected=legacy.copy();oracle=legacy.copy();fits=[]
 for c in range(k):
  donor_mass=masses[c].sum(axis=0);active=np.flatnonzero(donor_mass/sizes>1e-5)
  if len(active)<2:
   fits.append(dict(cluster=c,activeDonors=active.tolist(),status='skipped-fewer-than-two-active-donors'));continue
  mass=donor_mass[active];sums=totals[c].sum(axis=0)[active]
  center=(sums/(mass+1)[:,None]).sum(axis=0)/(mass/(mass+1)).sum()
  legacy[c,active]=(sums-mass[:,None]*center)/(mass+1)[:,None]
  beta,ref,check=fit(masses[c][:,active],totals[c][:,active])
  protected[c,active]=beta;oracle[c,active]=ref
  fits.append(dict(cluster=c,activeDonors=active.tolist(),status='fitted',**check))
 write(r/'models.json',dict(fits=fits,legacyEffects=legacy.tolist(),protectedEffects=protected.tolist(),svdEffects=oracle.tolist(),masses=masses.tolist(),totals=totals.tolist(),maximumMembershipRowSumError=max_row_error))
 original_error=0.;oracle_error=0.
 with (r/'legacy-reconstructed.bin').open('xb') as f,(r/'protected.bin').open('xb') as g:
  for (first,last,x),(i,j,w),(u,v,old) in zip(records(x_path,n,d),records(r_path,n,k),records(native/'scores.bin',n,d)):
   assert (first,last)==(i,j)==(u,v)
   legacy_y=x.copy();protected_y=x.copy();ref=x.copy()
   for b in range(B):
    mask=dc[first:last]==b
    legacy_y[mask]-=w[mask]@legacy[:,b];protected_y[mask]-=w[mask]@protected[:,b];ref[mask]-=w[mask]@oracle[:,b]
   original_error=max(original_error,float(np.max(np.abs(legacy_y-old))))
   oracle_error=max(oracle_error,float(np.max(np.abs(protected_y-ref))))
   assert np.allclose(legacy_y,old,rtol=1e-10,atol=1e-8)
   assert np.allclose(protected_y,ref,rtol=1e-9,atol=1e-8)
   append(f,first,legacy_y);append(g,first,protected_y)
 for name,path in inputs.items():assert sha(path)==hashes[name]
 write(r/'output-freeze.json',dict(status='passed',cells=n,clusters=k,maximumNativeReconstructionError=original_error,maximumProtectedSVDCoordinateError=oracle_error,files={name:sha(r/name) for name in ('models.json','legacy-reconstructed.bin','protected.bin')},executionFreezeSHA256=sha(r/'execution-freeze.json'),metricsRead=False,scope='Fixed-membership conditional regression only; no iterative native fit or biological qualification'))
 print(dict(status='passed',cells=n,nativeError=original_error,protectedSVDError=oracle_error))
if __name__=='__main__':main()
