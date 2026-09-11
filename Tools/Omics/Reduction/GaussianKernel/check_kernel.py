#!/usr/bin/env python3
"""Independent complete native Gaussian tile checks; no cohort labels used."""
import argparse,ctypes as ct,hashlib,json,math,time
from pathlib import Path
import numpy as np
D=ct.POINTER(ct.c_double);U=ct.c_uint32;L=ct.c_uint64;CANCEL=ct.CFUNCTYPE(ct.c_int32,ct.c_void_p)
def main():
 p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 lib=ct.CDLL(str(a.library));f=lib.nvivo_omics_gaussian_bias;f.restype=ct.c_int32
 f.argtypes=[D,D,U,D,U,U,ct.c_double,L,L,L,D,L,D,L,D,L,ct.POINTER(L),CANCEL,ct.c_void_p]
 scalar=lib.nvivo_omics_gaussian_scalar;scalar.restype=ct.c_int32
 scalar.argtypes=[D,D,U,D,U,ct.c_double,L,L,D,L,D,CANCEL,ct.c_void_p]
 def run(s,b,q,sigma=15.,work=None,capacity=None,cancel_after=None):
  s,b,q=[np.ascontiguousarray(x,dtype=np.float64) for x in (s,b,q)];na,d=s.shape;nq=len(q);delta=np.full((nq,d),999.);totals=np.full(nq,999.);w=np.empty(256*(d+nq));calls=0
  @CANCEL
  def cancel(_):
   nonlocal calls;calls+=1;return int(cancel_after is not None and calls>=cancel_after)
  ptr=lambda x:x.ctypes.data_as(D)
  evaluated=L()
  st=f(ptr(s),ptr(b),na,ptr(q),nq,d,sigma,s.size,q.size,2*na*nq*d if work is None else work,ptr(delta),delta.size,ptr(totals),totals.size,ptr(w),w.size if capacity is None else capacity,ct.byref(evaluated),cancel,None)
  if st==0:
   assert na*nq*d<=evaluated.value<=2*na*nq*d
   for i in range(nq):
    numerator=np.zeros(d);scalar_total=np.zeros(1)
    scalar_status=scalar(ptr(s),ptr(b),na,ptr(q[i]),d,sigma,s.size,d,ptr(numerator),d,ptr(scalar_total),CANCEL(lambda _:0),None)
    assert scalar_status==0
    expected=numerator/scalar_total[0] if scalar_total[0]>0 else np.zeros(d)
    np.testing.assert_allclose(delta[i],expected,rtol=1e-10,atol=1e-10)
    np.testing.assert_allclose(totals[i],scalar_total[0],rtol=1e-10,atol=1e-10)
  return st,delta,totals,calls
 def reference(s,b,q,sigma):
  out=np.zeros_like(q);tot=np.zeros(len(q))
  for i,v in enumerate(q):
   for anchor,bias in zip(s,b):
    dist=0.
    for x,y in zip(v,anchor):diff=float(x)-float(y);dist+=diff*diff
    w=math.exp(-.5*sigma*dist);tot[i]+=w
    for j in range(len(v)):out[i,j]+=w*float(bias[j])
   if tot[i]>0:out[i]/=tot[i]
  return out,tot
 rng=np.random.default_rng(17);cases=[];maximum=0.
 for d in (1,2,20,64):
  for na,nq in ((1,1),(255,7),(256,32),(257,31),(513,32)):
   s=rng.normal(size=(na,d))*.08;b=rng.normal(size=s.shape);q=rng.normal(size=(nq,d))*.08
   if na>1:s[1]=s[0];b[1]=-b[0]
   st,y,t,_=run(s,b,q);assert st==0;yr,tr=reference(s,b,q,15.)
   err=float(np.max(np.abs(y-yr)));maximum=max(maximum,err);assert err<=1e-10*(1+np.abs(yr).max());np.testing.assert_allclose(t,tr,rtol=1e-10,atol=1e-10)
   st2,y2,t2,_=run(s,b,q);assert st2==0;np.testing.assert_array_equal(y,y2);np.testing.assert_array_equal(t,t2)
   cases.append(dict(dimensions=d,anchors=na,queries=nq,maximumDeltaError=err,replayExact=True))
 for tag,s,b,q,sigma in [('large-shared-offset',1e12+rng.normal(size=(257,20))*.01,rng.normal(size=(257,20)),1e12+rng.normal(size=(32,20))*.01,15.),('underflow',np.zeros((257,2)),np.ones((257,2)),np.array([[0.,0.],[1e9,0.]]),1000.),('overflow-distance',np.array([[-1e308]]),np.array([[1.]]),np.array([[1e308]]),15.)]:
  st,y,t,_=run(s,b,q,sigma);assert st==0;yr,tr=reference(s,b,q,sigma);err=float(np.max(np.abs(y-yr)));assert err<=1e-10*(1+np.abs(yr).max());np.testing.assert_allclose(t,tr,rtol=1e-10,atol=1e-10);maximum=max(maximum,err);cases.append(dict(case=tag,maximumDeltaError=err))
 for exponent in (700.,708.,710.,735.,740.,744.,745.,746.):
  source=np.zeros((257,1));bias=np.linspace(-2.,4.,257)[:,None];query=np.array([[math.sqrt(exponent/7.5)]])
  st,y,t,_=run(source,bias,query);assert st==0;yr,tr=reference(source,bias,query,15.)
  err=float(np.max(np.abs(y-yr)));assert err<=1e-10*(1+np.abs(yr).max()),('near-underflow',exponent,err,y,yr,t,tr)
  assert np.array_equal(t==0,tr==0);maximum=max(maximum,err);cases.append(dict(case='near-underflow',exponent=exponent,maximumDeltaError=err))
 st,*_=run(np.zeros((1,1)),np.ones((1,1)),np.full((1,1),10.),work=1);assert st==2
 s=np.zeros((513,2));b=np.ones_like(s);q=np.zeros((2,2));reject=['fallback-work']
 for tag,kwargs,expected in [('work',dict(work=2051),2),('workspace',dict(capacity=1023),1),('cancel-before',dict(cancel_after=1),5),('cancel-later',dict(cancel_after=3),5),('nonfinite-sigma',dict(sigma=math.nan),1)]:
  st,*_=run(s,b,q,**kwargs);assert st==expected,(tag,st);reject.append(tag)
 for tag,which in [('source',0),('bias',1),('query',2)]:
  values=[s.copy(),b.copy(),q.copy()];values[which][-1,-1]=math.nan;st,*_=run(*values);assert st==4;reject.append('nonfinite-'+tag)
 record=dict(status='passed',cases=cases,rejections=reject,maximumDeltaError=maximum,librarySHA256=hashlib.sha256(a.library.read_bytes()).hexdigest(),checkerSHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),numpy=np.__version__,createdUnix=time.time(),scope='Every query and coordinate of independent numerical cases; not biological qualification')
 a.out.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k!='cases'}))
if __name__=='__main__':main()
