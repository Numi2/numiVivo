#!/usr/bin/env python3
"""Audit all fixed-mean reference grid alternatives with stable NB objectives."""
import argparse,gzip,hashlib,json,os,subprocess,time
from pathlib import Path
import numpy as np
import pandas as pd

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def adjusted_difference(y,mu,x,old,new):
    # Constants in y and mu cancel; each rising factorial is an exact integer sum.
    value=0.0
    for count,mean in zip(y,mu):
        k=np.arange(int(count),dtype=float)
        value+=float(np.sum(np.log1p(k*new)-np.log1p(k*old)))
        value-=(count+1/new)*np.log1p(new*mean)-(count+1/old)*np.log1p(old*mean)
    def logdet(alpha):
        w=np.sqrt(mu/(1+alpha*mu));r=np.linalg.qr(x*w[:,None],mode='r')
        return 2*np.log(np.abs(np.diag(r))).sum()
    return float(value-.5*(logdet(new)-logdet(old)))

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ['root','stages','out','r-library']:p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);owner=Path(__file__).resolve().parent
 write(a.out/'declaration.json',{f:digest(owner/f) for f in ['PROTOCOL.md','reference.py','reference.R','Main.swift']})
 selected=[];runs=[]
 for study in ['kang','hagai']:
  first=a.stages/study/'1'
  n=pd.DataFrame(json.loads(gzip.decompress((first/'native-stages.json.gz').read_bytes()))['features']).set_index('featureID')
  r=pd.read_csv(first/'stages.tsv.gz',sep='\t',index_col='featureID')
  joined=n.join(r,rsuffix='Reference')
  top=joined[(joined.geneWiseDispersion>=1e-6)&(joined.geneWiseDispersionReference<1e-6)].sort_values(['geneWiseDispersion','featureID'],ascending=[False,True]).head(16)
  for gene in top.index:selected.append(dict(study=study,seed=1,featureID=gene,reason='sixteen-largest-native-interior-reference-boundary-disagreements'))
 for hit in json.loads((a.root/'final-results/null-hit-diagnostics.json').read_text())['hits']:
  if not any((s['study'],s['seed'],s['featureID'])==(hit['study'],hit['seed'],hit['featureID']) for s in selected):
   selected.append(dict(study=hit['study'],seed=hit['seed'],featureID=hit['featureID'],reason='original-native-Hagai-sham-call'))
 write(a.out/'selection.json',selected)
 for study in ['kang','hagai']:
  for seed in range(1,11):
   d=a.root/study/str(seed);out=a.out/study/str(seed);out.parent.mkdir(parents=True,exist_ok=True)
   started=time.monotonic()
   with (out.parent/(str(seed)+'.log')).open('w') as log:
    result=subprocess.run(['/opt/homebrew/bin/Rscript',str(owner/'reference.R'),str(d/'r-input'),str(out)],stdout=log,stderr=subprocess.STDOUT,
      env={**os.environ,'R_LIBS_USER':str(a.r_library),'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
   record=dict(study=study,seed=seed,exitCode=result.returncode,seconds=time.monotonic()-started)
   runs.append(record);write(a.out/'runs.json',runs)
   if result.returncode:continue
   r=pd.read_csv(out/'reference.tsv.gz',sep='\t',index_col='featureID',float_precision='round_trip')
   old=pd.read_csv(a.stages/study/str(seed)/'stages.tsv.gz',sep='\t',index_col='featureID',float_precision='round_trip')
   assert list(r.index)==list(old.index) and np.array_equal(r.originalDispersion,old.geneWiseDispersion)
   means=pd.read_csv(out/'means.tsv.gz',sep='\t',index_col='featureID',float_precision='round_trip');assert list(means.index)==list(r.index)
   counts=pd.read_csv(d/'r-input/counts.tsv',sep='\t',index_col='featureID').loc[r.index]
   x=pd.read_csv(d/'r-input/design.tsv',sep='\t',index_col=0).to_numpy(float)
   assert list(counts.columns)==list(means.columns)
   differences=[adjusted_difference(y,mu,x,original,grid) for y,mu,original,grid in zip(counts.to_numpy(),means.to_numpy(),r.originalDispersion,r.gridDispersion)]
   r['gridObjectiveGain']=differences;r.to_csv(out/'checked.tsv.gz',sep='\t',compression={'method':'gzip','mtime':0})
   improved=r.gridObjectiveGain>1e-4;boundary=r.originalDispersion<1e-6
   record.update(status='checked-every-eligible-reference-profile',genes=len(r),originalBoundary=int(boundary.sum()),
     gridImproved=int(improved.sum()),boundaryImproved=int((boundary&improved).sum()),
     initialFloorBoundaryImproved=int((boundary&improved&(r.initialDispersion==1e-8)).sum()),
     maximumObjectiveGain=float(r.gridObjectiveGain.max()),minimumObjectiveGain=float(r.gridObjectiveGain.min()),
     referenceReproducedExactly=True,inputHashes={f:digest(d/'r-input'/f) for f in ['counts.tsv','design.tsv','samples.tsv','input.json']})
   write(a.out/'runs.json',runs);print(json.dumps(record),flush=True)
 assert len(runs)==20 and all(r.get('status')=='checked-every-eligible-reference-profile' for r in runs)
 requests=[]
 for s in selected:
  d=a.root/s['study']/str(s['seed'])/'r-input';out=a.out/s['study']/str(s['seed'])
  y=pd.read_csv(d/'counts.tsv',sep='\t',index_col=0).loc[s['featureID']].to_numpy(int)
  x=pd.read_csv(d/'design.tsv',sep='\t',index_col=0).to_numpy(float);samples=pd.read_csv(d/'samples.tsv',sep='\t',float_precision='round_trip')
  r=pd.read_csv(out/'reference.tsv.gz',sep='\t',index_col=0).loc[s['featureID']]
  n=pd.DataFrame(json.loads(gzip.decompress((a.stages/s['study']/str(s['seed'])/'native-stages.json.gz').read_bytes()))['features']).set_index('featureID').loc[s['featureID']]
  dispersions=[float(n.geneWiseDispersion),float(r.originalDispersion),float(r.gridDispersion)]+np.geomspace(1e-8,100,41).tolist()
  requests.append(dict(id=s['study']+'/'+str(s['seed'])+'/'+s['featureID'],design=x.tolist(),offsets=np.log(samples.sizeFactor).tolist(),
    contrast=[0,1]+[0]*(x.shape[1]-2),counts=y.tolist(),dispersions=sorted(set(dispersions))))
 write(a.out/'native-input.json',requests)
 write(a.out/'reference-complete.json',dict(status='completed-all-twenty-reference-profile-audits',nativeRequests=len(requests)))

if __name__=='__main__':main()
