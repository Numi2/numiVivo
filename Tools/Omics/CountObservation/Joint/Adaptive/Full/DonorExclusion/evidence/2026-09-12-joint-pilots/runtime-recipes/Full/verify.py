"""Stream full-gene independent cell likelihood and support-bound verification."""
from pathlib import Path
from fractions import Fraction
import collections,fcntl,gzip,hashlib,json,os,sys,time
import h5py,numpy as np
from scipy import optimize,special
root=Path(sys.argv[1]);origin=sys.argv[2];folder=Path(os.environ.get('NUMIVIVO_FULL_OUTPUT',str(root/origin)))
lock=(folder/'verifier.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(4<<20):h.update(b)
 return h.hexdigest()
manifest=read(root/(origin+'-manifest.json.gz'));prepared=read(root/(origin+'-prepare-state.json'));cache=root/(origin+'-cells.h5');assert sha(cache)==manifest['cacheSHA256']
state={'status':'running','origin':origin,'pid':os.getpid(),'startedUnix':time.time(),'genesVerified':0,'numericalValuesCompared':0,'maximumScaledError':0.,'leafBoundsChecked':0,'errors':[],'fitStates':{}}
def save():
 temp=folder/'verification-state.tmp';temp.write_text(json.dumps(state,indent=2)+'\n');temp.replace(folder/'verification-state.json')
save();gene_errors=[];comparisons=0;max_error=0.
def close(actual,expected,label,tol=2e-7):
 global comparisons,max_error
 a=np.asarray(actual,float);b=np.asarray(expected,float);assert a.shape==b.shape,label
 error=float(np.max(np.abs(a-b)/np.maximum(1,np.abs(b)))) if a.size else 0.;comparisons+=int(a.size);max_error=max(max_error,error)
 if not np.isfinite(error) or error>tol:gene_errors.append({'where':label,'maximumScaledError':error})
def likelihood(y,library,phi):
 e=np.asarray(library,float)/1e6;y=np.asarray(y,float);total=float(y.sum());exposure=float(e.sum())
 if total==0:mode=0.
 elif phi==0:mode=total/exposure
 else:
  def score(log_rate):
   q=special.expit(np.log(phi*e)+log_rate)
   return float(np.sum(y*(1-q)-q/phi))
  center=np.log(total/exposure);mode=float(np.exp(optimize.brentq(score,center-50,center+50,xtol=1e-13)))
 cache={}
 def ell(rate):
  rate=float(rate)
  if rate in cache:return cache[rate]
  if rate==0:value=0. if total==0 else -np.inf
  elif total==0:value=-exposure*rate if phi==0 else float(-np.log1p(phi*e*rate).sum()/phi)
  elif phi==0:value=total*np.log(rate/mode)-exposure*(rate-mode)
  else:value=total*np.log(rate/mode)-float(np.dot(y+1/phi,np.log1p(phi*e*(rate-mode)/(1+phi*e*mode))))
  cache[rate]=value;return value
 def slope(rate,scale):
  if rate==0:return -scale*exposure if total==0 else np.inf
  if phi==0:score=total-exposure*rate
  else:
   a=phi*e*rate;score=float(np.sum(y/(1+a)-(a/(1+a))/phi))
  return score*(1+scale/rate)
 return mode,ell,slope
def moments(cr,tr,w):
 c=np.repeat(cr,len(tr));t=np.tile(tr,len(cr));mc=float(w@c);mt=float(w@t);lc=float(w@np.log1p(c));lt=float(w@np.log1p(t))
 return {'controlMeanCPM':mc,'treatedMeanCPM':mt,'controlVarianceCPM2':float(w@((c-mc)**2)),'treatedVarianceCPM2':float(w@((t-mt)**2)),
  'covarianceCPM2':float(w@((c-mc)*(t-mt))),'responseMeanCPM':mt-mc,'responseVarianceCPM2':float(w@((t-c-(mt-mc))**2)),
  'controlMeanLog1pCPM':lc,'treatedMeanLog1pCPM':lt,'log1pResponseMean':lt-lc,'log1pResponseVariance':float(w@((np.log1p(t)-np.log1p(c)-(lt-lc))**2))}
def sampling(m,phi,libraries):
 e=np.array(libraries,float)/1e6;mean=float(e.sum()*m['treatedMeanCPM']);over=float(phi*(e@e)*(m['treatedVarianceCPM2']+m['treatedMeanCPM']**2));latent=float(e.sum()**2*m['treatedVarianceCPM2'])
 return dict(plannedCells=len(libraries),plannedLibraryCounts=sum(libraries),meanGeneCounts=mean,conditionalPoissonVariance=mean,conditionalCellOverdispersionVariance=over,latentRateVariance=latent,totalGeneCountVariance=mean+over+latent)
statuses=collections.Counter();receipts=sorted(folder.glob('*-receipt.json'));finished=[]
with h5py.File(cache,'r') as h:
 indices=[g['index'] for g in manifest['groups']];assert len(set(indices))==len(indices)
 if 'excludedDonorID' in manifest:assert all(g['donorID']!=manifest['excludedDonorID'] for g in manifest['groups'])
 for group in manifest['groups']:
  cached=h[str(group['index'])];assert cached.attrs['donorID']==group['donorID'] and cached.attrs['conditionID']==group['conditionID']
 library=[h[str(i)+'/libraryCounts'][:] for i in indices]
 pointers=[h[str(i)+'/indptr'][:] for i in indices]
 depth_maps=[np.unique(l,return_inverse=True,return_counts=True) for l in library]
 try:
  for receipt_path in receipts:
   receipt=read(receipt_path);tag=receipt_path.name.removesuffix('-receipt.json');input_path=folder/(tag+'-input.jsonl.gz');output_path=folder/(tag+'-output.jsonl.gz');result_path=folder/(tag+'-verification.json.gz')
   assert receipt['status']=='completed' and sha(input_path)==receipt['compressedInputSHA256'] and sha(output_path)==receipt['compressedOutputSHA256']
   if result_path.exists():
    result=read(result_path);assert result['nativeCompressedSHA256']==receipt['compressedOutputSHA256'] and result['status']=='passed-independent-numerics'
   else:
    records=[]
    with gzip.open(input_path,'rt') as fi,gzip.open(output_path,'rt') as fo:
     header_line=next(fi);header=json.loads(header_line);assert header['manifestSHA256']==prepared['manifestSHA256'] and header['cacheSHA256']==manifest['cacheSHA256']
     for position,line in enumerate(fo):
      inp_line=next(fi);inp=json.loads(inp_line);row=json.loads(line);j=receipt['first']+position;gene=manifest['genes'][j];feature=gene['featureID'];assert row['featureIndex']==inp['featureIndex']==j and row['featureID']==inp['featureID']==feature
      gene_errors=[];comparisons=0;max_error=0.;record={'featureIndex':j,'featureID':feature,'nativeStatus':row['status'],'leafBoundsChecked':0}
      cp=gene['controlCellDispersion'];tp=gene['treatedCellDispersion']
      if cp is None or tp is None:
       assert row['status']=='unavailableCellDispersion' and 'model' not in row and inp['groups']==[]
      else:
       # Reconstruct every source-cell vector; zeros are retained implicitly in
       # the CSC cache and explicitly in these one-gene reference vectors.
       likes={}
       for gi,group in enumerate(manifest['groups']):
        lo,hi=map(int,pointers[gi][j:j+2]);g=h[str(group['index'])];y=np.zeros(len(library[gi]),np.uint32);y[g['indices'][lo:hi]]=g['data'][lo:hi]
        d,inverse,n=depth_maps[gi];counts=np.bincount(inverse,weights=y).astype(np.uint64);ix=np.flatnonzero(counts)
        assert inp['groups'][gi]=={'bins':ix.tolist(),'counts':counts[ix].tolist()}
        assert header['groups'][gi]=={'donorID':group['donorID'],'conditionID':group['conditionID'],'libraryCounts':d.tolist(),'cellsPerLibrary':n.tolist()}
        phi=cp if group['conditionID']=='control' else tp;likes[(group['donorID'],group['conditionID'])]=likelihood(y,library[gi],phi)
       adaptive=row.get('model')
       if adaptive is None:
        assert row['status']=='fitError';record['nativeError']=row.get('detail')
       else:
        model=adaptive['model'];assert row['status']==adaptive['status'];donors=sorted({g['donorID'] for g in manifest['groups']});assert model['trainingDonorIDs']==donors
        cs=[likes[(d,'control')] for d in donors];ts=[likes[(d,'IFNB')] for d in donors];n=len(donors)
        assert model['featureID']==feature and model['controlCellDispersion']==cp and model['treatedCellDispersion']==tp
        expected_source=hashlib.sha256(header_line.rstrip('\n').encode()+inp_line.rstrip('\n').encode()).hexdigest();assert bytes(model['trainingSource']['bytes']).hex()==expected_source
        cr=np.array(model['controlRatesCPM']);tr=np.array(model['treatedRatesCPM']);w=np.array(model['probabilities']);assert np.all(w>=0) and abs(w.sum()-1)<1e-10
        a=np.exp(np.array([[l[1](r) for r in cr] for l in cs])[:,:,None]+np.array([[l[1](r) for r in tr] for l in ts])[:,None,:]).reshape(n,-1);fitted=a@w
        ll=float(np.log(fitted).sum());gap=max(0,float(np.max(a.T@(1/fitted)))/n-1)
        close(model['fittedDonorRelativeLikelihoods'],fitted,feature+'/likelihoods');close(model['relativeLogLikelihood'],ll,feature+'/objective');close(model['meanLogLikelihoodGap'],gap,feature+'/finiteGap')
        for k,v in moments(cr,tr,w).items():close(model['moments'][k],v,feature+'/'+k)
        steps=adaptive['steps'];assert all(b['relativeLogLikelihood']>=a['relativeLogLikelihood']-1e-9 for a,b in zip(steps,steps[1:]))
        cert=adaptive.get('certificate');record.update(finiteGridGap=gap,relativeLogLikelihood=ll,supportAdditions=len(steps)-1)
        if cert is not None:
         cscale=cert['controlScaleCPM'];tscale=cert['treatedScaleCPM'];close(cscale,1/(1+cp*max(max(library[i]) for i,g in enumerate(manifest['groups']) if g['conditionID']=='control')/1e6),feature+'/controlScale');close(tscale,1/(1+tp*max(max(library[i]) for i,g in enumerate(manifest['groups']) if g['conditionID']=='IFNB')/1e6),feature+'/treatedScale')
         cm=np.array(model['controlDonorMLERatesCPM']);tm=np.array(model['treatedDonorMLERatesCPM']);close(cm,[l[0] for l in cs],feature+'/controlModes');close(tm,[l[0] for l in ts],feature+'/treatedModes')
         bounds=[float(np.log1p(cm.min()/cscale)),float(np.log1p(cm.max()/cscale)),float(np.log1p(tm.min()/tscale)),float(np.log1p(tm.max()/tscale))];caches=[{},{}];inv=1/(fitted*n)
         def axis(x,side):
          cache=caches[side]
          if x not in cache:
           ls=cs if side==0 else ts;scale=cscale if side==0 else tscale;r=scale*np.expm1(x);cache[x]=(r,np.array([v[1](r) for v in ls]),np.array([v[2](r,scale) for v in ls]))
          return cache[x]
         paths=sorted(b['path'] for b in cert['boxes']);assert len(set(paths))==len(paths) and sum((Fraction(1,2**len(p)) for p in paths),Fraction(0))==1 and all(not b.startswith(a) for a,b in zip(paths,paths[1:]));upper=0.
         for box in cert['boxes']:
          b=bounds.copy()
          for bit in box['path']:
           assert bit in '01';cl,ch,tl,th=b;split=(bounds[1]>bounds[0]) and ((bounds[3]==bounds[2]) or (ch-cl)/(bounds[1]-bounds[0]) >= (th-tl)/(bounds[3]-bounds[2]));i=0 if split else 2;mid=(b[i]+b[i+1])/2;b[i if bit=='1' else i+1]=mid
          actual=[box[k] for k in ['controlLow','controlHigh','treatedLow','treatedHigh']];close(actual,b,feature+'/partition/'+box['path'],tol=1e-12)
          x0,x1,y0,y1=actual;x=(x0+x1)/2;y=(y0+y1)/2;ar0,a0,_=axis(x0,0);ar1,a1,_=axis(x1,0);br0,b0,_=axis(y0,1);br1,b1,_=axis(y1,1);_,am,sa=axis(x,0);_,bm,sb=axis(y,1)
          individual=float(np.sum(np.exp(np.where((cm>=ar0)&(cm<=ar1),0,np.maximum(a0,a1))+np.where((tm>=br0)&(tm<=br1),0,np.maximum(b0,b1)))*inv));tangent=0.
          for xx in [x0,x1]:
           for yy in [y0,y1]:
            dx=0 if xx==x else sa*(xx-x);dy=0 if yy==y else sb*(yy-y)
            with np.errstate(over='ignore',invalid='ignore'):value=float(np.sum(np.exp(am+bm+dx+dy)*inv))
            tangent=max(tangent,value) if not np.isnan(value) else np.inf
          raw=min(individual,tangent);bound=raw+1e-9*(1+raw);close(box['meanDirectionalUpperBound'],bound,feature+'/boxBound/'+box['path']);upper=max(upper,bound)
          sample=float(np.sum(np.exp(am+bm)*inv));assert sample<=bound+1e-8*max(1,bound)
         witness=float(sum(np.exp(cs[d][1](cert['bestControlRateCPM'])+ts[d][1](cert['bestTreatedRateCPM']))*inv[d] for d in range(n)));close(cert['maximumMeanDirectionalLowerBound'],witness,feature+'/witness');close(cert['maximumMeanDirectionalUpperBound'],upper,feature+'/upper');record.update(leafBoundsChecked=len(paths),independentUpper=upper,independentWitness=witness)
         if adaptive['status']=='boundedContinuousLikelihood':assert upper<=1+adaptive['plan']['meanLogLikelihoodGapTolerance']+2e-8
      record.update(numericalValuesCompared=comparisons,maximumScaledError=max_error,errors=gene_errors);records.append(record)
     assert len(records)==receipt['genes'] and fi.readline()==''
    result={'status':'passed-independent-numerics' if not any(r['errors'] for r in records) else 'failed-independent-numerics','nativeCompressedSHA256':receipt['compressedOutputSHA256'],'records':records}
    with result_path.open('xb') as raw,gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as gz:gz.write(json.dumps(result,sort_keys=True,separators=(',',':'),allow_nan=False).encode())
   for r in result['records']:
    state['genesVerified']+=1;state['numericalValuesCompared']+=r['numericalValuesCompared'];state['maximumScaledError']=max(state['maximumScaledError'],r['maximumScaledError']);state['leafBoundsChecked']+=r['leafBoundsChecked'];state['errors'].extend(r['errors']);statuses[r['nativeStatus']]+=1
   state['fitStates']=dict(statuses);state['lastShard']=tag;save();print(json.dumps({k:v for k,v in state.items() if k!='errors'}),flush=True)
  state['status']='verified-completed-shards' if not state['errors'] else 'failed-independent-numerics'
 except BaseException as e:state.update(status='failed',error=repr(e));raise
 finally:state.update(finishedUnix=time.time(),seconds=time.time()-state['startedUnix']);save()
sys.exit(0 if state['status']=='verified-completed-shards' and not state['errors'] else 1)
