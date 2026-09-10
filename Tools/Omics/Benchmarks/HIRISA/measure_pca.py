#!/usr/bin/env python3
"""Measure complete-source PCA moments and resource needs with bounded sparse chunks."""
import argparse,hashlib,importlib.metadata,inspect,json,time
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scanpy.preprocessing import _highly_variable_genes as hvg

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root;out=root/'pca'
 assert not (out/'measurement.json').exists();start=time.time();source=root/'hirisa.h5ad'
 source_audit=json.loads((root/'native-independent-verification.json').read_text())
 expected_source=source_audit['sourceSHA256'];assert source_audit['status']=='passed' and len(expected_source)==64
 assert sha(source)==expected_source
 plan=json.loads((out/'baseline-plan.json').read_text());options=plan['reduction'];settings=options['pca']
 with h5py.File(source,'r') as f:
  x=f['X'];n,m=map(int,x.attrs['shape']);assert (n,m)==(1612594,18082)
  genes=read_elem(f['var']).index.to_numpy(dtype=str);quality=f['obs/n_umis']
  seen=np.zeros(m,dtype=np.int64);mean=np.zeros(m);m2=np.zeros(m);logmean=np.zeros(m);logm2=np.zeros(m)
  total=0;umis=0;cells=0
  for first in range(0,n,1024):
   last=min(n,first+1024);ip=x['indptr'][first:last+1];left,right=int(ip[0]),int(ip[-1]);ip=ip-ip[0]
   col=x['indices'][left:right];raw=x['data'][left:right];assert np.all(raw>0)
   rows=np.repeat(np.arange(last-first),np.diff(ip));totals=np.bincount(rows,weights=raw,minlength=last-first)
   assert np.array_equal(totals,quality[first:last]);assert np.all(totals[rows]>0)
   values=np.log1p(raw.astype(np.float64)/totals[rows]*options['normalizationTarget']);linear=np.expm1(values)
   counts=np.bincount(col,minlength=m).astype(np.int64);positive=counts>0;combined=seen+counts
   def merge(v,old_mean,old_m2):
    chunk_mean=np.divide(np.bincount(col,weights=v,minlength=m),counts,out=np.zeros(m),where=positive)
    chunk_m2=np.bincount(col,weights=(v-chunk_mean[col])**2,minlength=m)
    delta=chunk_mean-old_mean
    old_m2[positive]+=chunk_m2[positive]+delta[positive]**2*seen[positive]*counts[positive]/combined[positive]
    old_mean[positive]+=delta[positive]*counts[positive]/combined[positive]
   merge(linear,mean,m2);merge(values,logmean,logm2);seen=combined
   total+=len(raw);umis+=int(totals.sum());cells+=last-first
   if last%131072==0 or last==n:print('cells',last,'entries',total,flush=True)
 assert total==3845991249 and umis==7700096227 and cells==n
 variance=(m2+mean**2*seen*(n-seen)/n)/(n-1);mean*=seen/n
 logvariance=(logm2+logmean**2*seen*(n-seen)/n)/(n-1);logmean*=seen/n
 assert np.isfinite(variance).all() and (variance>=0).all()
 protected=mean.copy();protected[protected==0]=1e-12
 dispersion=variance/protected;dispersion[dispersion==0]=np.nan;dispersion=np.log(dispersion)
 frame=pd.DataFrame({'means':np.log1p(protected),'dispersions':dispersion})
 # Use the installed, pinned Scanpy binning/dispersion/ranking implementation on
 # the complete measured feature moments. No synthetic count matrix is supplied.
 frame['mean_bin']=hvg._get_mean_bins(frame['means'],'seurat',settings['meanBins'])
 stats=hvg._get_disp_stats(frame,'seurat');normalized=(frame['dispersions']-stats['avg'])/stats['dev']
 cutoff=hvg._nth_highest(normalized.to_numpy().copy(),settings['highlyVariableFeatures'])
 selected=np.flatnonzero(np.nan_to_num(normalized.to_numpy(),nan=-np.inf)>=cutoff)
 entries=int(seen[selected].sum());passes=2*min(settings['maximumBasis'],len(selected))+3*settings['components']
 np.savez_compressed(out/'feature-moments.npz',featureIDs=genes,seen=seen,meanNormalized=mean,varianceNormalized=variance,
   meanLogNormalized=logmean,varianceLogNormalized=logvariance,meanBins=frame['mean_bin'].cat.codes.to_numpy(),
   logMeanForBinning=frame['means'].to_numpy(),logDispersion=dispersion,normalizedDispersion=normalized.to_numpy(),selectedFeatureIndices=selected)
 result=dict(schemaVersion=1,status='measured',sourceSHA256=sha(source),planSHA256=sha(out/'baseline-plan.json'),
  cells=n,features=m,canonicalNonzeros=total,UMIs=umis,selectedFeatures=len(selected),selectedEntries=entries,
  currentCOOCacheBytes=entries*16,plannedEntryVisits=entries*passes,currentCacheCap=2000000000,currentWorkCap=20000000000,
  currentPCAMatrixRowCap=1000000,cacheFitsCurrentCap=entries*16<=2000000000,workFitsCurrentCap=entries*passes<=20000000000,
  rowsFitCurrentWriter=n<=1000000,featureMomentsSHA256=sha(out/'feature-moments.npz'),measurementScriptSHA256=sha(Path(__file__)),
  scanpyHVGSourceSHA256=sha(Path(inspect.getsourcefile(hvg))),seconds=time.time()-start,
  versions={name:importlib.metadata.version(name) for name in ('scanpy','anndata','numpy','scipy','h5py','pandas')},
  scope='Independent complete-source chunked moments, original per-cell UMI checks and pinned Scanpy HVG selection. Not native PCA, a complete Scanpy pipeline, subspace agreement, biological or performance qualification.')
 (out/'measurement.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
