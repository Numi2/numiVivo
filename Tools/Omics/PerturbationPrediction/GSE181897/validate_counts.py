from pathlib import Path
import sys,json,gzip,hashlib,time,resource
r=Path(__file__).parent;sys.path.insert(0,str(r/'deps'))
import indexed_gzip,h5py,numpy as np
start=time.time();var=json.loads(gzip.decompress((r/'var-metadata.json.gz').read_bytes()));rna=np.array(var['genome'])=='GRCh38';adt=np.array(var['genome'])=='BD99AbSeq';assert np.all(rna|adt)
with indexed_gzip.IndexedGzipFile(str(r/'GSE181897_concat.4.raw.h5ad.gz'),index_file=str(r/'source.gzidx'),readbuf_size=1024*1024) as gz,h5py.File(gz,'r',rdcc_nbytes=64*1024*1024) as f:
 x=f['X'];n,m=map(int,x.attrs['shape']);ptr=x['indptr'][:].astype(np.int64);assert len(ptr)==n+1 and ptr[0]==0 and ptr[-1]==len(x['data'])==len(x['indices']) and np.all(np.diff(ptr)>=0)
 rowsums=np.zeros((n,2),dtype=np.uint64);maxvalue=0;storedzero=0;unsortedrows=0;duplicates=0
 for begin in range(0,n,1024):
  end=min(n,begin+1024);lo,hi=map(int,(ptr[begin],ptr[end]));v=x['data'][lo:hi];j=x['indices'][lo:hi];assert np.all(np.isfinite(v)) and np.all(v>=0) and np.all(v==np.floor(v)) and np.all((j>=0)&(j<m));maxvalue=max(maxvalue,float(v.max(initial=0)));assert maxvalue<=2**24
  storedzero+=int(np.count_nonzero(v==0));row=np.repeat(np.arange(end-begin),np.diff(ptr[begin:end+1]));assert len(row)==len(v)
  for c,mask in enumerate([rna[j],adt[j]]):
   sums=np.bincount(row[mask],weights=v[mask].astype(np.float64),minlength=end-begin);assert np.all(sums<2**53) and np.all(sums==np.floor(sums));rowsums[begin:end,c]=sums.astype(np.uint64)
  for a,b in zip(ptr[begin:end]-lo,ptr[begin+1:end+1]-lo):
   z=j[a:b];unsortedrows+=int(np.any(z[1:]<=z[:-1]));duplicates+=len(z)-len(np.unique(z))
  if begin%16384==0:print(json.dumps(dict(rows=end,records=hi,elapsedSeconds=time.time()-start)),flush=True)
 assert duplicates==0
 np.savez_compressed(r/'source-cell-library-counts.npz',RNA=rowsums[:,0],antibody=rowsums[:,1])
 # Author totals may refer to a prior feature-filtering stage; report every discrepancy.
 total_checks={}
 for label,col,c in [('RNA','mrna_n_counts',0),('antibody','adts_n_counts',1)]:
  author=f['obs/'+col][:].astype(np.float64);diff=rowsums[:,c].astype(np.float64)-author;total_checks[label]=dict(exactCells=int(np.count_nonzero(diff==0)),differentCells=int(np.count_nonzero(diff!=0)),maximumAbsoluteDifference=float(np.max(np.abs(diff))))
 result=dict(status='passed-count-domain-and-CSR',cells=n,features=m,storedRecords=int(ptr[-1]),RNAFeatures=int(rna.sum()),antibodyFeatures=int(adt.sum()),RNAUMIs=int(rowsums[:,0].sum()),antibodyUMIs=int(rowsums[:,1].sum()),maximumStoredCount=maxvalue,storedZeros=storedzero,unsortedRows=unsortedrows,duplicateFeaturesWithinRows=duplicates,authorLibraryComparisons=total_checks,elapsedSeconds=time.time()-start,peakRSSBytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,sourceSHA256=json.loads((r/'download.json').read_text())['compressedSHA256'],scriptSHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),aggregationStarted=False,predictionScoringStarted=False)
 (r/'source-count-validation.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result),flush=True)
