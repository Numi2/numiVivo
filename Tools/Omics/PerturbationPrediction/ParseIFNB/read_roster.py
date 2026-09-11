from remote_io import *
import numpy as np
f=Remote()
try:
 with h5py.File(f,'r') as h:
  result={};codes={};categories={}
  for key in ['cytokine','donor','sample','treatment']:
   g=h['obs/'+key];categories[key]=g['categories'].asstr()[:].tolist();codes[key]=g['codes'][:];print(key,categories[key] if key!='sample' else (len(categories[key]),categories[key][:3]),flush=True)
   assert len(codes[key])==9697974 and codes[key].min()>=0 and codes[key].max()<len(categories[key])
  relevant=[categories['cytokine'].index(x) for x in ['PBS','IFN-beta']];selected=np.flatnonzero(np.isin(codes['cytokine'],relevant));assert len(selected)>0
  counts={}
  for d,donor in enumerate(categories['donor']):
   counts[donor]={}
   for c,condition in zip(relevant,['PBS','IFN-beta']):counts[donor][condition]=int(np.count_nonzero((codes['cytokine']==c)&(codes['donor']==d)))
  runs=np.count_nonzero(np.diff(selected)!=1)+1;result=dict(sourceURL=URL,sourceETag=ETAG,sourceBytes=SIZE,totalCells=9697974,categories=categories,cohortCells=len(selected),contiguousRuns=int(runs),counts=counts,sourceFeatures=h['var/_index'].asstr()[:].tolist(),selectedFirstRows=selected[:20].tolist())
  (ROOT/'sources/roster.json').write_text(json.dumps(result,indent=2)+'\n');np.savez_compressed(ROOT/'sources/roster-codes.npz',**codes,selected=selected);print(json.dumps({k:v for k,v in result.items() if k not in ['categories','sourceFeatures']}),flush=True)
finally:(ROOT/'sources/roster-ranges.json').write_text(json.dumps(dict(url=URL,etag=ETAG,size=SIZE,requests=f.requests,bytesRead=sum(r['stop']-r['start']+1 for r in f.requests)),indent=2)+'\n')
