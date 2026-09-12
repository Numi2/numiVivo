from pathlib import Path
import json,hashlib,time,inspect
import scanpy as sc,anndata,numpy as np,scipy.sparse as sp
r=Path('/Users/n/numivivo-hvg-scverse-20260912')
(r/'protocol.json').write_text(json.dumps({'cohorts':['kang','hagai'],'flavor':'seurat','n_top_genes':2000,'n_bins':20,'normalize_total':10000,'selectionGate':'exact selected feature membership; differences retained, no native changes'},indent=2)+'\n')
results=[]
for cohort in ['kang','hagai']:
 p=Path('/Users/n/numivivo-pca-borrowed-20260912')/('pca-1' if cohort=='kang' else 'hagai-pca-1')
 m=json.loads((p/'model.json').read_text());a=anndata.read_h5ad(p/'original.h5ad');assert sp.issparse(a.X);a.X=a.X.astype(np.float64)
 sc.pp.normalize_total(a,target_sum=10000);sc.pp.log1p(a)
 start=time.monotonic();sc.pp.highly_variable_genes(a,flavor='seurat',n_top_genes=2000,n_bins=20,inplace=True);elapsed=time.monotonic()-start
 native={int(x) for x in m['selectedFeatureIndices']};ref=set(np.flatnonzero(a.var['highly_variable'].to_numpy()).tolist())
 assert a.var_names.tolist()==[x['featureID'] for x in m['features']]
 row={'cohort':cohort,'cells':a.n_obs,'features':a.n_vars,'nativeSelected':len(native),'scanpySelected':len(ref),'intersection':len(native&ref),'nativeOnly':[{'index':i,'featureID':str(a.var_names[i])} for i in sorted(native-ref)],'scanpyOnly':[{'index':i,'featureID':str(a.var_names[i])} for i in sorted(ref-native)],'exactMembershipPass':native==ref,'scanpyHVGSeconds':elapsed,'sourceModelSHA256':hashlib.sha256((p/'model.json').read_bytes()).hexdigest()}
 for key,nativekey in [('means','logMeanForBinning'),('dispersions','logDispersion'),('dispersions_norm','normalizedDispersion')]:
  x=np.asarray(a.var[key],dtype=float);y=np.array([f.get(nativekey,np.nan) for f in m['features']],dtype=float);mask=np.isfinite(x)&np.isfinite(y)
  row[key]={'maximumFiniteDifference':float(np.max(np.abs(x[mask]-y[mask]))),'finiteCompared':int(mask.sum()),'scanpyNonfinite':int((~np.isfinite(x)).sum()),'nativeNonfinite':int((~np.isfinite(y)).sum())}
 results.append(row);a.var.to_csv(r/(cohort+'-features.csv'));(r/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(row),flush=True)
