#!/usr/bin/env python3
"""Prepare full-source shared-panel inputs without changing count denominators."""
import argparse,collections,hashlib,json,subprocess,sys
from pathlib import Path
import anndata as ad,numpy as np,h5py
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'Logistic'))
from prepare import sha,write

def main():
 p=argparse.ArgumentParser();p.add_argument('--kang',type=Path,required=True);p.add_argument('--ding',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 kp=a.kang/'kang.h5ad';dp=a.ding/'ding.h5ad';k=ad.read_h5ad(kp);d=ad.read_h5ad(dp)
 assert k.shape==(24673,15706) and d.shape==(44031,33694)
 assert sha(dp)==json.loads((a.ding/'checks.json').read_text())['preparedSHA256']
 counts=collections.Counter(d.var['symbol'].astype(str));panel=sorted(x for x in k.var_names if counts[x]==1);assert len(panel)==14976
 ids=np.array([symbol if symbol in set(panel) else old for symbol,old in zip(d.var['symbol'].astype(str),d.var_names)],dtype=str)
 assert len(set(ids))==len(ids)
 d.var['reference_feature_id']=ids
 write(a.out/'feature-identities.json',dict(panel=panel,ding=[dict(originalID=str(old),symbol=str(symbol),nativeID=str(new)) for old,symbol,new in zip(d.var_names,d.var['symbol'],ids)],ambiguousKangSymbols=[x for x in k.var_names if counts[x]>1],absentKangSymbols=[x for x in k.var_names if counts[x]==0]))
 assigned=d.obs['annotationStatus'].astype(str).eq('source-assigned').to_numpy();assert int(assigned.sum())==29411
 # Native mapping chooses this added var column; original var index and metadata remain intact.
 d.write_h5ad(a.out/'ding-query.h5ad',compression='gzip');d[assigned].copy().write_h5ad(a.out/'ding-train.h5ad',compression='gzip')
 subprocess.run(['/bin/cp','-c',str(kp),str(a.out/'kang.h5ad')],check=True)
 records={}
 for name,original in [('kang',k),('ding-query',d),('ding-train',d[assigned])]:
  path=a.out/(name+'.h5ad');reloaded=ad.read_h5ad(path)
  assert (reloaded.X!=original.X).nnz==0;np.testing.assert_array_equal(reloaded.obs_names,original.obs_names);np.testing.assert_array_equal(reloaded.var_names,original.var_names)
  totals=np.asarray(reloaded.X.sum(axis=1)).ravel();assert np.all(totals>0)
  records[name]=dict(path=path.name,SHA256=sha(path),bytes=path.stat().st_size,shape=list(reloaded.shape),nonzeros=int(reloaded.X.nnz),totalCounts=int(totals.sum()),countsAndAxesExact=True)
  del reloaded
 mappingK=json.loads((a.kang/'kang-fit.json').read_text())['mapping'];mappingD=json.loads((a.ding/'fit.json').read_text())['mapping'];mappingD['featureIDColumn']='reference_feature_id';mappingD['groupColumn']='sourceCellType'
 for name,mapping in [('kang',mappingK),('ding',mappingD)]:
  mapping['sourceDescription']+=' Cross-study reference experiment: exact source-symbol panel; full source genes retained in normalization.'
  fit=dict(schemaVersion=1,mapping=mapping,featureNamespace='literal-source-gene-symbol-v1',labelProvenance='Original source labels in existing source-linked prepared Kang/Ding inputs; no cross-study label harmonization during fitting.',neighbors=15,logistic=dict(penalty=1,gradientTolerance=1e-8,maximumIterations=2000,maximumWork=100000000000),reduction=dict(normalizationTarget=10000,maximumCacheBytes=2000000000,maximumEntryVisits=2000000000,pca=dict(highlyVariableFeatures=2000,meanBins=20,components=20,maximumBasis=128,relativeResidualTolerance=1e-6,seed=7,featurePanel=panel)))
  write(a.out/(name+'-fit.json'),fit)
  queryMapping=dict(mapping);queryMapping.pop('groupColumn',None)
  write(a.out/(name+'-query.json'),dict(schemaVersion=1,mapping=queryMapping,featureNamespace='literal-source-gene-symbol-v1',maximumClassifierOperations=2000000000,maximumProjectionUpdates=2000000000))
 write(a.out/'input-freeze.json',dict(schemaVersion=1,protocolSHA256=sha(Path(__file__).with_name('PROTOCOL.md')),preparerSHA256=sha(__file__),sourceInputs={str(kp):sha(kp),str(dp):sha(dp)},records=records,files={p.name:sha(p) for p in sorted(a.out.iterdir()) if p.is_file()},queryCells=68704,fitStarted=False,scoringStarted=False))
 print(json.dumps(dict(status='prepared',panelFeatures=len(panel),records=records)),flush=True)
if __name__=='__main__':main()
