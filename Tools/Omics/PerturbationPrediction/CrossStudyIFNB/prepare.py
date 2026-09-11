#!/usr/bin/env python3
"""Freeze complete source-qualified pseudobulk inputs; no model fit or scoring."""
import argparse, collections, copy, gzip, hashlib, json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d): p.write_text(json.dumps(d,indent=2,sort_keys=True,allow_nan=False)+'\n')
def dense(bulk):
 m=bulk['matrix'];return sparse.csr_matrix((np.array(m['counts'],dtype=np.uint64),m['featureIndices'],m['rowOffsets']),shape=(m['cellCount'],m['featureCount'])).toarray()
def main():
 p=argparse.ArgumentParser();p.add_argument('--kang',type=Path,required=True);p.add_argument('--hirisa',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False,parents=True)
 k=json.loads((a.kang/'report.json').read_text());h=json.loads(gzip.decompress((a.hirisa/'inference-inputs/cohorts/02.json.gz').read_bytes()))
 assert sha(a.kang/'report.json')=='3f0389df9dd6b178f11c25608555daf60d66e895bc5b6b5ea9c50238fb5fd562'
 linkage=json.loads((a.hirisa/'native-independent-verification.json').read_text())
 assert linkage['status']=='passed' and linkage['nativeRawIngestionLinked']
 assert sha(a.hirisa/'inference-inputs/counts.npz')==linkage['inferenceCountsSHA256']=='ab60a81a7af8a0f4f3c22ec1090025c5f7433527e640a069dd5265b4c11b5722'
 hc=dense(h['pseudobulk'])
 with np.load(a.hirisa/'inference-inputs/counts.npz',allow_pickle=False) as f:
  assert f['featureIDs'].tolist()==h['pseudobulk']['featureIDs']
  for i,g in enumerate(h['pseudobulk']['groups']):
   assert len(g['sampleIDs'])==1
   np.testing.assert_array_equal(hc[i],f['counts'][f['sampleIDs'].tolist().index(g['sampleIDs'][0])])
 features={s:d['metadata']['features'] for s,d in [('Kang',k),('HIRISA',h)]}
 for s,d in [('Kang',k),('HIRISA',h)]:assert [x['id'] for x in features[s]]==d['pseudobulk']['featureIDs']
 names={'Kang':[x['id'] for x in features['Kang']],'HIRISA':[x['name'] for x in features['HIRISA']]}
 frequency={s:collections.Counter(x) for s,x in names.items()}
 common=sorted(x for x,n in frequency['Kang'].items() if n==1 and frequency['HIRISA'][x]==1)
 assert common
 panel=['symbol|'+x for x in common];datasets={};cohorts={};mapping_records=[]
 for study,d in [('Kang',k),('HIRISA',h)]:
  bulk=d['pseudobulk'];counts=dense(bulk);groups=bulk['groups'];ids=[]
  for j,f in enumerate(features[study]):
   name=names[study][j];shared=name in common
   fid='symbol|'+name if shared else study+'|'+f['id'];ids.append(fid)
   mapping_records.append(dict(study=study,sourceIndex=j,sourceID=f['id'],sourceName=name,transportID=fid,shared=shared,reason='exact-unique-source-symbol' if shared else ('ambiguous-source-symbol' if frequency[study][name]!=1 or frequency['HIRISA' if study=='Kang' else 'Kang'][name]>1 else 'unmatched-source-symbol')))
  assert len(ids)==len(set(ids));assert len(groups)==(16 if study=='Kang' else 10)
  samples=[]
  for i,g in enumerate(groups):
   assert g['organism']=='NCBITaxon:9606' and len(g['sampleIDs'])==1
   condition={'ctrl':'control','stim':'IFNB'}[g['condition']] if study=='Kang' else {'none':'control','IFNb':'IFNB'}[g['condition']]
   donor=study+':'+g['donorID'];samples.append(dict(id=g['sampleIDs'][0],donorID=donor,biologicalReplicateID=donor,condition=condition,organism=g['organism'],batchID='|'.join(g['batchIDs'])))
  obs=pd.DataFrame({'sample':pd.Series([x['id'] for x in samples],dtype=object),'comparison':pd.Series(['Bcell-context-transfer']*len(samples),dtype=object)})
  obs.index=pd.Index([study+'-aggregate-'+str(i) for i in range(len(samples))],dtype=object)
  obj=ad.AnnData(sparse.csr_matrix(counts),obs=obs,var=pd.DataFrame(index=pd.Index(ids,dtype=object)))
  datasets[study]=obj
  cohorts[study]=dict(samples=samples,sourceGroups=groups,sourceCellCount=len(d['metadata']['cells']),sourceFeatureCount=len(ids),libraryCounts=counts.sum(axis=1,dtype=np.uint64).tolist())
  np.savez_compressed(a.out/(study+'-source-counts.npz'),counts=counts,featureIDs=np.array(ids),sampleIDs=np.array([x['id'] for x in samples]))
  write(a.out/(study+'-cohort.json'),cohorts[study])
 write(a.out/'feature-mapping.json',mapping_records);write(a.out/'panel.json',panel)
 folds=[]
 for target,training in [('Kang','HIRISA'),('HIRISA','Kang')]:
  for donor in sorted({s['donorID'] for s in cohorts[target]['samples']}):
   query=[i for i,s in enumerate(cohorts[target]['samples']) if s['donorID']==donor and s['condition']=='control'];truth=[i for i,s in enumerate(cohorts[target]['samples']) if s['donorID']==donor and s['condition']=='IFNB'];assert len(query)==len(truth)==1
   for mode,origin in [('cross',training),('within',target)]:
    train=[i for i,s in enumerate(cohorts[origin]['samples']) if mode=='cross' or s['donorID']!=donor]
    name=mode+'-'+origin+'-to-'+target+'-'+donor.split(':')[1];root=a.out/name;root.mkdir()
    for role,study,indices in [('training',origin,train),('query',target,query)]:
     datasets[study][indices,:].copy().write_h5ad(root/(role+'.h5ad'),compression='gzip')
     mapping=dict(schemaVersion=1,id=name+'-'+role,evidence='measured',sourceDescription='Complete source-qualified donor pseudobulk rows, not individual cells. '+study+' '+role+'. Explicit Bcell context comparison; see frozen protocol and original source groups.',countUnit='umiCount',matrixPath='X',sampleColumn='sample',groupColumn='comparison',samples=[cohorts[study]['samples'][i] for i in indices])
     plan=dict(schemaVersion=1,mapping=mapping,featureNamespace='explicit-source-symbol-correspondence',perturbationID='IFNB-combined-study-context-transfer')
     if role=='training':plan.update(controlCondition='control',treatmentCondition='IFNB',provenance='Fixed Kang-HIRISA IFNB transfer protocol; full source denominators, exact shared symbols; different assay/timing/health/preparation. No unseen-study development or causal claim.',responseFeatureIDs=panel)
     write(root/(role+'.json'),plan)
    folds.append(dict(id=name,mode=mode,trainingStudy=origin,queryStudy=target,heldOutDonor=donor,trainingIndices=train,queryIndex=query[0],scoringIndex=truth[0]))
 write(a.out/'folds.json',folds)
 sources={str(p):sha(p) for p in [a.kang/'report.json',a.kang/'receipt.json',a.hirisa/'inference-inputs/cohorts/02.json.gz',a.hirisa/'native-independent-verification.json',a.hirisa/'inference-inputs/counts.npz',Path(__file__),Path(__file__).with_name('PROTOCOL.md')]}
 write(a.out/'source-bindings.json',sources)
 frozen={str(p.relative_to(a.out)):sha(p) for p in sorted(a.out.rglob('*')) if p.is_file()}
 write(a.out/'input-freeze.json',dict(schemaVersion=1,files=frozen,fitStarted=False,sharedFeatures=len(panel),sourceFeatureCounts={s:len(v) for s,v in features.items()},folds=len(folds),limitations='Previously inspected studies; exact source symbols only; biological replicate counts 8 and 5; no intervals.'))
 print(json.dumps(dict(sharedFeatures=len(panel),folds=len(folds),cells={s:c['sourceCellCount'] for s,c in cohorts.items()},freezeSHA256=sha(a.out/'input-freeze.json'))))
if __name__=='__main__':main()
