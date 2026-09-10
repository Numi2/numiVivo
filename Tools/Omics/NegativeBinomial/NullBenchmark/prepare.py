#!/usr/bin/env python3
"""Freeze balanced untreated-cell splits and independent donor/arm counts."""
import argparse,copy,hashlib,json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

SOURCES={'kang':'6e5609e5bd0635ac1bd2a8f9a2bf98ad4fa9a012e8e4e2ad926436d073d745f0','hagai':'2d0a69d121bcf6e9b63b82fafd06e597f6165a927b87b0cc51cb5c90c2caa5f5'}
PROTOCOL='4203a2d11df8e58afb973b814055e7d7b8aed128a1aefc5995c785cfc93f1f57'
def sha(path):
 with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(path,x):path.write_text(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n')
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--study',choices=SOURCES,required=True);p.add_argument('--source-bundle',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
 a=p.parse_args();assert sha(Path(__file__).with_name('PROTOCOL.md'))==PROTOCOL
 source=a.source_bundle/'original.h5ad';assert sha(source)==SOURCES[a.study]
 a.out.mkdir(parents=True,exist_ok=False)
 old=json.loads((a.source_bundle/'plan.json').read_text());mapping=old['mapping'];request=old['contrasts'][0]
 obj=ad.read_h5ad(source);raw=obj.X.tocsr()
 assert sparse.issparse(obj.X) and np.isfinite(raw.data).all() and (raw.data>=0).all() and (raw.data==np.floor(raw.data)).all()
 raw=raw.astype(np.int64)
 ids=obj.var[mapping['featureIDColumn']].astype(str).tolist() if mapping.get('featureIDColumn') else obj.var_names.tolist()
 barcodes=obj.obs[mapping['barcodeColumn']].astype(str).tolist() if mapping.get('barcodeColumn') else obj.obs_names.tolist()
 original_samples=obj.obs[mapping['sampleColumn']].astype(str).tolist()
 samples={s['id']:s for s in mapping['samples']}
 controls={s['id']:s for s in mapping['samples'] if s['condition']==request['controlCondition']}
 control_rows=[i for i,s in enumerate(original_samples) if s in controls]
 assert len(set((original_samples[i],barcodes[i]) for i in range(obj.n_obs)))==obj.n_obs
 assert len(controls)==(8 if a.study=='kang' else 3)
 assert all(s['donorID']==s['biologicalReplicateID'] for s in controls.values())
 assert len({s['donorID'] for s in controls.values()})==len(controls)
 edits=[];assignments={}
 for seed in range(1,11):
  directory=a.out/str(seed);directory.mkdir()
  labels=list(original_samples);new_samples=copy.deepcopy(mapping['samples']);by_arm={}
  for sample,spec in sorted(controls.items()):
   rows=[i for i,s in enumerate(original_samples) if s==sample]
   ordered=sorted(rows,key=lambda i:(hashlib.sha256(f'numivivo-null-v1|{a.study}|{seed}|{sample}|{barcodes[i]}'.encode()).digest(),i))
   for parity,arm in enumerate(['shamA','shamB']):
    members=ordered[parity::2];assert len(members)>=10
    sid=sample+'__'+arm;record=dict(spec,id=sid,condition=arm);new_samples.append(record)
    for i in members:labels[i]=sid
    by_arm[sid]=dict(sample=record,sourceRows=sorted(members))
  categories=sorted(set(labels));lookup={x:i for i,x in enumerate(categories)}
  column='numivivo_null_seed_'+str(seed)
  edits.append(dict(path='obs/'+column,mode='add',value=dict(categorical=dict(categories=categories,codes=[lookup[x] for x in labels],ordered=False))))
  sm=sorted(by_arm,key=lambda sid:(by_arm[sid]['sample']['biologicalReplicateID'],by_arm[sid]['sample']['condition']))
  rr=[];cc=[]
  for i,sid in enumerate(sm):rr.extend([i]*len(by_arm[sid]['sourceRows']));cc.extend(by_arm[sid]['sourceRows'])
  membership=sparse.csr_matrix((np.ones(len(rr),dtype=np.int64),(rr,cc)),shape=(len(sm),obj.n_obs))
  counts=(membership@raw).toarray();assert counts.sum()==raw[control_rows].sum()
  np.savez_compressed(directory/'reference.npz',counts=counts,sourceRows=control_rows)
  write(directory/'memberships.json',dict(observationSampleIDs=sm,groups=by_arm,sourceRows=control_rows,originalBarcodes=barcodes,originalSamples=original_samples,features=ids))
  plan=copy.deepcopy(old);plan['mapping']=dict(mapping,sampleColumn=column,samples=new_samples,
   sourceDescription=mapping['sourceDescription']+'; untreated-cell balanced sham split. Arms are paired subsamples, not additional biological replicates.')
  contrasts=[]
  for mode in ['default','active']:
   c=copy.deepcopy(request);c.update(id='sham-'+mode,controlCondition='shamA',treatmentCondition='shamB')
   assert c['negativeBinomialOptions']['trend']=='gammaParametric'
   if mode=='active':c['negativeBinomialOptions']['zeroTotalDonorPolicy']='activeDonorProfile'
   contrasts.append(c)
  plan['contrasts']=contrasts
  # Bind the actual annotated hash only after the native annotation command.
  write(directory/'unbound-plan.json',plan)
  # Exact independent design and all-positive median-ratio normalization.
  donor=[by_arm[s]['sample']['donorID'] for s in sm];donors=sorted(set(donor))
  conditions=[by_arm[s]['sample']['condition'] for s in sm]
  design=np.array([[1,int(arm=='shamB')]+[int(d==ref) for ref in donors[1:]] for d,arm in zip(donor,conditions)])
  assert np.linalg.matrix_rank(design)==len(donors)+1
  positive=(counts>0).all(axis=0);assert positive.sum()>=10
  logs=np.log(counts[:,positive]);ratios=np.exp(logs-logs.mean(axis=0));factors=np.median(ratios,axis=1);factors=np.exp(np.log(factors)-np.log(factors).mean())
  eligible=(counts.sum(axis=0)>=10)&((counts>0).sum(axis=0)>=3)
  ri=directory/'r-input';ri.mkdir()
  pd.DataFrame(counts.T,index=ids,columns=sm).to_csv(ri/'counts.tsv',sep='\t',index_label='featureID')
  pd.DataFrame(design,index=sm,columns=['intercept','treatment-minus-control']+['donor:'+x for x in donors[1:]]).to_csv(ri/'design.tsv',sep='\t',index_label='sampleID')
  pd.DataFrame(dict(sampleID=sm,donor=donor,condition=conditions,sizeFactor=factors,libraryCounts=counts.sum(axis=1))).to_csv(ri/'samples.tsv',sep='\t',index=False)
  write(ri/'input.json',dict(contrast=[0,1]+[0]*(len(donors)-1),minimumFeatureCounts=10,minimumExpressingPseudobulks=3,eligibleFeatures=int(eligible.sum()),expectedGenes=[],
   sourceH5ADSHA256=sha(source),protocolSHA256=PROTOCOL,study=a.study,seed=seed,design='Independent reconstruction of native paired design and all-positive median-ratio normalization; comparison must verify against available native report'))
  assignments[str(seed)]={s:len(by_arm[s]['sourceRows']) for s in sm}
 write(a.out/'annotation-plan.json',dict(schemaVersion=1,source=dict(bytes=list(bytes.fromhex(SOURCES[a.study]))),
  provenance='Untreated-cell sham benchmark, protocol SHA256='+PROTOCOL+'. Ten predeclared balanced SHA256 partitions. Original counts and all metadata retained.',edits=edits))
 write(a.out/'preparation.json',dict(study=a.study,sourceSHA256=sha(source),protocolSHA256=PROTOCOL,sourceCells=obj.n_obs,sourceGenes=obj.n_vars,untreatedCells=len(control_rows),untreatedUMIs=int(raw[control_rows].sum()),donors=len(controls),seeds=list(range(1,11)),assignments=assignments,predictionsOrTestsFitted=False))
 write(a.out/'native-input.json',dict(sourceSHA256=SOURCES[a.study],protocolSHA256=PROTOCOL,sourceRows=control_rows,
  plans={str(seed):json.loads((a.out/str(seed)/'unbound-plan.json').read_text()) for seed in range(1,11)}))
 print((a.out/'preparation.json').read_text())

if __name__=="__main__":main()
