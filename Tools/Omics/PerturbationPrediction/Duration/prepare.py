#!/usr/bin/env python3
"""Freeze all donor holdouts from the already verified complete source aggregates."""
import argparse,json,datetime
from pathlib import Path
import anndata as ad,numpy as np,pandas as pd
from scipy import sparse
from common import sha,write,verify_files,hours

def main():
 p=argparse.ArgumentParser();p.add_argument('--source-study',type=Path,required=True);p.add_argument('--study',type=Path,required=True);a=p.parse_args()
 src=a.source_study/'prediction-inputs';out=a.study/'inputs';out.mkdir(exist_ok=False)
 prior=json.loads((src/'input-freeze.json').read_text());assert prior['status']=='frozen-before-fit';verify_files(src,prior['files'])
 cohort=json.loads((src/'query-cohort.json').read_text());samples=cohort['samples'];panel=json.loads((src/'panel.json').read_text())
 with np.load(src/'query-source-counts.npz',allow_pickle=False) as z:counts=z['counts'];features=z['featureIDs'].tolist();sampleids=z['sampleIDs'].tolist()
 assert counts.shape==(21,36601) and sampleids==[s['id'] for s in samples] and features==cohort['featureIDs'];assert cohort['sourceCells']==126633 and len(panel)==12993
 donors=sorted({s['donorID'] for s in samples});assert len(donors)==3
 assert sorted(i for g in cohort['sourceGroups'] for i in g['sourceCellIndices'])==list(range(cohort['sourceCells']))
 # Exact source aggregates are small; copy them for a self-contained arithmetic rerun.
 for name in ['query-source-counts.npz','query-cohort.json','panel.json']:(out/name).write_bytes((src/name).read_bytes())
 folds=[]
 for donor in donors:
  fid=donor.split(':')[-1];fold=out/fid;fold.mkdir();training=[i for i,s in enumerate(samples) if s['donorID']!=donor];query=[i for i,s in enumerate(samples) if s['donorID']==donor and s['condition']=='control'];truth=[i for i,s in enumerate(samples) if s['donorID']==donor and s['condition']!='control'];assert len(training)==14 and len(query)==1 and len(truth)==6
  actual_hours=sorted(hours(samples[i]['condition']) for i in truth)
  for role,indices in [('training',training),('query',query)]:
   selected=[samples[i] for i in indices]
   obj=ad.AnnData(sparse.csr_matrix(counts[indices]),obs=pd.DataFrame({'sample':[s['id'] for s in selected],'population':['QC-PBMC']*len(indices)},index=pd.Index([s['id'] for s in selected],dtype=object)),var=pd.DataFrame(index=pd.Index(features,dtype=object)))
   obj.write_h5ad(fold/(role+'.h5ad'),compression='gzip')
   plan=dict(schemaVersion=1,mapping=dict(schemaVersion=1,id='GSE226572-duration-'+fid+'-'+role,evidence='measured',sourceDescription='Complete source-verified population pseudobulk transport rows, not individual cells; development holdout '+fid+' '+role,countUnit='umiCount',matrixPath='X',sampleColumn='sample',groupColumn='population',samples=selected),featureNamespace='explicit-source-symbol-correspondence',perturbationID='IFNB1a-1000U-per-mL-PBMC-36h-culture')
   if role=='training':plan.update(controlCondition='control',exposures=[dict(condition=c,hours=hours(c)) for c in sorted({s['condition'] for s in selected if s['condition']!='control'},key=hours)],responseFeatureIDs=panel,provenance='Frozen retrospective GSE226572 duration development protocol; two training donors only, six exposures each; no query treated rows.')
   else:plan['hours']=actual_hours
   write(fold/(role+'.json'),plan)
  folds.append(dict(id=fid,heldOutDonor=donor,trainingRows=training,queryRows=query,outcomeRows=truth,hours=actual_hours))
 write(out/'folds.json',folds)
 write(out/'input-freeze.json',dict(status='frozen-before-duration-fit',createdUTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),developmentNotExternalValidation=True,protocolSHA256=sha(Path(__file__).with_name('PROTOCOL.md')),preparerSHA256=sha(Path(__file__)),priorInputFreezeSHA256=sha(src/'input-freeze.json'),sourceStudy=str(a.source_study),sourceCells=cohort['sourceCells'],files={str(f.relative_to(out)):sha(f) for f in sorted(out.rglob('*')) if f.is_file()},fitStarted=False,scoringStarted=False))
 print(json.dumps(dict(folds=len(folds),outcomes=sum(len(f['outcomeRows']) for f in folds),sourceCells=cohort['sourceCells'],panel=len(panel),inputFreezeSHA256=sha(out/'input-freeze.json'))),flush=True)
if __name__=='__main__':main()
