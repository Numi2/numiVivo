#!/usr/bin/env python3
"""Freeze complete native count aggregates, exact feature panel and control-only queries."""
import argparse,collections,json
from pathlib import Path
import anndata as ad,numpy as np,pandas as pd
from scipy import sparse
from download import sha,write

def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);a=p.parse_args();root=a.study;native=root/'native-aggregation';out=root/'prediction-inputs';out.mkdir(exist_ok=False)
 execution=json.loads((native/'execution.json').read_text());assert execution['status']=='passed-native-aggregation-source-handoff-and-reconstruction'
 sources={cohort:json.loads((native/cohort/'report.json').read_text()) for cohort in ['kang','query']}
 names={c:[f['name'] if c=='query' else f['id'] for f in d['metadata']['features']] for c,d in sources.items()};freq={c:collections.Counter(v) for c,v in names.items()};panel_symbols=sorted(n for n,count in freq['kang'].items() if count==1 and freq['query'][n]==1);panel=['symbol|'+n for n in panel_symbols];assert panel
 cohorts={};mapping=[]
 for cohort,report in sources.items():
  for group in report['pseudobulk']['groups']:
   expected=[i for i,c in enumerate(report['metadata']['cells']) if c['sampleID'] in group['sampleIDs']];assert expected==group['sourceCellIndices'],'Native group membership differs from original sample identities'
  b=report['pseudobulk'];m=b['matrix'];counts=sparse.csr_matrix((np.array(m['counts'],dtype=np.uint64),m['featureIndices'],m['rowOffsets']),shape=(m['cellCount'],m['featureCount'])).toarray();ids=[]
  for i,f in enumerate(report['metadata']['features']):
   symbol=names[cohort][i];shared=symbol in panel_symbols;fid='symbol|'+symbol if shared else cohort+'|'+f['id'];ids.append(fid);mapping.append(dict(cohort=cohort,originalID=f['id'],symbol=symbol,transportID=fid,shared=shared,reason='exact-unique-symbol' if shared else 'unmatched-or-ambiguous-symbol'))
  assert len(set(ids))==len(ids);samples=[]
  for i,g in enumerate(b['groups']):
   condition={'ctrl':'control','stim':'IFNB'}[g['condition']] if cohort=='kang' else g['condition'];donor=cohort+':'+g['donorID'];samples.append(dict(id=cohort+'-aggregate-'+str(i),donorID=donor,biologicalReplicateID=donor,condition=condition,organism=g['organism'],batchID='|'.join(g['batchIDs'])))
  expected=16 if cohort=='kang' else 21;assert len(samples)==expected
  cohorts[cohort]=dict(samples=samples,sourceGroups=b['groups'],sourceCells=len(report['metadata']['cells']),featureIDs=ids,reportSHA256=sha(native/cohort/'report.json'))
  np.savez_compressed(out/(cohort+'-source-counts.npz'),counts=counts,featureIDs=np.array(ids),sampleIDs=np.array([s['id'] for s in samples]));write(out/(cohort+'-cohort.json'),cohorts[cohort])
  indices=list(range(len(samples))) if cohort=='kang' else [i for i,s in enumerate(samples) if s['condition']=='control'];assert len(indices)==(16 if cohort=='kang' else 3)
  selected=[samples[i] for i in indices];obs=pd.DataFrame(dict(sample=[s['id'] for s in selected],population=['QC-PBMC']*len(indices)),index=pd.Index([s['id'] for s in selected],dtype=object));obj=ad.AnnData(sparse.csr_matrix(counts[indices]),obs=obs,var=pd.DataFrame(index=pd.Index(ids,dtype=object)));role='training' if cohort=='kang' else 'query';obj.write_h5ad(out/(role+'.h5ad'),compression='gzip')
  plan=dict(schemaVersion=1,mapping=dict(schemaVersion=1,id='GSE226572-transfer-'+role,evidence='measured',sourceDescription='Complete source-verified whole-population pseudobulk rows, not individual cells. '+cohort+' '+role+'. Query contains only pooledzero-hourcontrols. See frozen protocol for QC and duration/context limitations.',countUnit='umiCount',matrixPath='X',sampleColumn='sample',groupColumn='population',samples=selected),featureNamespace='explicit-source-symbol-correspondence',perturbationID='IFNB-PBMC-fixed-Kang-response-transfer')
  if cohort=='kang':plan.update(controlCondition='control',treatmentCondition='IFNB',responseFeatureIDs=panel,donorResponseIntervalCoverage=.95,provenance='FrozenGSE226572protocol:Kang-six-hour population RNA response transferred unchanged to all query timepoints; no temporal fit or query outcome access.')
  write(out/(role+'.json'),plan)
 write(out/'panel.json',panel);write(out/'feature-mapping.json',mapping);write(out/'input-freeze.json',dict(status='frozen-before-fit',aggregationFreezeSHA256=sha(native/'aggregation-freeze.json'),protocolSHA256=sha(Path(__file__).with_name('PROTOCOL.md')),preparerSHA256=sha(Path(__file__)),files={str(f.relative_to(out)):sha(f) for f in sorted(out.iterdir()) if f.is_file()},trainingDonors=8,queryDonors=3,queryTreatedLibraries=18,sharedFeatures=len(panel),fitStarted=False,scoringStarted=False));print(json.dumps(dict(status='frozen-before-fit',sharedFeatures=len(panel),trainingDonors=8,queryDonors=3,treatedTimePoints=18)),flush=True)
if __name__=='__main__':main()
