#!/usr/bin/env python3
"""Verify native/reference cohorts and report every predeclared null comparison."""
import argparse,gzip,hashlib,json,subprocess
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import ndtr

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--allow-incomplete',action='store_true');p.add_argument('--remote-models',action='store_true',help='Read missing detailed reports from the retained macmini run; verify all payload hashes')
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);metrics=[];checks=[];missing=[]
def load(path):return json.loads(path.read_text())
def sha(b):return hashlib.sha256(b).hexdigest()
def bh(values):
 values=np.asarray(values,dtype=float);order=np.argsort(values,kind='stable');result=np.empty_like(values)
 if len(values):result[order]=np.minimum(1,np.minimum.accumulate((values[order]*len(values)/np.arange(1,len(values)+1))[::-1])[::-1])
 return result

def measure(study,seed,name,frame,eligible,scope='method-family',check_bh=True):
 valid=frame['pValue'].notna();pv=frame.loc[valid,'pValue'].to_numpy(float)
 assert np.isfinite(pv).all() and ((pv>=0)&(pv<=1)).all()
 if check_bh:assert np.allclose(bh(pv),frame.loc[valid,'adjustedPValue'].to_numpy(float),rtol=1e-10,atol=1e-12),(study,seed,name,'BH')
 q=frame['adjustedPValue'].dropna().to_numpy(float);assert np.isfinite(q).all() and ((q>=0)&(q<=1)).all()
 result=dict(study=study,seed=seed,method=name,scope=scope,eligibleGenes=eligible,reportedRows=len(frame),testedGenes=int(valid.sum()),withheldEligibleGenes=eligible-int(valid.sum()),
  rawPBelow001=int((pv<0.01).sum()),rawPBelow005=int((pv<0.05).sum()),rawPFractionBelow001=float((pv<0.01).mean()) if len(pv) else None,rawPFractionBelow005=float((pv<0.05).mean()) if len(pv) else None,
  bhBelow001=int((q<0.01).sum()),bhBelow005=int((q<0.05).sum()),bhBelow010=int((q<0.1).sum()),anyBHBelow005=bool((q<0.05).any()))
 metrics.append(result);return result

for study in ['kang','hagai']:
 for seed in range(1,11):
  d=a.root/study/str(seed);members=load(d/'memberships.json');ref=np.load(d/'reference.npz');y=ref['counts'];ids=members['features'];samples=members['observationSampleIDs']
  expected_x=pd.read_csv(d/'r-input/design.tsv',sep='\t',index_col=0).to_numpy(float)
  expected_s=pd.read_csv(d/'r-input/samples.tsv',sep='\t',float_precision='round_trip')
  eligible=(y.sum(axis=0)>=10)&((y>0).sum(axis=0)>=3);n_eligible=int(eligible.sum());tables={}
  for mode in ['default','active']:
   directory=d/mode;archive=directory/'archive.json';name='native-'+mode
   if not archive.exists():missing.append(dict(study=study,seed=seed,method=name));continue
   ar=load(archive);assert ar['status']=='native-publication-and-replay-verified-before-archival'
   payload={}
   for e in ar['entries']:
    assert e['logicalPath'] in ['plan.json','report.json','receipt.json'] and e['path']==e['logicalPath']+'.gz'
    f=directory/e['path']
    if f.exists():b=f.read_bytes()
    else:
     assert a.remote_models and e['logicalPath']=='report.json'
     remote=Path('/Users/n/numivivo-null-benchmark-20260910')/study/str(seed)/mode/e['path']
     b=subprocess.check_output(['ssh','macmini','cat '+str(remote)],timeout=60)
    assert sha(b)==e['sha256'];raw=gzip.decompress(b);assert sha(raw)==e['logicalSHA256'];payload[e['logicalPath']]=raw
   report=json.loads(payload['report.json']);receipt=json.loads(payload['receipt.json']);plan=json.loads(payload['plan.json'])
   assert bytes(receipt['source']['bytes']).hex()==ar['sourceSHA256']==bytes(plan['cellSelection']['source']['bytes']).hex()
   assert bytes(receipt['report']['bytes']).hex()==sha(payload['report.json']) and bytes(receipt['plan']['bytes']).hex()==sha(payload['plan.json'])
   assert report['sourceObservationIndices']==members['sourceRows']==ref['sourceRows'].tolist()
   assert report['sourceCellCount']==len(members['originalBarcodes'])
   assert report['pseudobulk']['featureIDs']==ids
   matrix=report['pseudobulk']['matrix'];native=sparse.csr_matrix((np.asarray(matrix['counts'],dtype=np.int64),matrix['featureIndices'],matrix['rowOffsets']),shape=(matrix['cellCount'],matrix['featureCount'])).toarray()
   assert np.array_equal(native,y)
   for j,group in enumerate(report['pseudobulk']['groups']):
    sid=samples[j];g=members['groups'][sid];assert group['sampleIDs']==[sid]
    assert [members['sourceRows'][i] for i in group['sourceCellIndices']]==g['sourceRows']
    for key in ['donorID','biologicalReplicateID','condition','organism']:assert group[key]==g['sample'][key]
   assert [v['barcode'] for v in report['metadata']['cells']]==[members['originalBarcodes'][i] for i in members['sourceRows']]
   assert sum(q['totalCounts'] for q in report['quality'])==int(y.sum())
   c=report['contrasts'][0];design=c['design'];assert len(report['contrasts'])==1
   assert np.array_equal(design['rows'],expected_x) and np.array_equal(design['libraryCounts'],y.sum(axis=1))
   assert np.allclose(design['sizeFactorValues'],expected_s['sizeFactor'],rtol=1e-12,atol=1e-12)
   assert [o['sampleIDs'][0] for o in design['observations']]==samples
   frame=pd.DataFrame(c['features']).set_index('featureID');assert list(frame.index)==ids
   assert np.array_equal(frame['status'].eq('filteredLowExpression').to_numpy(),~eligible)
   assert np.array_equal(frame['totalCounts'].to_numpy(),y.sum(axis=0))
   for col in ['pValue','adjustedPValue','log2FoldChange','standardError']:
    if col not in frame:frame[col]=np.nan
   assert int(frame['pValue'].notna().sum())==c['testedFeatures']
   tested=frame['pValue'].notna()
   z=frame.loc[tested,'log2FoldChange'].to_numpy(float)/frame.loc[tested,'standardError'].to_numpy(float)
   reference_p=2*ndtr(-np.abs(z))
   p_error=float(np.max(np.abs(reference_p-frame.loc[tested,'pValue'].to_numpy(float)))) if tested.any() else 0.0
   assert p_error<1e-10,(study,seed,mode,'Wald probability',p_error)
   result=measure(study,seed,name,frame,n_eligible)
   result['statuses']=dict(Counter(frame['status']));diag=c['negativeBinomial']['features']
   result['geneWiseLowerBoundaries']=sum(v.get('geneWiseLowerBoundary',False) for v in diag)
   result['geneWiseUpperBoundaries']=sum(v.get('geneWiseUpperBoundary',False) for v in diag)
   result['dispersionOutliers']=sum(v.get('dispersionOutlier',False) for v in diag)
   result['nonconvergedFinalFits']=sum(not v['finalFit']['converged'] for v in diag if 'finalFit' in v)
   compact=frame[['status','log2FoldChange','standardError','pValue','adjustedPValue']]
   compact.to_csv(a.out/f'{study}-{seed}-{name}.tsv.gz',sep='\t',compression={'method':'gzip','mtime':0})
   tables[name]=frame
   checks.append(dict(study=study,seed=seed,method=name,exactCountsMembershipDesign=True,reportSHA256=sha(payload['report.json']),maxWaldProbabilityError=p_error,maxSizeFactorError=float(np.max(np.abs(np.array(design['sizeFactorValues'])-expected_s['sizeFactor'])))))
  rdir=d/'r-reference'
  if not (rdir/'runs.json').exists():missing.append(dict(study=study,seed=seed,method='R'));continue
  statuses=load(rdir/'runs.json')
  for name,record in statuses.items():
   if record['status']!='completed':missing.append(dict(study=study,seed=seed,method=name,failure=record));continue
   frame=pd.read_csv(rdir/(name+'.tsv.gz'),sep='\t',index_col='featureID',float_precision='round_trip')
   assert list(frame.index)==np.asarray(ids)[eligible].tolist()
   result=measure(study,seed,name,frame,n_eligible)
   result['warnings']=record['warnings'];result['messages']=record['messages']
   if 'betaConverged' in frame:result['nonconvergedCoefficients']=int((~frame['betaConverged'].astype(bool)).sum())
   tables[name]=frame
  for mode in ['native_size_factors','package_normalization']:
   f=rdir/(mode+'-DESeq2-default-results.tsv.gz')
   if f.exists():
    frame=pd.read_csv(f,sep='\t',index_col='featureID',float_precision='round_trip').rename(columns={'pvalue':'pValue','padj':'adjustedPValue'})
    measure(study,seed,mode+'-DESeq2-default-policy',frame,n_eligible,scope='ordinary-DESeq2-policy',check_bh=False)
  if len(tables)==8:
   joint=sorted(set.intersection(*[set(f.index[f['pValue'].notna()]) for f in tables.values()]))
   for name,frame in tables.items():
    sub=frame.loc[joint].copy();sub['adjustedPValue']=bh(sub['pValue'].to_numpy(float));measure(study,seed,name,sub,len(joint),scope='joint-family')
native_runs=load(a.root/'native-runs.json')['runs'] if (a.root/'native-runs.json').exists() else []
for item in missing:
 if item['method'].startswith('native-'):
  suffix='/{}/{}/{}'.format(item['study'],item['seed'],item['method'][7:])
  attempts=[r for r in native_runs if r['command'][0]=='singlecell-h5ad-pseudobulk' and r['command'][-1].endswith(suffix)]
  if attempts and attempts[-1]['exitCode']!=0:item['failure']=attempts[-1]
pending=[m for m in missing if 'failure' not in m]
failures=[m for m in missing if 'failure' in m]
if not a.allow_incomplete:assert not pending,pending
summary=[]
for study in ['kang','hagai']:
 for method,scope in sorted({(m['method'],m['scope']) for m in metrics if m['study']==study}):
  rows=[m for m in metrics if (m['study'],m['method'],m['scope'])==(study,method,scope)]
  summary.append(dict(study=study,method=method,scope=scope,splits=len(rows),splitsWithAnyBH005=sum(r['anyBHBelow005'] for r in rows),totalBH005=sum(r['bhBelow005'] for r in rows),
   minimumTested=min(r['testedGenes'] for r in rows),maximumTested=max(r['testedGenes'] for r in rows),meanRawPFraction005=float(np.mean([r['rawPFractionBelow005'] for r in rows if r['rawPFractionBelow005'] is not None])) if any(r['rawPFractionBelow005'] is not None for r in rows) else None))
(a.out/'results.json').write_text(json.dumps(dict(status='partial' if pending else ('completed-with-fit-failures' if failures else 'completed-all-declared-comparisons'),missing=missing,pending=pending,failures=failures,perSplit=metrics,summary=summary,checks=checks,
 qualification='Conditional untreated-cell split null stress test. Overlapping splits and genes are not independent; no general FDR, independent-donor or biological-power qualification.'),indent=2,allow_nan=False)+'\n')
print(json.dumps(dict(status='partial' if pending else 'completed',nativeChecks=len(checks),failedMethods=len(failures),pending=len(pending))))
