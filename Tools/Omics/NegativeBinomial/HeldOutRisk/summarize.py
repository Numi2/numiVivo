#!/usr/bin/env python3
"""Animal-level summaries preserving overlap, unavailable folds and paired methods."""
import argparse,gzip,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
protocol=json.loads((r/'protocol.json').read_text());checks=[];observations=[];native_failures=[]
model_hashes=json.loads((r/'remote-model-hashes.json').read_text()) if (r/'remote-model-hashes.json').exists() else None
for rec in protocol['cases']:
 d=r/rec['id'];runs=json.loads((d/'runs.json').read_text())
 for run in runs:
  if run['stage']=='fit':
   digest=model_hashes[rec['id']+'/model.json.gz'] if model_hashes else sha(d/'model.json.gz')
   assert digest==run['outputSHA256']
  else:
   path=d/'scores.json.gz';assert sha(path)==run['outputSHA256'];assert hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest()==run['logicalSHA256']
 if any(x['exitCode'] for x in runs):native_failures.append(rec['id']);continue
 check=json.loads((d/'check.json').read_text());assert check['modelSHA256']==runs[0]['outputSHA256'] and check['scoresSHA256']==sha(d/'scores.json.gz')
 assert check['checkerSHA256']==sha(Path(__file__).with_name('check.py'));checks.append(check)
 for row in check['observations']:observations.append(dict(case=rec['id'],cellGroup=rec['cellGroup'],**row))
frame=pd.DataFrame(observations);frame.to_csv(r/'fold-scores.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
animal=frame.groupby(['cellGroup','sampleID','condition','method'],sort=True).agg(meanLogScore=('meanLogScore','mean'),meanSquaredLog1pError=('meanSquaredLog1pError','mean'),folds=('case','count')).reset_index()
animal.to_csv(r/'animal-scores.tsv',sep='\t',index=False)
secondary=animal.pivot(index=['cellGroup','sampleID'],columns='method',values='meanSquaredLog1pError')
groups=[];paired=[]
for group in sorted({c['cellGroup'] for c in protocol['cases']}):
 sub=animal[animal.cellGroup==group];parts=[c for c in checks if c['cellGroup']==group]
 expected=[c for c in protocol['cases'] if c['cellGroup']==group]
 complete=len(parts)==len(expected)==16 and all(c['primaryAvailable'] and c['status']=='passed' for c in parts) and len(sub)==24 and sub.folds.eq(4).all()
 wide=sub.pivot(index=['sampleID','condition'],columns='method',values='meanLogScore')
 differences={}
 for left,right in [('empirical','MLE'),('fixed','MLE'),('empirical','fixed')]:
  delta=wide[left]-wide[right];name=left+'Minus'+right
  differences[name]=dict(mean=float(delta.mean()),animalsImproved=int((delta>0).sum()),animals=len(delta),minimum=float(delta.min()),maximum=float(delta.max()))
  paired.extend(dict(cellGroup=group,sampleID=sample,condition=condition,comparison=name,delta=float(value)) for (sample,condition),value in delta.items())
 groups.append(dict(cellGroup=group,primaryAvailable=bool(complete),availableFolds=sum(c['primaryAvailable'] for c in parts),
  minimumTestedGenes=min(c['testedGenes'] for c in parts),maximumTestedGenes=max(c['testedGenes'] for c in parts),
  priorSDRange=[min(c['priorSDLog2'] for c in parts),max(c['priorSDLog2'] for c in parts)],comparisons=differences))
pd.DataFrame(paired).to_csv(r/'animal-differences.tsv',sep='\t',index=False)
all_complete=all(g['primaryAvailable'] for g in groups) and not native_failures
result=dict(status='completed-held-out-count-risk-evaluation',primaryAvailable=all_complete,groups=groups,
 primaryEmpiricalMinusMLE=float(np.mean([g['comparisons']['empiricalMinusMLE']['mean'] for g in groups])) if all_complete else None,
 fixedMinusMLE=float(np.mean([g['comparisons']['fixedMinusMLE']['mean'] for g in groups])) if all_complete else None,
 empiricalMinusFixed=float(np.mean([g['comparisons']['empiricalMinusfixed']['mean'] for g in groups])) if all_complete else None,
 secondarySquaredLog1pError=dict(empiricalMinusMLE=float((secondary['empirical']-secondary['MLE']).mean()),fixedMinusMLE=float((secondary['fixed']-secondary['MLE']).mean()),empiricalMinusFixed=float((secondary['empirical']-secondary['fixed']).mean()),empiricalWorseThanMLEAnimalPopulationSummaries=int((secondary['empirical']>secondary['MLE']).sum()),qualification='Secondary error; lower is better; overlapping observations, no independence-based uncertainty'),
 folds=len(checks),independentAnimals=8,testedGeneFoldPairs=sum(c['testedGenes'] for c in checks),conditionalFits=sum(c['conditionalFits'] for c in checks),
 numericalFailures=[dict(case=c['case'],errors=c['errors']) for c in checks if c['status']!='passed'],nativeFailures=native_failures,
 maximumErrors={key:max(v['maximum'] for c in checks for v in c['metrics'] if v['metric']==key) for key in ['scaledScore','standardDeviation','aggregateLogScore','meanSquaredLog1pError']},
 protocolSHA256=protocol['protocolSHA256'],checkerSHA256=sha(Path(__file__).with_name('check.py')),
 qualification=protocol['qualification'],scope='Plug-in NB score conditional on observed library depth; train-only filtering, normalization, dispersion and effect priors; no uncertainty integration, FDR or causal accuracy claim')
(r/'summary.json').write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['groups','numericalFailures']},indent=2))
