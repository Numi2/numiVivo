"""Read-only decision support over existing biological records; never reveals data.

Existing prediction/scoring owners stay byte-identical for v0.4 replay. This
module computes population coverage and decomposes numerical model differences.
It neither predicts a response nor changes a scientific recommendation threshold.
"""
from functools import lru_cache
from pathlib import Path
import copy
import numpy as np
from wetlab import read,sha,require

@lru_cache(maxsize=48)
def _population(config,digest,specimen,population,genes):
 from learned_spatial import load,source
 a,root=load(config);s=next(x for x in a['specimens'] if x['id']==specimen);data=source(s['source'],s['sourceSHA256']);
 if not s.get('inferenceSupported'):
  return {'modelID':a['id'],'specimen':specimen,'population':population,'populationAnchors':0,'reference':{'referenceCells':0,'safeHarbourCells':0,'barcodeNegativeCells':0,'sameProjectedTypeCells':0},'genes':[],'availableMarkers':[],'candidates':[],'recommendationSupported':False,'reason':'Geometry-only training specimen; no population prediction support. Select registered chip1 for experimental inference.','trainingRankingComparator':'UNAVAILABLE'}
 p=next(x for x in a['populations'] if x['id']==population);rows=[r for r in read(root/'query-rows.json') if r['cellType']+'|'+str(r['neighborhood'])==population];refs=sorted({j for r in rows for j in r['referenceRows']});features=read(root/'features.json');mask=np.load(root/'model-feature-mask.npy');index={g:i for i,g in enumerate(features)};labels=data.obs.guide_assignment.astype(str)
 def counts(ids):
  return {'referenceCells':len(ids),'safeHarbourCells':sum(labels.iloc[j] in ('mSafe','sgrna_mSafe') for j in ids),'barcodeNegativeCells':sum(labels.iloc[j]=='barcode-negative' for j in ids),'sameProjectedTypeCells':sum(str(data.obs.projected_cell_type.iloc[j])==p['cellType'] for j in ids)}
 coverage=[]
 for gene in genes:
  j=index.get(gene);measured=j is not None and bool(data.var.measured_in_source.iloc[j]);modeled=measured and bool(mask[j]);v=data.X[refs,j].toarray().ravel() if measured else []
  coverage.append({'gene':gene,'measured':bool(measured),'modeled':bool(modeled),'referenceCells':len(refs) if measured else 0,'detectedReferenceCells':int(np.count_nonzero(v)),'units':'raw UMI detection in selected population reference','scope':'selected population reference rows, not specimen totals','status':'MEASURED' if measured else 'UNAVAILABLE','predictionStatus':'MODEL INFERENCE' if modeled else 'UNAVAILABLE'})
 trainFile=Path(a['inputs'])/'chip2-rows.json';training=read(trainFile) if trainFile.exists() else None
 candidates=[]
 for target in [t['target'] for t in a['targets']]:
  for role in ('direct','neighbor'):
   matched=[r for r in rows if r['target']==target and r['role']==role];outcomes=sorted({j for r in matched for j in r['outcomeRows']});rr=sorted({j for r in matched for j in r['referenceRows']});c=counts(rr);tr=[r for r in training or [] if r['target']==target and r['role']==role and r['cellType']==p['cellType']];reasons=[]
   if len(outcomes)<5:reasons.append('Fewer than five outcome cells in this target / population')
   if c['safeHarbourCells']<5:reasons.append('Fewer than five verified safe-harbour reference cells in this population')
   if c['barcodeNegativeCells']:reasons.append('Barcode-negative reference exposure is unknown')
   reasons.append('Transfer uncertainty is uncalibrated; numerical ordering is not a recommendation')
   candidates.append({'target':target,'role':role,'population':population,'outcomeCells':len(outcomes),'outcomeQuantity':'metadata count only; RNA remains sealed until authorized','reference':c,'trainingTargetSeen':None if training is None else any(r['target']==target for r in training),'trainingPopulationGroups':None if training is None else len(tr),'trainingContext':'chip2 pooled mice 2 + 3; projected cell identity' if training is not None else 'UNAVAILABLE: retained training rows missing','specimenSeenInTraining':False,'canExecute':bool(outcomes and rr),'recommendationSupported':False,'reasons':reasons,'annotation':'control-projected identity, uncertain','biologicalUnits':s.get('animals',[]),'preservationSupported':False})
 # Suggestions are measurement availability, not recommended biological objectives.
 values=np.asarray((data.X[refs]>0).sum(0)).ravel() if refs else np.zeros(len(features));eligible=[j for j in range(len(features)) if mask[j] and bool(data.var.measured_in_source.iloc[j]) and values[j]>0];eligible=sorted(eligible,key=lambda j:(-values[j],features[j]))[:6]
 return {'modelID':a['id'],'modelVersion':a.get('modelVersion',a['id']),'specimen':specimen,'population':population,'populationAnchors':p['anchors'],'reference':counts(refs),'genes':coverage,'availableMarkers':[{'gene':features[j],'detectedReferenceCells':int(values[j]),'referenceCells':len(refs)} for j in eligible],'markerMeaning':'Measured marker suggestions only; researcher defines the objective','candidates':candidates,'recommendationSupported':False,'reason':'Insufficient verified controls and unqualified model transfer','trainingRankingComparator':'UNAVAILABLE for this spatial experiment; no-change is not a fitted training ranking','sourceSHA256':s['sourceSHA256'],'trainingMetadataSHA256':sha(trainFile) if trainFile.exists() else None}

def population_support(config,selection,genes):
 return copy.deepcopy(_population(str(Path(config).resolve()),sha(config),selection['specimen'],selection['population'],tuple(genes)))

def model_difference(card):
 """Per-gene utility decomposition, not an explanation of biological mechanism."""
 axes=[]
 for axis in card.get('axisResults',[]):
  ev=axis['result'].get('objectiveEvaluation');
  if not ev:continue
  arms=[a for a in ev['arms'] if a['role']=='direct' and a['predictedUtility'] is not None];utilities={a['target']:a['predictedUtility'] for a in arms};utilities['no-intervention']=0.;rank={t:i+1 for i,t in enumerate(sorted(utilities,key=lambda t:(-utilities[t],t)))};axes.append({'modelID':axis['modelID'],'conditionID':axis['conditionID'],'arms':arms,'rank':rank,'utilities':utilities})
 if len(axes)<2:return {'available':False,'reason':'Seal at least two compatible models to compare exact pooled-reference predictions; preview rankings use a different reference.','recommendationSupported':False}
 baseline=axes[0];comparisons=[]
 for other in axes[1:]:
  differences=[]
  for target in sorted(set(baseline['utilities'])&set(other['utilities'])):
   left=next((a for a in baseline['arms'] if a['target']==target),None);right=next((a for a in other['arms'] if a['target']==target),None);genes=[]
   if left and right:
    for l,r in zip(left['genes'],right['genes']):
     require(l['gene']==r['gene'],'Gene axis mismatch');genes.append({'gene':l['gene'],'leftChange':l['predicted']-l['control'],'rightChange':r['predicted']-r['control'],'utilityDifference':-((r['predicted']-r['control'])-(l['predicted']-l['control']))/len(left['genes']),'population':left['population']})
   differences.append({'target':target,'leftRank':baseline['rank'][target],'rightRank':other['rank'][target],'leftUtility':baseline['utilities'][target],'rightUtility':other['utilities'][target],'utilityDifference':other['utilities'][target]-baseline['utilities'][target],'genes':sorted(genes,key=lambda g:-abs(g['utilityDifference']))})
  comparisons.append({'leftModel':baseline['modelID'],'rightModel':other['modelID'],'leftCondition':baseline['conditionID'],'rightCondition':other['conditionID'],'candidates':differences,'rankingsChanged':any(d['leftRank']!=d['rightRank'] for d in differences),'meaning':'Numerical differences only; agreement is not calibrated confidence'})
 return {'available':True,'comparisons':comparisons,'recommendationSupported':False}

def gene_evidence(config,run,gene,target,role):
 from learned_spatial import LearnedSpatialResponseAdapter,check
 from safetensors.numpy import load_file
 run=Path(run);reg=check(run);owner=LearnedSpatialResponseAdapter();field=owner.feature(config,reg['specimen'],gene,run);require(field['modeled'],'Gene lacks measured/model feature support; select a supported marker');row=next(a for a in field['arms'] if a['target']==target and a['role']==role);features=read(run/'features.json');j=features.index(gene);base=load_file(str(run/'baselines.safetensors'));types=[tuple(x) for x in read(run/'baseline-types.json')];key=(target,role,reg['plan']['population'].split('|')[0]);tid=next(i for i,t in enumerate(reg['assay']['targets']) if t['target']==target);prior=base['typePrior'][types.index(key)] if key in types else base['prior'][tid];matched=float(row['control']+prior[j]);observed=row['observed'];return {'gene':gene,'target':target,'role':role,'population':reg['plan']['population'],'record':run.name,'units':'log1p(CPM)','control':row['control'],'predicted':row['predicted'],'observed':observed,'predictedResidual':row['residual'],'matchedBaseline':matched,'matchedResidual':None if observed is None else matched-observed,'noChangeResidual':None if observed is None else row['control']-observed,'outcomeCells':len(row['observedRows']) if observed is not None else None,'stage':'REVEALED OBSERVATION' if observed is not None else 'SEALED PREDICTION','mechanism':'UNAVAILABLE; numerical attribution only'}
