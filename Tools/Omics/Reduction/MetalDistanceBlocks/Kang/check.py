from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
import json,hashlib,numpy as np
from sklearn.metrics import adjusted_rand_score,balanced_accuracy_score
from scipy.optimize import linear_sum_assignment
r=Path(__file__).parent;p=json.loads((r/'protocol.json').read_text());terminal=json.loads((r/'terminal.json').read_text());assert len(terminal['results'])==15 and all(x['exitCode']==0 for x in terminal['results'])
meta=json.loads((r/'pca/metadata.json').read_text());assert hashlib.sha256((r/'pca/metadata.json').read_bytes()).hexdigest()==p['originalMetadataSHA256'];cells=meta['cells'];n=len(cells);assert n==p['cells']
samples={s['id']:s for s in meta['samples']};type_names,types=np.unique([c['group'] for c in cells],return_inverse=True);condition_names,conditions=np.unique([samples[c['sampleID']]['condition'] for c in cells],return_inverse=True);donor_names,donors=np.unique([samples[c['sampleID']]['donorID'] for c in cells],return_inverse=True)
rec=np.dtype([('row','<u4'),('col','<u4'),('value','<f8')]);neighbors={};edges={};annotation={}
for mode in ['cpu','metal']:
 a=np.fromfile(r/mode/'neighbors.bin',rec);assert len(a)==n*20 and np.array_equal(a['row'],np.repeat(np.arange(n),20));idx=a['col'].reshape(n,20);assert np.array_equal(idx[:,0],np.arange(n)) and np.all(idx<n) and all(len(set(x))==20 for x in idx);neighbors[mode]=idx
 e=np.fromfile(r/mode/'edges.bin',rec);assert np.all(np.isfinite(e['value'])) and np.all((e['value']>0)&(e['value']<=1));edges[mode]=e
 # Local annotation concordance only, not a donor-held-out predictor.
 type_votes=np.array([np.bincount(types[row[1:]],minlength=len(type_names)).argmax() for row in idx]);condition_votes=np.array([np.bincount(conditions[row[1:]],minlength=len(condition_names)).argmax() for row in idx])
 pertype=[]
 for i,name in enumerate(type_names):
  sel=types==i;present=np.unique(conditions[sel]);pertype.append(dict(type=str(name),cells=int(sel.sum()),neighborTypeRecall=float(np.mean(type_votes[sel]==i)),conditionBalancedAccuracy=float(balanced_accuracy_score(conditions[sel],condition_votes[sel])) if len(present)>1 else None,conditionClasses=len(present)))
 annotation[mode]=dict(perType=pertype,conditionBalancedAccuracy=float(balanced_accuracy_score(conditions,condition_votes)),typeHomophily=float(np.mean(types[idx[:,1:]]==types[:,None])),conditionHomophily=float(np.mean(conditions[idx[:,1:]]==conditions[:,None])))
missing=sum(len(set(a)-set(b)) for a,b in zip(neighbors['cpu'],neighbors['metal']));changedRows=int(np.count_nonzero(np.any(neighbors['cpu']!=neighbors['metal'],axis=1)))
keys={m:e['row'].astype(np.int64)*n+e['col'] for m,e in edges.items()};common,ic,im=np.intersect1d(keys['cpu'],keys['metal'],return_indices=True,assume_unique=True)
comparisons=[];g=p['gates']
for seed in p['seeds']:
 results={m:json.loads((r/(m+'-'+str(seed))/'result.json').read_text()) for m in neighbors}
 for result in results.values():assert result['cells']==[{k:c[k] for k in ['sampleID','barcode']} for c in cells]
 a=np.asarray(results['cpu']['labels']);b=np.asarray(results['metal']['labels']);table=np.zeros((int(a.max())+1,int(b.max())+1),dtype=np.int64);np.add.at(table,(a,b),1);ri,ci=linear_sum_assignment(table,maximize=True);agree=float(table[ri,ci].sum()/n);ari=float(adjusted_rand_score(a,b));typeari={m:float(adjusted_rand_score(types,res['labels'])) for m,res in results.items()}
 comparisons.append(dict(seed=seed,partitionARI=ari,matchedAssignmentAgreement=agree,matchedChangedCells=int(n-table[ri,ci].sum()),sourceTypeARI=typeari,clusters={m:len(res['clusterSizes']) for m,res in results.items()},passes=ari>=g['minimumPartitionARI'] and agree>=g['minimumMatchedAssignmentAgreement'] and typeari['cpu']-typeari['metal']<=g['maximumCellTypeARIRegression']))
recall=[]
for a,b in zip(annotation['cpu']['perType'],annotation['metal']['perType']):
 assert a['type']==b['type'];loss=a['neighborTypeRecall']-b['neighborTypeRecall'];conditionLoss=None if a['conditionBalancedAccuracy'] is None else a['conditionBalancedAccuracy']-b['conditionBalancedAccuracy'];recall.append(dict(type=a['type'],recallLoss=loss,conditionBalancedAccuracyLoss=conditionLoss,passes=loss<=g['maximumPerTypeNeighborRecallRegression'] and (conditionLoss is None or conditionLoss<=g['maximumConditionBalancedAccuracyRegression'])))
conditionloss=annotation['cpu']['conditionBalancedAccuracy']-annotation['metal']['conditionBalancedAccuracy'];result=dict(status='PASS' if all(x['passes'] for x in comparisons+recall) and conditionloss<=g['maximumConditionBalancedAccuracyRegression'] else 'FAIL',scope=p['scope'],protocolSHA256=hashlib.sha256((r/'protocol.json').read_bytes()).hexdigest(),cells=n,sourceTypes=type_names.tolist(),sourceConditions=condition_names.tolist(),sourceDonors=donor_names.tolist(),missingNeighbors=missing,neighborRecall=1-missing/(n*19),rowsWithOrderingOrMembershipChange=changedRows,cpuEdges=len(keys['cpu']),metalEdges=len(keys['metal']),commonEdges=len(common),maximumCommonEdgeWeightDifference=float(abs(edges['cpu']['value'][ic]-edges['metal']['value'][im]).max()),annotationConcordance=annotation,perTypeComparisons=recall,conditionBalancedAccuracyLoss=conditionloss,clusterComparisons=comparisons,biologicalGeneralizationQualified=False)
(r/'verification.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k not in ['annotationConcordance','perTypeComparisons']}))
