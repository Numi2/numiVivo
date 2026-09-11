#!/usr/bin/env python3
"""Trace original labels, complete native programs and frozen gene-space margins."""
import argparse,hashlib,json,sys,tarfile
from pathlib import Path
import anndata as ad,h5py,numpy as np
from scipy.special import softmax
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'Logistic'))
from prepare import sha,write
from bundles import read_json

PROGRAMS=[['HBB','HBA1','HBA2'],['PF4','PPBP','ITGA2B','GP9','TUBB1']]

def main():
 p=argparse.ArgumentParser();p.add_argument('--prior',type=Path,required=True);p.add_argument('--programs',type=Path,required=True);p.add_argument('--original-kang',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 freeze=json.loads((a.prior/'native-qualified/prediction-freeze.json').read_text());assert freeze['status']=='completed'
 programExecution=json.loads((a.programs/'execution.json').read_text());assert programExecution['status']=='completed' and all(c['exitCode']==0 for c in programExecution['commands'])
 cohorts={};programResults=[]
 for cohort,file in [('kang','kang.h5ad'),('ding','ding-query.h5ad')]:
  source=a.prior/'inputs'/file;data=ad.read_h5ad(source);ids=np.array(data.var_names if cohort=='kang' else data.var['reference_feature_id'],dtype=str);labels=np.asarray(data.obs['cell_type' if cohort=='kang' else 'sourceCellType'].astype(str));assert len(set(ids))==len(ids)
  x=data.X.astype(np.float64).tocsr();totals=np.asarray(x.sum(axis=1)).ravel();x=x.multiply((10000/totals)[:,None]).tocsr();x.data=np.log1p(x.data)
  root=a.programs/cohort;assert sha(root/'original.h5ad')==sha(source)
  metadata=json.loads((root/'metadata.json').read_text());np.testing.assert_array_equal(data.obs_names,[c['barcode'] for c in metadata['cells']]);np.testing.assert_array_equal(ids,[f['id'] for f in metadata['features']]);np.testing.assert_array_equal(data.obs['native_sample'].astype(str),[c['sampleID'] for c in metadata['cells']])
  def read(name,kind,columns):
   r=np.fromfile(root/name,dtype=np.dtype([('row','<u4'),('column','<u4'),('value',kind)]));assert len(r)==len(data)*columns
   np.testing.assert_array_equal(r['row'],np.repeat(np.arange(len(data)),columns));np.testing.assert_array_equal(r['column'],np.tile(np.arange(columns),len(data)));return r['value'].reshape(len(data),columns)
  scores=read('scores.bin','<f8',2);detected=read('detected-members.bin','<u8',2);nativeTotals=read('total-counts.bin','<u8',1).ravel();np.testing.assert_array_equal(nativeTotals,totals)
  direct=[];lookup={g:i for i,g in enumerate(ids)}
  for i,genes in enumerate(PROGRAMS):
   indices=[lookup[g] for g in genes];values=np.asarray(x[:,indices].mean(axis=1)).ravel();direct.append(values);np.testing.assert_allclose(scores[:,i],values,rtol=1e-12,atol=1e-12);np.testing.assert_array_equal(detected[:,i],np.asarray((data.X[:,indices]>0).sum(axis=1)).ravel())
  summaries=[]
  for label in sorted(set(labels)):
   for library in sorted(set(data.obs['native_sample'].astype(str))):
    mask=(labels==label)&np.asarray(data.obs['native_sample'].astype(str)==library)
    if not mask.any():continue
    summaries.append(dict(sourceLabel=label,library=library,cells=int(mask.sum()),libraryTotalQuantiles=np.quantile(totals[mask],[0,.25,.5,.75,1]).tolist(),scoreMeans=scores[mask].mean(axis=0).tolist(),scoreQuantiles=np.quantile(scores[mask],[0,.25,.5,.75,1],axis=0).tolist(),detectedAnyFractions=np.mean(detected[mask]>0,axis=0).tolist()))
  mega=labels==('Megakaryocytes' if cohort=='kang' else 'Megakaryocyte')
  record=dict(cohort=cohort,cells=len(data),sourceSHA256=sha(source),nativeScoreSHA256=sha(root/'scores.bin'),nativeScoresMaximumError=float(np.max(np.abs(scores-np.column_stack(direct)))),allDetectionsAndTotalsExact=True,sourceMegaCells=int(mega.sum()),sourceMegaProgramMeans=scores[mega].mean(axis=0).tolist(),sourceMegaDetectionFractions=np.mean(detected[mega]>0,axis=0).tolist(),allLabelLibrarySummaries=summaries)
  programResults.append(record);cohorts[cohort]=dict(data=data,ids=ids,labels=labels,x=x,totals=totals,mega=mega)
  np.savez_compressed(a.out/(cohort+'-all-cell-programs.npz'),barcodes=np.asarray(data.obs_names,dtype=str),sampleIDs=np.asarray(data.obs['native_sample'],dtype=str),sourceLabels=labels,scores=scores,detected=detected,totals=nativeTotals)
 # Decode the original compound dataframe directly, independently of AnnData.
 assert sha(a.original_kang)=='e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830'
 with h5py.File(a.original_kang,'r') as f:
  obs=f['obs'][()];codes=obs['cell_type'];categories=np.asarray(f['uns/cell_type_categories'].asstr()[()]);assert np.all((codes>=0)&(codes<len(categories)))
  barcodes=np.array([v.decode() if isinstance(v,bytes) else str(v) for v in obs['index']]);np.testing.assert_array_equal(barcodes,cohorts['kang']['data'].obs_names);np.testing.assert_array_equal(categories[codes],cohorts['kang']['labels'])
 sourceCheck=dict(originalSHA256=sha(a.original_kang),preparedSHA256=programResults[0]['sourceSHA256'],all24673LabelsAndBarcodesExact=True,categoryCodesInRange=True,originalCategories=categories.tolist(),originalMegaCode=int(np.flatnonzero(categories=='Megakaryocytes')[0]))
 marginResults=[]
 for item in freeze['archives']:
  direction=item['direction'];train,query=direction.split('-to-');archive=a.prior/'native-qualified'/item['path'];assert sha(archive)==item['SHA256']
  with tarfile.open(archive,'r:gz') as t:
   members=json.load(t.extractfile('members.json'));model=read_json(t,members,'reference/model.json');report=read_json(t,members,'mapped/report.json')
  classes=np.array(model['classes']);target='Megakaryocytes' if train=='kang' else 'Megakaryocyte';rival='CD4 T cells' if train=='kang' else 'CD4+ T cell';ti=int(np.flatnonzero(classes==target)[0]);ri=int(np.flatnonzero(classes==rival)[0]);params=np.array(model['logistic']['parameters']);scales=np.array(model['logistic']['scales']);means=np.array(model['logistic']['means']);v=np.array(model['reduction']['loadings']);centers=np.array(model['reduction']['projectionCenters']);selected=np.array(model['reduction']['selectedFeatureIndices']);features=np.array(model['featureIDs'])[selected]
  effective=(params[:,:-1]/scales)@v.T;offset=params[:,-1]-(means/scales)@params[:,:-1].T-centers@effective.T
  trainScores=np.array(model['reduction']['scores']);trainLogits=np.c_[(trainScores-means)/scales,np.ones(len(trainScores))]@params.T;labels=np.array(model['referenceLabels']);trainingMask=labels==target
  queryScores=np.array([c['scores'] for c in report['cells']]);queryLogits=np.c_[(queryScores-means)/scales,np.ones(len(queryScores))]@params.T
  raw=cohorts[query];lookup={g:i for i,g in enumerate(raw['ids'])};xq=raw['x'][:,[lookup[g] for g in features]];affine=np.asarray(xq@effective.T)+offset;np.testing.assert_allclose(affine,queryLogits,rtol=1e-11,atol=1e-11);np.testing.assert_allclose(softmax(affine,axis=1),[c['classProbabilities'] for c in report['cells']],rtol=1e-11,atol=1e-11)
  # Ding training uses every available assigned label; that row mask is explicit.
  tr=cohorts[train];tm=tr['mega'];assert int(tm.sum())==int(trainingMask.sum());trainLookup={g:i for i,g in enumerate(tr['ids'])};xt=tr['x'][:,[trainLookup[g] for g in features]]
  trainingMean=np.asarray(xt[tm].mean(axis=0)).ravel();queryMean=np.asarray(xq[raw['mega']].mean(axis=0)).ravel();weights=effective[ti]-effective[ri];contributions=(queryMean-trainingMean)*weights
  trainMargin=float((trainLogits[trainingMask,ti]-trainLogits[trainingMask,ri]).mean());queryMargin=float((queryLogits[raw['mega'],ti]-queryLogits[raw['mega'],ri]).mean());assert abs(contributions.sum()-(queryMargin-trainMargin))<1e-10
  record=dict(direction=direction,archiveSHA256=item['SHA256'],trainingSourceClass=target,fixedRival=rival,trainingSourceClassCells=int(trainingMask.sum()),trainingSourceClassRecall=float(np.mean(classes[trainLogits.argmax(axis=1)][trainingMask]==target)),querySourceClassCells=int(raw['mega'].sum()),geneAffineMaximumError=float(np.abs(affine-queryLogits).max()),trainingMeanMargin=trainMargin,queryMeanMargin=queryMargin,marginDifference=queryMargin-trainMargin,geneContributionSum=float(contributions.sum()),programGenesSelected={g:bool(g in features) for g in sum(PROGRAMS,[])},largestNegativeContributions=[dict(gene=str(features[j]),trainingMeanLog=float(trainingMean[j]),queryMeanLog=float(queryMean[j]),classContrastWeight=float(weights[j]),contribution=float(contributions[j])) for j in np.argsort(contributions)[:20]],largestPositiveContributions=[dict(gene=str(features[j]),contribution=float(contributions[j])) for j in np.argsort(contributions)[-20:][::-1]])
  marginResults.append(record);np.savez_compressed(a.out/(direction+'-gene-margin.npz'),featureIDs=features,geneWeights=weights,trainingSourceClassMeanLog=trainingMean,querySourceClassMeanLog=queryMean,contributions=contributions,queryCellMargins=queryLogits[:,ti]-queryLogits[:,ri],barcodes=np.asarray(raw['data'].obs_names,dtype=str))
 result=dict(status='passed-diagnostic-reconstruction',checkerSHA256=sha(__file__),priorPredictionFreezeSHA256=sha(a.prior/'native-qualified/prediction-freeze.json'),nativeProgramFreezeSHA256=sha(a.programs/'input-freeze.json'),sourceLabels=sourceCheck,programs=programResults,marginDecompositions=marginResults,scope='Post-result source-label and frozen-model diagnosis; not relabelling, a repaired predictor, source-class biological equivalence or a revised success gate.')
 write(a.out/'summary.json',result);print(json.dumps(dict(status=result['status'],cells=sum(r['cells'] for r in programResults),marginDirections=len(marginResults))),flush=True)
if __name__=='__main__':main()
