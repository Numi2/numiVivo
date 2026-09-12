from pathlib import Path
import json,hashlib
import numpy as np
s=Path(__file__).parent;source=Path('/Users/home/numivivo-cross-study-ifnb-20260911/inputs');freeze=json.loads((source/'input-freeze.json').read_text());pred=json.loads((s/'prediction-freeze.json').read_text());assert json.loads((s/'verification.json').read_text())['status']=='pass-all-models-all-penalties-all-predictions';panel=json.loads((source/'panel.json').read_text());folds=json.loads((source/'folds.json').read_text());cohorts={};counts={};ids={};logs={};order={}
for study in ['Kang','HIRISA']:
 for suffix in ['-cohort.json','-source-counts.npz']:
  name=study+suffix;assert hashlib.sha256((source/name).read_bytes()).hexdigest()==freeze['files'][name]
 cohorts[study]=json.loads((source/(study+'-cohort.json')).read_text())
 with np.load(source/(study+'-source-counts.npz')) as a:counts[study]=a['counts'].copy();ids[study]=a['featureIDs'].tolist()
 logs[study]=np.log1p(counts[study].astype('f8')/counts[study].sum(axis=1,dtype=np.uint64)[:,None]*1e6);order[study]=np.array([ids[study].index(p) for p in panel])
rows=[]
for f in folds:
 tag=f['id'];p=s/'native'/(tag+'-output.json');assert hashlib.sha256(p.read_bytes()).hexdigest()==pred['folds'][tag]['outputSHA256'];b=json.loads(p.read_text());a=json.loads((s/'native'/(tag+'-input.json')).read_text());tr=f['trainingStudy'];qu=f['queryStudy'];samples=cohorts[tr]['samples'];ds=sorted(set(samples[i]['donorID'] for i in f['trainingIndices']));assert ds==a['trainingDonorIDs']
 for condition,key in [('control','controls'),('IFNB','treated')]:
  positions=[]
  for d in ds:
   matches=[i for i in f['trainingIndices'] if samples[i]['donorID']==d and samples[i]['condition']==condition];assert len(matches)==1;positions+=matches
  np.testing.assert_allclose(a[key],logs[tr][positions][:,order[tr]],rtol=1e-12,atol=1e-12)
 q=logs[qu][f['queryIndex'],order[qu]];np.testing.assert_allclose(a['queryControl'],q,rtol=1e-12,atol=1e-12);truth=logs[qu][f['scoringIndex'],order[qu]];mean=np.mean(np.array(a['treated'])-a['controls'],axis=0);values={'candidate':b['predictedTreated'],'trainingMean':np.maximum(0,q+mean),'noChange':q};rmse={k:float(np.sqrt(np.mean((np.asarray(v)-truth)**2))) for k,v in values.items()};rows.append({'fold':tag,'mode':f['mode'],'trainingStudy':tr,'queryStudy':qu,'donor':f['heldOutDonor'],'features':len(panel),'rmse':rmse})
summaries=[]
for qu in ['Kang','HIRISA']:
 for mode in ['cross','within']:
  rr=[r for r in rows if r['queryStudy']==qu and r['mode']==mode];rmse={k:float(np.mean([r['rmse'][k] for r in rr])) for k in ['candidate','trainingMean','noChange']};gain=100*(1-rmse['candidate']/rmse['trainingMean']);summaries.append({'queryStudy':qu,'mode':mode,'donors':len(rr),'equalDonorMeanRMSE':rmse,'gainOverMeanPercent':gain,'strictlyBetterThanBoth':rmse['candidate']<min(rmse['trainingMean'],rmse['noChange']),'meanGainAtLeastFivePercent':gain>=5,'donorsWorseThanMean':[r['donor'] for r in rr if r['rmse']['candidate']>r['rmse']['trainingMean']]})
r={'scope':'Complete shared-panel, reused-study development evaluation. Native algorithm fit/prediction; independent normalization, all model/CV/prediction checks. No fresh validation or calibrated uncertainty.','predictionFreezeSHA256':hashlib.sha256((s/'prediction-freeze.json').read_bytes()).hexdigest(),'folds':rows,'summaries':summaries};(s/'scores.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(summaries,indent=2))
