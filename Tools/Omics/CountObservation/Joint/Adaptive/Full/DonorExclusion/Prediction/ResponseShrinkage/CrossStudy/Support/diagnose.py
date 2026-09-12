from pathlib import Path
import json,hashlib,gzip
import numpy as np
s=Path(__file__).parent;root=Path('/Users/home/numivivo-shrinkage-transfer-20260912');freeze=json.loads((root/'prediction-freeze.json').read_text());result={'scope':'Post-hoc diagnostic of frozen development predictions. Outside observed donor range is descriptive, not a calibrated OOD probability, prediction-error gate or causal attribution. No prediction is changed or excluded.','sourceFreezeSHA256':hashlib.sha256((root/'prediction-freeze.json').read_bytes()).hexdigest(),'folds':{}}
for tag,record in freeze['folds'].items():
 ps=[root/'native'/(tag+'-'+suffix+'.json') for suffix in ['input','output']]
 for p,key in zip(ps,['inputSHA256','outputSHA256']):assert hashlib.sha256(p.read_bytes()).hexdigest()==record[key]
 a,b=[json.loads(p.read_text()) for p in ps];m=b['model'];x=np.array(a['controls']);q=np.array(a['queryControl']);c=np.array(m['controlMean']);mu=np.array(m['responseMean']);slope=np.array(m['responseSlope']);pred=np.array(b['predictedTreated']);lo=x.min(0);hi=x.max(0);below=q<lo;above=q>hi;outside=below|above;correction=slope*(q-c);base=np.maximum(0,q+mu)
 np.testing.assert_allclose(pred,np.maximum(0,q+mu+correction),rtol=1e-10,atol=1e-10)
 stats={'features':len(q),'belowTrainingRange':int(below.sum()),'aboveTrainingRange':int(above.sum()),'outsideTrainingRangePercent':float(outside.mean()*100),'zeroWidthTrainingRange':int((hi==lo).sum()),'queryMinusTrainingMeanRMS':float(np.sqrt(np.mean((q-c)**2))),'slopeCorrectionRMS':float(np.sqrt(np.mean(correction**2))),'predictionMinusMeanBaselineRMS':float(np.sqrt(np.mean((pred-base)**2))),'slopeMedian':float(np.median(slope)),'slopeAbsoluteMedian':float(np.median(np.abs(slope))),'clippedPredictions':int((q+mu+correction<0).sum()),'strata':{}}
 for name,mask in [('inside',~outside),('outside',outside)]:
  stats['strata'][name]={'features':int(mask.sum()),'correctionSquaredSum':float(np.sum(correction[mask]**2))}
 result['folds'][tag]=stats
 with gzip.open(s/(tag+'-support.jsonl.gz'),'wt') as f:
  for i,fid in enumerate(a['featureIDs']):f.write(json.dumps({'featureID':fid,'trainingMinimum':float(lo[i]),'trainingMaximum':float(hi[i]),'queryControl':float(q[i]),'outsideObservedTrainingRange':bool(outside[i]),'slopeCorrection':float(correction[i]),'meanBaseline':float(base[i]),'frozenPrediction':float(pred[i])})+'\n')
 stats['supportSHA256']=hashlib.sha256((s/(tag+'-support.jsonl.gz')).read_bytes()).hexdigest()
result['status']='complete';(s/'control-only-diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
