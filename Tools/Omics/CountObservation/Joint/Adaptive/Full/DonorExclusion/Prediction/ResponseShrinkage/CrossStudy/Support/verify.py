from pathlib import Path
import json,gzip,hashlib,math
s=Path(__file__).parent;root=Path('/Users/home/numivivo-shrinkage-transfer-20260912');freeze=json.loads((root/'prediction-freeze.json').read_text());d=json.loads((s/'control-only-diagnostic.json').read_text());result={'folds':{}}
for tag,rec in freeze['folds'].items():
 inp=root/'native'/(tag+'-input.json');out=root/'native'/(tag+'-output.json');assert hashlib.sha256(inp.read_bytes()).hexdigest()==rec['inputSHA256'];assert hashlib.sha256(out.read_bytes()).hexdigest()==rec['outputSHA256'];a=json.loads(inp.read_text());b=json.loads(out.read_text());m=b['model'];path=s/(tag+'-support.jsonl.gz');assert hashlib.sha256(path.read_bytes()).hexdigest()==d['folds'][tag]['supportSHA256'];below=above=n=0;corrections=[];shift=[];seen=set()
 with gzip.open(path,'rt') as f:
  for i,line in enumerate(f):
   row=json.loads(line);assert row['featureID']==a['featureIDs'][i] and row['featureID'] not in seen;seen.add(row['featureID']);training=[r[i] for r in a['controls']];lo=min(training);hi=max(training);q=a['queryControl'][i];outside=not lo<=q<=hi
   assert row['trainingMinimum']==lo and row['trainingMaximum']==hi and row['queryControl']==q and row['outsideObservedTrainingRange']==outside
   correction=m['responseSlope'][i]*(q-m['controlMean'][i]);mean=max(0,q+m['responseMean'][i]);assert math.isclose(row['slopeCorrection'],correction,rel_tol=1e-12,abs_tol=1e-12);assert row['meanBaseline']==mean and row['frozenPrediction']==b['predictedTreated'][i]
   corrections.append(correction**2);shift.append((row['frozenPrediction']-mean)**2);below+=q<lo;above+=q>hi;n+=1
 assert n==len(a['featureIDs'])==d['folds'][tag]['features'];assert below==d['folds'][tag]['belowTrainingRange'] and above==d['folds'][tag]['aboveTrainingRange']
 assert math.isclose(math.sqrt(math.fsum(corrections)/n),d['folds'][tag]['slopeCorrectionRMS'],rel_tol=1e-12,abs_tol=1e-12)
 assert math.isclose(math.sqrt(math.fsum(shift)/n),d['folds'][tag]['predictionMinusMeanBaselineRMS'],rel_tol=1e-12,abs_tol=1e-12)
 result['folds'][tag]={'features':n,'below':below,'above':above}
result['status']='pass-scalar-reconstruction-all-support-records';result['records']=sum(v['features'] for v in result['folds'].values());(s/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(result['status'],result['records'])
