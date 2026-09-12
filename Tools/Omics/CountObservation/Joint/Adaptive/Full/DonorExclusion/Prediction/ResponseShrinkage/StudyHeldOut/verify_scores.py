from pathlib import Path
import json,math,hashlib
import numpy as np
s=Path(__file__).parent;r=json.loads((s/'scores.json').read_text());a=np.load(s/'cohorts.npz');meta=json.loads((s/'cohort-metadata.json').read_text());checked=0;maximum=0
for study,fold in r['folds'].items():
 pred=json.loads((s/'native'/(study+'.json')).read_text());truth=a[study+'_treated'];control=a[study+'_controls'];series={k:[] for k in ['candidate','equalStudyTrainingMean','noChange']}
 for i,donor in enumerate(meta['studies'][study]):
  for key,values in [('candidate',pred['predictedTreated'][i]),('equalStudyTrainingMean',pred['trainingMeanPredictions'][i]),('noChange',control[i])]:
   v=math.sqrt(math.fsum((float(p)-float(t))**2 for p,t in zip(values,truth[i]))/len(values));ref=fold['donorRMSE'][donor][key];assert math.isclose(v,ref,rel_tol=1e-11,abs_tol=1e-12);maximum=max(maximum,abs(v-ref));series[key].append(v);checked+=1
 for k,vs in series.items():assert math.isclose(math.fsum(vs)/len(vs),fold['equalDonorMeanRMSE'][k],rel_tol=1e-11,abs_tol=1e-12)
result={'status':'pass-independent-scalar-all-donor-scores','donorBaselineScores':checked,'maximumDifference':maximum,'scoresSHA256':hashlib.sha256((s/'scores.json').read_bytes()).hexdigest()};(s/'score-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
