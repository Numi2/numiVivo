"""Apply the frozen grid-initialization sensitivity criterion, without outcome scoring."""
from pathlib import Path
import gzip,json,sys,collections
root=Path(sys.argv[1]);protocol=json.loads((root/'initialization-protocol.json').read_bytes());result={}
for origin in ['Kang','HIRISA']:
 def read(directory):
  fits={};queries={}
  for line in gzip.open(directory/(origin+'-adaptive.jsonl.gz'),'rt'):
   row=json.loads(line);feature=row['modelID'].split(':grid=')[0]
   if row['kind']=='fit':fits[feature]=row
   else:queries[(feature,row['queryID'])]=row
  return fits,queries
 pf,pq=read(root);af,aq=read(root/'initialization17');records=[]
 assert set(pf)==set(af) and set(pq)==set(aq)
 for feature in sorted(pf):
  ps=pf[feature]['status'];ass=af[feature]['status']
  if ps!='boundedContinuousLikelihood' or ass!='boundedContinuousLikelihood':records.append({'feature':feature,'status':'unavailableOrUnconverged','primaryStatus':ps,'alternateStatus':ass});continue
  worst=None;maximum=0;count=0
  for key,p in pq.items():
   if key[0]!=feature:continue
   a=aq[key];assert p['status']==a['status']=='conditionalPrediction';count+=1
   for field in protocol['fields']:
    section='moments' if field in p['prediction']['moments'] else 'plannedTreatedCountMoments';x=p['prediction'][section][field];y=a['prediction'][section][field];e=abs(x-y)/max(1,abs(x))
    if e>maximum:maximum=e;worst={'queryID':key[1],'field':field,'initialGrid9':x,'initialGrid17':y}
  pm=pf[feature]['model']['model'];am=af[feature]['model']['model'];ll=(am['relativeLogLikelihood']-pm['relativeLogLikelihood'])/len(pm['trainingDonorIDs'])
  records.append({'feature':feature,'status':'passed' if maximum<=protocol['primaryToleranceMaximumScaledQueryMomentChange'] else 'failedInitializationSensitivity','queries':count,'maximumScaledQueryMomentChange':maximum,'meanLogLikelihoodChange':ll,'worst':worst})
 result[origin]={'states':dict(collections.Counter(r['status'] for r in records)),'records':records}
result['qualification']='Numerical sensitivity at fixed fitted cell dispersions and on the same known development data; this is not statistical parameter uncertainty or biological validation.'
(root/'initialization-comparison.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result))
