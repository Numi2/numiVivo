"""Compare declared endpoints with the previous grid65; this is not outcome scoring."""
from pathlib import Path
import gzip,json,sys,collections
root=Path(sys.argv[1]);parent=Path(sys.argv[2]);result={}
for origin in ['Kang','HIRISA']:
 previous={}
 with gzip.open(parent/(origin+'-native.jsonl.gz'),'rt') as f:
  for line in f:
   row=json.loads(line)
   if row['kind']=='query' and row['modelID'].endswith('grid=65') and 'prediction' in row:previous[(row['modelID'].split(':grid=')[0],row['queryID'])]=row['prediction']
 changes={};states=collections.Counter()
 with gzip.open(root/(origin+'-adaptive.jsonl.gz'),'rt') as f:
  for line in f:
   row=json.loads(line)
   if row['kind']!='query':continue
   states[row['status']]+=1
   if 'prediction' not in row:continue
   feature=row['modelID'].split(':grid=')[0];old=previous[(feature,row['queryID'])];new=row['prediction'];entry=changes.setdefault(feature,{'queries':0,'maximumScaledChange':0.0});entry['queries']+=1
   for section in ['moments','plannedTreatedCountMoments']:
    for field,a in new[section].items():
     b=old[section][field];error=abs(a-b)/max(1,abs(a))
     if error>entry['maximumScaledChange']:entry.update(maximumScaledChange=error,worst={'queryID':row['queryID'],'field':field,'adaptive':a,'grid65':b})
 result[origin]={'states':dict(states),'changes':changes}
result['qualification']='Comparison of fitted likelihood models on known development controls, not held-out treated outcomes or proof of uniquely identified/stable posterior moments.'
(root/'grid65-comparison.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result))
