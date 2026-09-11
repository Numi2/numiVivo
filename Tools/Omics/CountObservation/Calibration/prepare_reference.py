"""Retain unavailable queries and construct independent posterior check inputs."""
from pathlib import Path
import collections,gzip,hashlib,json,sys
root=Path(sys.argv[1]);q=json.loads(gzip.decompress((root/'query-input.json.gz').read_bytes()));rows=[json.loads(line) for line in gzip.decompress((root/'applied.jsonl.gz').read_bytes()).splitlines()]
queries={x['id']:x for x in q['queries']};assert len(rows)==len(q['queries'])*2
cases=[];available=[];seen=set();summary=[]
for row in rows:
 key=(row['model'],row['id']);assert key not in seen;seen.add(key)
 assert row['status']!='error',row
 query=queries[row['id']];model=q['models'][row['model']];feature=next(f for f in model['features'] if f['featureID']==query['featureID'])
 if row['status']=='unavailableCalibration':
  assert row['failure']==feature['status'] and feature['status']!='availableConditionalMomentCalibration';continue
 assert row['status']=='availableConditionalPosterior' and feature['status']=='availableConditionalMomentCalibration'
 cases.append(dict(id=row['id'],counts=query['counts'],libraryCounts=query['libraryCounts'],plannedLibraryCounts=query['plannedLibraryCounts'],cellDispersion=feature['cellDispersion'],gammaPriorShape=feature['gammaPriorShape'],gammaPriorRatePerCPM=feature['gammaPriorRatePerCPM']))
 available.append(row)
for i,origin in enumerate(['Kang','HIRISA']):
 selected=[r for r in rows if r['model']==i];supported=[r for r in selected if r['status']=='availableConditionalPosterior']
 summary.append({'origin':origin,'queries':len(selected),'states':dict(collections.Counter(r['status'] for r in selected)),'unavailableReasons':dict(collections.Counter(r['failure'] for r in selected if 'failure' in r)),'availableZeroCountQueries':sum(r['posterior']['geneCounts']==0 for r in supported),'allAvailableZeroCountVariancesPositive':all(r['posterior']['varianceCPM']>0 for r in supported if r['posterior']['geneCounts']==0)})
(root/'posterior-cases.json').write_text(json.dumps(cases,separators=(',',':'))+'\n')
(root/'posterior-native.jsonl').write_text(''.join(json.dumps(r,separators=(',',':'))+'\n' for r in available))
(root/'application-summary.json').write_text(json.dumps({'origins':summary,'totalQueries':len(rows),'supportedQueries':len(available),'queryTreatedCountsUsed':False,'perturbationResponseApplied':False},indent=2)+'\n')
print(json.dumps(summary))
