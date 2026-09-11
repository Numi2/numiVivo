"""Apply frozen training calibration to existing original control-cell vectors."""
from pathlib import Path
import gzip,hashlib,json,sys
root=Path(sys.argv[1]);old=Path('/Users/home/numivivo-count-observation-20260911/inputs')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
freeze=json.loads((old/'freeze.json').read_bytes());assert sha(old/'cases.json')==freeze['files']['cases.json']
queries=[]
for c in json.loads((old/'cases.json').read_bytes()):
 if c['cellDispersion']!=0:continue
 donor,feature=c['id'].split(':symbol|');feature='symbol|'+feature.split(':phi=')[0]
 queries.append(dict(id=donor+':'+feature,featureID=feature,donorID='GSE181897:exp_id:'+donor,counts=c['counts'],libraryCounts=c['libraryCounts'],plannedLibraryCounts=c['plannedLibraryCounts']))
assert len(queries)==62*16
models=[];sources={}
for origin in ['Kang','HIRISA']:
 p=root/(origin+'-native.json.gz');sources[origin]=sha(p)
 assert json.loads((root/(origin+'-verification.json')).read_bytes())['status']=='passed'
 model=next(m for m in json.loads(gzip.decompress(p.read_bytes()))['models'] if m['conditionID']=='control');models.append(model)
raw=json.dumps({'models':models,'queries':queries},sort_keys=True,separators=(',',':')).encode()
with (root/'query-input.json.gz').open('wb') as f,gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as g:g.write(raw)
(root/'query-freeze.json').write_text(json.dumps({'trainingReportsSHA256':sources,'originalControlCaseSHA256':sha(old/'cases.json'),'queryInputSHA256':hashlib.sha256(raw).hexdigest(),'compressedQueryInputSHA256':sha(root/'query-input.json.gz'),'queriesPerTrainingOrigin':len(queries),'preparedBeforeExecution':True,'queryTreatedCountValuesUsed':False,'endpoint':'Latent control rate and future SAME-condition count moments; no treated response','plannedSampling':'Original frozen count-observation protocol: 20 cells, 1000 RNA UMIs per cell'},indent=2)+'\n')
