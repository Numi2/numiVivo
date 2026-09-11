"""Retain full diagnostic arrays, execution records and frozen model moments."""
import argparse,gzip,hashlib,json,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
a.out.mkdir(parents=True,exist_ok=False)
def sha(b):return hashlib.sha256(b).hexdigest()
for name in ['protocol.json','diagnose.py','environment.json','repeat-check.json','run.log','repeat.log']:
 shutil.copyfile(a.study/name,a.out/name)
shutil.copytree(a.study/'results',a.out/'results')
shutil.copyfile(a.study/'repeat/results/pipeline.json',a.out/'repeat-pipeline.json')
parent=Path('/Users/home/numivivo-gse181897-20260911')
index=json.loads((parent/'prediction-execution/retention.json').read_bytes())
models={}
for origin in ['Kang','HIRISA']:
 def read(bundle,name):
  h=index['bundles'][bundle][name]['SHA256'];p=parent/'prediction-execution/objects'/(h+'.gz')
  assert sha(p.read_bytes())==index['objects'][h]['compressedSHA256']
  b=gzip.decompress(p.read_bytes());assert sha(b)==h
  return h,json.loads(b)
 model_sha,model=read('model-'+origin,'model.json')
 report_sha,report=read('prediction-'+origin+'-00','report.json')
 critical=report['predictions'][0]['meanResponsePredictiveInterval']['studentCriticalValue']
 models[origin]=dict(modelSHA256=model_sha,predictionSHA256=report_sha,featureIDs=model['featureIDs'],
  meanResponse=model['meanResponse'],donorResponseVariances=model['donorResponseVariances'],
  donors=len(model['trainingDonors']),studentCriticalValue=critical)
raw=(json.dumps(models,sort_keys=True,separators=(',',':'))+'\n').encode()
(a.out/'model-moments.json.gz').write_bytes(gzip.compress(raw,compresslevel=9,mtime=0))
files={str(p.relative_to(a.out)):dict(bytes=p.stat().st_size,SHA256=sha(p.read_bytes())) for p in sorted(a.out.rglob('*')) if p.is_file()}
manifest=dict(schemaVersion=1,files=files,retrospective=True,predictionModelChanged=False,
 parentPredictionCommit='2a294637b986f7be52a08ec4a1d19efe27538715',parentStudy=str(parent),
 parentPredictionFreezeSHA256=sha((parent/'prediction-execution/prediction-freeze.json').read_bytes()),
 scope='Full diagnostic arrays and model moments are included. Rerunning cell bootstrap requires the original qualified parent source H5AD and count/role artifacts; no new independent prediction or uncertainty calibration is claimed.')
(a.out/'manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
print(json.dumps(dict(files=len(files),bytes=sum(x['bytes'] for x in files.values()))))
