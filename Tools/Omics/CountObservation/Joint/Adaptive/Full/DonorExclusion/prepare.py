"""Freeze training-only inputs from already qualified original-cell moments."""
from pathlib import Path
import gzip,hashlib,json,sys,time
root=Path(sys.argv[1]);parent=Path(sys.argv[2]);root.mkdir(exist_ok=True);(root/'inputs').mkdir(exist_ok=True)
protocol={'baseCommit':'217614c515e0666cf795ffb9d94b31a49789129e','phase':'Train-only cell-dispersion extraction for subsequent donor-excluded joint fitting','allGenes':True,'queryTreatedCountsExcluded':True,'jointLikelihoodFits':'not yet run','uncertainty':'Plug-in dispersions only; no parameter uncertainty integrated','folds':[],'createdUnix':time.time()}
for origin in ['Kang','HIRISA']:
 p=parent/(origin+'-input.json.gz');compressed=p.read_bytes();raw=gzip.decompress(compressed);source=json.loads(raw);donors=sorted({g['donorID'] for g in source['groups']});features=source['featureIDs'];assert len(donors)==(8 if origin=='Kang' else 5)
 (root/(origin+'-parent-input.json.gz')).write_bytes(compressed)
 for index,excluded in enumerate(donors):
  groups=[g for g in source['groups'] if g['donorID']!=excluded];assert len(groups)==2*(len(donors)-1)
  payload={k:source[k] for k in ['origin','featureIDs','controlConditionID','treatedConditionID']};payload.update(excludedDonorID=excluded,groups=groups,parentInputSHA256=hashlib.sha256(raw).hexdigest());data=json.dumps(payload,sort_keys=True,separators=(',',':'),allow_nan=False).encode();tag=f'{origin}-{index:02d}';path=root/'inputs'/(tag+'.json.gz');assert not path.exists();path.write_bytes(gzip.compress(data,mtime=0))
  protocol['folds'].append({'tag':tag,'origin':origin,'excludedDonorID':excluded,'trainingDonorIDs':[d for d in donors if d!=excluded],'featureCount':len(features),'trainingCells':sum(g['moments']['cells'] for g in groups),'inputSHA256':hashlib.sha256(data).hexdigest(),'compressedInputSHA256':hashlib.sha256(path.read_bytes()).hexdigest(),'parentInputSHA256':hashlib.sha256(raw).hexdigest()})
(root/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');print('Frozen',len(protocol['folds']),'training-only full-gene calibration inputs')
