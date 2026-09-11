"""Bind each verified training-only calibration to its original sparse groups.

Cache files are hard-linked read-only inputs, not copied or rewritten. Only the
explicit training-group indices are admitted to the joint fitting stream.
"""
from pathlib import Path
import copy,gzip,hashlib,json,os,sys,time
import h5py
root,cal,full=map(Path,sys.argv[1:4]);root.mkdir(exist_ok=False)
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1<<20):h.update(b)
 return h.hexdigest()
protocol=read(cal/'protocol.json');checked=read(cal/'verification.json');assert checked['status']=='passed-all-training-only-calibrations'
parents={}
for origin in ['Kang','HIRISA']:
 m=read(full/(origin+'-manifest.json.gz'));assert sha(full/(origin+'-cells.h5'))==m['cacheSHA256'];parents[origin]=m
plan={'status':'prepared-all-training-only-joint-folds','phase':'Joint fitting only; predictions and held-out outcome scoring remain subsequent stages','ownerCommit':'fa01e5f5b85e761ce3006da442912f53c630826d','createdUnix':time.time(),'folds':[],'noFullCohortWarmStarts':True,'noHeldOutTreatedInputs':True}
for fold in protocol['folds']:
 tag=fold['tag'];origin=fold['origin'];parent=parents[origin];native=cal/'native'/(tag+'.json.gz');receipt=read(cal/'native'/(tag+'-receipt.json'));assert sha(native)==receipt['compressedOutputSHA256'];model=read(native);assert model['inputSHA256']==fold['inputSHA256'];directory=root/tag;directory.mkdir();m=copy.deepcopy(parent);m['excludedDonorID']=fold['excludedDonorID'];m['trainingCalibrationSHA256']=sha(native);m['parentManifestSHA256']=sha(full/(origin+'-manifest.json.gz'));m['groups']=[g for g in parent['groups'] if g['donorID']!=fold['excludedDonorID']];m['qualification']='Training-only donor-excluded joint fit inputs. Parent source metadata and sparse cache retain the original source cohort; only the explicit groups below enter fitting. No full-cohort phi or fitted weights are reused.'
 assert sorted({g['donorID'] for g in m['groups']})==fold['trainingDonorIDs']==model['control']['trainingDonorIDs']==model['treated']['trainingDonorIDs']
 m['genes']=[]
 for j,(c,t) in enumerate(zip(model['control']['features'],model['treated']['features'])):
  assert parent['genes'][j]['featureID']==c['featureID']==t['featureID'];m['genes'].append({'featureIndex':j,'featureID':c['featureID'],'controlCellDispersion':c.get('cellDispersion'),'treatedCellDispersion':t.get('cellDispersion'),'controlCalibrationStatus':c['status'],'treatedCalibrationStatus':t['status']})
 assert len(m['genes'])==fold['featureCount'];cache=directory/(origin+'-cells.h5');os.link(full/(origin+'-cells.h5'),cache)
 with h5py.File(cache,'r') as h:
  for g in m['groups']:
   x=h[str(g['index'])];assert x.attrs['donorID']==g['donorID'] and x.attrs['conditionID']==g['conditionID']
   assert len(x['sourceRows'])==g['cells'] and sha(native)==m['trainingCalibrationSHA256']
 raw=json.dumps(m,sort_keys=True,separators=(',',':'),allow_nan=False).encode();(directory/(origin+'-manifest.json.gz')).write_bytes(gzip.compress(raw,mtime=0));available=sum(g['controlCellDispersion'] is not None and g['treatedCellDispersion'] is not None for g in m['genes']);prepared={'status':'completed','origin':origin,'genes':len(m['genes']),'sourceCells':sum(g['cells'] for g in m['groups']),'availablePairedGenes':available,'cacheSHA256':m['cacheSHA256'],'manifestSHA256':hashlib.sha256(raw).hexdigest(),'excludedDonorID':m['excludedDonorID'],'cacheReusedWithoutCopy':True};(directory/(origin+'-prepare-state.json')).write_text(json.dumps(prepared,indent=2)+'\n');(directory/'runtime-freeze.json').write_bytes((full/'gap-repair/runtime-freeze.json').read_bytes());plan['folds'].append(dict(fold,root=str(directory),availablePairedGenes=available,manifestSHA256=prepared['manifestSHA256'],cacheIndices=[g['index'] for g in m['groups']]))
(root/'protocol.json').write_text(json.dumps(plan,indent=2)+'\n');print(json.dumps({'folds':len(plan['folds']),'availableOriginGeneFolds':sum(f['availablePairedGenes'] for f in plan['folds']),'copiedSparseCountBytes':0}))
