#!/usr/bin/env python3
"""Actual CLI panel failure checks and unchanged full-universe historical replay."""
import argparse,copy,hashlib,json,os,shutil,subprocess
from pathlib import Path

def write(p,d):p.write_text(json.dumps(d,sort_keys=True,separators=(',',':')))
def main():
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--historical-fold',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
 src=a.inputs/'within-HIRISA-to-HIRISA-2616BW';plan=json.loads((src/'training.json').read_text());original=copy.deepcopy(plan);del plan['responseFeatureIDs'];write(a.out/'default.json',plan);commands=[]
 def run(args,name,error=None):
  r=subprocess.run([str(a.binary),*map(str,args)],capture_output=True,text=True);(a.out/(name+'.log')).write_text(r.stdout+r.stderr)
  commands.append(dict(arguments=list(map(str,args)),exitCode=r.returncode,expectedError=error));write(a.out/'commands.json',commands)
  if error:assert r.returncode!=0 and error in r.stdout+r.stderr,(name,r.stdout,r.stderr)
  else:assert r.returncode==0,(name,r.stdout,r.stderr)
 run(['singlecell-perturbation-fit',src/'training.h5ad','--plan',a.out/'default.json','--output',a.out/'model'],'default-fit')
 run(['singlecell-perturbation-verify',a.out/'model'],'default-model-verify')
 run(['singlecell-perturbation-predict',src/'query.h5ad','--plan',src/'query.json','--reference',a.out/'model','--output',a.out/'prediction'],'default-predict')
 run(['singlecell-perturbation-prediction-verify',a.out/'prediction'],'default-prediction-verify')
 old=json.loads((a.historical_fold/'model.json').read_text());new=json.loads((a.out/'model/model.json').read_text())
 keys=['method','selectedFeatureIndices','contextCenters','contextScales','contexts','meanResponse','medianResponse','dualCoefficients','maximumSolveResidual','qualification']
 for key in keys:assert old[key]==new[key],key
 old_prediction=json.loads((a.historical_fold/'prediction.json').read_text())['predictions'][0];new_prediction=json.loads((a.out/'prediction/report.json').read_text())['predictions'][0]
 for key in ['control','libraryCounts','estimates']:assert old_prediction[key]==new_prediction[key],key
 for name,ids,error in [('empty',[],'response feature panel'),('duplicate',[original['responseFeatureIDs'][0]]*2,'response feature panel'),('unmeasured',['missing-source-feature'],'unmeasured training features')]:
  bad=copy.deepcopy(original);bad['responseFeatureIDs']=ids;write(a.out/(name+'.json'),bad)
  run(['singlecell-perturbation-fit',src/'training.h5ad','--plan',a.out/(name+'.json'),'--output',a.out/name],name,error)
  assert not (a.out/name).exists()
 oversized=a.out/'oversized.json';oversized.write_bytes(b' '*(2_097_152+1))
 run(['singlecell-perturbation-fit',src/'training.h5ad','--plan',oversized,'--output',a.out/'oversized'],'oversized','I/O bound exceeded');assert not (a.out/'oversized').exists()
 assert not list(a.out.glob('.numivivo-perturbation-*'))
 write(a.out/'checks.json',dict(status='passed',commands=len(commands),fullUniverseFeatures=len(new['featureIDs']),historicalNumericalModelExact=True,historicalAllFourPredictionsExact=True,failedOutputsAbsent=True,historicalModelSHA256=hashlib.sha256((a.historical_fold/'model.json').read_bytes()).hexdigest(),historicalPredictionSHA256=hashlib.sha256((a.historical_fold/'prediction.json').read_bytes()).hexdigest(),binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest()))
 print((a.out/'checks.json').read_text())
if __name__=='__main__':main()
