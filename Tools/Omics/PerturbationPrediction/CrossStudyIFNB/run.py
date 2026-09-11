#!/usr/bin/env python3
"""Execute all frozen native fits and predictions, then freeze outputs before scoring."""
import argparse,hashlib,json,os,platform,subprocess,time
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 freeze=json.loads((a.inputs/'input-freeze.json').read_text())
 for name,digest in freeze['files'].items():assert sha(a.inputs/name)==digest,name
 write(a.out/'execution.json',dict(binary=str(a.binary.resolve()),binarySHA256=sha(a.binary),inputFreezeSHA256=sha(a.inputs/'input-freeze.json'),runnerSHA256=sha(__file__),platform=platform.platform(),hdf5Library=os.environ.get('NUMIVIVO_HDF5_LIBRARY'),startedUnix=time.time()))
 commands=[]
 def run(args,log):
  start=time.time();r=subprocess.run([str(a.binary.resolve()),*map(str,args)],capture_output=True,text=True);log.write_text(r.stdout+r.stderr)
  commands.append(dict(arguments=list(map(str,args)),exitCode=r.returncode,elapsedSeconds=time.time()-start));write(a.out/'commands.json',commands)
  assert r.returncode==0,(args,r.stdout,r.stderr)
 for index,f in enumerate(json.loads((a.inputs/'folds.json').read_text())):
  src=a.inputs/f['id'];dest=a.out/f['id'];dest.mkdir()
  run(['singlecell-perturbation-fit',src/'training.h5ad','--plan',src/'training.json','--output',dest/'model'],dest/'fit.log')
  run(['singlecell-perturbation-verify',dest/'model'],dest/'verify-model.log')
  args=['singlecell-perturbation-predict',src/'query.h5ad','--plan',src/'query.json','--reference',dest/'model','--output']
  run(args+[dest/'prediction'],dest/'predict.log')
  run(['singlecell-perturbation-prediction-verify',dest/'prediction'],dest/'verify-prediction.log')
  if index==0:
   run(args+[dest/'repeat'],dest/'repeat.log');assert (dest/'repeat/report.json').read_bytes()==(dest/'prediction/report.json').read_bytes()
  print(json.dumps(dict(completed=index+1,fold=f['id'])),flush=True)
 files={str(p.relative_to(a.out)):sha(p) for p in sorted(a.out.rglob('*')) if p.is_file()}
 write(a.out/'prediction-freeze.json',dict(schemaVersion=1,files=files,commands=len(commands),completedFolds=index+1,scoringStarted=False,completedUnix=time.time()))
 print(json.dumps(dict(completedFolds=index+1,commands=len(commands),freezeSHA256=sha(a.out/'prediction-freeze.json'))),flush=True)
if __name__=='__main__':main()
