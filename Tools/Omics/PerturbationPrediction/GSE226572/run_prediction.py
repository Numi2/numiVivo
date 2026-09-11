#!/usr/bin/env python3
"""Execute and freeze all native predictions before opening treated outcomes for scoring."""
import argparse,json,os,re,subprocess,time
from pathlib import Path
from download import sha,write

def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);a=p.parse_args();src=a.study/'prediction-inputs';out=a.study/'native-prediction';out.mkdir(exist_ok=False)
 f=json.loads((src/'input-freeze.json').read_text());assert f['status']=='frozen-before-fit'
 for name,h in f['files'].items():assert sha(src/name)==h,name
 execution=dict(status='running',binarySHA256=sha(a.binary),hdf5SHA256=sha(Path(os.environ['NUMIVIVO_HDF5_LIBRARY'])),inputFreezeSHA256=sha(src/'input-freeze.json'),runnerSHA256=sha(Path(__file__)),commands=[]);write(out/'execution.json',execution)
 commands=[('fit',['singlecell-perturbation-fit',src/'training.h5ad','--plan',src/'training.json','--output',out/'model']),('verify-model',['singlecell-perturbation-verify',out/'model']),('predict',['singlecell-perturbation-predict',src/'query.h5ad','--plan',src/'query.json','--reference',out/'model','--output',out/'prediction']),('verify-prediction',['singlecell-perturbation-prediction-verify',out/'prediction']),('repeat',['singlecell-perturbation-predict',src/'query.h5ad','--plan',src/'query.json','--reference',out/'model','--output',out/'repeat'])]
 for name,args in commands:
  start=time.time();r=subprocess.run(['/usr/bin/time','-l',str(a.binary),*map(str,args)],capture_output=True,text=True);(out/(name+'.log')).write_text(r.stdout+r.stderr);peak=re.search(r'(\d+)\s+maximum resident set size',r.stderr);execution['commands'].append(dict(name=name,args=list(map(str,args)),exitCode=r.returncode,seconds=time.time()-start,maximumResidentBytes=int(peak[1]) if peak else None));write(out/'execution.json',execution);assert r.returncode==0,(name,r.stderr);print(name+' passed',flush=True)
 assert (out/'prediction/report.json').read_bytes()==(out/'repeat/report.json').read_bytes();execution.update(status='completed',repeatedPredictionBytesExact=True);write(out/'execution.json',execution)
 write(out/'prediction-freeze.json',dict(status='completed-before-scoring',executionSHA256=sha(out/'execution.json'),files={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()},scoringStarted=False));print('Complete predictions frozen',flush=True)
if __name__=='__main__':main()
