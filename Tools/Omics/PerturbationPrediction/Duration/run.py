#!/usr/bin/env python3
"""Fit, reconstruct and freeze the actual native duration-model predictions."""
import argparse,json,os,re,subprocess,time
from pathlib import Path
from common import sha,write,verify_files

def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);a=p.parse_args();src=a.study/'inputs';out=a.study/'native';out.mkdir(exist_ok=False)
 frozen=json.loads((src/'input-freeze.json').read_text());assert frozen['status']=='frozen-before-duration-fit';verify_files(src,frozen['files']);assert frozen['protocolSHA256']==sha(Path(__file__).with_name('PROTOCOL.md'))
 execution=dict(status='running',binarySHA256=sha(a.binary),hdf5SHA256=sha(Path(os.environ['NUMIVIVO_HDF5_LIBRARY'])),inputFreezeSHA256=sha(src/'input-freeze.json'),runnerSHA256=sha(Path(__file__)),commands=[]);write(out/'execution.json',execution)
 for fold in json.loads((src/'folds.json').read_text()):
  fid=fold['id'];dest=out/fid;dest.mkdir();f=src/fid
  commands=[('fit',['singlecell-duration-fit',f/'training.h5ad','--plan',f/'training.json','--output',dest/'model']),('verify-model',['singlecell-duration-verify',dest/'model']),('predict',['singlecell-duration-predict',f/'query.h5ad','--plan',f/'query.json','--reference',dest/'model','--output',dest/'prediction']),('verify-prediction',['singlecell-duration-prediction-verify',dest/'prediction']),('repeat',['singlecell-duration-predict',f/'query.h5ad','--plan',f/'query.json','--reference',dest/'model','--output',dest/'repeat'])]
  for name,args in commands:
   start=time.time();r=subprocess.run(['/usr/bin/time','-l',str(a.binary),*map(str,args)],capture_output=True,text=True);(dest/(name+'.log')).write_text(r.stdout+r.stderr);peak=re.search(r'(\d+)\s+maximum resident set size',r.stderr)
   execution['commands'].append(dict(fold=fid,name=name,args=list(map(str,args)),exitCode=r.returncode,seconds=time.time()-start,maximumResidentBytes=int(peak[1]) if peak else None));write(out/'execution.json',execution);assert r.returncode==0,(fid,name,r.stderr);print(fid+' '+name+' passed',flush=True)
  assert (dest/'prediction/report.json').read_bytes()==(dest/'repeat/report.json').read_bytes()
 execution.update(status='completed',allThreeRepeatedPredictionBytesExact=True);write(out/'execution.json',execution)
 write(out/'prediction-freeze.json',dict(status='completed-before-duration-scoring',developmentNotExternalValidation=True,executionSHA256=sha(out/'execution.json'),files={str(f.relative_to(out)):sha(f) for f in sorted(out.rglob('*')) if f.is_file()},scoringStarted=False));print('All 18 native outcomes frozen',flush=True)
if __name__=='__main__':main()
