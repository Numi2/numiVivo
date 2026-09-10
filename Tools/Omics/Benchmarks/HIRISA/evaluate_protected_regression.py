#!/usr/bin/env python3
"""Evaluate the frozen regression alternative with unchanged full-cohort diagnostics."""
import argparse,copy,os,shutil,subprocess,time
from pathlib import Path
from prepare_annotation_retention import read,write
from check_integration_response import sha,validate_design
from check_program_calibration import evaluate
from check_integration_programs import compare


def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ('study','root','python'):p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args();s=a.study;r=a.root;owner=Path(__file__).parent;outputs=read(r/'output-freeze.json')
 assert outputs['status']=='passed' and outputs['metricsRead'] is False
 for name,h in outputs['files'].items():assert sha(r/name)==h
 e=r/'evaluation';e.mkdir();old=s/'annotation-retention';f=copy.deepcopy(read(old/'freeze.json'))
 for name in ('protocol.md','ledger.json','rows.npz'):shutil.copy2(old/name,e/name)
 f.update(parentAnnotationFreezeSHA256=sha(old/'freeze.json'),candidateOutputFreezeSHA256=sha(r/'output-freeze.json'),developmentProtocolSHA256=sha(r/'protocol.md'),declaredAtUnix=time.time())
 for name,file in [('protected','protected.bin'),('legacy-reconstructed','legacy-reconstructed.bin')]:f['matrices'][name]=dict(path=str((r/file).relative_to(s)),SHA256=sha(r/file))
 write(e/'freeze.json',f)
 write(r/'evaluation-execution-freeze.json',dict(createdUnix=time.time(),outputFreezeSHA256=sha(r/'output-freeze.json'),annotationFreezeSHA256=sha(e/'freeze.json'),sourceFiles={p.name:sha(p) for p in owner.glob('*.py')},pid=os.getpid()))
 env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'};commands=[]
 def run(name,args):
  command=[str(a.python),*map(str,args)];start=time.time()
  with (r/(name+'.log')).open('x') as stream:
   process=subprocess.Popen(command,env=env,stdout=stream,stderr=subprocess.STDOUT)
   write(r/(name+'-process.json'),dict(pid=process.pid,driverPID=os.getpid(),command=command,startedUnix=start));code=process.wait()
  status=dict(name=name,returnCode=code,seconds=time.time()-start);write(r/(name+'-status.json'),status);commands.append(status);print(status,flush=True);assert code==0,name
 run('tests',[owner/'test_protected_regression.py'])
 for name in ('baseline','legacy-reconstructed','protected'):
  args=[owner/'check_annotation_retention.py','--study',s,'--root',e,'--matrix',name,'--out',e/(name+'.json')]
  if name!='baseline':args+=['--baseline',e/'baseline.json']
  run(name+'-annotation',args)
  run(name+'-annotation-oracle',[owner/'check_annotation_retention_oracle.py','--study',s,'--root',e,'--result',e/(name+'.json'),'--out',e/(name+'-oracle.json')])
 for first,second in [('baseline','baseline'),('legacy-reconstructed','native')]:
  new,prior=read(e/(first+'.json')),read(old/(second+'.json'))
  for x,y in zip(new['folds'],prior['folds']):
   for k in ('id','confusion','metrics','accuracy','balancedRecall'):assert x[k]==y[k],(first,k)
 common=['--source',s/'native-release/original.h5ad','--metadata',s/'integration-full/native/metadata.json','--design',s/'design.json','--scores',r/'protected.bin']
 run('response',[owner/'check_integration_response.py',*common,'--protocol',s/'integration-response-v2/protocol.json','--baseline',s/'integration-response-v2/baseline.json','--out',e/'response.json'])
 run('mixed-programs',[owner/'check_integration_programs.py',*common,'--protocol',s/'integration-programs/protocol.json','--program-reference',s/'integration-programs/reference','--baseline',s/'integration-programs/baseline.json','--out',e/'mixed-programs.json'])
 mixed=read(e/'mixed-programs.json');protocol=read(s/'program-calibration/protocol.json');_,ids=validate_design(read(s/'design.json'))
 calibrated=evaluate(mixed,ids,protocol)
 baseline=read(s/'program-calibration/results/baseline.json')
 calibrated.update(status='measured',protocol=protocol,scope='Same frozen within-library decoder on one protected fixed-membership regression development candidate.',bindings={k:v for k,v in mixed['bindings'].items() if k!='baseline'},previousResultSHA256=sha(e/'mixed-programs.json'),calibrationProtocolSHA256=sha(s/'program-calibration/protocol.json'),evaluatorSHA256=sha(owner/'check_program_calibration.py'),driverSHA256=sha(Path(__file__)))
 calibrated['comparisons']=compare(baseline,calibrated,protocol)
 calibrated['allProgramGradientGatesPassed']=all(all(x['gates'].values()) for x in calibrated['comparisons'])
 write(e/'within-programs.json',calibrated)
 run('within-programs-oracle',[owner/'check_calibration_dense_oracle.py',*common,'--protocol',s/'integration-programs/protocol.json','--program-reference',s/'integration-programs/reference','--calibration-protocol',s/'program-calibration/protocol.json','--result',e/'within-programs.json','--out',e/'within-programs-oracle.json'])
 for name,h in outputs['files'].items():assert sha(r/name)==h
 write(r/'complete.json',dict(status='completed-development-diagnostic',commands=commands,allOriginalCellsPreserved=True,legacyConfusionAndMetricsExact=True,allNumericalChecksPassed=True,outputFreezeSHA256=sha(r/'output-freeze.json'),results={p.name:sha(p) for p in e.glob('*.json')}))
if __name__=='__main__':main()
