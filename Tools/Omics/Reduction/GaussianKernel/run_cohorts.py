#!/usr/bin/env python3
"""Run original complete cohorts under scan and tiled matching, retaining all attempts."""
import argparse,hashlib,json,os,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for n in ('binary','source','prior-manifest','out'):p.add_argument('--'+n,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def write(p,x):p.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
os.environ['VECLIB_MAXIMUM_THREADS']='1'
old=json.loads(a.prior_manifest.read_text());inputs={}
for c in ('hagai','kang','ding'):
 entries=next(v['files'] for v in old if v['cohort']==c)
 for name in ('pca/scores.bin','pca/metadata.json','mnn/plan.json','mnn/report.json','mnn/anchors.bin','mnn/scores.bin'):
  entry=next(v for v in entries if v['path']==name);q=a.source/c/name;assert sha(q)==entry['sha256'];inputs[c+'/'+name]=entry['sha256']
write(a.out/'execution-freeze.json',dict(pid=os.getpid(),createdUnix=time.time(),binarySHA256=sha(a.binary),inputSHA256=inputs,priorManifestSHA256=sha(a.prior_manifest),driverSHA256=sha(Path(__file__)),protocolSHA256=sha(Path(__file__).with_name('PROTOCOL.md')),nativeSourceBindings=(a.binary.parent/'runtime/sources.sha256').read_text(),scope='Complete original cohorts, all coordinates and anchors; scan/tiled change only, unchanged biological evaluation required after numerical checks. No million-cell claim.'))
commands=[]
for c in ('hagai','kang','ding'):
 for mode in ('scan','tiled','tiled-replay'):
  label=c+'-'+mode;command=['/usr/bin/time','-l',str(a.binary),str(a.source/c),str(a.out/label),'scan' if mode=='scan' else 'tiled'];start=time.time()
  with (a.out/(label+'.log')).open('x') as log:
   process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT);write(a.out/(label+'-process.json'),dict(pid=process.pid,driverPID=os.getpid(),command=command,startedUnix=start));code=process.wait()
  commands.append(dict(label=label,returnCode=code,seconds=time.time()-start));write(a.out/'commands.json',commands);print(commands[-1],flush=True);assert code==0,label
 for name in ('scores.bin','anchors.bin','report.json'):
  assert sha(a.out/(c+'-tiled')/name)==sha(a.out/(c+'-tiled-replay')/name), 'tiled replay differs: '+name
for n,h in inputs.items():assert sha(a.source/n)==h
write(a.out/'complete.json',dict(status='passed',commands=commands,allOriginalAnchorBytesAndAssemblyStructuresExact=True,allOriginalScoresAndStepMetricsWithinTolerance=True,metricsRead=False,allTiledReplaysExact=True,binarySHA256=sha(a.binary),results={str(q.relative_to(a.out)):sha(q) for q in a.out.rglob('*') if q.is_file()},qualification='Complete original numerical comparison. Biological evaluation pending; no million-cell or Metal claim.'))
