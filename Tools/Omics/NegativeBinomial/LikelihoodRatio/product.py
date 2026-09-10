#!/usr/bin/env python3
"""Real H5AD publication/replay and exact unchanged-Wald compatibility."""
import argparse,gzip,hashlib,json,os,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for name in ['root','baseline','binary','hdf5']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.root.mkdir(parents=True,exist_ok=False)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2)+'\n')
archive=json.loads((a.baseline/'archive.json').read_text());source=Path(archive['source']);assert sha(source)==archive['sourceSHA256']
for e in archive['entries']:
 assert sha(a.baseline/e['path'])==e['sha256'];assert hashlib.sha256(gzip.decompress((a.baseline/e['path']).read_bytes())).hexdigest()==e['logicalSHA256']
baseline=gzip.decompress((a.baseline/'report.json.gz').read_bytes());plan=gzip.decompress((a.baseline/'plan.json.gz').read_bytes())
env={**os.environ,'NUMIVIVO_HDF5_LIBRARY':str(a.hdf5)};runs=[];binary_sha=sha(a.binary);hdf5_sha=sha(a.hdf5)
for method in ['wald','likelihoodRatio']:
 authored=plan
 if method=='likelihoodRatio':
  obj=json.loads(plan)
  for c in obj['contrasts']:c['negativeBinomialOptions']['testMethod']='likelihoodRatio'
  authored=json.dumps(obj,sort_keys=True,indent=2).encode()+b'\n'
 input_path=a.root/(method+'-plan.json');input_path.write_bytes(authored);output=a.root/method
 for stage,command in [('publish',['singlecell-h5ad-pseudobulk',source,'--plan',input_path,'--output',output]),('replay',['singlecell-h5ad-pseudobulk-verify',output])]:
  assert sha(a.binary)==binary_sha;start=time.monotonic()
  with (a.root/(method+'-'+stage+'.log')).open('w') as log:run=subprocess.run(['/usr/bin/time','-l',str(a.binary)]+list(map(str,command)),stdout=log,stderr=subprocess.STDOUT,env=env)
  runs.append(dict(method=method,stage=stage,command=list(map(str,command)),exitCode=run.returncode,seconds=time.monotonic()-start));write(a.root/'runs.json',dict(binarySHA256=binary_sha,hdf5SHA256=hdf5_sha,runs=runs));assert run.returncode==0
 assert sha(output/'original.h5ad')==archive['sourceSHA256']
 if method=='wald':assert (output/'report.json').read_bytes()==baseline
 entries=[]
 for name in ['plan.json','report.json','receipt.json']:
  raw=(output/name).read_bytes();packed=gzip.compress(raw,mtime=0);target=output/(name+'.gz');target.write_bytes(packed);assert gzip.decompress(packed)==raw
  entries.append(dict(path=target.name,sha256=sha(target),logicalPath=name,logicalSHA256=hashlib.sha256(raw).hexdigest()))
 # Remove only redundant publication snapshots after replay and exact checks.
 handles=subprocess.run(['lsof','+D',str(output)],capture_output=True);assert handles.returncode==1 and not handles.stdout
 write(output/'archive.json',dict(status='published-and-replay-verified-before-archival',source=str(source),sourceSHA256=sha(source),entries=entries,completeSourceRetained=True))
 for name in ['original.h5ad','plan.json','report.json','receipt.json']:(output/name).unlink()
write(a.root/'complete.json',dict(status='real-H5AD-Wald-and-LRT-published-and-replayed',originalWaldReportExactlyEqual=True,binarySHA256=binary_sha,sourceSHA256=archive['sourceSHA256'],baselineReportLogicalSHA256=hashlib.sha256(baseline).hexdigest()))
