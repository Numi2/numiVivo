#!/usr/bin/env python3
"""Sample the unchanged owner; every profiling output retains exact prior values."""
import hashlib,json,os,shutil,subprocess,time
from pathlib import Path
root=Path('/Users/n/numivivo-normalization-profile-20260911');previous=Path('/Users/n/numivivo-metal-normalization-20260911');store='/Users/n/numivivo-legacy-count-route-20260911/count-store'
assert shutil.disk_usage(root).free>4_500_000_000
assert hashlib.sha256((previous/'runtime/h5ad-check').read_bytes()).hexdigest()=='11155a9405926245f55cdd34d242834d7bb0608b69acac17523426376972326e'
env=dict(os.environ,NUMIVIVO_HDF5_LIBRARY='/Users/n/numivivo-multiassay-hdf5-20260909/libhdf5.dylib')
results=[]
for backend,label in [('cpu-fp64','cpu'),('metal-fp32','metal')]:
 with (root/(label+'.stdout')).open('wb') as out,(root/(label+'.log')).open('wb') as err:
  started=time.monotonic();p=subprocess.Popen([str(root/'profile-counts'),store,backend,str(previous/(label+'-0')),str(root/label)],stdout=out,stderr=err,env=env)
  (root/(label+'-process.json')).write_text(json.dumps(dict(pid=p.pid,backend=backend,timeUnix=time.time()))+'\n')
  sample=subprocess.run(['/usr/bin/sample',str(p.pid),'5','1','-file',str(root/(label+'-sample.txt'))],capture_output=True)
  (root/(label+'-sample.stdout')).write_bytes(sample.stdout);(root/(label+'-sample.log')).write_bytes(sample.stderr)
  code=p.wait();row=dict(backend=backend,returnCode=code,sampleReturnCode=sample.returncode,seconds=time.monotonic()-started,pid=p.pid)
  (root/(label+'-status.json')).write_text(json.dumps(row,indent=2)+'\n');results.append(row)
  assert code==sample.returncode==0,(row,(root/(label+'.log')).read_text())
  expected=json.loads((previous/(label+'-0')/'receipt.json').read_text());checks=json.loads((root/label/'checks.json').read_text());assert len(checks)==8 and all(x==expected for x in checks)
(root/'checks.json').write_text(json.dumps(dict(status='passed',calls=16,entriesPerCall=14184532,allReceiptsExact=True,rows=results,scope='Stack sampling with profiler overhead; not a replacement speed benchmark',binarySHA256=hashlib.sha256((root/'profile-counts').read_bytes()).hexdigest()),indent=2)+'\n')
print('passed: 16 exact full-cohort calls sampled')
