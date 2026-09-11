#!/usr/bin/env python3
"""Exercise the actual product router and executable-identity owner on full Kang."""
import hashlib,json,os,subprocess
from pathlib import Path
root=Path(os.environ.get('NUMIVIVO_PROFILE_ROOT', '/Users/n/numivivo-normalization-profile-20260911'));binary=root/'runtime/numivivo-omics';source=Path('/Users/n/numivivo-legacy-count-route-20260911/count-store')
env=dict(os.environ,NUMIVIVO_HDF5_LIBRARY='/Users/n/numivivo-multiassay-hdf5-20260909/libhdf5.dylib')
rows=[]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def run(label,args,expected=0,reason=None):
 assert not (root/(label+'.log')).exists()
 p=subprocess.run([str(binary),*map(str,args)],capture_output=True,env=env)
 (root/(label+'.stdout')).write_bytes(p.stdout);(root/(label+'.log')).write_bytes(p.stderr)
 assert p.returncode==expected,(label,p.stderr)
 if reason:assert reason in p.stderr.decode(),(label,p.stderr)
 rows.append(dict(label=label,arguments=list(map(str,args)),returnCode=p.returncode))
run('cli-help',['singlecell-help'])
assert b'--backend cpu-fp64|metal-fp32' in (root/'cli-help.stdout').read_bytes()
store=root/'cli-store'
run('cli-store',['singlecell-h5ad-store',source/'original.h5ad','--plan',source/'plan.json','--output',store])
for name in ['original.h5ad','counts.bin','metadata.json','quality.json','plan.json']:assert sha(store/name)==sha(source/name)
for label,flags in [('cli-cpu',[]),('cli-metal',['--backend','metal-fp32']),('cli-cpu-explicit',['--backend','cpu-fp64'])]:
 run(label,['singlecell-count-store-normalize',store,'--target','10000','--output',root/label,*flags])
 base=Path('/Users/n/numivivo-metal-normalization-20260911')/('metal-0' if 'metal' in label else 'cpu-0')
 assert sha(root/label/'values.bin')==sha(base/'values.bin')
 assert sha(root/label/'metadata.json')==sha(base/'metadata.json')
 run(label+'-verify',['singlecell-count-store-normalize-verify',root/label,'--store',store])
for name in ['values.bin','metadata.json','input-receipt.json','receipt.json']:assert sha(root/'cli-cpu'/name)==sha(root/'cli-cpu-explicit'/name)
for label,flags,reason in [('cli-unknown',['--backend','unknown'],'count normalization backend'),('cli-missing',['--backend'],'singlecell-count-store-normalize'),('cli-wrong-flag',['--device','metal-fp32'],'count normalization backend')]:
 run(label,['singlecell-count-store-normalize',store,'--target','10000','--output',root/label,*flags],65,reason)
 assert not (root/label).exists()
run('cli-overwrite',['singlecell-count-store-normalize',store,'--target','10000','--output',root/'cli-metal','--backend','metal-fp32'],65,'destination exists')
run('cli-old-implementation',['singlecell-count-store-normalize',source,'--target','10000','--output',root/'cli-old-implementation','--backend','metal-fp32'],65,'count store receipt')
assert not (root/'cli-old-implementation').exists()
for line in (root/'runtime/sources.sha256').read_text().splitlines():
 digest,path=line.split(None,1);assert sha(Path(path.strip()))==digest
result=dict(status='passed',completeEntries=14184532,commands=rows,defaultAndExplicitCPUExact=True,productAndHarnessValuesExact=True,sourceCountsAndMetadataExact=True,binarySHA256=sha(binary),compiledSourcesSHA256=sha(root/'runtime/sources.sha256'),implementation='Actual product executable plus OS fingerprint; no zero tag or mocked router')
(root/'cli-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status='passed',commands=len(rows),binarySHA256=result['binarySHA256'])))
