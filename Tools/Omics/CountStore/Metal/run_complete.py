#!/usr/bin/env python3
"""Run full original Kang native CPU/Metal publications and exact replay."""
import hashlib,json,os,platform,re,shutil,subprocess,time
from pathlib import Path
root=Path('/Users/n/numivivo-metal-normalization-20260911')
store=Path('/Users/n/numivivo-legacy-count-route-20260911/count-store')
repo=Path('/Users/n/numivivo-h5ad-20260909')
binary=root/'runtime/h5ad-check'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
assert shutil.disk_usage(root).free>3_500_000_000
assert sha(store/'original.h5ad')=='89c2c4d6aea0163676b45574aaafd3593525dce269ee315911680d8f31c014a0'
assert sha(store/'counts.bin')=='c4ef502f8c1874e851381de121712f9f1c45cc560a50d1a0576e415c796203c7'
receipt=json.loads((store/'receipt.json').read_text());assert receipt['entries']==14184532
for line in (root/'runtime/sources.sha256').read_text().splitlines():
 digest,path=line.split(None,1);assert sha(Path(path.strip()))==digest
protocol=repo/'Tools/Omics/CountStore/Metal/PROTOCOL.md'
freeze=dict(base='9cb829b8a800e3e89e0437d912d4e500234a2148',protocolSHA256=sha(protocol),binarySHA256=sha(binary),baselineBinarySHA256=sha(root/'cpu-baseline'),timeUnix=time.time(),platform=platform.platform(),countReceiptSHA256=sha(store/'receipt.json'))
(root/'run-freeze.json').write_text(json.dumps(freeze,indent=2)+'\n')
# Record competing processes, but do not terminate any external workload.
ps=subprocess.check_output(['ps','-axo','pid,ppid,pcpu,command'],text=True)
(root/'workload-before.txt').write_text(ps)
assert not any('xctest -XCTest' in line or 'NumiBrainPackageTests.xctest' in line for line in ps.splitlines()), 'External test workload still active'
env=dict(os.environ,NUMIVIVO_HDF5_LIBRARY='/Users/n/numivivo-multiassay-hdf5-20260909/libhdf5.dylib')
rows=[]
def run(label,args):
 assert not (root/(label+'.status.json')).exists()
 t=time.monotonic()
 with (root/(label+'.stdout')).open('wb') as out,(root/(label+'.log')).open('wb') as err:
  p=subprocess.run(['/usr/bin/time','-l',*map(str,args)],stdout=out,stderr=err,env=env)
 log=(root/(label+'.log')).read_text();rss=re.search(r'(\d+)\s+maximum resident set size',log)
 row=dict(label=label,arguments=list(map(str,args)),returnCode=p.returncode,seconds=time.monotonic()-t,peakRSSBytes=int(rss.group(1)) if rss else None)
 (root/(label+'.status.json')).write_text(json.dumps(row,indent=2)+'\n');rows.append(row)
 assert p.returncode==0,(label,log[-3000:])
 for name in ['receipt.json','metadata.json','input-receipt.json','values.bin']:
  if (root/label/name).exists():row.setdefault('files',{})[name]=dict(bytes=(root/label/name).stat().st_size,SHA256=sha(root/label/name))
 (root/'run-progress.json').write_text(json.dumps(rows,indent=2)+'\n')
run('cpu-before',[root/'cpu-baseline',store,root/'cpu-before'])
for i in range(3):
 for backend,label in [('cpu-fp64','cpu'),('metal-fp32','metal')]:
  name=f'{label}-{i}';run(name,[binary,'normalize-count-store',store,'10000',backend,root/name])
for label in ['cpu-0','metal-0']:
 run(label+'-verify',[binary,'verify-normalized-count-store',root/label,store])
for label in ['cpu-before','cpu-1','cpu-2']:
 assert next(r for r in rows if r['label']==label)['files']==next(r for r in rows if r['label']=='cpu-0')['files']
for label in ['metal-1','metal-2']:
 assert next(r for r in rows if r['label']==label)['files']==next(r for r in rows if r['label']=='metal-0')['files']
assert 'execution' not in json.loads((root/'cpu-0/receipt.json').read_text())
assert json.loads((root/'metal-0/receipt.json').read_text())['execution']['backend']=='metal-fp32'
(root/'native-checks.json').write_text(json.dumps(dict(status='passed',entries=14184532,rows=rows,cpuDefaultOldBytesExact=True,metalRepeatBytesExact=True,cpuAndMetalNativeReplay=True,freeze=freeze),indent=2)+'\n')
print(json.dumps(dict(status='passed',runs=len(rows),binarySHA256=freeze['binarySHA256'])))
