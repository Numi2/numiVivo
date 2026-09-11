#!/usr/bin/env python3
"""Compare exact full-cohort bundles before timing a shared-writer promotion."""
import hashlib,json,os,re,shutil,subprocess,time
from pathlib import Path
r=Path(os.environ.get('NUMIVIVO_PROFILE_ROOT', '/Users/n/numivivo-normalization-profile-20260911'));old=Path('/Users/n/numivivo-metal-normalization-20260911');store=Path('/Users/n/numivivo-legacy-count-route-20260911/count-store')
binaries={'old':old/'runtime/h5ad-check','new':r/'runtime-matched/h5ad-check'}
env=dict(os.environ,NUMIVIVO_HDF5_LIBRARY='/Users/n/numivivo-multiassay-hdf5-20260909/libhdf5.dylib')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
assert shutil.disk_usage(r).free>1_000_000_000
freeze=dict(base='128c86801f630912e24c213418cb9c8ad7ae0bde',binaries={k:sha(v) for k,v in binaries.items()},timeUnix=time.time())
(r/'comparison-freeze.json').write_text(json.dumps(freeze,indent=2)+'\n')
orders=[['old-cpu','new-cpu','old-metal','new-metal'],['new-metal','old-metal','new-cpu','old-cpu'],['old-metal','new-metal','old-cpu','new-cpu']]
rows=[];restoration=[];verified=set()
for round_index,order in enumerate(orders):
 for arm in order:
  owner,backend=arm.split('-');label=f'compare-{round_index}-{arm}';dest=r/label
  args=[str(binaries[owner]),'normalize-count-store',str(store),'10000','cpu-fp64' if backend=='cpu' else 'metal-fp32',str(dest)]
  with (r/(label+'.stdout')).open('wb') as out,(r/(label+'.log')).open('wb') as err:
   start=time.monotonic();p=subprocess.run(['/usr/bin/time','-l',*args],stdout=out,stderr=err,env=env);seconds=time.monotonic()-start
  text=(r/(label+'.log')).read_text();row=dict(label=label,owner=owner,backend=backend,round=round_index,seconds=seconds,returnCode=p.returncode,peakRSSBytes=int(re.search(r'(\d+)\s+maximum resident set size',text).group(1)))
  rows.append(row);(r/'comparison-progress.json').write_text(json.dumps(rows,indent=2)+'\n');assert p.returncode==0,text
  reference=old/(backend+'-0');identities={}
  for name in ['values.bin','metadata.json','input-receipt.json','receipt.json']:
   digest=sha(dest/name);assert digest==sha(reference/name),(label,name)
   identities[name]=digest
  row['completeBundleBytesExact']=True;row['files']=identities
  if owner=='new' and backend not in verified:
   p=subprocess.run([str(binaries[owner]),'verify-normalized-count-store',str(dest),str(store)],capture_output=True,env=env)
   (r/(label+'-verify.log')).write_bytes(p.stderr);(r/(label+'-verify.stdout')).write_bytes(p.stdout);assert p.returncode==0,p.stderr;verified.add(backend)
  handles=subprocess.run(['lsof','-nP','--',str(dest/'values.bin')],capture_output=True);assert handles.returncode==1 and not handles.stdout and not handles.stderr
  restoration.append(dict(path=str(dest/'values.bin'),retainedPath=str(reference/'values.bin'),SHA256=identities['values.bin'],bytes=(dest/'values.bin').stat().st_size));(dest/'values.bin').unlink()
assert verified=={'cpu','metal'}
(r/'comparison.json').write_text(json.dumps(dict(status='passed',rows=rows,restoration=restoration,freeze=freeze,completeEntries=14184532,nativeReplaysPassed=True,allBundleBytesExact=True),indent=2)+'\n');print('passed: 12 complete bundle comparisons and both new-owner replays')
