from pathlib import Path
import json,subprocess,time,hashlib,os
r=Path.home()/'numivivo-gse181897-20260911';src=r/'handoff';out=r/'native-handoff-qualified';out.mkdir(exist_ok=False)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
freeze=json.loads((src/'input-freeze.json').read_text())
for name,digest in freeze['files'].items():assert sha(src/name)==digest,name
binary=Path.home()/'numivivo-donor-intervals-20260911/runtime/numivivo-omics'
assert sha(binary)=='0f08cc423a10a2910962b25039c90731f2ceff013ec32e9a6df09d547a87ae4b'
library=Path.home()/'numivivo-multiassay-hdf5-20260909/libhdf5.dylib'
assert library.is_file()
env=dict(os.environ,NUMIVIVO_HDF5_LIBRARY=str(library))
commands=[]
for name,args in [
 ('aggregate',['singlecell-h5ad-pseudobulk',src/'B-source-codes-RNA.h5ad','--plan',src/'aggregate-plan.json','--output',out/'aggregate']),
 ('verify',['singlecell-h5ad-pseudobulk-verify',out/'aggregate']),
 ('repeat',['singlecell-h5ad-pseudobulk',src/'B-source-codes-RNA.h5ad','--plan',src/'aggregate-plan.json','--output',out/'repeat'])]:
 start=time.time()
 with (out/(name+'.log')).open('w') as f:result=subprocess.run([str(binary),*map(str,args)],stdout=f,stderr=subprocess.STDOUT,env=env)
 commands.append(dict(name=name,arguments=list(map(str,args)),exitCode=result.returncode,elapsedSeconds=time.time()-start));(out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n');assert result.returncode==0,(name,(out/(name+'.log')).read_text())
assert sha(out/'aggregate/report.json')==sha(out/'repeat/report.json')
(out/'execution.json').write_text(json.dumps(dict(status='passed-native-aggregation-replay-and-repeat',binarySHA256=sha(binary),hdf5LibrarySHA256=sha(library),hdf5LibraryPath=str(library),inputFreezeSHA256=sha(src/'input-freeze.json'),reportSHA256=sha(out/'aggregate/report.json'),commands=commands,conditionMeaningQualified=False,fitStarted=False),indent=2)+'\n')
print((out/'execution.json').read_text())
