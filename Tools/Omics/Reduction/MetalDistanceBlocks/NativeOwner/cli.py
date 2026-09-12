from pathlib import Path
import subprocess,json,os,time,hashlib
r=Path(__file__).parent;binary=r/'build/numivivo-omics';old=Path('/Users/n/numivivo-mnn-cli-lifecycle-20260912')
env=dict(os.environ,NUMIVIVO_HDF5_LIBRARY='/Users/n/numivivo-integration-reference-py/lib/python3.13/site-packages/h5py/.dylibs/libhdf5.320.0.0.dylib')
plan=dict(schemaVersion=1,inputKind='fitted',storage='binary',neighbors=dict(neighbors=20,maximumDistancePairs=100000000,representation='pca'),execution=dict(workers=1,queryBlockRows=128,candidateBlockRows=8192,backend='metalFP32'))
(r/'metal-plan.json').write_text(json.dumps(plan));plan['execution'].pop('backend');(r/'cpu-plan.json').write_text(json.dumps(plan));plan['execution']['backend']='metalFP32';plan['storage']='json';(r/'metal-json-plan.json').write_text(json.dumps(plan))
results=[]
def run(name,args,expected=0):
 t=time.time()
 with (r/(name+'.log')).open('w') as f:code=subprocess.call(['/usr/bin/time','-l',str(binary)]+list(map(str,args)),stdout=f,stderr=subprocess.STDOUT,env=env)
 results.append(dict(step=name,exitCode=code,seconds=time.time()-t));(r/'cli-status.json').write_text(json.dumps(results,indent=2));assert code==expected,name
run('fresh-pca',['singlecell-h5ad-pca','/Users/n/numivivo-native-mnn-20260910/hagai/pca/original.h5ad','--plan',old/'pca-plan.json','--output',r/'pca'])
for mode in ['cpu','metal','metal-json']:
 run(mode,['singlecell-pca-neighbors',r/'pca','--plan',r/(mode+'-plan.json'),'--output',r/mode])
 run(mode+'-verify',['singlecell-pca-neighbors-verify',r/mode])
for name,change in [('workers',dict(workers=2)),('backend',dict(backend='unknown'))]:
 p=json.loads((r/'metal-plan.json').read_text());p['execution'].update(change);f=r/('invalid-'+name+'.json');f.write_text(json.dumps(p));dest=r/('rejected-'+name)
 run('reject-'+name,['singlecell-pca-neighbors',r/'pca','--plan',f,'--output',dest],65);assert not dest.exists()
p=json.loads((r/'metal-plan.json').read_text());p['approximation']={};f=r/'invalid-hnsw.json';f.write_text(json.dumps(p))
run('reject-hnsw',['singlecell-pca-neighbors',r/'pca','--plan',f,'--output',r/'rejected-hnsw'],65);assert not (r/'rejected-hnsw').exists()
(r/'cli-terminal.json').write_text(json.dumps(dict(results=results,binarySHA256=hashlib.sha256(binary.read_bytes()).hexdigest(),completedUnix=time.time()),indent=2))
