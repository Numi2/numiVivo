from pathlib import Path
import json,hashlib,subprocess,os,time
r=Path(__file__).parent;parent=Path('/Users/n/numivivo-native-metal-knn-20260912');binary=parent/'build/numivivo-omics'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((parent/'manifest.json').read_text());expected={x['path']:x['sha256'] for x in manifest['files']}
for mode in ['cpu','metal']:
 for name in ['neighbors.bin','edges.bin','offsets.bin','bandwidths.bin','graph.json','receipt.json']:
  assert sha(parent/mode/name)==expected[mode+'/'+name]
assert sha(binary)==expected['build/numivivo-omics']
protocol=dict(sourceCommit='e6248e770730419d9dc94847e38407dce04aff95',binarySHA256=sha(binary),inputManifestSHA256=sha(parent/'manifest.json'),seeds=[7,19,42],resolution=1,maximumSweeps=100,maximumLevels=32,levelTolerance=1e-7,cells=13863,scope='Complete real-Hagai CPU FP64 versus Metal FP32 graph clustering stability; no author cell labels and no biological accuracy claim',acceptance=dict(minimumARI=0.99,minimumMatchedAssignmentAgreement=0.99,requireAllSeeds=True),createdUnix=time.time())
(r/'protocol.json').write_text(json.dumps(protocol,indent=2));results=[]
env=dict(os.environ,NUMIVIVO_HDF5_LIBRARY='/Users/n/numivivo-integration-reference-py/lib/python3.13/site-packages/h5py/.dylibs/libhdf5.320.0.0.dylib')
for seed in protocol['seeds']:
 plan=dict(schemaVersion=1,clustering=dict(resolution=1,seed=seed,maximumSweeps=100,maximumLevels=32,levelTolerance=1e-7),maximumEdgeVisits=1000000000)
 path=r/('plan-'+str(seed)+'.json');path.write_text(json.dumps(plan))
 for mode in ['cpu','metal']:
  name=mode+'-'+str(seed);dest=r/name
  for step,args in [('fit',['singlecell-graph-cluster',parent/mode,'--plan',path,'--output',dest]),('verify',['singlecell-graph-cluster-verify',dest])]:
   start=time.time()
   with (r/(name+'-'+step+'.log')).open('w') as f:code=subprocess.call(['/usr/bin/time','-l',str(binary)]+list(map(str,args)),stdout=f,stderr=subprocess.STDOUT,env=env)
   results.append(dict(mode=mode,seed=seed,step=step,exitCode=code,seconds=time.time()-start));(r/'status.json').write_text(json.dumps(results,indent=2));assert code==0
(r/'terminal.json').write_text(json.dumps(dict(results=results,completedUnix=time.time()),indent=2))
