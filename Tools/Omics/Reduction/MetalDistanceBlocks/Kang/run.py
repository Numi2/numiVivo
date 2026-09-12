from pathlib import Path
import json,hashlib,subprocess,os,time
r=Path(__file__).parent;source=Path('/Users/n/numivivo-native-mnn-20260910/kang');owner=Path('/Users/n/numivivo-native-metal-knn-20260912');binary=owner/'build/numivivo-omics'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
prior=json.loads(Path('/Users/n/numivivo-local-mnn-20260911/prior-artifacts.json').read_text());item=next(x for x in prior if x['cohort']=='kang')
for name in ['pca/original.h5ad','pca/scores.bin','pca/metadata.json']:
 assert sha(source/name)==next(x['sha256'] for x in item['files'] if x['path']==name)
assert sha(binary)==json.loads((owner/'cli-terminal.json').read_text())['binarySHA256']
protocol=dict(ownerCommit='e6248e770730419d9dc94847e38407dce04aff95',binarySHA256=sha(binary),sourceSHA256=sha(source/'pca/original.h5ad'),originalPCASHA256=sha(source/'pca/scores.bin'),originalMetadataSHA256=sha(source/'pca/metadata.json'),pcaPlanSHA256=sha(source/'fit.json'),cells=24673,seeds=[7,19,42],neighborsIncludingSelf=20,graphPrecision=['CPU FP64','Metal FP32'],input='Original full fitted PCA, not integrated MNN',scope='Kang backend numerical and source-annotation preservation; development cohort, no new biological prediction accuracy',gates=dict(minimumPartitionARI=0.99,minimumMatchedAssignmentAgreement=0.99,maximumCellTypeARIRegression=0.01,maximumPerTypeNeighborRecallRegression=0.01,maximumConditionBalancedAccuracyRegression=0.01),annotations='Original cell group, sample donorID and condition; never inputs to neighbor distance or clustering',createdUnix=time.time())
(r/'protocol.json').write_text(json.dumps(protocol,indent=2));results=[]
env=dict(os.environ,NUMIVIVO_HDF5_LIBRARY='/Users/n/numivivo-integration-reference-py/lib/python3.13/site-packages/h5py/.dylibs/libhdf5.320.0.0.dylib')
def run(label,args):
 start=time.time()
 with (r/(label+'.log')).open('w') as f:code=subprocess.call(['/usr/bin/time','-l',str(binary)]+list(map(str,args)),stdout=f,stderr=subprocess.STDOUT,env=env)
 results.append(dict(step=label,exitCode=code,seconds=time.time()-start));(r/'status.json').write_text(json.dumps(results,indent=2));assert code==0,label
run('pca',['singlecell-h5ad-pca',source/'pca/original.h5ad','--plan',source/'fit.json','--output',r/'pca'])
assert sha(r/'pca/scores.bin')==protocol['originalPCASHA256'];assert sha(r/'pca/metadata.json')==protocol['originalMetadataSHA256']
for mode in ['cpu','metal']:
 execution=dict(workers=1,queryBlockRows=128,candidateBlockRows=8192)
 if mode=='metal':execution['backend']='metalFP32'
 plan=dict(schemaVersion=1,inputKind='fitted',storage='binary',neighbors=dict(neighbors=20,maximumDistancePairs=500000000,representation='pca'),execution=execution)
 path=r/(mode+'-plan.json');path.write_text(json.dumps(plan))
 run(mode,['singlecell-pca-neighbors',r/'pca','--plan',path,'--output',r/mode])
for seed in protocol['seeds']:
 path=r/('cluster-'+str(seed)+'.json');path.write_text(json.dumps(dict(schemaVersion=1,clustering=dict(resolution=1,seed=seed,maximumSweeps=100,maximumLevels=32,levelTolerance=1e-7),maximumEdgeVisits=1000000000)))
 for mode in ['cpu','metal']:
  name=mode+'-'+str(seed)
  run(name,['singlecell-graph-cluster',r/mode,'--plan',path,'--output',r/name])
  run(name+'-verify',['singlecell-graph-cluster-verify',r/name])
(r/'terminal.json').write_text(json.dumps(dict(results=results,completedUnix=time.time()),indent=2))
