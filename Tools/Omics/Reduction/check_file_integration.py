#!/usr/bin/env python3
"""Fitted/query integration lifecycle, exact/HNSW downstream consumers and rehashed corruption controls."""
import argparse,copy,hashlib,json,shutil,struct,subprocess
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
p=argparse.ArgumentParser(description=__doc__)
for name in ['binary','oracle','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);commands=[]
def write(path,v):path.write_text(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False))
def run(label,args,ok=True,binary=None):
 r=subprocess.run([str(binary or a.binary),*map(str,args)],capture_output=True,text=True)
 (a.out/(label+'.log')).write_text(r.stdout+r.stderr);commands.append(dict(label=label,arguments=list(map(str,args)),exitCode=r.returncode,expectedSuccess=ok));write(a.out/'commands.json',commands)
 assert r.returncode==(0 if ok else 65),(label,r.stderr);return r
samples=[dict(id='s'+str(i),biologicalReplicateID='d'+str(i),donorID='d'+str(i),condition='shared',batchID='b'+str(i),organism='fixture') for i in range(2)]
mapping=dict(schemaVersion=1,id='file-integration-fixture',evidence='synthetic',countUnit='umiCount',matrixPath='X',sourceDescription='Numerical storage fixture only; no biological evidence',sampleColumn='sample',samples=samples,mitochondrialFeatureIDs=[])
rng=np.random.default_rng(20260909);counts=sparse.csr_matrix(rng.poisson(4,size=(64,12)).astype(np.int64));obs=pd.DataFrame({'sample':['s'+str(i%2) for i in range(64)]},index=['c'+str(i) for i in range(64)])
source=a.out/'source.h5ad';ad.AnnData(X=counts,obs=obs,var=pd.DataFrame(index=['g'+str(i) for i in range(12)])).write_h5ad(source)
fit=dict(schemaVersion=1,featureNamespace='fixture-genes',mapping=mapping,reduction=dict(pca=dict(components=3,highlyVariableFeatures=12,maximumBasis=12,meanBins=2)))
write(a.out/'fit.json',fit);run('fit',['singlecell-h5ad-pca',source,'--plan',a.out/'fit.json','--output',a.out/'fitted'])
query_source=a.out/'query-source.h5ad';query_data=ad.read_h5ad(source);query_data.obs_names=['query-'+v for v in query_data.obs_names];query_data.write_h5ad(query_source)
write(a.out/'query.json',dict(schemaVersion=1,featureNamespace='fixture-genes',mapping=mapping));run('query',['singlecell-h5ad-pca-query',query_source,'--plan',a.out/'query.json','--reference',a.out/'fitted','--output',a.out/'query'])
options=dict(clusters=3,maximumIterations=3,seed=7);write(a.out/'options.json',options)
for kind in ['fitted','query']:
 plan=dict(schemaVersion=1,inputKind=kind,integration=options);path=a.out/(kind+'-plan.json');write(path,plan);root=a.out/(kind+'-integrated')
 run(kind+'-publish',['singlecell-pca-integrate',a.out/kind,'--plan',path,'--output',root]);run(kind+'-verify',['singlecell-pca-integrate-verify',root]);run(kind+'-repeat',['singlecell-pca-integrate',a.out/kind,'--plan',path,'--output',a.out/(kind+'-repeat')])
 assert (root/'receipt.json').read_bytes()==(a.out/(kind+'-repeat/receipt.json')).read_bytes()
 oracle=a.out/(kind+'-legacy.json');run(kind+'-legacy',[a.out/kind,3,a.out/'options.json',oracle],binary=a.oracle);legacy=json.loads(oracle.read_text());report=json.loads((root/'report.json').read_text())
 for file,key,width in [('scores.bin','scores',3),('memberships.bin','memberships',3),('assignment-scores.bin','assignmentScores',3)]:
  expected=b''.join(struct.pack('<IId',i,j,v) for i,row in enumerate(legacy[key]) for j,v in enumerate(row));assert (root/file).read_bytes()==expected
 for key in ['objectives','relativeImprovements','assignmentCenters','levels','cellLevels','stoppingReason','maximumRidgeResidual']:assert report[key]==legacy[key]
 batch_plan=a.out/(kind+'-batch-plan.json');write(batch_plan,dict(schemaVersion=1,inputKind=kind,integration=dict(options,covariate='batch')));batch_root=a.out/(kind+'-batch')
 run(kind+'-batch-publish',['singlecell-pca-integrate',a.out/kind,'--plan',batch_plan,'--output',batch_root]);run(kind+'-batch-verify',['singlecell-pca-integrate-verify',batch_root])
 assert json.loads((batch_root/'report.json').read_text())['levels']==['b0','b1']
 for filename in ['scores.bin','memberships.bin','assignment-scores.bin']:assert (batch_root/filename).read_bytes()==(root/filename).read_bytes()
 adaptive_plan=a.out/(kind+'-adaptive-plan.json');write(adaptive_plan,dict(schemaVersion=1,inputKind=kind,integration=dict(options,ridge=.2,ridgeScaling='expectedClusterBatchMass')));adaptive=a.out/(kind+'-adaptive')
 run(kind+'-adaptive-publish',['singlecell-pca-integrate',a.out/kind,'--plan',adaptive_plan,'--output',adaptive]);run(kind+'-adaptive-verify',['singlecell-pca-integrate-verify',adaptive])
 report=json.loads((adaptive/'report.json').read_text());assert report['method']=='diversity-soft-clustering-expected-mass-ridge-Double-v1' and len(report['ridgePenalties'])==3
 tampered=a.out/(kind+'-adaptive-penalties-tamper');shutil.copytree(adaptive,tampered);report['ridgePenalties'][0][0]+=1;write(tampered/'report.json',report)
 receipt=json.loads((tampered/'receipt.json').read_text());receipt['report']=dict(bytes=list(hashlib.sha256((tampered/'report.json').read_bytes()).digest()));write(tampered/'receipt.json',receipt)
 run(kind+'-adaptive-rehashed-penalties',['singlecell-pca-integrate-verify',tampered],False)
 for mode in ['exact','hnsw']:
  gp=dict(schemaVersion=1,inputKind='integrated',neighbors=dict(neighbors=5,representation='integrated'),execution=dict(workers=1),storage='binary')
  if mode=='hnsw':gp['approximation']=dict(connections=8,constructionWidth=64,searchWidth=64)
  path=a.out/(kind+'-'+mode+'-graph.json');write(path,gp);graph=a.out/(kind+'-'+mode+'-graph');run(kind+'-'+mode+'-graph',['singlecell-pca-neighbors',root,'--plan',path,'--output',graph]);run(kind+'-'+mode+'-graph-verify',['singlecell-pca-neighbors-verify',graph])
  cp=a.out/'cluster-plan.json';write(cp,dict(schemaVersion=1,clustering=dict()));ep=a.out/'embed-plan.json';write(ep,dict(schemaVersion=1,embedding=dict(epochs=20)))
  for action,plan in [('cluster',cp),('embed',ep)]:
   target=a.out/(kind+'-'+mode+'-'+action);run(kind+'-'+mode+'-'+action,['singlecell-graph-'+action,graph,'--plan',plan,'--output',target]);run(kind+'-'+mode+'-'+action+'-verify',['singlecell-graph-'+action+'-verify',target])
base=a.out/'fitted-integrated'
for file,key in [('scores.bin','scores'),('memberships.bin','memberships'),('assignment-scores.bin','assignmentScores'),('report.json','report'),('metadata.json','metadata')]:
 root=a.out/('tamper-'+key);shutil.copytree(base,root);path=root/file
 if file.endswith('.bin'):
  data=bytearray(path.read_bytes());data[8]^=1;path.write_bytes(data)
 else:
  value=json.loads(path.read_text())
  if key=='report':value['objectives'][0]+=1
  else:value['cells'][0]['sampleID']='s1'
  write(path,value)
 receipt=json.loads((root/'receipt.json').read_text());receipt[key]=dict(bytes=list(hashlib.sha256(path.read_bytes()).digest()));write(root/'receipt.json',receipt)
 run('rehashed-'+key,['singlecell-pca-integrate-verify',root],False)
root=a.out/'tamper-parent';shutil.copytree(base,root);path=root/'input/scores.bin';data=bytearray(path.read_bytes());data[8]^=1;path.write_bytes(data)
receipt=json.loads((root/'input/receipt.json').read_text());receipt['scores']=dict(bytes=list(hashlib.sha256(data).digest()));write(root/'input/receipt.json',receipt)
receipt=json.loads((root/'receipt.json').read_text());receipt['input']=dict(bytes=list(hashlib.sha256((root/'input/receipt.json').read_bytes()).digest()));write(root/'receipt.json',receipt)
run('rehashed-parent',['singlecell-pca-integrate-verify',root],False)
for label,plan in [('budget',dict(schemaVersion=1,integration=dict(options,maximumWork=1))),('unknown',dict(schemaVersion=1,doner=1)),('recursive',dict(schemaVersion=1,inputKind='integrated',integration=options))]:
 path=a.out/(label+'.json');write(path,plan);target=a.out/('rejected-'+label);run(label,['singlecell-pca-integrate',a.out/'fitted','--plan',path,'--output',target],False);assert not target.exists()
run('overwrite',['singlecell-pca-integrate',a.out/'fitted','--plan',a.out/'fitted-plan.json','--output',base],False)
for label,edit in [('confounded',lambda m:[s.update(condition=s['donorID']) for s in m['samples']]),('unknown-donor',lambda m:[s.pop('donorID') for s in m['samples']])]:
 spec=copy.deepcopy(fit);edit(spec['mapping']);path=a.out/(label+'-fit.json');write(path,spec);target=a.out/(label+'-pca');run(label+'-fit',['singlecell-h5ad-pca',source,'--plan',path,'--output',target]);r=run(label,['singlecell-pca-integrate',target,'--plan',a.out/'fitted-plan.json','--output',a.out/('rejected-'+label)],False)
 assert ('confounded with condition' if label=='confounded' else 'known selected covariate') in r.stderr
for kind,representation in [('fitted','integrated'),('integrated','pca')]:
 plan=dict(schemaVersion=1,inputKind=kind,neighbors=dict(neighbors=5,representation=representation));path=a.out/(kind+'-mismatch.json');write(path,plan);run(kind+'-mismatch',['singlecell-pca-neighbors',base,'--plan',path,'--output',a.out/(kind+'-rejected-graph')],False)
assert not any(p.name.startswith(('.integration-matrix-','.numivivo-integration-')) for p in a.out.rglob('*'))
summary=dict(status='passed',commands=len(commands),expectedRejections=sum(not c['expectedSuccess'] for c in commands),allFrozenMatrixBitsExact=True,repeatedReceiptsExact=True,fittedAndQuery=True,donorAndBatch=True,adaptiveFittedAndQuery=True,rehashedAdaptivePenaltiesRejected=True,exactAndHNSWDownstreamClusteringEmbedding=True,rehashedOutputsAndParentRejected=True,confoundedAndUnknownCovariateRejected=True,recursiveAndRepresentationMismatchRejected=True,scratchRemoved=True,binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest())
write(a.out/'checks.json',summary);print(json.dumps(summary),flush=True)
