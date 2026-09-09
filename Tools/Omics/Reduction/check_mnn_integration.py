#!/usr/bin/env python3
"""Native MNN fitted/query lifecycle, independent reference, and corruption controls."""
import argparse,copy,hashlib,json,shutil,subprocess
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import scanorama
from check_full_integration import matrix
p=argparse.ArgumentParser(description=__doc__)
for name in ['binary','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);commands=[]
def write(path,value):path.write_text(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False))
def run(label,args,ok=True):
 r=subprocess.run([str(a.binary),*map(str,args)],capture_output=True,text=True);(a.out/(label+'.log')).write_text(r.stdout+r.stderr)
 commands.append(dict(label=label,arguments=list(map(str,args)),exitCode=r.returncode,expectedSuccess=ok));write(a.out/'commands.json',commands)
 assert r.returncode==(0 if ok else 65),(label,r.stderr);return r
samples=[dict(id='s'+str(i),biologicalReplicateID='d'+str(i),donorID='d'+str(i),condition='shared',batchID='b'+str(i),organism='fixture') for i in range(3)]
mapping=dict(schemaVersion=1,id='mnn-numerical-fixture',evidence='synthetic',countUnit='umiCount',matrixPath='X',sourceDescription='Numerical lifecycle fixture; not biological validation',sampleColumn='sample',samples=samples,mitochondrialFeatureIDs=[])
rng=np.random.default_rng(20260910);counts=sparse.csr_matrix(rng.poisson(4,size=(72,12)).astype(np.int64));obs=pd.DataFrame({'sample':['s'+str(i%3) for i in range(72)]},index=['c'+str(i) for i in range(72)])
source=a.out/'source.h5ad';ad.AnnData(X=counts,obs=obs,var=pd.DataFrame(index=['g'+str(i) for i in range(12)])).write_h5ad(source)
fit=dict(schemaVersion=1,featureNamespace='fixture-genes',mapping=mapping,reduction=dict(pca=dict(components=3,highlyVariableFeatures=12,maximumBasis=12,meanBins=2)))
write(a.out/'fit.json',fit);run('fit',['singlecell-h5ad-pca',source,'--plan',a.out/'fit.json','--output',a.out/'fitted'])
query=ad.read_h5ad(source);query.obs_names=['query-'+x for x in query.obs_names];query.write_h5ad(a.out/'query-source.h5ad');write(a.out/'query-plan.json',dict(schemaVersion=1,featureNamespace='fixture-genes',mapping=mapping))
run('query',['singlecell-h5ad-pca-query',a.out/'query-source.h5ad','--plan',a.out/'query-plan.json','--reference',a.out/'fitted','--output',a.out/'query'])
reconstructions=[]
for kind in ['fitted','query']:
 plan=dict(schemaVersion=1,inputKind=kind,mnn=dict(neighbors=5,minimumAlignment=0));path=a.out/(kind+'-mnn-plan.json');write(path,plan);root=a.out/(kind+'-mnn')
 run(kind+'-publish',['singlecell-pca-integrate',a.out/kind,'--plan',path,'--output',root]);run(kind+'-verify',['singlecell-pca-integrate-verify',root]);run(kind+'-repeat',['singlecell-pca-integrate',a.out/kind,'--plan',path,'--output',a.out/(kind+'-repeat')])
 assert (root/'receipt.json').read_bytes()==(a.out/(kind+'-repeat/receipt.json')).read_bytes()
 receipt=json.loads((root/'receipt.json').read_text());assert 'anchors' in receipt and 'memberships' not in receipt and 'assignmentScores' not in receipt
 assert not (root/'memberships.bin').exists() and not (root/'assignment-scores.bin').exists()
 metadata=json.loads((root/'metadata.json').read_text());assert metadata==json.loads((a.out/kind/'metadata.json').read_text());sm={s['id']:s for s in metadata['samples']};batch=np.array([sm[c['sampleID']]['donorID'] for c in metadata['cells']]);x=matrix(a.out/kind/'scores.bin',72,3);indices=[np.flatnonzero(batch==b) for b in np.unique(batch)];scale=np.median(np.linalg.norm(x,axis=1));ds=[x[ix]/scale for ix in indices]
 alignments,matches=scanorama.find_alignments(ds,knn=5,approx=False,alpha=0,verbose=0);reference=scanorama.assemble(ds,knn=5,sigma=15,approx=False,alpha=0,batch_size=256,verbose=0,alignments=alignments,matches=matches);expected=np.empty_like(x)
 for ix,y in zip(indices,reference):expected[ix]=y*scale
 actual=matrix(root/'scores.bin',72,3);np.testing.assert_allclose(actual,expected,atol=1e-9,rtol=1e-10)
 report=json.loads((root/'report.json').read_text());assert abs(report['medianRowNorm']-scale)<1e-12
 records=np.fromfile(root/'anchors.bin',dtype=[('first','<u4'),('second','<u4'),('distance','<f8')]);assert np.isfinite(records['distance']).all()
 for pair in report['alignments']:
  i,j=pair['firstLevel'],pair['secondLevel'];got=records[pair['anchorOffset']:pair['anchorOffset']+pair['anchors']]
  assert list(zip(map(int,got['first']),map(int,got['second'])))==sorted((int(indices[i][u]),int(indices[j][v])) for u,v in matches[i,j])
 np.savez_compressed(a.out/(kind+'-reference.npz'),scores=expected);reconstructions.append(dict(kind=kind,maximumCoordinateError=float(np.max(np.abs(actual-expected))),anchorsExact=True))
 batch_plan=copy.deepcopy(plan);batch_plan['mnn']['covariate']='batch';write(a.out/(kind+'-batch-plan.json'),batch_plan)
 run(kind+'-batch',['singlecell-pca-integrate',a.out/kind,'--plan',a.out/(kind+'-batch-plan.json'),'--output',a.out/(kind+'-batch')]);assert (root/'scores.bin').read_bytes()==(a.out/(kind+'-batch/scores.bin')).read_bytes()
 for mode in ['exact','hnsw']:
  gp=dict(schemaVersion=1,inputKind='integrated',neighbors=dict(neighbors=5,representation='integrated'),execution=dict(workers=1),storage='binary')
  if mode=='hnsw':gp['approximation']=dict(connections=8,constructionWidth=64,searchWidth=64)
  gplan=a.out/(kind+'-'+mode+'-plan.json');write(gplan,gp);graph=a.out/(kind+'-'+mode+'-graph');run(kind+'-'+mode+'-graph',['singlecell-pca-neighbors',root,'--plan',gplan,'--output',graph]);run(kind+'-'+mode+'-graph-verify',['singlecell-pca-neighbors-verify',graph])
  for action,options in [('cluster',dict(clustering={})),('embed',dict(embedding=dict(epochs=20)))]:
   ep=a.out/(action+'-plan.json');write(ep,dict(schemaVersion=1,**options));target=a.out/(kind+'-'+mode+'-'+action);run(kind+'-'+mode+'-'+action,['singlecell-graph-'+action,graph,'--plan',ep,'--output',target]);run(kind+'-'+mode+'-'+action+'-verify',['singlecell-graph-'+action+'-verify',target])
base=a.out/'fitted-mnn'
for filename,key in [('scores.bin','scores'),('anchors.bin','anchors'),('report.json','report')]:
 root=a.out/('tamper-'+key);shutil.copytree(base,root);path=root/filename
 if filename.endswith('.bin'):
  raw=bytearray(path.read_bytes());raw[0]^=1;path.write_bytes(raw)
 else:
  value=json.loads(path.read_text());value['medianRowNorm']*=2;write(path,value)
 receipt=json.loads((root/'receipt.json').read_text());receipt[key]=dict(bytes=list(hashlib.sha256(path.read_bytes()).digest()));write(root/'receipt.json',receipt);run('rehashed-'+key,['singlecell-pca-integrate-verify',root],False)
for label,plan in [('work',dict(schemaVersion=1,mnn=dict(maximumWork=1))),('memory',dict(schemaVersion=1,mnn=dict(maximumResidentBytes=1))),('kernel-budget',dict(schemaVersion=1,mnn=dict(neighbors=5,minimumAlignment=0,maximumWork=5184))),('assembly-memory',dict(schemaVersion=1,mnn=dict(neighbors=5,minimumAlignment=0,maximumResidentBytes=36864))),('neighbors',dict(schemaVersion=1,mnn=dict(neighbors=100))),('two-methods',dict(schemaVersion=1,integration={},mnn={})),('recursive',dict(schemaVersion=1,inputKind='integrated',mnn={}))]:
 path=a.out/(label+'-reject-plan.json');write(path,plan);target=a.out/('rejected-'+label);run(label,['singlecell-pca-integrate',a.out/'fitted','--plan',path,'--output',target],False);assert not target.exists()
run('overwrite',['singlecell-pca-integrate',a.out/'fitted','--plan',a.out/'fitted-mnn-plan.json','--output',base],False)
assert not any(p.name.startswith('.numivivo-integration-') for p in a.out.rglob('*'))
write(a.out/'checks.json',dict(status='passed',commands=len(commands),expectedRejections=sum(not c['expectedSuccess'] for c in commands),fittedAndQuery=True,donorAndBatch=True,independentReconstruction=reconstructions,exactAndHNSWGraphClusteringEmbedding=True,rehashedWitnessesRejected=True,methodSpecificWitnesses=True,scratchRemoved=True,binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),qualification='Numerical and artifact lifecycle fixtures only; full-cohort biological evaluation is separate.'))
