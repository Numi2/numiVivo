#!/usr/bin/env python3
"""Independent complete-graph Louvain references and streamed native validation."""
import argparse,hashlib,importlib.metadata,itertools,json,random,time
from pathlib import Path
import igraph as ig
import ijson
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import adjusted_rand_score
RECORD=np.dtype([('row','<u4'),('column','<u4'),('value','<f8')])
INTEGER=np.dtype([('row','<u4'),('column','<u4'),('bits','<u8')])
GRAPH_FILES=['graph.json','plan.json','execution.json','neighbors.bin','bandwidths.bin','offsets.bin','edges.bin']
PCA_FILES=['scores.bin','loadings.bin','metadata.json','model.json','quality.json','plan.json']
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,v):
 with p.open('x') as f:json.dump(v,f,sort_keys=True,indent=2,allow_nan=False);f.write('\n')
def field(p,key):
 with p.open('rb') as f:return next(ijson.items(f,key,use_float=True))
def bound_source(root,protocol):
 report=read(root/'graph.json');receipt=read(root/'receipt.json')
 assert report['cells']==protocol['cells'] and report['method'].startswith('approximate-HNSW')
 for key,name in [('graph','graph.json'),('plan','plan.json'),('input','input/receipt.json'),('executionReport','execution.json')]:
  assert bytes(receipt[key]['bytes']).hex()==sha(root/name),name
 for key in ['neighbors','bandwidths','offsets','edges']:
  assert bytes(report[key]['bytes']).hex()==sha(root/(key+'.bin')),key
 assert (root/'edges.bin').stat().st_size==report['connectivityEntries']*16
 return report,np.memmap(root/'edges.bin',mode='r',dtype=RECORD)
def source(root,protocol,graph_check):
 assert graph_check['status']=='passed' and graph_check['receiptSHA256']==sha(root/'receipt.json')
 assert graph_check['protocolSHA256']==protocol['graphProtocolSHA256']
 return bound_source(root,protocol)
def equivalent_source(root,reference_root,protocol,graph_check):
 """Reuse original numerical evidence only after checking all fitted payloads.

 Receipts keep their original implementation bindings. This verifies a relation
 between two artifacts, rather than synthesizing replacement receipts/checks.
 """
 old_report,old_edges=source(reference_root,protocol,graph_check);del old_edges
 report,edges=bound_source(root,protocol)
 assert read(root/'plan.json')['inputKind']==read(reference_root/'plan.json')['inputKind']=='fitted'
 graph_names=GRAPH_FILES;pca_names=PCA_FILES
 bindings={}
 for name in graph_names+['input/'+name for name in pca_names]:
  current=root/name;original=reference_root/name
  assert current.is_file() and original.is_file() and not current.is_symlink() and not original.is_symlink(),name
  assert current.stat().st_size==original.stat().st_size,('Reference reuse size',name)
  digest=sha(current);assert digest==sha(original),('Reference reuse content',name)
  bindings[name]=dict(bytes=current.stat().st_size,SHA256=digest)
 current_pca=read(root/'input/receipt.json');original_pca=read(reference_root/'input/receipt.json')
 assert current_pca['source']==original_pca['source'],'Reference reuse source identity'
 assert current_pca['implementation']==read(root/'receipt.json')['implementation'],'Current PCA/graph implementation mismatch'
 assert original_pca['implementation']==read(reference_root/'receipt.json')['implementation'],'Original PCA/graph implementation mismatch'
 for name in pca_names:
  key={'scores.bin':'scores','loadings.bin':'loadings','metadata.json':'metadata','model.json':'model','quality.json':'quality','plan.json':'plan'}[name]
  assert bytes(current_pca[key]['bytes']).hex()==bindings['input/'+name]['SHA256'],name
  assert bytes(original_pca[key]['bytes']).hex()==bindings['input/'+name]['SHA256'],name
 proof=dict(mode='complete-fitted-graph-payload-equivalence',originalGraphReceiptSHA256=sha(reference_root/'receipt.json'),
  currentGraphReceiptSHA256=sha(root/'receipt.json'),originalPCAReceiptSHA256=sha(reference_root/'input/receipt.json'),
  currentPCAReceiptSHA256=sha(root/'input/receipt.json'),originalSourceSHA256=bytes(current_pca['source']['bytes']).hex(),
  equivalentPayloads=bindings,originalGraphCheckReceiptSHA256=graph_check['receiptSHA256'],
  scope='Every graph payload and all six fitted PCA payloads have identical bytes; original count-source fingerprints agree. Original graph check and igraph reference receipts remain unchanged. This establishes numerical-reference applicability, not validation of a different native executable or new biological qualification.')
 return report,edges,proof
def prepare(a,p):
 start=time.monotonic();report,edges=source(a.graph,p,read(a.graph_check));n=report['cells']
 assert ig.__version__==p['referenceVersion']
 # Take exactly one copy of each already-verified symmetric edge.
 upper=edges['row']<edges['column']
 pairs=np.column_stack([edges['row'][upper],edges['column'][upper]]).astype(np.int64)
 weights=edges['value'][upper].copy();del upper
 graph=ig.Graph(n=n,edges=pairs,directed=False);del pairs
 assert graph.vcount()==n and graph.ecount()==len(weights)
 references=[];arrays={}
 for seed in p['referenceSeeds']:
  ig.set_random_number_generator(random.Random(seed));began=time.monotonic()
  partition=graph.community_multilevel(weights=weights,resolution=p['resolution'])
  labels=np.asarray(partition.membership,dtype=np.int32);arrays['seed'+str(seed)]=labels
  objective=graph.modularity(labels,weights=weights,resolution=p['resolution'],directed=False)
  references.append(dict(seed=seed,clusters=len(partition),modularity=objective,seconds=time.monotonic()-began))
  write(a.out/('seed-'+str(seed)+'.json'),references[-1]);print('Reference seed',seed,references[-1],flush=True)
 np.savez_compressed(a.out/'partitions.npz',**arrays)
 write(a.out/'reference.json',dict(status='passed',cells=n,undirectedEdges=len(weights),graphReceiptSHA256=sha(a.graph/'receipt.json'),graphSHA256=sha(a.graph/'graph.json'),
  graphCheckSHA256=sha(a.graph_check),protocolSHA256=sha(a.protocol),partitionsSHA256=sha(a.out/'partitions.npz'),
  references=references,seconds=time.monotonic()-start,versions={v:importlib.metadata.version(v) for v in ['igraph','numpy','scipy','scikit-learn']},
  scope='Independent complete-graph weighted Louvain partitions; matching resolution, separate heuristic visitation and stopping rules. No native result was used to choose these partitions.'))
def check(a,p):
 start=time.monotonic();root=a.bundle;graphroot=root/'input'
 # Before/after hashes bind one immutable observation, including reference reuse.
 paths=[root/name for name in ['receipt.json','plan.json','result.json','execution.json']]
 for parent in {graphroot,getattr(a,'reference_graph',None)}-{None}:
  paths += [parent/name for name in GRAPH_FILES+['receipt.json','input/receipt.json']+['input/'+name for name in PCA_FILES]]
 paths += [a.protocol,a.graph_check,a.reference/'reference.json',a.reference/'partitions.npz']
 assert all(path.is_file() and not path.is_symlink() for path in paths)
 before={path:sha(path) for path in paths}
 receipt=read(root/'receipt.json')
 for key,name in [('input','input/receipt.json'),('plan','plan.json'),('result','result.json'),('executionReport','execution.json')]:
  assert bytes(receipt[key]['bytes']).hex()==sha(root/name),name
 gr=read(graphroot/'receipt.json');assert gr['implementation']==receipt['implementation']
 reference_graph=getattr(a,'reference_graph',None);reuse=None
 if reference_graph is None:report,edges=source(graphroot,p,read(a.graph_check))
 else:report,edges,reuse=equivalent_source(graphroot,reference_graph,p,read(a.graph_check))
 n=report['cells']
 reference=read(a.reference/'reference.json')
 assert reference['status']=='passed' and reference['protocolSHA256']==sha(a.protocol)
 assert reference['graphReceiptSHA256']==sha((reference_graph or graphroot)/'receipt.json')
 assert reference['graphSHA256']==sha(graphroot/'graph.json')
 assert reference['graphCheckSHA256']==sha(a.graph_check)
 assert reference['partitionsSHA256']==sha(a.reference/'partitions.npz')
 plan=read(root/'plan.json');options=field(root/'result.json','options')
 assert options==plan['clustering'] and options['seed']==p['nativeSeed'] and options['resolution']==p['resolution']
 execution=read(root/'execution.json');assert 0<execution['edgeVisits']<=plan['maximumEdgeVisits']
 native=root/'result.json';metadata=graphroot/'input/metadata.json'
 ir=read(graphroot/'input/receipt.json');assert sha(metadata)==bytes(ir['metadata']['bytes']).hex()
 sentinel=object();count=0
 with native.open('rb') as f,metadata.open('rb') as m:
  for left,right in itertools.zip_longest(ijson.items(f,'cells.item'),ijson.items(m,'cells.item'),fillvalue=sentinel):
   assert left is not sentinel and right is not sentinel
   assert left==dict(sampleID=right['sampleID'],barcode=right['barcode'])
   count+=1
 assert count==n
 with native.open('rb') as f:
  iterator=ijson.items(f,'labels.item');labels=np.fromiter(iterator,dtype=np.int64,count=n)
  assert next(iterator,None) is None
 sizes=np.asarray(field(native,'clusterSizes'),dtype=np.int64);k=len(sizes)
 assert 0<k<=n and (sizes>0).all() and labels.min()==0 and labels.max()==k-1
 np.testing.assert_array_equal(np.bincount(labels,minlength=k),sizes)
 first=np.full(k,n,dtype=np.int64);np.minimum.at(first,labels,np.arange(n))
 assert (np.diff(first)>0).all()
 degrees=np.bincount(edges['row'],weights=edges['value'],minlength=n)
 volume=np.bincount(labels,weights=degrees,minlength=k);total=float(degrees.sum())
 internal=labels[edges['row']]==labels[edges['column']]
 objective=float(edges['value'][internal].sum()/total-p['resolution']*np.sum((volume/total)**2)) if total else 0.0
 native_objective=float(field(native,'modularity'))
 np.testing.assert_allclose(native_objective,objective,rtol=1e-10,atol=1e-10)
 offsets_file=np.memmap(graphroot/'offsets.bin',mode='r',dtype=INTEGER)
 offsets=offsets_file['bits'].astype(np.int64)
 assert len(offsets)==n+1 and offsets[0]==0 and offsets[-1]==len(edges)
 within=sparse.csr_matrix((internal.astype(np.uint8),edges['column'],offsets),shape=(n,n));within.eliminate_zeros()
 _,components=connected_components(within,directed=False);del within
 _,representatives=np.unique(components,return_index=True)
 per_cluster=np.bincount(labels[representatives],minlength=k)
 disconnected=int(np.count_nonzero(per_cluster>1))
 assert disconnected==field(native,'disconnectedCommunities')
 levels=field(native,'levels')
 assert execution['aggregatedLevels']==len(levels)-1
 assert all(b['modularity']+1e-10>=a['modularity'] for a,b in zip(levels,levels[1:]))
 assert levels[-1]['communities']==k
 assert all(0<v['sweeps']<=options['maximumSweeps'] for v in levels)
 partitions=np.load(a.reference/'partitions.npz',allow_pickle=False);references=[]
 for value in reference['references']:
  membership=partitions['seed'+str(value['seed'])];assert membership.shape==(n,)
  references.append(dict(**value,adjustedRand=float(adjusted_rand_score(labels,membership))))
 gap=max(v['modularity'] for v in references)-objective;gate=gap<=p['maximumReferenceModularityDeficit']
 pairs=[dict(first=x,second=y,adjustedRand=float(adjusted_rand_score(partitions['seed'+str(x)],partitions['seed'+str(y)]))) for i,x in enumerate(p['referenceSeeds']) for y in p['referenceSeeds'][i+1:]]
 result=dict(status='passed' if gate else 'reference-quality-gap',cells=n,clusters=k,nativeModularity=native_objective,independentModularity=objective,
  disconnectedCommunities=disconnected,referenceMaximumModularityGap=gap,references=references,referencePairwiseARI=pairs,
  allOriginalCellIdentitiesChecked=True,allCellLabelsChecked=True,allGraphEdgesInObjective=True,allCommunityConnectivityChecked=True,
  edgeVisits=execution['edgeVisits'],resultSHA256=sha(native),graphReceiptSHA256=sha(graphroot/'receipt.json'),protocolSHA256=sha(a.protocol),
  seconds=time.monotonic()-start,scope=p['scope'],referenceReuse=reuse)
 assert all(sha(path)==digest for path,digest in before.items()),'Inputs changed during independent check'
 result['inputBindingsStableThroughoutCheck']=True
 write(a.out/'checks.json',result);print(json.dumps(result),flush=True)
 assert gate,'Frozen reference objective deficit exceeded; retain result before changing settings'
if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=['prepare','check'])
 for key in ['protocol','out','graph-check']:parser.add_argument('--'+key,type=Path,required=True)
 for key in ['graph','bundle','reference','reference-graph']:parser.add_argument('--'+key,type=Path)
 a=parser.parse_args();p=read(a.protocol);a.out.mkdir(parents=True,exist_ok=False)
 if a.mode=='prepare':assert a.graph is not None and a.reference_graph is None;prepare(a,p)
 else:assert a.bundle is not None and a.reference is not None;check(a,p)
