#!/usr/bin/env python3
"""Full UMAP schedule/optimizer comparison and bounded query-to-full-cohort quality."""
import argparse,hashlib,importlib.metadata,json
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from umap.umap_ import find_ab_params,make_epochs_per_sample
from umap.layouts import optimize_layout_euclidean
p=argparse.ArgumentParser(description=__doc__)
for name in ['bundle','oracle-result','protocol','out']:p.add_argument('--'+name,type=Path,required=True)
p.add_argument('--case',choices=['baron','hagai','norman'],required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
def sha(path):
 with path.open('rb') as f:return hashlib.file_digest(f,'sha256').digest()
def document(path):return json.loads(path.read_text())
r=a.bundle;protocol=document(a.protocol);receipt=document(r/'receipt.json');e=document(r/'result.json');plan=document(r/'plan.json');execution=document(r/'execution.json');root=r/'input';g=document(root/'graph.json');gr=document(root/'receipt.json');ir=document(root/'input/receipt.json')
for key,name in [('input','input/receipt.json'),('plan','plan.json'),('result','result.json'),('executionReport','execution.json')]:assert sha(r/name)==bytes(receipt[key]['bytes'])
for key,name in [('input','input/receipt.json'),('plan','plan.json'),('graph','graph.json'),('executionReport','execution.json')]:assert sha(root/name)==bytes(gr[key]['bytes'])
for key in ['neighbors','bandwidths','offsets','edges']:assert sha(root/(key+'.bin'))==bytes(g[key]['bytes'])
for key in ['metadata','scores']:assert sha(root/('input/'+key+('.json' if key=='metadata' else '.bin')))==bytes(ir[key]['bytes'])
assert (r/'result.json').read_bytes()==a.oracle_result.read_bytes(),'Frozen pre-refactor trajectory differs'
o=e['options'];assert o==plan['embedding']==protocol['embedding'];metadata=document(root/'input/metadata.json')
assert e['cells']==[dict(sampleID=v['sampleID'],barcode=v['barcode']) for v in metadata['cells']]
n=len(e['cells']);assert n==g['cells']==protocol['cohorts'][a.case]
record=np.dtype([('row','<u4'),('column','<u4'),('bits','<u8')]);scores=np.fromfile(root/'input/scores.bin',dtype=record)
assert len(scores)==n*g['dimensions'];np.testing.assert_array_equal(scores['row'],np.repeat(np.arange(n),g['dimensions']));np.testing.assert_array_equal(scores['column'],np.tile(np.arange(g['dimensions']),n))
x=scores['bits'].view('<f8').reshape(n,g['dimensions']);y=np.asarray(e['coordinates']);initial=np.asarray(e['initialCoordinates']);assert y.shape==initial.shape==(n,o['dimensions']) and np.isfinite(y).all()
curve_a,curve_b=find_ab_params(o['spread'],o['minimumDistance']);np.testing.assert_allclose([e['curveA'],e['curveB']],[curve_a,curve_b],rtol=2e-5,atol=2e-5)
xx=np.linspace(0,3*o['spread'],300);target=np.where(xx<o['minimumDistance'],1,np.exp(-(xx-o['minimumDistance'])/o['spread']));curve=1/(1+e['curveA']*xx**(2*e['curveB']))
np.testing.assert_allclose(e['curveSquaredError'],np.sum((curve-target)**2),rtol=1e-10,atol=1e-12)
edges=np.fromfile(root/'edges.bin',dtype=record);weights=edges['bits'].view('<f8');keep=weights>=weights.max()/o['epochs'];maximum=weights.max();head=edges['row'][keep].astype(np.int64);tail=edges['column'][keep].astype(np.int64);weights=weights[keep]
assert len(weights)==e['retainedDirectedEdges'] and len(edges)-len(weights)==e['discardedDirectedEdges']
assert execution['scheduleBytes']==len(weights)*32 and execution['maximumScheduleMappedBytes']==16*1024*1024
assert not any(v.name.startswith('.embedding-schedule-') for v in r.iterdir())
assert e['edgeVisits']==len(weights)*o['epochs'] and e['edgeVisits']<=o['maximumUpdates']
# Native scheduling has a declared max/weight period; compare it with UMAP's
# algebraically equivalent formula but keep their floating-point increments distinct.
intervals=maximum/weights;reference_intervals=make_epochs_per_sample(weights,o['epochs']);np.testing.assert_allclose(intervals,reference_intervals,rtol=2e-14,atol=2e-14)
next_positive=intervals.copy();ni=intervals/o['negativeSampleRate'];next_negative=ni.copy();positive_count=negative_count=0
for epoch in range(o['epochs']):
 active=np.flatnonzero(next_positive<=epoch);samples=np.maximum(0,((epoch-next_negative[active])/ni[active]).astype(np.int64))
 positive_count+=len(active);negative_count+=int(samples.sum());next_positive[active]+=intervals[active];next_negative[active]+=samples*ni[active]
assert positive_count==e['attractiveUpdates'] and negative_count==e['negativeSamples'] and positive_count+negative_count<=o['maximumUpdates']
print('Full schedule and frozen trajectory checks passed',flush=True)
coordinates=dict(initial=initial,native=y)
for seed in protocol['referenceSeeds']:
 values=np.ascontiguousarray(initial,dtype=np.float32);rng=np.random.RandomState(seed).randint(-(2**31)+1,2**31-1,3).astype(np.int64)
 values=optimize_layout_euclidean(values,values,head,tail,o['epochs'],n,reference_intervals,e['curveA'],e['curveB'],rng,gamma=o['repulsionStrength'],initial_alpha=o['learningRate'],negative_sample_rate=o['negativeSampleRate'],parallel=False,move_other=True)
 assert np.isfinite(values).all();coordinates['seed'+str(seed)]=values
 np.savez_compressed(a.out/('reference-'+str(seed)+'.npz'),coordinates=values)
 print('Full optimizer seed '+str(seed)+' passed',flush=True)
query=protocol['qualityQueries'][a.case]
rows=np.arange(n) if query=='all' else np.sort(np.random.default_rng(query['seed']).choice(n,size=query['count'],replace=False))
np.save(a.out/'quality-rows.npy',rows);k=protocol['evaluationNeighbors'];indices=np.arange(n);sums={name:[0,0,0] for name in coordinates}
# For each query, rank against all cohort cells. Only 64-row distance blocks
# and cell-scale rank vectors are live; never build a full n-by-n matrix.
for start in range(0,len(rows),64):
 selected=rows[start:start+64];high=cdist(x[selected],x);low={name:cdist(values[selected],values) for name,values in coordinates.items()}
 for slot,i in enumerate(selected):
  hi=np.lexsort((indices,high[slot]));hi=hi[hi!=i];hr=np.empty(n,dtype=np.int64);hr[hi]=np.arange(1,n)
  for name,distances in low.items():
   lo=np.lexsort((indices,distances[slot]));lo=lo[lo!=i];lr=np.empty(n,dtype=np.int64);lr[lo]=np.arange(1,n)
   sums[name][0]+=int(np.maximum(hr[lo[:k]]-k,0).sum());sums[name][1]+=int(np.maximum(lr[hi[:k]]-k,0).sum());sums[name][2]+=len(np.intersect1d(hi[:k],lo[:k]))
 if start%1024==0:print('Quality queries '+str(min(start+64,len(rows)))+'/'+str(len(rows)),flush=True)
denominator=len(rows)*k*(2*n-3*k-1)
metrics={name:dict(trustworthiness=1-2*v[0]/denominator,continuity=1-2*v[1]/denominator,neighborRecall=v[2]/(len(rows)*k)) for name,v in sums.items()}
references=[dict(seed=seed,**metrics['seed'+str(seed)]) for seed in protocol['referenceSeeds']];native=metrics['native']
passed=(native['trustworthiness']>=protocol['minimumTrustworthiness'] and native['trustworthiness']>=min(v['trustworthiness'] for v in references)-protocol['maximumTrustworthinessDeficitFromWorstReference'] and native['neighborRecall']>=min(v['neighborRecall'] for v in references)-protocol['maximumRecallDeficitFromWorstReference'])
result=dict(status='passed' if passed else 'reference-quality-gap',cells=n,dimensions=o['dimensions'],epochs=o['epochs'],allFrozenLegacyResultBytesExact=True,resultSHA256=sha(r/'result.json').hex(),graphSHA256=sha(root/'graph.json').hex(),protocolSHA256=sha(a.protocol).hex(),curveA=e['curveA'],curveB=e['curveB'],curveSquaredError=e['curveSquaredError'],edgeVisits=e['edgeVisits'],attractiveUpdates=positive_count,negativeSamples=negative_count,qualityQueries=len(rows),candidateCells=n,qualityScope='complete' if len(rows)==n else 'fixed sampled queries against all cells',initial=metrics['initial'],native=native,references=references,versions={name:importlib.metadata.version(name) for name in ['numpy','scipy','umap-learn']},qualification='Full-cohort trajectories and exact schedule counts; neighborhood preservation uses the declared query scope. Not optimizer convergence, density preservation, cell labels, donor integration or biological qualification')
(a.out/'checks.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result),flush=True)
assert passed,'Predeclared reference quality gate failed; retain and investigate'
