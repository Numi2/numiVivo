#!/usr/bin/env python3
"""UMAP curve/optimizer reference and chunked neighborhood preservation checks."""
import argparse,hashlib,importlib.metadata,json
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from umap.umap_ import find_ab_params,make_epochs_per_sample
from umap.layouts import optimize_layout_euclidean
p=argparse.ArgumentParser()
p.add_argument('--report',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args()
r=json.loads(a.report.read_text()); e=r['embedding']; g=r['neighbors']; o=e['options']
x=np.asarray(r['integration']['scores'] if g['options'].get('representation') == 'integrated' else r['reduction']['scores']); y=np.asarray(e['coordinates']); initial=np.asarray(e['initialCoordinates'])
n=len(x); k=15
assert e['cells']==g['cells']==r['reduction']['cells']
assert y.shape==initial.shape==(n,o['dimensions']) and np.isfinite(y).all()
curve_a,curve_b=find_ab_params(o['spread'],o['minimumDistance'])
np.testing.assert_allclose([e['curveA'],e['curveB']],[curve_a,curve_b],rtol=2e-5,atol=2e-5)
xx=np.linspace(0,3*o['spread'],300);target=np.where(xx<o['minimumDistance'],1,np.exp(-(xx-o['minimumDistance'])/o['spread']))
curve=1/(1+e['curveA']*xx**(2*e['curveB']))
np.testing.assert_allclose(e['curveSquaredError'],np.sum((curve-target)**2),rtol=1e-10,atol=1e-12)
weights=np.asarray(g['weights']);head=np.repeat(np.arange(n),np.diff(g['rowOffsets']));tail=np.asarray(g['columnIndices'])
keep=weights>=weights.max()/o['epochs'];weights=weights[keep];head=head[keep];tail=tail[keep]
assert len(weights)==e['retainedDirectedEdges']
assert len(g['weights'])-len(weights)==e['discardedDirectedEdges']
assert e['edgeVisits']==len(weights)*o['epochs'] and e['edgeVisits']<=o['maximumUpdates']
intervals=make_epochs_per_sample(weights,o['epochs'])
# Reconstruct the schedule independently, including floating-point increments.
next_positive=intervals.copy(); ni=intervals/o['negativeSampleRate'];next_negative=ni.copy();positive_count=negative_count=0
for epoch in range(o['epochs']):
    active=np.flatnonzero(next_positive<=epoch)
    samples=((epoch-next_negative[active])/ni[active]).astype(np.int64)
    positive_count+=len(active);negative_count+=int(samples.sum())
    next_positive[active]+=intervals[active];next_negative[active]+=samples*ni[active]
assert positive_count==e['attractiveUpdates'] and negative_count==e['negativeSamples']
assert positive_count+negative_count<=o['maximumUpdates']

def quality(embedding):
    trust_penalty=continuity_penalty=overlap=0
    for start in range(0,n,64):
        high=cdist(x[start:start+64],x);low=cdist(embedding[start:start+64],embedding)
        for slot in range(len(high)):
            i=start+slot
            hi=np.lexsort((np.arange(n),high[slot]));lo=np.lexsort((np.arange(n),low[slot]))
            hi=hi[hi!=i];lo=lo[lo!=i]
            hr=np.empty(n,dtype=int);lr=np.empty(n,dtype=int)
            hr[hi]=np.arange(1,n);lr[lo]=np.arange(1,n)
            trust_penalty+=int(np.maximum(hr[lo[:k]]-k,0).sum())
            continuity_penalty+=int(np.maximum(lr[hi[:k]]-k,0).sum())
            overlap+=len(np.intersect1d(hi[:k],lo[:k]))
    denominator=n*k*(2*n-3*k-1)
    return dict(trustworthiness=1-2*trust_penalty/denominator,continuity=1-2*continuity_penalty/denominator,
                neighborRecall=overlap/(n*k))
native_quality=quality(y);initial_quality=quality(initial);references=[]
for seed in [7,19,41]:
    coordinates=np.ascontiguousarray(initial,dtype=np.float32)
    rng=np.random.RandomState(seed).randint(-(2**31)+1,2**31-1,3).astype(np.int64)
    result=optimize_layout_euclidean(coordinates,coordinates,head,tail,o['epochs'],n,intervals,e['curveA'],e['curveB'],rng,
         gamma=o['repulsionStrength'],initial_alpha=o['learningRate'],negative_sample_rate=o['negativeSampleRate'],
         parallel=False,move_other=True)
    references.append(dict(seed=seed,**quality(result)))
# Predeclared engineering gates on these datasets, not biological acceptance.
passed=(native_quality['trustworthiness']>=0.90 and
        native_quality['trustworthiness']>=min(v['trustworthiness'] for v in references)-0.02 and
        native_quality['neighborRecall']>=min(v['neighborRecall'] for v in references)-0.05)
summary=dict(status='passed' if passed else 'reference-quality-gap',reportSHA256=hashlib.sha256(a.report.read_bytes()).hexdigest(),
    cells=n,dimensions=o['dimensions'],epochs=o['epochs'],evaluationNeighbors=k,curveA=e['curveA'],curveB=e['curveB'],
    curveSquaredError=e['curveSquaredError'],attractiveUpdates=positive_count,negativeSamples=negative_count,
    initial=initial_quality,native=native_quality,references=references,
    qualification='Numerical and descriptive neighborhood preservation checks; not optimizer convergence, density preservation, biological labels or donor integration',
    versions={name:importlib.metadata.version(name) for name in ['numpy','scipy','umap-learn']})
a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
print(json.dumps(summary,indent=2));assert passed,'Native neighborhood preservation fails the predeclared reference gate'
