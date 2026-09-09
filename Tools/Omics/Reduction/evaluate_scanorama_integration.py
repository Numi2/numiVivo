#!/usr/bin/env python3
"""Independent Scanpy/Scanorama reference on frozen full-cohort native PCA."""
import argparse,hashlib,json,time
from importlib.metadata import version
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import scanorama
import scanpy.external as sce
from check_full_integration import matrix,metrics as donor_metrics
from evaluate_ding_integration import metrics as ding_metrics

def save(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cohort',choices=['kang','hagai','ding'],required=True);p.add_argument('--scale',choices=['raw','median-norm'],default='raw')
for name in ['inputs','native','protocol','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);protocol=json.loads(a.protocol.read_text());assert protocol['status'] in ['declared-before-scanorama-results','declared-before-scale-audit-results']
assert version('scanorama')==protocol['reference']['scanorama']
assert (a.scale=='raw')==(protocol['status']=='declared-before-scanorama-results')
if a.cohort=='ding':
 inp=np.load(a.inputs/'inputs.npz',allow_pickle=False);x=inp['scores'];batch=inp['methods'];baseline=json.loads((a.inputs/'baseline.json').read_text());types=inp['types'];known=inp['known'];input_path=a.inputs/'inputs.npz'
else:
 inp=np.load(a.inputs/'evaluation-inputs.npz',allow_pickle=False);x=inp['scores'];batch=inp['donors'];baseline=json.loads((a.inputs/'checks.json').read_text())['baseline'];types=inp['cellTypes'];known=np.ones(len(x),dtype=bool);input_path=a.inputs/'evaluation-inputs.npz'
np.testing.assert_array_equal(x,matrix(a.native/'pca/scores.bin',*x.shape))
metadata=json.loads((a.native/'pca/metadata.json').read_text());sample={s['id']:s for s in metadata['samples']};assert [sample[c['sampleID']]['batchID' if a.cohort=='ding' else 'donorID'] for c in metadata['cells']]==batch.tolist()
levels=np.unique(batch);rows=[np.flatnonzero(batch==b) for b in levels];order=np.concatenate(rows);assert len(np.unique(order))==len(x);scale=float(np.median(np.linalg.norm(x,axis=1))) if a.scale=='median-norm' else 1.0;assert np.isfinite(scale) and scale>0
normalized=x/scale;datasets=[normalized[ix].copy() for ix in rows]
start=time.perf_counter();alignments,matches=scanorama.find_alignments(datasets,knn=20,approx=False,alpha=.1,verbose=0)
# Preserve exact directed global-index pairs before any label-based evaluation.
np.savez_compressed(a.out/'anchors.npz',**{f'{i}-{j}':np.array([(rows[i][s],rows[j][t]) for s,t in sorted(v)],dtype=np.int64).reshape(-1,2) for (i,j),v in matches.items()})
save(a.out/'alignment-order.json',[dict(source=int(i),reference=int(j),sourceBatch=str(levels[i]),referenceBatch=str(levels[j])) for i,j in alignments])
# The official Scanpy wrapper operates on PCA in obsm. A zero-column sparse X
# avoids duplicating counts; no counts or source labels are used for alignment.
adata=ad.AnnData(X=csr_matrix((len(x),0)),obs=pd.DataFrame({'batch':batch[order]},index=[str(i) for i in order]))
adata.obsm['X_pca']=normalized[order].copy()
sce.pp.scanorama_integrate(adata,'batch',basis='X_pca',adjusted_basis='X_scanorama',knn=20,sigma=15,approx=False,alpha=.1,batch_size=256,verbose=1,alignments=alignments,matches=matches)
y=np.empty_like(x);y[order]=adata.obsm['X_scanorama']*scale;assert np.isfinite(y).all() and y.shape==x.shape
seconds=time.perf_counter()-start;np.savez_compressed(a.out/'scores.npz',scores=y,sourceRow=np.arange(len(x)),groupedSourceRow=order)
print(a.cohort+' Scanorama corrected',flush=True)
# Diagnostics use source labels only after correction; missing labels are not filled.
anchor_diagnostics=[]
for (i,j),v in matches.items():
 pairs=np.array([(rows[i][s],rows[j][t]) for s,t in sorted(v)],dtype=np.int64).reshape(-1,2)
 available=known[pairs[:,0]]&known[pairs[:,1]]
 anchor_diagnostics.append(dict(sourceBatch=str(levels[i]),referenceBatch=str(levels[j]),pairs=len(pairs),sourceAssignedPairs=int(available.sum()),sameSourceTypeFraction=float(np.mean(types[pairs[available,0]]==types[pairs[available,1]])) if available.any() else None,selectedForAssembly=[i,j] in [list(p) for p in alignments]))
save(a.out/'anchor-diagnostics.json',anchor_diagnostics)
margin=protocol['gateMargins']
if a.cohort=='ding':
 m=ding_metrics(y,inp,30,4,a.out,'scanorama')
 g=dict(mixingImproved=m['sameMethodExcess']<baseline['sameMethodExcess'],cellTypesPreserved=m['cellTypeBalancedAccuracy']>=baseline['cellTypeBalancedAccuracy']-margin['maximumCellTypeBalancedAccuracyLoss'],everyTypeRecallPreserved=all(m['cellTypeRecall'][t]>=v-margin['maximumPerTypeRecallLoss'] for t,v in baseline['cellTypeRecall'].items()),programsPreserved=all(m['programs'][p]['spearman']>=v['spearman']-margin['maximumProgramSpearmanLoss'] for p,v in baseline['programs'].items()),withinStratumProgramsPreserved=all(m['programs'][p]['withinStratumSpearman']>=v['withinStratumSpearman']-margin['maximumProgramSpearmanLoss'] for p,v in baseline['programs'].items()))
else:
 m=donor_metrics(y,inp['donors'],inp['conditions'],inp['cellTypes'],inp['program'],30,a.out/'scanorama-neighbors.npz')
 g=dict(mixingImproved=m['sameDonorExcess']<baseline['sameDonorExcess'],conditionPreserved=m['conditionBalancedAccuracy']>=baseline['conditionBalancedAccuracy']-margin['maximumConditionBalancedAccuracyLoss'],programPreserved=m['programSpearman']>=baseline['programSpearman']-margin['maximumProgramSpearmanLoss'],withinStratumProgramPreserved=m['withinStratumProgramSpearman']>=baseline['withinStratumProgramSpearman']-margin['maximumProgramSpearmanLoss'],completeClassifierStrata=not m['missingClassifierStrata'])
 if baseline['cellTypeBalancedAccuracy'] is not None:g.update(cellTypesPreserved=m['cellTypeBalancedAccuracy']>=baseline['cellTypeBalancedAccuracy']-margin['maximumCellTypeBalancedAccuracyLoss'],everyTypeRecallPreserved=all(m['cellTypeRecall'][t]>=v-margin['maximumPerTypeRecallLoss'] for t,v in baseline['cellTypeRecall'].items()))
save(a.out/'checks.json',dict(status='full-cohort-measured',cohort=a.cohort,protocolSHA256=sha(a.protocol),inputSHA256=sha(input_path),metadataSHA256=sha(a.native/'pca/metadata.json'),scoresSHA256=sha(a.out/'scores.npz'),cells=len(x),components=x.shape[1],versions={n:version(n) for n in ['scanorama','scanpy','numpy','scipy','anndata','annoy','fbpca','geosketch','intervaltree','sortedcontainers']},integrationSeconds=seconds,inputScale=a.scale,scaleValue=scale,baseline=baseline,metrics=m,gates=g,allMeasuredGatesPassed=all(g.values()),qualification=protocol['qualification']))
print(a.cohort+' measured '+json.dumps(g),flush=True)
