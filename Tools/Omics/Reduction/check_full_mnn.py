#!/usr/bin/env python3
"""Verify native full-cohort MNN against frozen independent Scanpy/Scanorama outputs."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from check_full_integration import matrix,metrics as donor_metrics
from evaluate_ding_integration import metrics as ding_metrics
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cohort',choices=['kang','hagai','ding'],required=True)
for name in ['native','inputs','reference','protocol','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,v):(a.out/name).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
prior=json.loads((a.reference/'checks.json').read_text());assert prior['status']=='full-cohort-measured' and prior['inputScale']=='median-norm'
path=a.inputs/('inputs.npz' if a.cohort=='ding' else 'evaluation-inputs.npz');assert sha(path)==prior['inputSHA256'];inputs=np.load(path,allow_pickle=False);x=inputs['scores'];batch=inputs['methods' if a.cohort=='ding' else 'donors'];baseline=prior['baseline'];n,d=x.shape
np.testing.assert_array_equal(x,matrix(a.native/'pca/scores.bin',n,d));assert sha(a.native/'pca/metadata.json')==prior['metadataSHA256']
root=a.native/'mnn';metadata=json.loads((root/'metadata.json').read_text());assert metadata==json.loads((a.native/'pca/metadata.json').read_text());sm={s['id']:s for s in metadata['samples']};assert [sm[c['sampleID']]['batchID' if a.cohort=='ding' else 'donorID'] for c in metadata['cells']]==batch.tolist()
plan=json.loads((root/'plan.json').read_text());assert 'integration' not in plan and plan['inputKind']=='fitted';o=plan['mnn'];assert o['neighbors']==20 and o['sigma']==15 and o['minimumAlignment']==.1
report=json.loads((root/'report.json').read_text());assert report['method']=='mutual-nearest-neighbor-panorama-median-norm-Double-v1'
levels=np.unique(batch);assert report['levels']==levels.tolist();np.testing.assert_array_equal(report['cellLevels'],np.searchsorted(levels,batch));assert abs(report['medianRowNorm']-prior['scaleValue'])<1e-12
frozen=np.load(a.reference/'scores.npz');assert sha(a.reference/'scores.npz')==prior['scoresSHA256'];y=matrix(root/'scores.bin',n,d);np.testing.assert_allclose(y,frozen['scores'],rtol=1e-9,atol=1e-8)
anchors=np.fromfile(root/'anchors.bin',dtype=[('first','<u4'),('second','<u4'),('distance','<f8')]);reference_anchors=np.load(a.reference/'anchors.npz');assert len(anchors)==sum(v['anchors'] for v in report['alignments']);assert np.isfinite(anchors['distance']).all()
norms=np.linalg.norm(x,axis=1)[:,None];unit=np.divide(x,norms,out=np.zeros_like(x),where=norms>0);maximum_distance_error=0.;offset=0
assert len(anchors)==sum(len(reference_anchors[key]) for key in reference_anchors.files)
assert {f"{v['firstLevel']}-{v['secondLevel']}" for v in report['alignments']}=={key for key in reference_anchors.files if len(reference_anchors[key])}
for pair in report['alignments']:
 assert pair['anchorOffset']==offset;offset+=pair['anchors'];i,j=pair['firstLevel'],pair['secondLevel'];a1=anchors[pair['anchorOffset']:offset];actual=np.stack([a1['first'],a1['second']],axis=1).astype(np.int64);expected=reference_anchors[f'{i}-{j}'];np.testing.assert_array_equal(actual,expected)
 first,second=actual[:,0],actual[:,1];assert (batch[first]==levels[i]).all() and (batch[second]==levels[j]).all()
 distances=np.sum(np.abs(unit[first]-unit[second]),axis=1);error=float(np.max(np.abs(distances-a1['distance'])));maximum_distance_error=max(maximum_distance_error,error);assert error<1e-12
 assert pair['firstMatchedCells']==len(np.unique(first)) and pair['secondMatchedCells']==len(np.unique(second));score=max(pair['firstMatchedCells']/int(np.sum(batch==levels[i])),pair['secondMatchedCells']/int(np.sum(batch==levels[j])));assert score==pair['score']
expected_order=[(x['source'],x['reference']) for x in json.loads((a.reference/'alignment-order.json').read_text())];actual_order=[(report['alignments'][i]['firstLevel'],report['alignments'][i]['secondLevel']) for i in report['assemblyOrder']];assert actual_order==expected_order
assert report['distanceScalarTerms']==sum(int(np.sum(batch==i))*int(np.sum(batch==j))*d for bi,i in enumerate(levels) for j in levels[bi+1:]);assert report['kernelScalarTerms']==sum(step['correctedCells']*step['anchors']*d for step in report['steps'])
assert report['distanceScalarTerms']+report['kernelScalarTerms']<=o['maximumWork'] and report['estimatedMaximumLatentResidentBytes']<=o['maximumResidentBytes']
math=dict(originalPCABitsAndMetadataExact=True,allAnchorsAndAlignmentOrderExact=True,maximumAnchorDistanceError=maximum_distance_error,maximumCoordinateError=float(np.max(np.abs(y-frozen['scores']))),medianRowNorm=report['medianRowNorm'],workAccountingExact=True,zeroWeightCellsByStep=[s['zeroWeightCells'] for s in report['steps']],referenceChecksSHA256=sha(a.reference/'checks.json'))
save('numerical-checks.json',dict(status='passed',**math));print(a.cohort+' native/reference numerical checks passed',flush=True)
margin=json.loads(a.protocol.read_text())['gateMargins'];start=time.perf_counter()
if a.cohort=='ding':
 m=ding_metrics(y,inputs,30,4,a.out,'native');g=dict(mixingImproved=m['sameMethodExcess']<baseline['sameMethodExcess'],cellTypesPreserved=m['cellTypeBalancedAccuracy']>=baseline['cellTypeBalancedAccuracy']-margin['maximumCellTypeBalancedAccuracyLoss'],everyTypeRecallPreserved=all(m['cellTypeRecall'][t]>=v-margin['maximumPerTypeRecallLoss'] for t,v in baseline['cellTypeRecall'].items()),programsPreserved=all(m['programs'][p]['spearman']>=v['spearman']-margin['maximumProgramSpearmanLoss'] for p,v in baseline['programs'].items()),withinStratumProgramsPreserved=all(m['programs'][p]['withinStratumSpearman']>=v['withinStratumSpearman']-margin['maximumProgramSpearmanLoss'] for p,v in baseline['programs'].items()))
else:
 m=donor_metrics(y,inputs['donors'],inputs['conditions'],inputs['cellTypes'],inputs['program'],30,a.out/'native-neighbors.npz');g=dict(mixingImproved=m['sameDonorExcess']<baseline['sameDonorExcess'],conditionPreserved=m['conditionBalancedAccuracy']>=baseline['conditionBalancedAccuracy']-margin['maximumConditionBalancedAccuracyLoss'],programPreserved=m['programSpearman']>=baseline['programSpearman']-margin['maximumProgramSpearmanLoss'],withinStratumProgramPreserved=m['withinStratumProgramSpearman']>=baseline['withinStratumProgramSpearman']-margin['maximumProgramSpearmanLoss'],completeClassifierStrata=not m['missingClassifierStrata'])
 if baseline['cellTypeBalancedAccuracy'] is not None:g.update(cellTypesPreserved=m['cellTypeBalancedAccuracy']>=baseline['cellTypeBalancedAccuracy']-margin['maximumCellTypeBalancedAccuracyLoss'],everyTypeRecallPreserved=all(m['cellTypeRecall'][t]>=v-margin['maximumPerTypeRecallLoss'] for t,v in baseline['cellTypeRecall'].items()))
save('checks.json',dict(status='full-cohort-native-measured',cohort=a.cohort,cells=n,components=d,numerical=math,metrics=m,gates=g,allMeasuredGatesPassed=all(g.values()),evaluationSeconds=time.perf_counter()-start,qualification='Native implementation and frozen independently invoked Scanpy/Scanorama reference on the same already inspected full cohorts. Original margins retained; incomplete Kang classifier strata and partial Ding source annotations remain. Transductive correction and source-label metrics are not unseen-study or prospective biological qualification.'))
print(a.cohort+' native preservation '+json.dumps(g),flush=True)
